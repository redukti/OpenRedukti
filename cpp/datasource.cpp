/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017-2019 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#include <datasource.h>
#include <logger.h>

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <stdexcept>

namespace redukti
{

bool parse_delimited(const char *input_start, const char *input_end, std::vector<const char *> &out_tokens,
		     buffer_ptr<char> &buf, const char *delimiters) noexcept
{
	out_tokens.clear();

	if (input_end <= input_start) {
		return false;
	}
	auto input_size = input_end - input_start;
	if (input_size > (int)buf.size()) {
		error("Input size %d is bigger than buffer size %d\n", (int)input_size, (int)buf.size());
		return false;
	}
	const char *input_ptr = input_start;
	char *wordp = buf.begin();

	while (*input_ptr && input_ptr != input_end) {
		char *word = wordp;
		*wordp = 0;

		bool inquote = false;
		while (*input_ptr && input_ptr != input_end) {
			if (word == wordp) {
				// we are at the beginning for a word, so look
				// for potential quote
				if (*input_ptr == '"' && !inquote) {
					// We are in a quoted word
					inquote = true;
					input_ptr++;
					continue;
				}
			}
			if (inquote) {
				// We are in a quoted word
				if (*input_ptr == '"') {
					// Check if it is an escape - i.e.
					// double quote
					if (input_ptr + 1 < input_end && *(input_ptr + 1) == '"') {
						// escape so we add the quote
						// character
						*wordp++ = '"';
						input_ptr += 2;
						continue;
					} else {
						// not escape so the quoted word
						// ends here
						inquote = false;
						*wordp++ = 0;
						input_ptr++;
						if (input_ptr < input_end &&
						    (*input_ptr == ',' || *input_ptr == '\t' ||
						     (delimiters && strchr(delimiters, *input_ptr)))) {
							// Skip delimiter
							// following quote
							input_ptr++;
						}
						break;
					}
				} else {
					// still in quoted word
					*wordp++ = *input_ptr++;
					continue;
				}
			} else {
				// Not in quoted word
				if (*input_ptr == ',' || *input_ptr == '\t' ||
				    (delimiters && strchr(delimiters, *input_ptr))) {
					// word ends due to delimiter
					*wordp++ = 0;
					input_ptr++;
					break;
				} else if (*input_ptr == '\r' || *input_ptr == '\n') {
					// skip line feed or CRLF
					*wordp++ = 0;
					if (*input_ptr == '\r' && input_ptr + 1 < input_end &&
					    *(input_ptr + 1) == '\n') {
						input_ptr++;
					}
					input_ptr++;
					break;
				} else {
					*wordp++ = *input_ptr++;
				}
			}
		}
		out_tokens.push_back(word);
	}
	return true;
}

class CSVDataSourceImpl
{
	private:
	FILE *fp;
	std::string _name;
	bool _hasHeadings;
	bool _headingsRead;
	buffer_ptr<char> _line_buf;

	public:
	CSVDataSourceImpl(const std::string &name, bool hasHeadings) noexcept;
	~CSVDataSourceImpl() noexcept;
	bool next(Line *line) noexcept;
	bool read_next_line() noexcept;
	bool is_valid() const noexcept { return fp != nullptr; }

	private:
	CSVDataSourceImpl(const CSVDataSourceImpl &) = delete;
	CSVDataSourceImpl &operator=(const CSVDataSourceImpl &) = delete;
};

CSVDataSourceImpl::CSVDataSourceImpl(const std::string &name, bool hasHeadings) noexcept
    : _name(name), _hasHeadings(hasHeadings), _headingsRead(false), _line_buf(64 * 1024)
{
	fp = fopen(name.c_str(), "r");
	if (fp == nullptr) {
		error("Unable to open file %s\n", name.c_str());
	}
}

CSVDataSourceImpl::~CSVDataSourceImpl() noexcept
{
	if (fp != nullptr)
		fclose(fp);
	fp = nullptr;
}

bool CSVDataSourceImpl::read_next_line() noexcept
{
	char *p = fgets(_line_buf.begin(), _line_buf.size(), fp);
	if (p) {
		auto len = strlen(p);
		if (len == 0) {
			error("Line is empty\n");
			return false;
		}
		if (p[len - 1] != '\n') {
			error("Line exceeds the size of the buffer or is missing "
			      "newline\n");
			return false;
		}
	}
	if (p == nullptr) {
		return false;
	}
	return true;
}

bool CSVDataSourceImpl::next(Line *line) noexcept
{
	if (fp == nullptr)
		return false;
	if (_hasHeadings) {
		if (!_headingsRead) {
			if (!read_next_line())
				return false;
			parse_delimited(_line_buf.begin(), _line_buf.begin() + strlen(_line_buf.begin()),
					line->headings, line->heading_buf);
			_headingsRead = true;
		}
	} else {
		line->headings.clear();
	}
	if (!read_next_line())
		return false;
	parse_delimited(_line_buf.begin(), _line_buf.begin() + strlen(_line_buf.begin()) + 1, line->fields,
			line->line_buf);
	return true;
}

CSVDataSource::CSVDataSource(std::string filename, bool hasHeadings) noexcept
{
	impl = std::unique_ptr<CSVDataSourceImpl>(new CSVDataSourceImpl(filename, hasHeadings));
}

CSVDataSource::~CSVDataSource() noexcept {}

bool CSVDataSource::next(Line *line) noexcept { return impl->next(line); }

bool CSVDataSource::is_valid() const noexcept { return impl->is_valid(); }

//////////////////////////////// TESTS

int test_buffer_ptr()
{
	int failure_count = 0;
	buffer_ptr<char> buf(100);
	if (buf.size() != 100) {
		failure_count++;
	}
	auto start = std::begin(buf);
	auto end = std::end(buf);
	if ((end - start) != (int)(buf.size() * sizeof(char))) {
		failure_count++;
	}
	if (buf.ref_count() != 1) {
		failure_count++;
	}
	{
		auto buf_copy = buf;
		if (buf.ref_count() != 2) {
			std::cerr << "unexpected reference count - expected 2 - got " << buf.ref_count() << "\n";
			failure_count++;
		}
		buf_copy = buf;
		if (buf.ref_count() != 2) {
			std::cerr << "unexpected reference count - expected 2 - got " << buf.ref_count() << "\n";
			failure_count++;
		}
	}
	if (buf.ref_count() != 1) {
		failure_count++;
	}
	buf[4] = 'a';
	if (buf[4] != 'a') {
		std::cerr << "buf[4] = " << buf[4] << "\n";
		failure_count++;
	}
	std::fill(std::begin(buf), std::end(buf), 'z');
	auto n = std::count(std::begin(buf), std::end(buf), 'z');
	if (n != 100) {
		failure_count++;
	}
	return failure_count;
}

int test_ds()
{
	int failure_count = 0;
	char data[] = {"This,is,test,,data,\"embedded,\"\"data\",final\n"};
	const char *expected[] = {"This", "is", "test", "", "data", "embedded,\"data", "final"};

	buffer_ptr<char> buf(1024);
	std::vector<const char *> tokens;
	parse_delimited(std::begin(data), std::end(data), tokens, buf);
	int n = 0;
	for (auto word : tokens) {
		if (strcmp(word, expected[n]) != 0) {
			fprintf(stderr, "tok = [%s], expected = [%s]\n", word, expected[n]);
			failure_count++;
		}
		n++;
	}
	return failure_count;
}

int test_datasource()
{
	int failure_count = 0;
	failure_count += test_buffer_ptr();
	failure_count += test_ds();
	if (failure_count == 0)
		printf("DataSource Tests OK\n");
	else
		printf("DataSource Tests FAILED\n");
	return failure_count;
}

} // namespace redukti
