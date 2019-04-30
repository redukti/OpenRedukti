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

#ifndef _REDUKTI_DATASOURCE_H_
#define _REDUKTI_DATASOURCE_H_

#include <buffer.h>

#include <memory>
#include <string>
#include <vector>

namespace redukti
{

/**
 * Tokenizer that takes as input a string, and returns
 * tokens. Token are identified using delimiters COMMA,
 * TAB, LF, or CRLF. Tokens can be surrounded in double quotes
 * which would then allow anything in quotes to be treated
 * as a token, the double quote character itself can be
 * escaped by using two consecutive double quotes.
 * Initialise the tokenizer with new input. Allows
 * reuse of the tokenizer. The start and end must point
 * to beginning and end of the input as per STL rules.
 */
extern bool parse_delimited(const char *start, const char *end, std::vector<const char *> &out_tokens,
			    buffer_ptr<char> &buf, const char *delims = nullptr) noexcept;

struct Line {
	/**
	 * Buffer that holds the headings line
	 */
	buffer_ptr<char> heading_buf;
	/**
	 * Buffer that holds the current line
	 */
	buffer_ptr<char> line_buf;
	/**
	 * The headings point to offsets in the heading_buf.
	 */
	std::vector<const char *> headings;
	/**
	 * The fields point to offsets in the line_buf
	 */
	std::vector<const char *> fields;

	Line(size_t n) : heading_buf(n), line_buf(n) {}
};

class CSVDataSourceImpl;
class CSVDataSource
{
	private:
	std::unique_ptr<CSVDataSourceImpl> impl;

	public:
	/**
	 * Constructor takes a file name and boolean flag to indicate
	 * that a heading line is expected to be the first line.
	 */
	CSVDataSource(std::string filename, bool hasHeadings) noexcept;
	~CSVDataSource() noexcept;
	bool is_valid() const noexcept;
	/**
	 * Retrieve the next line, and return false if EOF.
	 */
	bool next(Line *line) noexcept;

	template <typename EvalFunc> void eval(EvalFunc &func) noexcept
	{
		auto line = std::make_unique<Line>(1024);
		int linenum = 0;
		while (next(line.get())) {
			linenum++;
			func(line.get());
		}
	}
};

extern int test_datasource();

} // namespace redukti

#endif
