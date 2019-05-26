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
#include <logger.h>

#include <assert.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <sstream>
#include <thread>

unsigned volatile Redukti_log_mask = LOG_INFO | LOG_WARN | LOG_ERROR | LOG_FATAL_ERROR;

void redukti_log_message(int log_level, const char *filename, int line_number, const char *function, FILE *file,
			 const char *format, ...)
{
	assert(file);

	if (!(log_level & Redukti_log_mask))
		return;

	switch (log_level) {
	case LOG_TRACE:
		fputs("TRACE ", file);
		break;
	case LOG_DEBUG:
		fputs("DEBUG ", file);
		break;
	case LOG_WARN:
		fputs("WARN ", file);
		break;
	case LOG_INFO:
		fputs("INFO ", file);
		break;
	case LOG_FATAL_ERROR:
		fputs("FATAL ERROR ", file);
		break;
	default:
		fputs("ERROR ", file);
		break;
	}

	std::stringbuf buf;
	std::ostream os(&buf);

	os << std::this_thread::get_id();
	fprintf(file, "tid(%s) ", buf.str().c_str());

	if (filename && function) {
		fprintf(file, "%s:%d (%s) ", filename, line_number, function);
	}
	va_list args;
	va_start(args, format);
	vfprintf(file, format, args);
	va_end(args);
	if (log_level == LOG_FATAL_ERROR)
		exit(1);
}
