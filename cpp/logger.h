/**
 * DO NOT REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Contributor(s):
 *
 * The Original Software is OpenRedukti (https://github.com/redukti/OpenRedukti).
 * The Initial Developer of the Original Software is REDUKTI LIMITED (http://redukti.com).
 * Authors: Dibyendu Majumdar
 *
 * Copyright 2017 REDUKTI LIMITED. All Rights Reserved.
 *
 * The contents of this file are subject to the the GNU General Public License
 * Version 3 (https://www.gnu.org/licenses/gpl.txt).
 */

#ifndef _REDUKTI_LOGGER_H
#define _REDUKTI_LOGGER_H

#include <stdarg.h>
#include <stdio.h>

namespace redukti
{

#define REDUKTI_DEBUG_ENABLED 1

extern unsigned log_mask;

enum LogLevel { LOG_DEBUG = 1, LOG_WARN = 2, LOG_INFO = 4, LOG_ERROR = 8, LOG_FATAL_ERROR = 16 };

extern void log_message(int LogLevel, const char *filename, int line_number, const char *function, FILE *file,
			const char *format, ...);

#if REDUKTI_DEBUG_ENABLED
#if _MSC_VER
#define debug(format, ...)                                                                                             \
	((void)((log_mask & LOG_DEBUG) &&                                                                              \
		(log_message(LOG_DEBUG, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__), 0)))
#define warn(format, ...) log_message(LOG_WARN, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#define inform(format, ...) log_message(LOG_INFO, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#define error(format, ...) log_message(LOG_ERROR, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#define die(format, ...) log_message(LOG_FATAL_ERROR, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#else
#define debug(format, ...)                                                                                             \
	((void)((log_mask & LOG_DEBUG) &&                                                                              \
		(log_message(LOG_DEBUG, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__), 0)))
#define warn(format, ...) log_message(LOG_WARN, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#define inform(format, ...) log_message(LOG_INFO, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#define error(format, ...) log_message(LOG_ERROR, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#define die(format, ...) log_message(LOG_FATAL_ERROR, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#endif
#else
#define debug(format, ...) ((void)0)
#define warn(format, ...) log_message(LOG_WARN, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#define inform(format, ...) log_message(LOG_INFO, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#define error(format, ...) log_message(LOG_ERROR, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#define die(format, ...) log_message(LOG_FATAL_ERROR, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#endif
} // namespace redukti

#endif