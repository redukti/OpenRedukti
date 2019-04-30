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

#ifndef _REDUKTI_LOGGER_H
#define _REDUKTI_LOGGER_H

#include <stdarg.h>
#include <stdio.h>

#define REDUKTI_DEBUG_ENABLED 1

#ifdef __cplusplus
extern "C"
#else
extern
#endif
    volatile unsigned Redukti_log_mask;

enum LogLevel { LOG_TRACE = 1, LOG_DEBUG = 2, LOG_WARN = 4, LOG_INFO = 8, LOG_ERROR = 16, LOG_FATAL_ERROR = 32 };

// The logger implementation is very simple
// and not suited for high performance.
// Recommended usage is only for errors, warnings and debug
// messages. Info level messages should not be in performance
// critical sections.

#ifdef __cplusplus
extern "C"
#else
extern
#endif
    void
    redukti_log_message(int LogLevel, const char *filename, int line_number, const char *function, FILE *file,
			const char *format, ...);

#define is_debug_enabled() (Redukti_log_mask & LOG_DEBUG)
#define is_trace_enabled() (Redukti_log_mask & LOG_TRACE)

#if REDUKTI_DEBUG_ENABLED
#if _MSC_VER
#define trace(format, ...)                                                                                             \
	((void)(is_trace_enabled() &&                                                                                  \
		(redukti_log_message(LOG_TRACE, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__), 0)))
#define debug(format, ...)                                                                                             \
	((void)(is_debug_enabled() &&                                                                                  \
		(redukti_log_message(LOG_DEBUG, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__), 0)))
#define warn(format, ...) redukti_log_message(LOG_WARN, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#define inform(format, ...) redukti_log_message(LOG_INFO, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#define error(format, ...) redukti_log_message(LOG_ERROR, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#define die(format, ...) redukti_log_message(LOG_FATAL_ERROR, __FILE__, __LINE__, __func__, stderr, format, __VA_ARGS__)
#else
#define trace(format, ...)                                                                                             \
	((void)(is_trace_enabled() &&                                                                                  \
		(redukti_log_message(LOG_TRACE, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__), 0)))
#define debug(format, ...)                                                                                             \
	((void)(is_debug_enabled() &&                                                                                  \
		(redukti_log_message(LOG_DEBUG, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__), 0)))
#define warn(format, ...) redukti_log_message(LOG_WARN, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#define inform(format, ...) redukti_log_message(LOG_INFO, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#define error(format, ...) redukti_log_message(LOG_ERROR, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#define die(format, ...)                                                                                               \
	redukti_log_message(LOG_FATAL_ERROR, __FILE__, __LINE__, __func__, stderr, format, ##__VA_ARGS__)
#endif
#else
#define trace(format, ...) ((void)0)
#define debug(format, ...) ((void)0)
#define warn(format, ...) redukti_log_message(LOG_WARN, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#define inform(format, ...) redukti_log_message(LOG_INFO, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#define error(format, ...) redukti_log_message(LOG_ERROR, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#define die(format, ...) redukti_log_message(LOG_FATAL_ERROR, nullptr, __LINE__, nullptr, stderr, format, __VA_ARGS__)
#endif

#endif