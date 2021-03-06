/*
 * Copyright (c) 2012 David Rodrigues
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/*
 * Taken from https://github.com/dmcrodrigues/macro-logger
 * To use, compile with LOG_LEVEL=level_number, where allowed level_numbers are:
 * 0 NO_LOGS
 * 1 ERROR
 * 2 INFO
 * 3 DEBUG
 */

#ifndef __MACROLOGGER_H__
#define __MACROLOGGER_H__

#ifdef __OBJC__
#import <Foundation/Foundation.h>
#else
#include <time.h>
#include <string.h>
#endif

// === auxiliar functions
static inline char *timenow();

#define _FILE strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__

#define NO_LOG          0x00
#define ERROR_LEVEL     0x01
#define INFO_LEVEL      0x02
#define DEBUG_LEVEL     0x03

#ifndef LOG_LEVEL
#define LOG_LEVEL   DEBUG_LEVEL
#endif

#ifdef __OBJC__

#if __has_feature(objc_arc)
#define AUTORELEASEPOOL_BEGIN   @autoreleasepool {
#define AUTORELEASEPOOL_END     }
#define RELEASE(OBJ)            OBJ = nil
#else
#define AUTORELEASEPOOL_BEGIN   NSAutoreleasePool *_pool = [[NSAutoreleasePool alloc] init];
#define AUTORELEASEPOOL_END     [_pool release];
#define RELEASE(OBJ)            [OBJ release];
#endif

#define PRINTFUNCTION(format, ...)      objc_print(@format, __VA_ARGS__)
#define STDOUTPRINTFUNCTION(format, ...)      objc_print(@format, __VA_ARGS__)
#else
#define PRINTFUNCTION(format, ...)      fprintf(stderr, format, __VA_ARGS__)
#define STDOUTPRINTFUNCTION(format, ...)      printf(format, __VA_ARGS__)

#endif

/* #define LOG_FMT             "%s | %s | %s | %s:%d | " */
#define LOG_ARGS(LOG_TAG)   timenow(), LOG_TAG, _FILE, __FUNCTION__, __LINE__
#define LOG_FMT ""

#define NEWLINE     "\n"

#define ERROR_TAG   "ERROR"
#define INFO_TAG    "INFO"
#define DEBUG_TAG   "DEBUG"
#define STDOUT_TAG   "STDOUT"

#if LOG_LEVEL >= DEBUG_LEVEL
/* #define LOG_DEBUG(message, args...)     PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(DEBUG_TAG), ## args) */
#define LOG_DEBUG(message, args...)     PRINTFUNCTION(LOG_FMT message NEWLINE, ## args)
#define LOG_DEBUG_MASTER(message, args...) if (rank == 0) PRINTFUNCTION(LOG_FMT message NEWLINE, ## args)
#else
#define LOG_DEBUG(message, args...)
#define LOG_DEBUG_MASTER(message, args...)
#endif

#if LOG_LEVEL >= INFO_LEVEL
#define LOG_INFO(message, args...)      PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(INFO_TAG), ## args)
#define LOG_INFO_MASTER(message, args...) if (rank == 0) PRINTFUNCTION(LOG_FMT message NEWLINE, ## args)
#else
#define LOG_INFO(message, args...)
#define LOG_DEBUG_MASTER(message, args...)
#endif

#if LOG_LEVEL >= ERROR_LEVEL
#define LOG_ERROR(message, args...)     PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(ERROR_TAG), ## args)
#else
#define LOG_ERROR(message, args...)
#endif

#if LOG_LEVEL >= NO_LOGS
#define LOG_IF_ERROR(condition, message, args...) if (condition) PRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(ERROR_TAG), ## args)
#else
#define LOG_IF_ERROR(condition, message, args...)
#endif

/* Always log to stdout using LOG */
#define LOG_STDOUT(message, args...) STDOUTPRINTFUNCTION(LOG_FMT message NEWLINE, LOG_ARGS(STDOUT_TAG), ## args)

static inline char *timenow() {
    static char buffer[64];
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 64, "%Y-%m-%d %H:%M:%S", timeinfo);

    return buffer;
}

#ifdef __OBJC__

static inline void objc_print(NSString *format, ...) {
    AUTORELEASEPOOL_BEGIN
    va_list args;
    va_start(args, format);
    NSString *logStr = [[NSString alloc] initWithFormat:format arguments:args];
    fprintf(stderr, "%s", [logStr UTF8String]);
    RELEASE(logStr);
    va_end(args);
    AUTORELEASEPOOL_END
}

#endif

#endif
