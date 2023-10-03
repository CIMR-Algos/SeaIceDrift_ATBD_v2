#ifndef _SRC_FMUTIL_CONFIG_H
#define _SRC_FMUTIL_CONFIG_H 1
 
/* src/fmutil_config.h. Generated automatically at end of configure. */
/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to 1 if you have the <inttypes.h> header file. */
#ifndef FMUTIL_HAVE_INTTYPES_H 
#define FMUTIL_HAVE_INTTYPES_H  1 
#endif

/* Define if proj is available */
#ifndef FMUTIL_HAVE_LIBPROJ 
#define FMUTIL_HAVE_LIBPROJ  1 
#endif

/* Define to 1 if you have the <math.h> header file. */
#ifndef FMUTIL_HAVE_MATH_H 
#define FMUTIL_HAVE_MATH_H  1 
#endif

/* Define to 1 if you have the <memory.h> header file. */
#ifndef FMUTIL_HAVE_MEMORY_H 
#define FMUTIL_HAVE_MEMORY_H  1 
#endif

/* Define to 1 if you have the <stdarg.h> header file. */
#ifndef FMUTIL_HAVE_STDARG_H 
#define FMUTIL_HAVE_STDARG_H  1 
#endif

/* Define to 1 if you have the <stdint.h> header file. */
#ifndef FMUTIL_HAVE_STDINT_H 
#define FMUTIL_HAVE_STDINT_H  1 
#endif

/* Define to 1 if you have the <stdio.h> header file. */
#ifndef FMUTIL_HAVE_STDIO_H 
#define FMUTIL_HAVE_STDIO_H  1 
#endif

/* Define to 1 if you have the <stdlib.h> header file. */
#ifndef FMUTIL_HAVE_STDLIB_H 
#define FMUTIL_HAVE_STDLIB_H  1 
#endif

/* Define to 1 if you have the <strings.h> header file. */
#ifndef FMUTIL_HAVE_STRINGS_H 
#define FMUTIL_HAVE_STRINGS_H  1 
#endif

/* Define to 1 if you have the <string.h> header file. */
#ifndef FMUTIL_HAVE_STRING_H 
#define FMUTIL_HAVE_STRING_H  1 
#endif

/* Define to 1 if you have the `strtof' function. */
#ifndef FMUTIL_HAVE_STRTOF 
#define FMUTIL_HAVE_STRTOF  1 
#endif

/* Define to 1 if you have the `strtol' function. */
#ifndef FMUTIL_HAVE_STRTOL 
#define FMUTIL_HAVE_STRTOL  1 
#endif

/* Define to 1 if you have the <sys/stat.h> header file. */
#ifndef FMUTIL_HAVE_SYS_STAT_H 
#define FMUTIL_HAVE_SYS_STAT_H  1 
#endif

/* Define to 1 if you have the <sys/types.h> header file. */
#ifndef FMUTIL_HAVE_SYS_TYPES_H 
#define FMUTIL_HAVE_SYS_TYPES_H  1 
#endif

/* Define to 1 if you have the <time.h> header file. */
#ifndef FMUTIL_HAVE_TIME_H 
#define FMUTIL_HAVE_TIME_H  1 
#endif

/* Define to 1 if you have the <unistd.h> header file. */
#ifndef FMUTIL_HAVE_UNISTD_H 
#define FMUTIL_HAVE_UNISTD_H  1 
#endif

/* Name of package */
#ifndef FMUTIL_PACKAGE 
#define FMUTIL_PACKAGE  "libfmutil" 
#endif

/* Define to the address where bug reports for this package should be sent. */
#ifndef FMUTIL_PACKAGE_BUGREPORT 
#define FMUTIL_PACKAGE_BUGREPORT  "o.godoy@met.no" 
#endif

/* Define to the full name of this package. */
#ifndef FMUTIL_PACKAGE_NAME 
#define FMUTIL_PACKAGE_NAME  "METNO Remote Sensing Applications Utility Library" 
#endif

/* Define to the full name and version of this package. */
#ifndef FMUTIL_PACKAGE_STRING 
#define FMUTIL_PACKAGE_STRING  "METNO Remote Sensing Applications Utility Library 1.0" 
#endif

/* Define to the one symbol short name of this package. */
#ifndef FMUTIL_PACKAGE_TARNAME 
#define FMUTIL_PACKAGE_TARNAME  "libfmutil" 
#endif

/* Define to the home page for this package. */
#ifndef FMUTIL_PACKAGE_URL 
#define FMUTIL_PACKAGE_URL  "" 
#endif

/* Define to the version of this package. */
#ifndef FMUTIL_PACKAGE_VERSION 
#define FMUTIL_PACKAGE_VERSION  "1.0" 
#endif

/* Define to 1 if you have the ANSI C header files. */
#ifndef FMUTIL_STDC_HEADERS 
#define FMUTIL_STDC_HEADERS  1 
#endif

/* Version number of package */
#ifndef FMUTIL_VERSION 
#define FMUTIL_VERSION  "1.0" 
#endif

/* Trying to add stuff from fmutil.h */
#ifndef DEG_TO_RAD
#define DEG_TO_RAD 0.0174532925199432958
#endif
#ifndef RAD_TO_DEG
#define RAD_TO_DEG 57.29577951308232
#endif

#ifndef FMRETURNCODES
typedef enum {
    FM_OK,
    FM_SYNTAX_ERR,
    FM_IO_ERR,
    FM_MEMALL_ERR,
    FM_MEMFREE_ERR,
    FM_VAROUTOFSCOPE_ERR,
    FM_OTHER_ERR
} fmreturncode;
#endif

/* once: _SRC_FMUTIL_CONFIG_H */
#endif
