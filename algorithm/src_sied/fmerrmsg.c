/*
 * NAME:
 * fmerrmsg
 * fmlogmsg
 *
 * PURPOSE:
 * Write error and log messages to screen or wherever needed...
 *
 * REQUIREMENTS:
 * NA
 *
 * INPUT:
 * NA
 *
 * OUTPUT:
 * NA
 *
 * NOTES:
 * Time of logging event is added to output. This time is taken from the
 * time function to retain compability backwards.
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 * Øystein Godøy, DNMI/FOU, 21/12/2001
 *
 * MODIFIED:
 * Øystein Godøy, met.no/FOU, 23.02.2004
 * Added LOGMSG and changed interface...
 * Thomas Lavergne, met.no, 03.12.2007: add the optional use of
 * printf-like format and parameters to fmlogmsg() and fmerrmsg()
 *
 * ID:
 * $Id: fmerrmsg.c 391 2011-11-08 10:07:39Z thomasl $
 */

#include "fmutil_config.h"
#include <stdio.h>
#ifdef FMUTIL_HAVE_STDARG_H
#include <stdarg.h>
#endif
#include <time.h>

#ifdef FMUTIL_HAVE_STDARG_H 
void fmerrmsg(char *where, char *what,...) 
#else
void fmerrmsg(char *where, char *what) 
#endif
{
#ifdef FMUTIL_HAVE_STDARG_H 
    va_list ap;
#endif
    time_t wtime;
    struct tm *mytime;
    char ftime[20];

    wtime = time(NULL);
    mytime = localtime(&wtime);
    strftime(ftime,19,"%F %R",mytime);

    fprintf(stderr,"\n");
    fprintf(stderr," ERRMSG [%s %s]:\n", where, ftime);
    fprintf(stderr,"\t");
#ifdef FMUTIL_HAVE_STDARG_H 
    va_start(ap,what);
    vfprintf(stderr, what, ap);
    va_end (ap);
#else
    fprintf(stderr,"%s", what);
#endif
    fprintf(stderr,"\n");
    fprintf(stderr," \n");
}

#ifdef FMUTIL_HAVE_STDARG_H 
void fmlogmsg(char *where, char *what,...) 
#else
void fmlogmsg(char *where, char *what) 
#endif
{

#ifdef FMUTIL_HAVE_STDARG_H 
    va_list ap;
#endif
    time_t wtime;
    struct tm *mytime;
    char ftime[20];

    wtime = time(NULL);
    mytime = localtime(&wtime);
    strftime(ftime,19,"%F %R",mytime);

    fprintf(stdout," \n");
    fprintf(stdout," LOGMSG [%s %s]:\n", where, ftime);
    fprintf(stdout,"\t");
#ifdef FMUTIL_HAVE_STDARG_H 
    va_start(ap,what);
    vfprintf(stdout, what, ap);
    va_end (ap);
#else
    fprintf(stdout,"%s", what);
#endif
    fprintf(stdout,"\n");
    fprintf(stdout," \n");
}


