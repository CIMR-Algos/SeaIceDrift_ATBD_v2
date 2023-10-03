/*
 * NAME:
 * NA
 *
 * PURPOSE:
 * NA
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
 * NA
 *
 * BUGS:
 * NA
 *
 * AUTHOR:
 * NA
 *
 * MODIFIED:
 * Thomas Lavergne, met.no, 03.12.2007: add the optional use of
 * printf-like format and parameters to fmlogmsg() and fmerrmsg()
 *
 * ID:
 * $Id: fmerrmsg.h 391 2011-11-08 10:07:39Z thomasl $
 */ 

#ifndef FMERRMSG_H
#define FMERRMSG_H

#ifdef FMUTIL_HAVE_STDARG_H
extern void fmlogmsg(char *where, char *what,...); 
extern void fmerrmsg(char *where, char *what,...); 
#else
extern void fmlogmsg(char *where, char *what); 
extern void fmerrmsg(char *where, char *what); 
#endif


#endif /* FMERRMSG */
