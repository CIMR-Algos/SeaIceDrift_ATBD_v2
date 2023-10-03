/* *****************************************************************
 * COPYRIGHT:
 * EUMETSAT
 *
 * PRODUCED BY:
 * Norwegian Meteorological Institute (met.no)
 * Research and Development Department
 * P.O.BOX 43 - Blindern, N-0313 OSLO, NORWAY
 *
 * This SW was developed by met.no and the Danish Meteorological
 * Institute (DMI) within the context of the Co-operation Agreement
 * for the development of a pilot SAF on Ocean and Sea Ice.
 * *****************************************************************/ 

/*
 * NAME:
 *    errorcodes.h
 *
 * PURPOSE:
 *    To implement the standard exit/error codes
 *    for OSI-SAF applications.
 *
 * REQUIREMENTS:
 *    C preprocessor (might not be straightforward to use within Fortran programs)
 *
 * INPUT:
 *    NA
 *
 * OUTPUT:
 *    NA
 *
 * NOTES:
 *    Developers are welcome to add some new error codes to implement
 *    their specific needs. They should take care of using ranges of 
 *    values that are not reserved for other needs.
 *
 * BUGS:
 *    NA
 *
 * AUTHOR:
 *    Thomas Lavergne (met.no) - 20th June 2007
 *
 * MODIFIED:
 * NA
 *
 * CVS_ID:
 * $Id$
 */ 
#ifndef ERRORCODES_H
#define ERRORCODES_H

#define OSISAF_EXIT_CORRECT      0
#define OSISAF_ERROR_CMDSYNTAX   1
#define OSISAF_ERROR_SYSIO       2
#define OSISAF_ERROR_SYSMEM      3

#define OSISAF_ERROR_CMDVALUE   10
#define OSISAF_ERROR_OTHER      10

#endif 
