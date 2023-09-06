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
 *   icedrift_uncertianties.c 
 *
 * PURPOSE:
 *   Store the uncertainties in each single-sensor product for use in the multi-sensor analysis
 * 
 * HISTORY:
 *   Thomas Lavergne, met.no/FoU, 03.11.2009    :    Upated with values from the 3 winters validation (JGR paper)
 */
#include <stdio.h>
#include "icedrift_instruments.h"

int get_Cobs_instrument (int instrumentType, float *sigX, float *sigY, float *corrXY) {
   int ret = 1;
   switch(instrumentType) {
      case INSTRUMENT_AMSR:
         *sigX = 2.50; *sigY = 2.50; *corrXY = 0.0;
         break;
      case INSTRUMENT_AMSR2:
         *sigX = 2.50; *sigY = 2.50; *corrXY = 0.0;
         break;
      case INSTRUMENT_SSMI:
         *sigX = 3.50; *sigY = 3.50; *corrXY = 0.0;
         break;
      case INSTRUMENT_SSMIS:
         *sigX = 3.50; *sigY = 3.50; *corrXY = 0.0;
         break;
      case INSTRUMENT_ASCAT:
         *sigX = 4.50; *sigY = 4.50; *corrXY = 0.00;
         break;
      default:
         fprintf(stderr,"ERROR (%s) Do not know this instrument (%d)\n",__func__,instrumentType);
         ret = 0;
   }
   return ret;

}
