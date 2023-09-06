

/*
 * PURPOSE:
 * To write ice-drift satellite products as processed by OSISAF chain.
 * 
 * BUGS:
 * None known.
 * 
 * AUTHOR: 
 * Thomas Lavergne, 10.10.2008 (extracted from icedrift_io.c)
 *
 * MODIFIED:
 * Thomas Lavergne, 10.10.2008 : remove some un-needed fields and have most of 
 *                               the remaining 'optionals'.
 * Thomas Lavergne, 11.11.2008 : introduce time_bounds in the netcdf file.
 * Thomas Lavergne, 09.03.2009 : add the uncertainty estimates
 * CVS:
 * $Id: write_icedrift_product_file.c 12095 2017-07-01 18:56:35Z thomasl $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <netcdf.h>
#include <projects.h>
#include <fmutil.h>
#include <rsprod.h>
#define USE_CF_ROUTINES
#include <osisaf_ice_netcdf.h>
#include "icedrift_common.h"
#include "icedrift_filenames.h"

int write_icedriftproduct_tonetcdf(char outdir[],char outfname[],size_t out_dims[],
      float *driftX, float *driftY, fmsec1970 *t0, fmsec1970 *t1, short *flag,
      float *sigdX,  float *sigdY,  float *corrdXdY, short *uflag,
      float *length, float *dir, float *outlatB, float *outlonB, float *outlatE, float *outlonE,
      float *corr,short *navg,float *avgX, float *avgY, float *length_avg,float *length_diff, float *stdX, float *stdY, float fillvalf,
      float **driftX_exp, float **driftY_exp, short *flag_exp, short nbflag_exp,
      short *pattern, short fillvals, int fillvali,
      char Instrument[], char Platform[], char isostartdate[], char isoenddate[], 
      char out_gridf[], char out_area[], char out_projstr[], float out_A[], float out_B[], PJ *out_pj) {

   int ret;
   int ncid;


   rsprod_file *file;
   if (rsprod_file_create(&file,NULL,NULL)) {
      fprintf(stderr,"ERROR (%s) Did not manage to create the file object.\n",__func__);
      return 1;
   }

   /* create the CF1.0 remapping information */
   rsprod_dimensions *outdims;
   char *mapInfo,*coords;
   size_t lengths[2]; lengths[0] = out_dims[0]; lengths[1] = out_dims[1];
   fmsec1970 timeVals[3];
   timeVals[1] = isodatetime2fmsec1970_alt(isostartdate);
   timeVals[2] = isodatetime2fmsec1970_alt(isoenddate);
   timeVals[0] = timeVals[1];
   float ll_bbox[4];
   ret = prepare_CF1_struct(file,lengths,out_A,out_B,out_projstr,&outdims,&mapInfo,&coords,ll_bbox,2,
         timeVals,"seconds since 1970-01-01 00:00:00");
   if (ret) {
      fprintf(stderr,"ERROR (%s) Impossible to create the CF1 remapping/time information.\n",__func__);
      return 1;
   }
   ret = rsprod_dims_addOne(outdims,0,"time",1,1);
   if (ret) {
      fprintf(stderr,"ERROR (%s) could not add the 'time' dimension.\n",progname);
      return 1;
   }

   float fmin,fmax;
   short smin,smax;
   int   imin,imax;


   /* transform the time we got (secs since 1970-01-01) into dtime to time_bounds values */
   int *dt0 = fmMalloc(out_dims[TDIM]*sizeof(*dt0));
   int *dt1 = fmMalloc(out_dims[TDIM]*sizeof(*dt1));
   for (size_t e = 0 ; e < out_dims[TDIM] ; e++) {
      if (driftX[e] != fillvalf) {
         //printf("E: %u, t0: %ld, timeVals[1]: %ld\n",e,t0[e],timeVals[1]);
         dt0[e] = (t0[e]   - timeVals[1]);
         dt1[e] = (t1[e]   - timeVals[2]);
      } else {
         dt0[e] = dt1[e] = fillvali;
      }
   }


   /* ++++++++++++++++++++++++++++++++++++++++++++ */
   /*    MANDATORY FIELD                           */
   /* ++++++++++++++++++++++++++++++++++++++++++++ */
   rsprod_field *driftX_out;
   if (rsprod_field_createStandard(
	    &driftX_out,"driftX",RSPROD_FLOAT,out_dims[TDIM],outdims,"drift along X axis of grid","drift_x","km",
	    &fillvalf,NULL,NULL,coords,mapInfo,driftX)) {
      fprintf(stderr,"ERROR (%s) could not create Field object <driftX>.\n",__func__);
      return 1;
   }
   rsprod_field *driftY_out;
   if (rsprod_field_createStandard(
	    &driftY_out,"driftY",RSPROD_FLOAT,out_dims[TDIM],outdims,"drift along Y axis of grid","drift_y","km",
	    &fillvalf,NULL,NULL,coords,mapInfo,driftY)) {
      fprintf(stderr,"ERROR (%s) could not create Field object <driftY>.\n",__func__);
      return 1;
   }
   char tunit[1024];
   sprintf(tunit,"seconds since %s",isostartdate);
   rsprod_field *t0_out;
   if (rsprod_field_createStandard(
	    &t0_out,"t0",RSPROD_INT,out_dims[TDIM],outdims,"time of start pixel","time_start",tunit,
	    &fillvali,NULL,NULL,coords,mapInfo,dt0)) {
      fprintf(stderr,"ERROR (%s) could not create Field object <t0>.\n",__func__);
      return 1;
   }
   rsprod_field *t1_out;
   sprintf(tunit,"seconds since %s",isoenddate);
   if (rsprod_field_createStandard(
	    &t1_out,"t1",RSPROD_INT,out_dims[TDIM],outdims,"time of end pixel","time_end",tunit,
	    &fillvali,NULL,NULL,coords,mapInfo,dt1)) {
      fprintf(stderr,"ERROR (%s) could not create Field object <t1>.\n",__func__);
      return 1;
   }
   rsprod_field *flag_out;
   if (rsprod_field_createStandard(
	    &flag_out,"flag",RSPROD_SHORT,out_dims[TDIM],outdims,"Quality/Processing Flag","qual_flag","1",&fillvals,NULL,NULL,
	    coords,mapInfo,flag)) {
      fprintf(stderr,"ERROR (%s) could not create Field object <flag>.\n",__func__);
      return 1;
   }

   /* ++++++++++++++++++++++++++++++++++++++++++++ */
   /*     OPTIONAL FIELD                           */
   /* ++++++++++++++++++++++++++++++++++++++++++++ */
   rsprod_field **driftXexp_out,**driftYexp_out;
   if (driftX_exp && driftY_exp) {
      if ( (nbflag_exp <= 0) || (flag_exp == NULL) ) {
         fprintf(stderr,"ERROR (%s) Unable to write the driftX_exp variable since flag values are wrong.\n",progname);
         return 1;
      }
      driftXexp_out = fmMalloc(nbflag_exp*sizeof(*driftXexp_out));
      driftYexp_out = fmMalloc(nbflag_exp*sizeof(*driftYexp_out));
      char fieldName[32];
      char stdName[32];
      for (int ff = 0 ; ff < nbflag_exp ; ff++) {
         sprintf(fieldName,"driftX_%02d",flag_exp[ff]);
         sprintf(stdName,"drift along X (flag %02d)",flag_exp[ff]);
         if (rsprod_field_createStandard(
               &(driftXexp_out[ff]),fieldName,RSPROD_FLOAT,out_dims[TDIM],outdims,stdName,"drift_x","km",
	           &fillvalf,NULL,NULL,coords,mapInfo,driftX_exp[ff])) {
            fprintf(stderr,"ERROR (%s) could not create Field object <%s>.\n",__func__,fieldName);
            return 1;
         }
         sprintf(fieldName,"driftY_%02d",flag_exp[ff]);
         sprintf(stdName,"drift along Y (flag %02d)",flag_exp[ff]);
         if (rsprod_field_createStandard(
               &(driftYexp_out[ff]),fieldName,RSPROD_FLOAT,out_dims[TDIM],outdims,stdName,"drift_y","km",
	           &fillvalf,NULL,NULL,coords,mapInfo,driftY_exp[ff])) {
            fprintf(stderr,"ERROR (%s) could not create Field object <%s>.\n",__func__,fieldName);
            return 1;
         }
      }
   }

   rsprod_field *sigdX_out;
   if (sigdX) {
      if (rsprod_field_createStandard(
	       &sigdX_out,"sigdX",RSPROD_FLOAT,out_dims[TDIM],outdims,"standard deviation on driftX","stddevX","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,sigdX)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <sigdX>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *sigdY_out;
   if (sigdY) {
      if (rsprod_field_createStandard(
	       &sigdY_out,"sigdY",RSPROD_FLOAT,out_dims[TDIM],outdims,"standard deviation on driftY","stddevY","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,sigdY)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <sigdY>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *corrdXdY_out;
   if (corrdXdY) {
      if (rsprod_field_createStandard(
	       &corrdXdY_out,"corrdXdY",RSPROD_FLOAT,out_dims[TDIM],outdims,
               "correlation of errors in driftX and driftY","corrdXdY","1",
	       &fillvalf,NULL,NULL,coords,mapInfo,corrdXdY)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <corrdXdY>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *uflag_out;
   if (uflag) {
      if (rsprod_field_createStandard(
	    &uflag_out,"uflag",RSPROD_SHORT,out_dims[TDIM],outdims,"Uncertainty Flag","uncert_flag","1",
            &fillvals,NULL,NULL,coords,mapInfo,uflag)) {
         fprintf(stderr,"ERROR (%s) could not create Field object <uflag>.\n",__func__);
         return 1;
      }
   }
   rsprod_field *length_out;
   if (length) {
      if (rsprod_field_createStandard(
	       &length_out,"length",RSPROD_FLOAT,out_dims[TDIM],outdims,"length of the drift","length","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,length)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <length>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *dir_out;
   if (dir) {
      if (rsprod_field_createStandard(
	       &dir_out,"dir",RSPROD_FLOAT,out_dims[TDIM],outdims,"direction of the drift","direction","degrees from North",
	       &fillvalf,NULL,NULL,coords,mapInfo,dir)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <dir>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *latB_out;
   if (outlatB) {
      if (rsprod_field_createStandard(
	       &latB_out,"latStart",RSPROD_FLOAT,out_dims[TDIM],outdims,"latitude of start of drift point","lat_start","degrees",
	       &fillvalf,NULL,NULL,coords,mapInfo,outlatB)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <latStart>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *lonB_out;
   if (outlonB) {
      if (rsprod_field_createStandard(
	       &lonB_out,"lonStart",RSPROD_FLOAT,out_dims[TDIM],outdims,"longitude of start of drift point","lon_start","degrees",
	       &fillvalf,NULL,NULL,coords,mapInfo,outlonB)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <lonStart>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *latE_out;
   if (outlatE) {
      if (rsprod_field_createStandard(
	       &latE_out,"latEnd",RSPROD_FLOAT,out_dims[TDIM],outdims,"latitude of end of drift point","lat_end","degrees",
	       &fillvalf,NULL,NULL,coords,mapInfo,outlatE)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <latEnd>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *lonE_out;
   if (outlonE) {
      if (rsprod_field_createStandard(
	       &lonE_out,"lonEnd",RSPROD_FLOAT,out_dims[TDIM],outdims,"longitude of end of drift point","lon_end","degrees",
	       &fillvalf,NULL,NULL,coords,mapInfo,outlonE)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <lonEnd>.\n",__func__);
	 return 1;
      }
   }

   /* ++++++++++++++++++++++++++++++++++++++++++++ */
   /*     DEBUGGING FIELD                          */
   /* ++++++++++++++++++++++++++++++++++++++++++++ */
   rsprod_field *corr_out;
   if (corr) {
      if (rsprod_field_createStandard(
	       &corr_out,"corr",RSPROD_FLOAT,out_dims[TDIM],outdims,"cross-correlation at best estimate","cross-correlation","1",
	       &fillvalf,NULL,NULL,coords,mapInfo,corr)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <corr>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *avgN_out;
   if (navg) {
      if (rsprod_field_createStandard(
	       &avgN_out,"navg",RSPROD_SHORT,out_dims[TDIM],outdims,"number of neighbour vectors to compute average","n_avg_drift","1",
	       &fillvals,NULL,NULL,coords,mapInfo,navg)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <navg>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *avgX_out;
   if (avgX) {
      if (rsprod_field_createStandard(
	       &avgX_out,"avgX",RSPROD_FLOAT,out_dims[TDIM],outdims,"local average drift along X axis of grid","avg_drift_x","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,avgX)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <avgX>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *avgY_out;
   if (avgY) {
      if (rsprod_field_createStandard(
	       &avgY_out,"avgY",RSPROD_FLOAT,out_dims[TDIM],outdims,"local average drift along Y axis of grid","avg_drift_y","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,avgY)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <avgY>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *avglen_out;
   if (length_avg) {
      if (rsprod_field_createStandard(
	       &avglen_out,"avglen",RSPROD_FLOAT,out_dims[TDIM],outdims,"length of local average drift","avg_len","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,length_avg)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <avglen>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *difflen_out;
   if (length_diff) {
      if (rsprod_field_createStandard(
	       &difflen_out,"difflen",RSPROD_FLOAT,out_dims[TDIM],outdims,"length of difference vector AVG-DRIFT","diff_len","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,length_diff)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <difflen>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *stdX_out;
   if (stdX) {
      if (rsprod_field_createStandard(
	       &stdX_out,"stdX",RSPROD_FLOAT,out_dims[TDIM],outdims,"local variability of drift along X axis of grid","std_drift_x","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,stdX)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <stdX>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *stdY_out;
   if (stdY) {
      if (rsprod_field_createStandard(
	       &stdY_out,"stdY",RSPROD_FLOAT,out_dims[TDIM],outdims,"local variability drift along Y axis of grid","std_drift_y","km",
	       &fillvalf,NULL,NULL,coords,mapInfo,stdY)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <stdY>.\n",__func__);
	 return 1;
      }
   }
   rsprod_field *pattern_out;
   if (pattern) {
      if (rsprod_field_createStandard(
	       &pattern_out,"pattern",RSPROD_SHORT,out_dims[TDIM],outdims,"Size of the correlation pattern","pattern_index","1",
	       &fillvals,NULL,NULL,
	       coords,mapInfo,pattern)) {
	 fprintf(stderr,"ERROR (%s) could not create Field object <pattern>.\n",__func__);
	 return 1;
      }
   }

   /* create the list of Global attributes */
   rsprod_attributes *globalAttributes;
   rsprod_attributes_create(&globalAttributes);
   ret = 0;
   rsprod_attr *attr;
   ret += rsprod_attr_createWithCopyValues(&attr,"filename",RSPROD_CHAR,strlen(outfname),outfname) + 
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"Conventions",RSPROD_CHAR,strlen("CF-1.0"),"CF-1.0") + 
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"instrument",RSPROD_CHAR,strlen(Instrument),Instrument) + 
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"platform",RSPROD_CHAR,strlen(Platform),Platform) + 
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"processing_software",RSPROD_CHAR,strlen(__func__),__func__) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"start_date_and_time",RSPROD_CHAR,strlen(isostartdate),isostartdate) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"end_date_and_time",RSPROD_CHAR,strlen(isoenddate),isoenddate) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"remapping_gridfile",RSPROD_CHAR,strlen(out_gridf),out_gridf) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"remapping_gridname",RSPROD_CHAR,strlen(out_area),out_area) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   ret += rsprod_attr_createWithCopyValues(&attr,"remapping_projstr",RSPROD_CHAR,strlen(out_projstr),out_projstr) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   float gridparams[4]; 
   gridparams[0]=out_A[YDIM]; gridparams[1]=out_B[YDIM]; gridparams[2]=out_A[XDIM]; gridparams[3]=out_B[XDIM]; 
   ret += rsprod_attr_createWithCopyValues(&attr,"remapping_gridparams",RSPROD_FLOAT,4,gridparams) +
      rsprod_attributes_addCopyAttr(globalAttributes,attr);
   rsprod_attr_delete(attr);
   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot create list of  global attributes.\n",__func__);
      return 1;
   }


   /* add the global attributes to the rsprod_file object */
   file->glob_attr = globalAttributes;

   /* add the data fields to the rsprod_file object */
   ret = 0;
   ret += rsprod_file_addDataset(file,driftX_out);
   ret += rsprod_file_addDataset(file,driftY_out);
   ret += rsprod_file_addDataset(file,t0_out);
   ret += rsprod_file_addDataset(file,t1_out);
   ret += rsprod_file_addDataset(file,flag_out);
   if (driftX_exp && driftY_exp) {
      for (int f = 0 ; f < nbflag_exp ; f++) {
         ret += rsprod_file_addDataset(file,driftXexp_out[f]);
         ret += rsprod_file_addDataset(file,driftYexp_out[f]);
      }
   }
   if (sigdX)
      ret += rsprod_file_addDataset(file,sigdX_out);
   if (sigdY)
      ret += rsprod_file_addDataset(file,sigdY_out);
   if (corrdXdY)
      ret += rsprod_file_addDataset(file,corrdXdY_out);
   if (uflag)
      ret += rsprod_file_addDataset(file,uflag_out);
   if (length)
      ret += rsprod_file_addDataset(file,length_out);
   if (dir)
      ret += rsprod_file_addDataset(file,dir_out); 
   if (outlatB)
      ret += rsprod_file_addDataset(file,latB_out);
   if (outlonB)
      ret += rsprod_file_addDataset(file,lonB_out); 
   if (outlatE)
      ret += rsprod_file_addDataset(file,latE_out);
   if (outlonE)
      ret += rsprod_file_addDataset(file,lonE_out); 
   if (corr)
      ret += rsprod_file_addDataset(file,corr_out);
   if (navg)
      ret += rsprod_file_addDataset(file,avgN_out);
   if (avgX)
      ret += rsprod_file_addDataset(file,avgX_out);
   if (avgY)
      ret += rsprod_file_addDataset(file,avgY_out); 
   if (length_avg)
      ret += rsprod_file_addDataset(file,avglen_out); 
   if (length_diff)
      ret += rsprod_file_addDataset(file,difflen_out); 
   if (stdX) 
      ret += rsprod_file_addDataset(file,stdX_out); 
   if (stdY) 
      ret += rsprod_file_addDataset(file,stdY_out); 
   if (pattern)
      ret += rsprod_file_addDataset(file,pattern_out);
   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot add the fields to rsprod file.\n",__func__);
      return 1;
   }

   /* output filename */
   char *fulloutfile = fmMalloc((strlen(outfname)+strlen(outdir)+1+1)*sizeof(char));
   sprintf(fulloutfile,"%s/%s",outdir,outfname);
   printf("Full Name to ouput file is <%s>\n",fulloutfile);
   
   /* open the netcdf file */
   ret = nc_create(fulloutfile,NC_WRITE,&ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Cannot nc_open(NC_WRITE) netcdf file:\n\t%s\n",__func__,fulloutfile);
      return 1;
   }

   /* write the rsprod_file object to netCDF */
   ret = rsprod_file_writeToNetCDF(file,ncid);
   if (ret) {
      fprintf(stderr,"ERROR (%s) Cannot write to netCDF file <%s>\n",__func__,fulloutfile);
      return 1;
   }
   
   /* Close NetCDF file */
   ret = nc_close(ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Cannot nc_close() on file.\n",__func__);
      return 1;
   }

   return 0;
}

