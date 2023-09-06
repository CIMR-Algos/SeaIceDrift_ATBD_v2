

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <netcdf.h>
#define USE_CF_ROUTINES
#include <osisaf_ice_netcdf.h>
#include <projects.h>
#include <fmutil.h>
#include <rsprod.h>
#include "icedrift_common.h"
#include "icedrift_filenames.h"
#include "icedrift_model.h"
#include "icedrift_io.h"

int load_nwp_field(char *filename,size_t *ne, size_t *nx, size_t *ny, float **DegHours) {

   int ok;
   double deghours_T;
   ok = nc_file_read(filename,"DegHours::f DegHours::# DegHours:threshold:*d",
         DegHours,ne,&deghours_T);
   if (!ok) {
      fprintf(stderr,"ERROR (%s) Did not find required NWP fields in file <%s>.\n",__func__,filename);
      return 0;
   }
   if (roundf(deghours_T*100) != 27315) {
      fprintf(stderr,"ERROR (%s) Threshold for DegHours is not 273.15! (is %f).\n",__func__,deghours_T);
      return 0;
   }

   return 1;
}

int load_image_from_netcdf_new(char *filename, short nbWaveBands, char **WaveBands,
      size_t *ne, size_t *nx, size_t *ny, 
      float **lat, float **lon, short loadLatLon,
      float ***avgT,float ***obs, short ***flag,
      short **mask, float *fillvalue,
      float *map_Ax, float *map_Ay, float *map_Bx, float *map_By, char *map_projstr, fmsec1970 *time) {

   int ret;

   /* build the request strings for reading with rsprod_file_loadFromNetCDF */
   /* (xc,yc) + (lat+lon) + (Img + Timg + Fimg) x nbWB + ice_edge */
   size_t nb_needed_datasets = 2 + (loadLatLon?2:0) + nbWaveBands*3 + 1; 
   char **needed_datasets = fmMalloc(nb_needed_datasets*sizeof(*needed_datasets));
   size_t d = 0;
   needed_datasets[d] = fmMalloc(strlen("xc")+1); strcpy(needed_datasets[d],"xc"); d++;
   needed_datasets[d] = fmMalloc(strlen("yc")+1); strcpy(needed_datasets[d],"yc"); d++;
   if (loadLatLon) {
   needed_datasets[d] = fmMalloc(strlen("lat")+1); strcpy(needed_datasets[d],"lat"); d++;
   needed_datasets[d] = fmMalloc(strlen("lon")+1); strcpy(needed_datasets[d],"lon"); d++;
   }
   for (short b = 0 ; b < nbWaveBands ; b++) {

      /* imageName */
      char *ImgName = WaveBands[b];

      /* channel name */
      unsigned int nbtokens;
      char **toks;
      fmstrsplit(ImgName,"_",&nbtokens,&toks);
      char *ChnlName = toks[0];

      /* time name */
      char TimeName[32];
      sprintf(TimeName,"%s_avgTime",ChnlName); /* for bt37v_Lap we use bt37v_avgTime */

      /* flag name */
      char FlagName[32];
      sprintf(FlagName,"%s_flag",ImgName);    /* for bt37v_Lap we use bt37v_Lap_flag */

      /* add those names to the list of datasets to be read */
      needed_datasets[d] = fmMalloc(strlen(ImgName)+1); strcpy(needed_datasets[d],ImgName); d++;
      needed_datasets[d] = fmMalloc(strlen(TimeName)+1); strcpy(needed_datasets[d],TimeName); d++;
      needed_datasets[d] = fmMalloc(strlen(FlagName)+1); strcpy(needed_datasets[d],FlagName); d++;

   }
   /* load the ice_edge */
   needed_datasets[d] = fmMalloc(strlen("ice_edge")+1); strcpy(needed_datasets[d],"ice_edge"); d++;

   /* open the netCDF file and read from it */
   int ncid;
   ret = nc_open(filename,NC_NOCLOBBER,&ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Could not open netCDF file <%s> for read access.\n",__func__,filename);
      return 1;
   }
   
   rsprod_file *img_file = NULL;
   ret = rsprod_file_loadFromNetCDF(&img_file,nb_needed_datasets,needed_datasets,ncid);
   if (ret != nb_needed_datasets) {
      fprintf(stderr,"ERROR (%s) Cannot read all datasets from the image file using rsprod library (ret = %d).\n",__func__,ret);
      for (size_t n = 0 ; n < nb_needed_datasets ; n++) { printf("%2d:<%s> ",n+1,needed_datasets[n]); }printf("\n");
      return 1;
   }
   printf("Loaded the following datasets:\n");
   for (size_t n = 0 ; n < nb_needed_datasets ; n++) { printf("%2d:<%s> ",n+1,needed_datasets[n]); }printf("\n");

   /* close the netCDF file */
   ret = nc_close(ncid);
   if (ret != NC_NOERR) {
      fprintf(stderr,"ERROR (%s) Cannot close <%s>\n",__func__,filename);
      return 1;
   } 

   return 0;
}

int load_image_from_netcdf(char *filename, short nbWaveBands, char **WaveBands,
      size_t *ne, size_t *nx, size_t *ny, 
      float **lat, float **lon, short loadLatLon,
      float ***avgT,float ***obs, short ***flag,
      short **mask, float *fillvalue,
      float *map_Ax, float *map_Ay, float *map_Bx, float *map_By, char *map_projstr, fmsec1970 *time) {

   int ret;
   int ncid, nxid, nyid;
   int *obsid, latid, lonid, iceid, *avgTid,*flagid;

   char gridfilename[1028];
   char map_area[32], areanc[32], infilenc[1028];
   
   unsigned int nbToks;
   char **tokens, **newtokens;

   char errmsg[] = "load_image_from_netcdf (error) ";
   char ncerrmsg[1024];
 
   /* prepare local variables */
   memset(map_projstr,0,1028); /* netcdf char strings are not 0-terminated */

   if (!is_tcimage_filename(filename)) {
      fprintf(stderr,"ERROR (%s) Not a valid tcimage file:\n\t%s\n",__func__,filename);
      return 1;
   }

   /* Open the netCDF formatted daily map */
   ret = nc_open(filename, NC_NOWRITE, &ncid); 
   check_ncstatus(ncid, ret, errmsg);

   /* Image mapping and dimension info */
   size_t imgDims[2];
   float  A[2],B[2];
   short  wTime = 1;
   fmsec1970 timeVal[1];
   char   timeUnit[64];
   ret = decode_CF1_struct(ncid,imgDims,A,B,map_projstr,wTime,timeVal,timeUnit,"crs");
   if (ret) {
      fprintf(stderr,"ERROR (%s) cannot decode the projection information.\n",__func__);
      return 1;
   }
   *map_Ax = A[0]; *map_Bx = B[0];
   *map_Ay = A[1]; *map_By = B[1];
   *nx = imgDims[0];
   *ny = imgDims[1];
   *ne = *nx * *ny;
   double time_value = timeVal[0];
   ret = correct_time_from_unit(time_value,timeUnit,time);
   if (ret) {
      fprintf(stderr,"ERROR (%s) cannot use the timeUnit information <%s>.\n",__func__,timeUnit);
      return 1;
   }

   /* get the varid for the variables we are interested in */
   obsid  = fmMalloc(nbWaveBands*sizeof(int));
   avgTid = fmMalloc(nbWaveBands*sizeof(int));
   flagid = fmMalloc(nbWaveBands*sizeof(int));
   for (short b = 0 ; b < nbWaveBands ; b++) {
      /* observation (image) */
      sprintf(ncerrmsg,"%s: search for %s",__func__,WaveBands[b]);
      //printf("%s\n",ncerrmsg);
      ret = nc_inq_varid(ncid,WaveBands[b],&(obsid[b])); check_ncstatus(ncid, ret, ncerrmsg);

      /* average sensing time */
      unsigned int nbtokens;
      char **toks;
      fmstrsplit(WaveBands[b],"_",&nbtokens,&toks);
      char timeName[] = "dtime";
      sprintf(ncerrmsg,"%s: search for %s",__func__,timeName);
      ret = nc_inq_varid(ncid,timeName,&(avgTid[b])); check_ncstatus(ncid, ret, ncerrmsg);
      
      /* observation flag */
      sprintf(timeName,"%s_flag",WaveBands[b]); /* for bt37v_Lap we use bt37v_Lap_flag */
      sprintf(ncerrmsg,"%s: search for %s",__func__,timeName);
      ret = nc_inq_varid(ncid,timeName,&(flagid[b]));
      if (ret != NC_NOERR) {
         /* missing flag : we will make it */
         flagid[b] = -1;
      }
   }
   ret = nc_inq_varid(ncid,"ice_edge",&iceid); check_ncstatus(ncid, ret, errmsg);
   if (loadLatLon) {
      ret = nc_inq_varid(ncid,"lon",&lonid); check_ncstatus(ncid, ret, errmsg);
      ret = nc_inq_varid(ncid,"lat",&latid); check_ncstatus(ncid, ret, errmsg);
   }

   /* Now it should be ok to read the datasets */
   if (loadLatLon) {
      *lon = fmMalloc (*ne * sizeof(float));
      ret = nc_get_var_float(ncid,lonid,*lon); check_ncstatus(ncid, ret, errmsg);
      *lat = fmMalloc (*ne * sizeof(float));
      ret = nc_get_var_float(ncid,latid,*lat); check_ncstatus(ncid, ret, errmsg); 
   }
	 
   /* fill value for image */
   ret=nc_get_att_float(ncid,obsid[0],"_FillValue",fillvalue); check_ncstatus(ncid, ret, errmsg);
     
   /* image information (obs + icemask + time) */
   *obs  = fmMalloc (nbWaveBands * sizeof(float *));
   *avgT = fmMalloc (nbWaveBands * sizeof(float *));
   *flag = fmMalloc (nbWaveBands * sizeof(short *));
   for (short b = 0 ; b < nbWaveBands ; b++) {
      (*obs)[b]  = fmMalloc (*ne * sizeof(float));
      ret   = nc_get_var_float(ncid,obsid[b],(*obs)[b]); check_ncstatus(ncid, ret, errmsg);
      (*avgT)[b] = fmMalloc (*ne * sizeof(float));
      ret   = nc_get_var_float(ncid,avgTid[b],(*avgT)[b]); check_ncstatus(ncid, ret, errmsg);
      (*flag)[b] = fmMalloc (*ne * sizeof(short));
      if (flagid[b] != -1) {
         ret   = nc_get_var_short(ncid,flagid[b],(*flag)[b]); check_ncstatus(ncid, ret, errmsg);
      } else {
         for (size_t te = 0 ; te < *ne ; te++) {
            if ( (*obs)[b][te] == *fillvalue ) {
               (*flag)[b][te] = -1.;
            } else {
               (*flag)[b][te] = 0;
            }
         }
      }
   }
   *mask = fmMalloc (*ne * sizeof(short));
   ret   = nc_get_var_short(ncid,iceid,*mask); check_ncstatus(ncid, ret, errmsg); 

   /* Close the file */
   ret = nc_close(ncid); check_ncstatus(ncid, ret, errmsg);

   return 0;
}

int loadOutPositions_ascii(int npoints, char posFile[],double *olat,double *olon) {

   int ret;
   FILE *fp;
   double llat,llon;
   unsigned int cptl = 1;
   
   if (!(fp = fopen(posFile,"r"))) {
      fprintf(stderr,"ERROR (loadOutPositions_ascii) cannot read from file <%s>\n",posFile);
      return 0;
   }

   /* first line: number of drift positions and position file containind lat, lon for those */
   while ( ((ret = fscanf(fp,"%lf%lf\n",&llat,&llon)) == 2) && (cptl <= npoints) ) {
      olat[cptl-1] = llat;
      olon[cptl-1] = llon;
      printf("Read number %u in <%s> is {%f,%f}\n",cptl,posFile,llat,llon);
      cptl++;
   }

   printf("Done reading in <%s>\n",posFile);

   fclose(fp);

   return cptl-1;

}


