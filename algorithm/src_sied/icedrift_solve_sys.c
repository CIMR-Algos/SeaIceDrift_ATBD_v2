

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <projects.h>
#include <fmutil.h>
#include <errorcodes.h>
#include "icedrift_flags.h"
#include "icedrift_common.h"
#include "icedrift_prepost.h"
#include "icedrift_model.h"
#include "icedrift_solve_filter.h"
#include "icedrift_solve_common.h"

/* global variables for communications with the I/O in icedrift_prepost.c */
char OptimMetric[5];
extern char *IceDriftProductDescription;

/* global variables for communications with the 'wrappers' in icedrift_solve_common.c */
int     xdim;
double *x;
size_t  currentVector;
short   currentPattern;

void setProcessingMethodString(char *meth);
void setFilteringLevelString(char *lev);
void setProductDesctiptionString(void);

char progname[] = "icedrift_solve_sys";

int findBestX_sys(size_t p, short pat, size_t xnbbins, double xbins[], size_t ynbbins, double ybins[], 
      double x[], double xbest[], double bestval[], short processingflag[], int plotCorrMap);

/* if defined, the correlation map of the first vector
 * will be saved and plot as contour levels using gnuplot */
#undef GNUPLOT_CORRELATION_MAP

int main(int argc, char *argv[]) {

   int ret;
   int n;
   double *xbest;
   short  *processingflag;
   double *bestval;
   


   if ( (argc < 2) || (argc > 3)) {
      fmerrmsg("icedrift_solve_sys","Expecting one (or to) parameters on command line.\n\tThe 'divfactor' and optionnaly an 'outdir'.");
      exit(OSISAF_ERROR_CMDSYNTAX);
   }

   float divfactor = atof(argv[1]);
   if (argc == 3)
      strcpy(parametersFile,argv[2]);
   
   int numberUnprocReported=0;
   numbmod(&xdim);

   /* open log/report file. */
   FILE *logf = fopen(reportFile,"w");
   if (!logf) {
      fmerrmsg(progname,"Cannot open the log file for writing.");
      exit(OSISAF_ERROR_OTHER);
   }

   x     = fmMalloc(sizeof(double)*xdim);
   xbest = fmMalloc(sizeof(double)*xdim);
   processingflag = fmMalloc(sizeof(double)*xdim);

   initmod(&xdim,x);
   size_t ndim = NUNKNOWNS;

   double xlpix     = divfactor * img_Ax; 
   int   xnbbins    = rint(maxdriftdistance/xlpix);
   double xrmaxdist = xnbbins*xlpix;
   
   double ylpix     = divfactor * img_Ay; 
   int    ynbbins   = rint(maxdriftdistance/ylpix);
   double yrmaxdist = ynbbins*ylpix;

   printf("New rmaxdist: (%f,%f)\n",xrmaxdist,yrmaxdist);

   xnbbins = 2 * xnbbins + 1;
   ynbbins = 2 * ynbbins + 1;
   size_t anbbins = 1;
   size_t bnbbins = 21;

   printf("SOLVE_SYS: divfactor = %f. OptimMetric is <%s>\n",divfactor,OptimMetric);
   printf("\tnxbin = %u, nybin = %u\n",xnbbins,ynbbins);

   /* initialize the discretized LUT on which we will
    * evaluate the function */
   double step;
   double *xbins = fmMalloc(xnbbins * sizeof(double));
   double *ybins = fmMalloc(ynbbins * sizeof(double));
   double *abins = fmMalloc(anbbins * sizeof(double));
   double *bbins = fmMalloc(bnbbins * sizeof(double));
      
   step = xlpix;
   for (size_t e = 0 ; e < xnbbins ; e++) {
      double bins  = e  * step;
      xbins[e] = -xrmaxdist + bins;
   }
   step = ylpix;
   for (size_t e = 0 ; e < ynbbins ; e++) {
      double bins  = e * step;
      ybins[e] = -yrmaxdist + bins;
   }
   printf("\tstep in X: %f km\n",xbins[1]-xbins[0]);
   printf("\tstep in Y: %f km\n",ybins[1]-ybins[0]);

   char IceDriftProcessingMeth[20];
   sprintf(IceDriftProcessingMeth,"sys-%05.2fkm",xbins[1]-xbins[0]);
   setProcessingMethodString(IceDriftProcessingMeth);
   setProductDescriptionString();


#ifdef GNUPLOT_CORRELATION_MAP
   /* allocate memory space to store the correlation map */
   float **corrMap = fmMalloc(sizeof(*corrMap)*xnbbins);
   for ( size_t ex = 0 ; ex < xnbbins ; ex++ ) {
      corrMap[ex] = fmMalloc(sizeof(**corrMap)*ynbbins);
   }
#endif
   

   bestval = fmMalloc(sizeof(double) * NDRIFTPIXELS);

   /* prepare and initialize the variables for the filtering steps */
   double correlation_limit = 0.5;
   float  radiusNeighbours  = 100;
   size_t nbNeighbours,neighboursCenter;
   short *neighboursMask;
   long  *neighboursIndexes[3];
   ret = setNeighbourhoodPattern(radiusNeighbours,&nbNeighbours,&neighboursCenter,&neighboursMask,neighboursIndexes);
   if (ret) {
      fmerrmsg(progname,"Cannot set the neighborhood pattern of radius %f km. \n",radiusNeighbours);
      exit(OSISAF_ERROR_OTHER);
   }
   double *neighboursX  = fmMalloc(nbNeighbours*sizeof(double));
   double *neighboursY  = fmMalloc(nbNeighbours*sizeof(double));

   /* initialize all the drift position to 'unprocessed' (to be processed) */
   size_t *tmppointer;
   size_t *tobeprocessed = fmMalloc(NDRIFTPIXELS*sizeof(size_t));
   short  *pattern       = fmMalloc(NDRIFTPIXELS*sizeof(short));
   size_t nb_tobeprocessed = 0;
   for (size_t p = 0 ; p < NDRIFTPIXELS ; p++) {
      processingflag[p] = ICEDRIFT_UNPROCESSED;
      tobeprocessed[p]  = p;
      pattern[p]        = -1;
      nb_tobeprocessed++;
   }
   short atleast_one_unprocessed = (nb_tobeprocessed!=0);
   if (!atleast_one_unprocessed) goto finish_processing;
  
  
   size_t *new_tobeprocessed = fmMalloc(NDRIFTPIXELS*sizeof(size_t));
   size_t new_nbtobeprocessed;
   /* go through all the locations and flag those which are too close to the borders of the image grid (or outside) */
   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      short pat = 0;
      size_t p = tobeprocessed[q];
      long xi,yi;

      if (iwcs[p] == img_dims[TDIM]) {
	 /* the center of the pattern is outside the image grid. flag and to not use */
	 processingflag[p] = ICEDRIFT_OUTSIDE_IMGBORDER;
	 continue;
      }
      locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
      center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];
      pattern[p] = pat;

      load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
	 center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

      long cptOut=0;
      for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
	 if ( pattern_mask[pat][o] && (pattern_isvalid[pat][0][o] == TCIMAGE_OUTSIDE_GRID) ) {
	    /* found one point in the pattern that is outside the image grid */
	    cptOut++;
	    break;
	 }
      }
      if (cptOut != 0) {
	 processingflag[p] = ICEDRIFT_CLOSETO_IMGBORDER;
	 continue;
      }

      /* store the points who passed all the tests */
      new_tobeprocessed[new_nbtobeprocessed] = p;
      new_nbtobeprocessed++;
   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;
   

   /* go through all the locations and flag those whose central point is over land */
   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      size_t p = tobeprocessed[q];
	 
      if (pointIsOnLand(iwcs[p],icelandmask[BEG])) {
	 processingflag[p] = ICEDRIFT_CENTER_OVER_LAND;
      } else {
	 new_tobeprocessed[new_nbtobeprocessed] = p;
	 new_nbtobeprocessed++;
      }

   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;


   /* go through all the remaining locations and flag those for which
    * all pixels in the subimage are 1) nodata or 2) over the ocean or land */
   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {

      size_t p = tobeprocessed[q];
      long xi,yi;
      locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
      center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];

      /* First work with the largest pattern (index 0). */
      short pat = 0;
      pattern[p] = pat;

      load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
	 center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

      load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
	 center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);	    
      long cptUnproc=0;
      long cptMissing=0;
      long cptIce=0;
      long cptPoints=0;
      for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
	 if ( pattern_mask[pat][o] ) {
	    cptPoints++;
	    if ( (pattern_isvalid[pat][0][o] == TCIMAGE_NODATA) || (pattern_isvalid2[pat][0][o] == TCIMAGE_NODATA) ) {
	       cptMissing++;
	    }
	    if ((pattern_icemask[pat][o] == 3) && (pattern_icemask2[pat][o] == 3)) {
	       cptIce++;
	    }
	    if ( (pattern_isvalid[pat][0][o] == TCIMAGE_UNPROCESSED) || (pattern_isvalid2[pat][0][o] == TCIMAGE_UNPROCESSED) ) {
	       cptUnproc++;
	    }
	 }
      }
      if ( cptIce == 0 ) {
	 processingflag[p] = ICEDRIFT_NOICE;
	 continue;
      } 

      /* no everything is ice. */
      if ( cptIce < cptPoints ) { 

	 /* try with the second pattern */
	 pat = 1;
	 pattern[p] = pat;
	 load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
	    center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

	 load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
	    center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);	    
	 cptMissing=0;
	 cptUnproc=0;
	 cptIce=0;
	 cptPoints=0;
	 for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
	    if ( pattern_mask[pat][o] ) {
	       cptPoints++;
	       if ( (pattern_isvalid[pat][0][o] == TCIMAGE_NODATA) || (pattern_isvalid2[pat][0][o] == TCIMAGE_NODATA) ) {
		  cptMissing++;
	       }
	       if ((pattern_icemask[pat][o] == 3) && (pattern_icemask2[pat][o] == 3)) {
		  cptIce++;
	       }
	       if ( (pattern_isvalid[pat][0][o] == TCIMAGE_UNPROCESSED) || (pattern_isvalid2[pat][0][o] == TCIMAGE_UNPROCESSED) ) {
	          cptUnproc++;
	       }
	    }
	 }
      
      } 
      /* check for ice */

      /* if not enough ice is found, flag and skip this drift location */
      if (cptIce != cptPoints) {
	 processingflag[p] = ICEDRIFT_CLOSETO_COAST_OR_EDGE;
	 continue;
      } 
      /* => if we make it here, all the pixels in the pattern (index 'pat') are supposedly ice */
      /*    we can now check if there are missing data */
      if ( cptUnproc != 0 ) {
	 /* this is bad: it shows that we avoided to process (tc or laplacian) an area which would 
	  * be usefull. Use a special flag and report to log file
	  */
	 if ( numberUnprocReported < 10 ) {
	    fprintf(stderr,"WARNING (%s) Found an interesting UNPROCESSED area. Check the flags to TC daily images\n",progname);
	 } else if (numberUnprocReported == 10) {
	    fprintf(stderr,"WARNING (%s) Stop reporting 'Unprocessed area' (%d).\n",progname,numberUnprocReported);
	 } else {
	    ;
	 }
	 numberUnprocReported++;
	 processingflag[p] = ICEDRIFT_CLOSETO_UNPROCESSED;
	 continue;
      }

      if ( cptMissing != 0 ) {

	 /* there are some missing data. Can we try reducing the pattern? */
	 if (pat != (NBPATTERNS-1)) {
	    /* load the smaller pattern and check again */
	    pat = 1;
	    pattern[p] = pat;

	    load_subimage(&ret,obs[BEG],TCflag[BEG],icelandmask[BEG],img_dims,
	       center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], pattern_img[pat], pattern_isvalid[pat], pattern_icemask[pat]);

	    load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
	       center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);	    
	    cptMissing=0;
	    cptPoints=0;
	    for (size_t o = 0 ; o < pattern_size[pat] ; o++) {
	       if ( pattern_mask[pat][o] ) {
		  cptPoints++;
		  if ( (pattern_isvalid[pat][0][o] == TCIMAGE_NODATA) || (pattern_isvalid2[pat][0][o] == TCIMAGE_NODATA) ) {
		     cptMissing++;
		  }
	       }
	    }
	 }
	 if (cptMissing != 0) {
	    /* still some missing points. Cannot process. */
	    processingflag[p] = ICEDRIFT_CLOSETO_MISSING;
	    continue;
	 } 
	 /* else, we try to keep this location */
      } 


      /* keep those who passed the previous tests (and record their pattern number) */
      new_tobeprocessed[new_nbtobeprocessed] = p;
      new_nbtobeprocessed++;

   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;

   new_nbtobeprocessed = 0;
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      size_t p = tobeprocessed[q];

      if (!(p%500))
         printf("p: %d\n",p);

      short pat = pattern[p];

      if (!( (pat >= 0) && (pat <NBPATTERNS) )) {
	 printf("ERROR with pattern. (p: %u q:%u)\n",p,q);
	 processingflag[p] = ICEDRIFT_FAILS;
	 continue;
      } 
      
      currentVector = p;
   
      /* load the subimage in end_drift image */
      long xi,yi;
      locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
      center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];

      load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
         center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
	 pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);
   
      /* set the correlation function parameters (sigmoid) */
      setBestKnowledge(img_lat[iwcs[p]],img_lon[iwcs[p]]);

      /* perform the optimization */
      int oret = findBestX_sys(p,pat,xnbbins,xbins,ynbbins,ybins, 
	 x,xbest,bestval, processingflag, 0);
      if (oret) {
	 processingflag[p] = ICEDRIFT_FAILS;
	 continue;
      } else {
	 processingflag[p] = ICEDRIFT_OK;
      }
      new_tobeprocessed[new_nbtobeprocessed] = p;
      new_nbtobeprocessed++;

   }
   nb_tobeprocessed = new_nbtobeprocessed;
   /* swap the old and new memory space */
   tmppointer = tobeprocessed;
   tobeprocessed = new_tobeprocessed;
   new_tobeprocessed = tmppointer;
   atleast_one_unprocessed = (nb_tobeprocessed!=0);
   /* test if all is done */
   if (!atleast_one_unprocessed) goto finish_processing;

finish_processing:
   fmlogmsg(progname,"Finished processing, write the netcdf file.");

   /* post-process and transform the drift result */
   float *latB = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *lonB = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *latE = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *lonE = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *leng = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *dire = fmMalloc(NDRIFTPIXELS * sizeof(float));

   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      size_t p = tobeprocessed[q];
      long xcoordp,ycoordp;
      locfmijmap(p,out_dims[XDIM],&xcoordp,&ycoordp);

      double xdrift = xbest[p*NUNKNOWNS + XDRIFT_IDX];
      double ydrift = xbest[p*NUNKNOWNS + YDRIFT_IDX];	 

      /* compute the latitude and longitude of start point */
      latB[p] = img_lat[iwcs[p]]; lonB[p] = img_lon[iwcs[p]];
      /* compute the latitude and longitude of final point */ 
      double lat,lon;
      long xcoordi,ycoordi;
      locfmijmap(iwcs[p],img_dims[XDIM],&xcoordi,&ycoordi);
      remap_xy2ll(-(xdrift/img_Ax)+(xcoordi+0.5),-(ydrift/img_Ay)+(ycoordi+0.5),img_pj,img_Ax,img_Bx,img_Ay,img_By,&lat,&lon);
      latE[p] = lat; lonE[p] = lon;

      /*
      printf("ICEDRIFT is x:%f, y:%f. (%f,%f)(%f,%f) -> (%f,%f)(%f %f)\n",
	    xdrift,ydrift,
	    (float)xcoordi+0.5,(float)ycoordi+0.5,
	    latB[p],lonB[p],
	    -(xdrift/img_Ax)+(xcoordi+0.5),-(ydrift/img_Ay)+(ycoordi+0.5),
	    latE[p],lonE[p]);
	    */

      /* compute the length and direction of the arrows */ 
      double dlen,ddir;
      compute_distance(latB[p],lonB[p],latE[p],lonE[p],&dlen);
      compute_directionToNorth(latB[p],lonB[p],latE[p],lonE[p],&ddir);
      leng[p] = dlen; dire[p] = ddir;
   }

   /* compute the zonal average drift vector */ 
   float  *xDrift = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *yDrift = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *xavg = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *yavg = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *avgLat = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float  *avgLon = fmMalloc(NDRIFTPIXELS * sizeof(float));
   size_t *navg = fmMalloc(NDRIFTPIXELS * sizeof(size_t));

   float *len_avg  = fmMalloc(NDRIFTPIXELS * sizeof(float));
   float *len_diff = fmMalloc(NDRIFTPIXELS * sizeof(float));
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {

      size_t p = tobeprocessed[q];
      long xcoordp,ycoordp;
      locfmijmap(p,out_dims[XDIM],&xcoordp,&ycoordp);

      size_t nbValidNeighbours   = 0;
      size_t nbInvalidNeighbours = 0;
      for (size_t n = 0 ; n < nbNeighbours ; n++) {
	 if (neighboursMask[n]) {

	    long xneighbour = xcoordp + neighboursIndexes[XDIM][n];
	    long yneighbour = ycoordp + neighboursIndexes[YDIM][n];
	    
	    if ( !pointIsInGrid((double)(xneighbour),(double)(yneighbour),out_dims) ) continue;

	    long windexNeighbour;
	    locfmivec(&windexNeighbour,xneighbour,yneighbour,out_dims[XDIM]);

	    if ( (processingflag[windexNeighbour] == ICEDRIFT_OK) && 
		  (bestval[windexNeighbour] >= correlation_limit) ) {
	       neighboursX[nbValidNeighbours] = xbest[windexNeighbour*NUNKNOWNS + XDRIFT_IDX];
	       neighboursY[nbValidNeighbours] = xbest[windexNeighbour*NUNKNOWNS + YDRIFT_IDX];
	       nbValidNeighbours++;
	    } else if (processingflag[windexNeighbour] == ICEDRIFT_CLOSETO_MISSING) {
	       ;
	    } else {
	       nbInvalidNeighbours++;
	    }
	 }
      }

      navg[p] = nbValidNeighbours;
      xavg[p] = 1000; yavg[p] = 1000;
      if (nbValidNeighbours >= 5) { /* TODO: replace by a statistical test: does this average have a meaning */
      //if (nbInvalidNeighbours != nbNeighbours) { /* TODO: replace by a statistical test: does this average have a meaning */
	 double xavgP,yavgP,xstdP,ystdP;
	 ret = compute_mean_vector(&xavgP,&yavgP,&xstdP,&ystdP,nbValidNeighbours,neighboursX,neighboursY);
	 xavg[p] = xavgP;yavg[p] = yavgP;
	 navg[p] = nbValidNeighbours;

	 /* get the lat and lon associated with the average drift */
	 long xcoordi,ycoordi;
	 double latAvg,lonAvg;
	 locfmijmap(iwcs[p],img_dims[XDIM],&xcoordi,&ycoordi);
	 remap_xy2ll(-(xavgP/img_Ax)+(xcoordi+0.5),-(yavgP/img_Ay)+(ycoordi+0.5),img_pj,img_Ax,img_Bx,img_Ay,img_By,&latAvg,&lonAvg);
	 avgLat[p] = latAvg;
	 avgLon[p] = lonAvg;

	 /* compute the length of the average drift */
	 double lenavg;
	 compute_distance(latB[p],lonB[p],latAvg,lonAvg,&lenavg);

	 /* compute the length of the difference vector: avg - drift */
	 double lendiff;
	 compute_distance(latAvg,lonAvg,latE[p],lonE[p],&lendiff);

	 /* store those values for writing in the file */
	 len_avg[p]  = lenavg;
	 len_diff[p] = lendiff;

	 /*
	 if (p == 8072) {
	    printf("p:%u AVG:(%f,%f) VEC:(%f %f) and lendiff[p] = %f\n",
		  p,avgLat[p],avgLon[p],latE[p],lonE[p],len_diff[p]);
	 }
	 */

      } else {
	 len_diff[p] = -1.; 
      }
   }
  
   setFilteringLevelString("lev0");
   setProductDescriptionString();
   /* write the unfiltered data . */
   ret = prepareAndWrite_icedriftProduct(&xdim,xbest,bestval,processingflag,pattern,leng,dire,latB,lonB,latE,lonE,
         navg,xavg,yavg,len_avg,len_diff,NULL,NULL,
         NULL,NULL,NULL,NULL);
   if (ret) {
      printf("error in writing unfiltered product to netcdf file\n");
   }
   fprintf(logf,"OUTF:<%s/%s>\n",outdir,outfname);

   /* filtering with neighbours. */
   double difflenLimit = 10;
   int invalid_exists = 1;
   size_t nbcorr = 0;
   while (invalid_exists) {
      /* sort by order of diff_length */
      ret = sort_icedrift_bylength(nb_tobeprocessed,tobeprocessed,len_diff);
      /* test the first (and thus worst) candidate */
      size_t p = tobeprocessed[0];
      if (len_diff[p] > difflenLimit) { 
	 
	 if (processingflag[p] == ICEDRIFT_REFUSED_BY_NEIGHBOURS) {
	    len_diff[p] = -3;
	    continue;
	 }

	 nbcorr++;
	 //printf("Start a new inversion for location p:%u (len_diff[p]=%f)\n",p,len_diff[p]);
	 /* correct (modify or delete this vector) */
	 processingflag[p] = ICEDRIFT_REFUSED_BY_NEIGHBOURS; 
	 double currentBestX = xbest[p*NUNKNOWNS+XDRIFT_IDX];
	 double currentBestY = xbest[p*NUNKNOWNS+YDRIFT_IDX];
	 double currentBestC = bestval[p];
	 /*
	 printf("\tDRIFT: (%f,%f){%f} AVG: {%f,%f}\n",
	       xbest[p*NUNKNOWNS+XDRIFT_IDX],xbest[p*NUNKNOWNS+YDRIFT_IDX],bestval[p],
	       xavg[p],yavg[p]);
	       */

	 /* load the subimage in end_drift image */
	 short pat = pattern[p];
	 long xi,yi;
	 locfmijmap(iwcs[p],img_dims[XDIM],&xi,&yi);
	 center_coord[XDIM] = xi; center_coord[YDIM] = yi ; center_coord[TDIM] = iwcs[p];
	 load_subimage(&ret,obs[END],TCflag[END],icelandmask[END],img_dims,
	    center_coord,pattern_size[pat], pattern_mask[pat], pattern_windexes[pat], 
	    pattern_img2[pat], pattern_isvalid2[pat], pattern_icemask2[pat]);

	 /* modify the sigmoid  */
	 set_sigmoidlength(difflenLimit);
	 /* set the correlation function parameters (sigmoid center) */
	 setBestKnowledge(avgLat[p],avgLon[p]);

	 /* perform the optimization */
	 int oret = findBestX_sys(p,pat,xnbbins,xbins,ynbbins,ybins, 
	    x,xbest,bestval, processingflag, 0);
	 if ((processingflag[p] != ICEDRIFT_OK) || (bestval[p] <= 0.)) {
	 //   printf("Optimization failed.\n");
	    processingflag[p] = ICEDRIFT_REFUSED_BY_NEIGHBOURS;
	 } else {
	    processingflag[p] = ICEDRIFT_CORRECT_BY_NEIGHBOURS;
	 }

	 /* recompute the average for all the neighbouring pixels (which do not include 'p')*/
         long xcoordp,ycoordp;
         locfmijmap(p,out_dims[XDIM],&xcoordp,&ycoordp);
         for (size_t n = 0 ; n < nbNeighbours ; n++) {

	    if (neighboursMask[n]) {

	       long xneighbour = xcoordp + neighboursIndexes[XDIM][n];
	       long yneighbour = ycoordp + neighboursIndexes[YDIM][n];
	       
	       if ( !pointIsInGrid((double)(xneighbour),(double)(yneighbour),out_dims) ) continue;

	       long windexNeighbour;
	       locfmivec(&windexNeighbour,xneighbour,yneighbour,out_dims[XDIM]);
	       /* windexNeighbour is the absolute 1D index for the n-th neighbour of 'p' */
	       size_t nbValidNeighbours   = 0;
	       size_t nbInvalidNeighbours = 0;
	       /* recompute the average around windexNeighbour */
	       for (size_t np = 0 ; np < nbNeighbours ; np++) {
		  if (neighboursMask[np]) {

		     long xneighbour_np_n = xneighbour + neighboursIndexes[XDIM][np];
		     long yneighbour_np_n = yneighbour + neighboursIndexes[YDIM][np];
	    
		     if ( !pointIsInGrid((double)(xneighbour_np_n),(double)(yneighbour_np_n),out_dims) ) continue;

		     long windexNeighbour_np;
		     locfmivec(&windexNeighbour_np,xneighbour_np_n,yneighbour_np_n,out_dims[XDIM]);

		     if (  ((processingflag[windexNeighbour_np] == ICEDRIFT_OK) ||
			   (processingflag[windexNeighbour_np] == ICEDRIFT_CORRECT_BY_NEIGHBOURS))
			   /* && (bestval[windexNeighbour_np] >= correlation_limit) */ ) {
			neighboursX[nbValidNeighbours] = xbest[windexNeighbour_np*NUNKNOWNS + XDRIFT_IDX];
			neighboursY[nbValidNeighbours] = xbest[windexNeighbour_np*NUNKNOWNS + YDRIFT_IDX];
			nbValidNeighbours++;
		     } else if (processingflag[windexNeighbour_np] == ICEDRIFT_CLOSETO_MISSING) {
			;
		     } else {
			nbInvalidNeighbours++;
		     }
		  }
	       }
	       navg[windexNeighbour] = nbValidNeighbours;
	       xavg[windexNeighbour] = 1000; yavg[windexNeighbour] = 1000;
	       if (nbValidNeighbours >= 5) { /* TODO: replace by a statistical test: does this average have a meaning */
	 
		  double xavgP,yavgP,xstdP,ystdP;
          ret = compute_mean_vector(&xavgP,&yavgP,&xstdP,&ystdP,nbValidNeighbours,neighboursX,neighboursY);
          xavg[windexNeighbour] = xavgP;yavg[windexNeighbour] = yavgP;
		  navg[windexNeighbour] = nbValidNeighbours;

		  /* get the lat and lon associated with the average drift */
		  long xcoordi,ycoordi;
		  double latAvg,lonAvg;
		  locfmijmap(iwcs[windexNeighbour],img_dims[XDIM],&xcoordi,&ycoordi);
		  //remap_xy2ll(-(xavgP/img_Ax)+(xcoordi+0.5),-(yavgP/img_Ay)+(ycoordi+0.5),img_pj,img_Ax,img_Bx,img_Ay,img_By,&latAvg,&lonAvg);
		  remap_xy2ll(-(xavgP/img_Ax)+xcoordi,-(yavgP/img_Ay)+ycoordi,img_pj,img_Ax,img_Bx,img_Ay,img_By,&latAvg,&lonAvg);
		  avgLat[windexNeighbour] = latAvg;
		  avgLon[windexNeighbour] = lonAvg;

		  /* compute the length of the average drift */
		  double lenavg;
		  compute_distance(latB[windexNeighbour],lonB[windexNeighbour],latAvg,lonAvg,&lenavg);

		  /* compute the length of the difference vector: avg - drift */
		  double lendiff;
		  compute_distance(latAvg,lonAvg,latE[windexNeighbour],lonE[windexNeighbour],&lendiff);

		  /* store those values for writing in the file */
		  len_avg[windexNeighbour]  = lenavg;
		  len_diff[windexNeighbour] = lendiff;

	       } else {
		  len_diff[windexNeighbour] = -2.; 
	       }

	    }
	 }

	 /* transform the new vector (in 'p') to len,dir, lat and lon and recompute its difference to the average in 'p' */
	 if (processingflag[p] == ICEDRIFT_CORRECT_BY_NEIGHBOURS) {
	    double xdrift = xbest[p*NUNKNOWNS+XDRIFT_IDX];
	    double ydrift = xbest[p*NUNKNOWNS+YDRIFT_IDX];;
	    double lat,lon;
	    long xcoordi,ycoordi;
	    locfmijmap(iwcs[p],img_dims[XDIM],&xcoordi,&ycoordi);
	    //remap_xy2ll(-(xdrift/img_Ax)+(xcoordi+0.5),-(ydrift/img_Ay)+(ycoordi+0.5),img_pj,img_Ax,img_Bx,img_Ay,img_By,&lat,&lon);
	    remap_xy2ll(-(xdrift/img_Ax)+xcoordi,-(ydrift/img_Ay)+ycoordi,img_pj,img_Ax,img_Bx,img_Ay,img_By,&lat,&lon);
	    latE[p] = lat; lonE[p] = lon;

	    /* compute the length and direction of the arrow */ 
	    double dlen,ddir;
	    compute_distance(latB[p],lonB[p],latE[p],lonE[p],&dlen);
	    compute_directionToNorth(latB[p],lonB[p],latE[p],lonE[p],&ddir);
	    leng[p] = dlen; dire[p] = ddir;

	    /* compute length of the difference vector to average (the average vector in 'p' is not modified) */
	    double lendiff;
	    compute_distance(avgLat[p],avgLon[p],latE[p],lonE[p],&lendiff);

	    /* store those values for writing in the file */
	    len_diff[p] = lendiff;
	    //printf("p:%u. New distance to average is %f\n",p,len_diff[p]);
	 } else {
	    len_diff[p] = -3;
	    //printf("p:%u. Set len_diff to %f\n",p,len_diff[p]);
	 }

      } else {
	 /* the worst vector passes the test. So we can stop filtering. */
	 invalid_exists=0;
      }
      /*
      if (nbcorr >= 5)
	 invalid_exists=0;
	 */
   }
   
   for (size_t p = 0 ; p < NDRIFTPIXELS ; p++) {
      if ( (processingflag[p] != ICEDRIFT_OK) && (processingflag[p] != ICEDRIFT_CORRECT_BY_NEIGHBOURS) ) continue;
      if ( (navg[p] < 5) || (xavg[p] == 1000)) {
//	 printf("\tp:%u => navg[p]=%d, xavg[p]=%f, processingflag[p]=%d\n",p,navg[p],xavg[p],processingflag[p]);
	 processingflag[p] = ICEDRIFT_NOAVERAGE;
      }
   }

   fmlogmsg(progname,"Finished processing, write the netcdf file.");
   setFilteringLevelString("lev1");
   setProductDescriptionString();
   ret = prepareAndWrite_icedriftProduct(&xdim,xbest,bestval,processingflag,pattern,leng,dire,latB,lonB,latE,lonE,
         navg,xavg,yavg,len_avg,len_diff,NULL,NULL,
         NULL,NULL,NULL,NULL);
   if (ret) {
      printf("error in writing netcdf file\n");
   }
   fprintf(logf,"OUTF:<%s/%s>\n",outdir,outfname);
   
   /* last level of filtering: remove the vectors for which we still have no average vector as
    * well as those with too bad a correlation. */
   for (size_t q = 0 ; q < nb_tobeprocessed ; q++) {
      size_t p = tobeprocessed[q];
      if ((processingflag[p] == ICEDRIFT_OK) && (processingflag[p] != ICEDRIFT_CORRECT_BY_NEIGHBOURS)) continue;
      if (bestval[p] < 0.3) {
	 processingflag[p] = ICEDRIFT_LOWCORRELATION;
      }
   }
   setFilteringLevelString("lev2");
   setProductDescriptionString();
   ret = prepareAndWrite_icedriftProduct(&xdim,xbest,bestval,processingflag,pattern,leng,dire,latB,lonB,latE,lonE,
         navg,xavg,yavg,len_avg,len_diff,NULL,NULL,
         NULL,NULL,NULL,NULL);
   if (ret) {
      printf("error in writing netcdf file\n");
   }
   fprintf(logf,"OUTF:<%s/%s>\n",outdir,outfname);

   fclose(logf);

   for (size_t q = 0 ; q < NDRIFTPIXELS ; q++) {
      if ((processingflag[q] == ICEDRIFT_OK) ||
	    (processingflag[q] == ICEDRIFT_CORRECT_BY_NEIGHBOURS)) {
	 xDrift[q] = xbest[q*NUNKNOWNS+XDRIFT_IDX];
	 yDrift[q] = xbest[q*NUNKNOWNS+YDRIFT_IDX];
      }
   }
   
   /* write an ascii observation file to get a different symbol when the drift is zero */
   ret = prepareAndWrite_zeroDriftObsFile(NDRIFTPIXELS,latB,lonB,xDrift,yDrift,processingflag);

   /* As a last step, compute the differential fluxes (divergence, curl, ...) of the drift field and add them in the file. */
   /*
   float *div, *curl;
   short *def_flg;
   ret = compute_differential_fluxes(NDRIFTPIXELS,out_dims,xDrift,yDrift,processingflag,&div,&curl,&def_flg);
   if (ret) {
      fprintf(stderr,"ERROR in computing the differential fluxes (deformation).\n");
   }
   printf("Write deformation result to netcdf files...\n");
   float As[2],Bs[2];
   As[0] = 62.5; As[1] = 62.5; 
   Bs[0] = -3750; Bs[1] = 5750.0; 
   char oprojstr[] = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45";
   PJ *opj = NULL;
   ret = initproj_from_str(oprojstr, &opj);
   ret = prepareAndWrite_deformationFile(NDRIFTPIXELS,out_dims,def_flg,div,curl,As,Bs,opj,oprojstr);
   */

   free(x);
   free(xbest);

   return EXIT_SUCCESS;
}

int findBestX_sys(size_t p, short pat, size_t xnbbins, double xbins[], size_t ynbbins, double ybins[], 
      double x[], double xbest[], double bestval[], short processingflag[], int plotCorrMap) {
      
   size_t sta = p*NUNKNOWNS;
   double *xloc = &(x[sta]);
   size_t ea = 0;
   bestval[p] = -1;
   for (size_t ex = 0 ; ex < xnbbins ; ex ++) {
      //printf("ex: %d\n",ex);
      for (size_t ey = 0 ; ey < ynbbins ; ey ++) {
	    //printf("\tey: %d\n",ey);

	 double val = 0;
	 double xtry[NUNKNOWNS];
	 x[sta+XDRIFT_IDX] = xbins[ex];
	 x[sta+YDRIFT_IDX] = ybins[ey];
	 model_onevector_corr(&xdim,x,&val,p,pat);

#define PPLOT 999999
	 if (p == PPLOT) {
	    printf("%f %f %f\n",xbins[ex],ybins[ey],val);
	 }

#ifdef GNUPLOT_CORRELATION_MAP
	 /* store the correlation map value */
	 if (plotCorrMap) {
	    corrMap[ex][ey] = val;
	 }
#endif 

	 /* make the selection for this x,y pair */
	 if (val > bestval[p]) {
	    xbest[sta + XDRIFT_IDX] = xloc[XDRIFT_IDX];
	    xbest[sta + YDRIFT_IDX] = xloc[YDRIFT_IDX];
	    bestval[p] = val;
	 }
      }
   }

   if (p == PPLOT) {
      printf("%f %f %f\n",xbest[sta + XDRIFT_IDX],xbest[sta + YDRIFT_IDX],bestval[p]);
   }
#ifdef  GNUPLOT_CORRELATION_MAP
   /* Only for the first vector (if successfully processed) */
   if ( plotCorrMap ) {

      if (processingflag[p]) {
	 printf("WARNING: UNABLE TO PLOT the correlation map because the processing was not successfull.\n");
      } else {

	 printf("Description string is <%s>\n",IceDriftProductDescription);

	 char corrMapDir[] = "/home/thomasl/Documents/MERSEA/data/corrMaps";
	 char corrMapFileName[1028];
	 sprintf(corrMapFileName,"%s/corrmap-%s.txt",corrMapDir,IceDriftProductDescription);

	 printf("Write the correlation map in <%s>\n",corrMapFileName);

	 /* write the correlation map to ascii file */
	 float halfXstep = 0.5 * (xbins[1]-xbins[0]);
	 float halfYstep = 0.5 * (ybins[1]-ybins[0]);
	 FILE *corrMapF = fopen(corrMapFileName,"w");
	 if (!corrMapF) {
	    printf("ERROR Unable to open file <%s>\n",corrMapFileName);
	    exit(1);
	 }
	 fprintf(corrMapF,"# %f %f %f %f\n",xbest[sta + XDRIFT_IDX],-xbest[sta + YDRIFT_IDX],olat[p],olon[p]);
	 for (size_t ex = 0 ; ex < xnbbins ; ex ++) {
	    float valx = xbins[ex] ;
	    for (size_t ey = 0 ; ey < ynbbins ; ey++) {
	       float valy = ybins[ynbbins-ey-1] ;
	       float valcorr = corrMap[ex][ey]/bestval[p];
	       fprintf(corrMapF,"%f %f %+07.4f\n",
		     valx,valy,valcorr);
	    }
	    fprintf(corrMapF,"\n");
	 }
	 fclose(corrMapF);

//	    /* write a gnuplot command file */
//	    float bestx = xbest[sta + XDRIFT_IDX];
//	    float besty = -xbest[sta + YDRIFT_IDX];
//	    float lim = 30;
//	    float xAlow = bestx - lim;
//	    float xAhig = bestx + lim;
//	    float yAlow = besty - lim;
//	    float yAhig = besty + lim;
//
//	    char epsFname[1028];
//	    char epsDir[] = "/tmp";
//	    sprintf(epsFname,"%s/corrmap-%s.eps",epsDir,
//		  IceDriftProductDescription);
//	    
//	    FILE *cmdFile = fopen("/tmp/extractedCorrMap.gnuplot","w");
//            fprintf(cmdFile,"set terminal postscript eps enhanced \"Helvetica\" 20\n");
//            fprintf(cmdFile,"set output '%s'\n",epsFname);
//	    fprintf(cmdFile,"set xrange [%f:%f]\n",xAlow,xAhig);
//	    fprintf(cmdFile,"set yrange [%f:%f]\n",yAlow,yAhig);
//	    fprintf(cmdFile,"set xlabel 'drift X [km]'\n");
//	    fprintf(cmdFile,"set ylabel 'drift Y [km]'\n");
//	    fprintf(cmdFile,"set xzeroaxis lt 0 lw 0.3\n");
//	    fprintf(cmdFile,"set yzeroaxis lt 0 lw 0.3\n");
//	    fprintf(cmdFile,"set nokey\n");
//	    //fprintf(cmdFile,"set grid lw 3.0\n");
//	    fprintf(cmdFile,"set contour\n");
//	    fprintf(cmdFile,"unset surface\n");
//	    fprintf(cmdFile,"set style data lines\n");
//	    fprintf(cmdFile,"set contour base\nset cntrparam bspline\nset cntrparam levels incremental 0.6,0.05,1.0\n");
//	    fprintf(cmdFile,"set size square\n");
//	    fprintf(cmdFile,"set view map\n");
//	    fprintf(cmdFile,"set arrow  from 0,0,0 to %f,%f,0 front\n",bestx,besty);
//	    fprintf(cmdFile,"set multiplot\n");
//	    fprintf(cmdFile,"splot '/tmp/extractedCorrMap.dat' using 1:2:3\n");
//	    //fprintf(cmdFile,"set parametric\n");
//	    //fprintf(cmdFile,"plot [t=-pi:pi] %f*sin(t),%f*cos(t)\n",maxdriftdistance,maxdriftdistance);
//
//	    fclose(cmdFile);
//
//	    /* launch gunplot */
//	    ret = system("gnuplot /tmp/extractedCorrMap.gnuplot");
//	    printf("Image %s is ready.\n",epsFname); 

      }
	   
   }

#endif  /* GNUPLOT_CORRELATION_MAP */

   return 0;

}
