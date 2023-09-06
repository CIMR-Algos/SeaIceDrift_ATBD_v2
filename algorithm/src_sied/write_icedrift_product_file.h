#ifndef WRITE_ICEDRIFT_PRODUCT_FILE_H
#define WRITE_ICEDRIFT_PRODUCT_FILE_H

int write_icedriftproduct_tonetcdf(char outdir[],char outfname[],size_t out_dims[],
      float *driftX, float *driftY, fmsec1970 *t0, fmsec1970 *t1, short *flag,
      float *sigdX,  float *sigdY,  float *corrdXdY, short *uflag,
      float *length, float *dir, float *outlatB, float *outlonB, float *outlatE, float *outlonE,
      float *corr,short *navg,float *avgX, float *avgY, float *length_avg,float *length_diff, float *stdX, float *stdY, float fillvalf,
      float **driftX_exp, float **driftY_exp, short *flag_exp, short nbflag_exp,
      short *pattern, short fillvals, int fillvali,
      char Instrument[], char Platform[], char isostartdate[], char isoenddate[], 
      char out_gridf[], char out_area[], char out_projstr[], float out_A[], float out_B[], PJ *out_pj);

#endif /* WRITE_ICEDRIFT_PRODUCT_FILE_H */
