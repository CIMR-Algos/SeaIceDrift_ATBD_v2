"""CFFI building code to wrap the ice drift core C code into a library
accessible from a Python wrapper function."""

import os
from cffi import FFI
ffibuilder = FFI()

# Code paths
icedrift_path = os.getcwd()

ffibuilder.set_source("_idcore", # name of the output C extension
"""
    # include <stdio.h>
    # include <stdlib.h>
    # include <string.h>
    # include <float.h>
    # include <math.h>
    # include <proj.h>
    # include <assert.h>
    # include "fmextra.h"
    # include "fmerrmsg.h"
    # include "fmutil_config.h"
    # include "icedrift_flags.h"
    # include "icedrift_common.h"
    # include "icedrift_model.h"
    # include "icedrift_solve_filter.h"
    # include "optimization_simplex.h"
    # include "vector_anydim.h"
    # include "memory.h"
    # include "icedrift_solve_common.h"
    # include "icedrift_common.h"
    # include "icedrift_prepost.h"
    # include "icedrift_solve_core.h"

""",
    include_dirs=[os.path.join(icedrift_path, 'src_simplex'), '/usr/include'],

    sources=['icedrift_solve_core.c',
             'fmerrmsg.c',
             'fmextra.c',
             'icedrift_instruments.c',
             'icedrift_model.c',
             'icedrift_common.c',
             'icedrift_prepost.c',
             'icedrift_solve_common.c',
             'icedrift_solve_filter.c',
             os.path.join(icedrift_path, 'src_simplex/optimization_simplex.c'),
             os.path.join(icedrift_path, 'src_simplex/vector_anydim.c'),
             os.path.join(icedrift_path, 'src_simplex/memory.c')
            ],
    define_macros=[('ACCEPT_USE_OF_DEPRECATED_PROJ_API_H', '1')],
    libraries=['c', 'proj'])


# cdef() expects a single string declaring the C types, functions and
# globals, in valid C syntax.
ffibuilder.cdef("""
                extern double **obs[2];
                extern short **TCflag[2];
                extern short *icelandmask[2];
                extern double *img_lat;
                extern double *img_lon;
                extern double img_Ax;
                extern double img_Bx;
                extern double img_Ay;
                extern double img_By;
                extern size_t img_dims[3];
                extern char *img_projstr;
                extern short nbWaveBands;
                extern short twghtStart;
                extern short twghtEnd;
                extern double maxdriftdistance;
                extern double sigmoid_length;
                extern char *OptimMetric;
                extern double pattern_radius[2];
                extern double radiusNeighbours;
                extern int NDRIFTPIXELS;
                extern char *out_area;
                extern char *out_projstr;
                extern double out_Ax;
                extern double out_Bx;
                extern double out_Ay;
                extern double out_By;
                extern size_t out_dims[3];
                extern double *olat;
                extern double *olon;
                extern unsigned int *owcs;
                extern unsigned int *iwcs;
                extern float *driftX;
                extern float *driftY;
                extern short *pflag;
                extern float *sigdX;
                extern float *sigdY;
                extern float *corrdXdY;
                extern short *uflag;
                extern float *length;
                extern float *dir;
                extern float *outlonB;
                extern float *outlatB;
                extern float *outlonE;
                extern float *outlatE;
                extern float *outfc;
                extern short *snavg;
                extern float *avgX;
                extern float *avgY;
                extern float *length_avg;
                extern float *length_diff;
                extern float *stdX;
                extern float *stdY;
                extern short *patternIndex;
                extern char *reportFile;

                int core(void);
                void *malloc(size_t size);
                void free(void *ptr);
                """)


if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
