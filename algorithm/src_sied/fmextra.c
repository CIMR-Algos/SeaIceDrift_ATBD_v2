
#include <stdlib.h>
#include <stdio.h>
#include "fmerrmsg.h"

/* Copied from fmstorage.c by Øystein Godøy */
long fmivec(long x, long y, unsigned long nx) {
    long i;

    i = y*nx+x;
    return(i);
}
