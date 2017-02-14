/* Vitali Morozov <morozov@anl.gov>, Argonne Leadership Computing Facility */
#include <math.h>
#include <stdio.h>
#include "../../../../source/flashUtilities/general/mangle_names.h"

// Version for sizeof( unsigned long long ) = 8
double FTOC(sqrt3)( double *x )
{
    double d, d3, dr;
    unsigned long long *px = (unsigned long long *)x, *py = (unsigned long long *) &d;

    /* initial approximation, 5 bits of precision */
    *py = *px / 3 + 0x2aa0000000000000;

    /* 1 Halley's iteration -> 15 bits of precision */
    d3 = d * d * d;
    dr = d3 + *x;
    d = d * ( dr + *x ) / ( d3 + dr);

    /* 2 Halley's iterations -> 45 bits of precision */
    d3 = d * d * d;
    dr = d3 + *x;
    d = d * ( dr + *x ) / ( d3 + dr);

    return d;
}
