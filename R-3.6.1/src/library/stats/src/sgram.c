/* sgram.f -- translated by f2c (version 20190311).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__0 = 0;
static integer c__4 = 4;
static integer c__3 = 3;

/* Output from Public domain Ratfor, version 1.0 */
/* PURPOSE */
/*       Calculation of the cubic B-spline smoothness prior */
/*       for "usual" interior knot setup. */
/*       Uses BSPVD and INTRV in the CMLIB */
/*       sgm[0-3](nb)    Symmetric matrix 'SIGMA' */
/*                       whose (i,j)'th element contains the integral of */
/*                       B''(i,.) B''(j,.) , i=1,2 ... nb and j=i,...nb. */
/*                       Only the upper four diagonals are computed. */
/* Subroutine */ int sgram_(doublereal *sg0, doublereal *sg1, doublereal *sg2,
	 doublereal *sg3, doublereal *tb, integer *nb)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, ii, jj;
    static doublereal yw1[4], yw2[4], wpt, work[16];
    static integer mflag, ileft, lentb;
    static doublereal vnikx[12]	/* was [4][3] */;
    extern /* Subroutine */ int bsplvd_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern integer interv_(doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);

/*      implicit none */
/* indices */
/*     ------------- */

/* in ../../../appl/interv.c */
    /* Parameter adjustments */
    --tb;
    --sg3;
    --sg2;
    --sg1;
    --sg0;

    /* Function Body */
    lentb = *nb + 4;
/* Initialise the sigma vectors */
    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sg0[i__] = 0.;
	sg1[i__] = 0.;
	sg2[i__] = 0.;
	sg3[i__] = 0.;
    }
    ileft = 1;
    i__1 = *nb;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*        Calculate a linear approximation to the second derivative of the */
/*        non-zero B-splines over the interval [tb(i),tb(i+1)]. */
	i__2 = *nb + 1;
	ileft = interv_(&tb[1], &i__2, &tb[i__], &c__0, &c__0, &ileft, &mflag)
		;
/*        Left end second derivatives */
	bsplvd_(&tb[1], &lentb, &c__4, &tb[i__], &ileft, work, vnikx, &c__3);
/*        Put values into yw1 */
	for (ii = 1; ii <= 4; ++ii) {
	    yw1[ii - 1] = vnikx[ii + 7];
	}
/*        Right end second derivatives */
	bsplvd_(&tb[1], &lentb, &c__4, &tb[i__ + 1], &ileft, work, vnikx, &
		c__3);
/*        Slope*(length of interval) in Linear Approximation to B'' */
	for (ii = 1; ii <= 4; ++ii) {
	    yw2[ii - 1] = vnikx[ii + 7] - yw1[ii - 1];
	}
/*        Calculate Contributions to the sigma vectors */
	wpt = tb[i__ + 1] - tb[i__];
	if (ileft >= 4) {
	    for (ii = 1; ii <= 4; ++ii) {
		jj = ii;
		sg0[ileft - 4 + ii] += wpt * (yw1[ii - 1] * yw1[jj - 1] + (
			yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1])
			 * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		jj = ii + 1;
		if (jj <= 4) {
		    sg1[ileft + ii - 4] += wpt * (yw1[ii - 1] * yw1[jj - 1] + 
			    (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii 
			    - 1]) * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		}
		jj = ii + 2;
		if (jj <= 4) {
		    sg2[ileft + ii - 4] += wpt * (yw1[ii - 1] * yw1[jj - 1] + 
			    (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii 
			    - 1]) * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		}
		jj = ii + 3;
		if (jj <= 4) {
		    sg3[ileft + ii - 4] += wpt * (yw1[ii - 1] * yw1[jj - 1] + 
			    (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii 
			    - 1]) * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		}
	    }
	} else if (ileft == 3) {
	    for (ii = 1; ii <= 3; ++ii) {
		jj = ii;
		sg0[ileft - 3 + ii] += wpt * (yw1[ii - 1] * yw1[jj - 1] + (
			yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1])
			 * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		jj = ii + 1;
		if (jj <= 3) {
		    sg1[ileft + ii - 3] += wpt * (yw1[ii - 1] * yw1[jj - 1] + 
			    (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii 
			    - 1]) * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		}
		jj = ii + 2;
		if (jj <= 3) {
		    sg2[ileft + ii - 3] += wpt * (yw1[ii - 1] * yw1[jj - 1] + 
			    (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii 
			    - 1]) * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		}
	    }
	} else if (ileft == 2) {
	    for (ii = 1; ii <= 2; ++ii) {
		jj = ii;
		sg0[ileft - 2 + ii] += wpt * (yw1[ii - 1] * yw1[jj - 1] + (
			yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1])
			 * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		jj = ii + 1;
		if (jj <= 2) {
		    sg1[ileft + ii - 2] += wpt * (yw1[ii - 1] * yw1[jj - 1] + 
			    (yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii 
			    - 1]) * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
		}
	    }
	} else if (ileft == 1) {
	    for (ii = 1; ii <= 1; ++ii) {
		jj = ii;
		sg0[ileft - 1 + ii] += wpt * (yw1[ii - 1] * yw1[jj - 1] + (
			yw2[ii - 1] * yw1[jj - 1] + yw2[jj - 1] * yw1[ii - 1])
			 * .5 + yw2[ii - 1] * yw2[jj - 1] * .333);
	    }
	}
    }
    return 0;
} /* sgram_ */

