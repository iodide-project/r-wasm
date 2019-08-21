/* stxwx.f -- translated by f2c (version 20190311).
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
static integer c__1 = 1;

/* Output from Public domain Ratfor, version 1.0 */
/* Subroutine */ int stxwx_(doublereal *x, doublereal *z__, doublereal *w, 
	integer *k, doublereal *xknot, integer *n, doublereal *y, doublereal *
	hs0, doublereal *hs1, doublereal *hs2, doublereal *hs3)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal eps, work[16];
    static integer mflag, ileft, lenxk;
    static doublereal vnikx[4]	/* was [4][1] */;
    extern /* Subroutine */ int bsplvd_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    extern integer interv_(doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);

/*      implicit none */
/* local */

/* in ../../../appl/interv.c */
    /* Parameter adjustments */
    --w;
    --z__;
    --x;
    --hs3;
    --hs2;
    --hs1;
    --hs0;
    --y;
    --xknot;

    /* Function Body */
    lenxk = *n + 4;
/*     Initialise the output vectors */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.;
	hs0[i__] = 0.;
	hs1[i__] = 0.;
	hs2[i__] = 0.;
	hs3[i__] = 0.;
    }
/* Compute X' W^2 X -> hs0,hs1,hs2,hs3  and X' W^2 Z -> y */
/* Note that here the weights w(i) == sqrt(wt[i])  where wt[] where original weights */
    ileft = 1;
    eps = 1e-10;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n + 1;
	ileft = interv_(&xknot[1], &i__2, &x[i__], &c__0, &c__0, &ileft, &
		mflag);
/*        if(mflag==-1) {write(6,'("Error in hess ",i2)')mflag;stop} */
/*        if(mflag==-1) {return} */
	if (mflag == 1) {
	    if (x[i__] <= xknot[ileft] + eps) {
		--ileft;
	    } else {
		return 0;
	    }
/*        else{write(6,'("Error in hess ",i2)')mflag;stop}} */
	}
	bsplvd_(&xknot[1], &lenxk, &c__4, &x[i__], &ileft, work, vnikx, &c__1)
		;
	j = ileft - 3;
/* Computing 2nd power */
	d__1 = w[i__];
	y[j] += d__1 * d__1 * z__[i__] * vnikx[0];
/* Computing 2nd power */
	d__1 = w[i__];
/* Computing 2nd power */
	d__2 = vnikx[0];
	hs0[j] += d__1 * d__1 * (d__2 * d__2);
/* Computing 2nd power */
	d__1 = w[i__];
	hs1[j] += d__1 * d__1 * vnikx[0] * vnikx[1];
/* Computing 2nd power */
	d__1 = w[i__];
	hs2[j] += d__1 * d__1 * vnikx[0] * vnikx[2];
/* Computing 2nd power */
	d__1 = w[i__];
	hs3[j] += d__1 * d__1 * vnikx[0] * vnikx[3];
	j = ileft - 2;
/* Computing 2nd power */
	d__1 = w[i__];
	y[j] += d__1 * d__1 * z__[i__] * vnikx[1];
/* Computing 2nd power */
	d__1 = w[i__];
/* Computing 2nd power */
	d__2 = vnikx[1];
	hs0[j] += d__1 * d__1 * (d__2 * d__2);
/* Computing 2nd power */
	d__1 = w[i__];
	hs1[j] += d__1 * d__1 * vnikx[1] * vnikx[2];
/* Computing 2nd power */
	d__1 = w[i__];
	hs2[j] += d__1 * d__1 * vnikx[1] * vnikx[3];
	j = ileft - 1;
/* Computing 2nd power */
	d__1 = w[i__];
	y[j] += d__1 * d__1 * z__[i__] * vnikx[2];
/* Computing 2nd power */
	d__1 = w[i__];
/* Computing 2nd power */
	d__2 = vnikx[2];
	hs0[j] += d__1 * d__1 * (d__2 * d__2);
/* Computing 2nd power */
	d__1 = w[i__];
	hs1[j] += d__1 * d__1 * vnikx[2] * vnikx[3];
	j = ileft;
/* Computing 2nd power */
	d__1 = w[i__];
	y[j] += d__1 * d__1 * z__[i__] * vnikx[3];
/* Computing 2nd power */
	d__1 = w[i__];
/* Computing 2nd power */
	d__2 = vnikx[3];
	hs0[j] += d__1 * d__1 * (d__2 * d__2);
    }
    return 0;
} /* stxwx_ */

