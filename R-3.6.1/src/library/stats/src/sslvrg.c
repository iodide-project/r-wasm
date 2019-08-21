/* sslvrg.f -- translated by f2c (version 20190311).
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

static integer c__3 = 3;
static integer c__4 = 4;
static integer c__0 = 0;
static integer c__1 = 1;

/* Output from Public domain Ratfor, version 1.0 */
/* Smoothing Spline LeVeRaGes = SSLVRG */
/* ----------------------------------- leverages = H_ii = diagonal entries of Hat matrix */
/* Subroutine */ int sslvrg_(doublereal *penalt, doublereal *dofoff, 
	doublereal *x, doublereal *y, doublereal *w, doublereal *ssw, integer 
	*n, doublereal *knot, integer *nk, doublereal *coef, doublereal *sz, 
	doublereal *lev, doublereal *crit, integer *icrit, doublereal *lambda,
	 doublereal *xwy, doublereal *hs0, doublereal *hs1, doublereal *hs2, 
	doublereal *hs3, doublereal *sg0, doublereal *sg1, doublereal *sg2, 
	doublereal *sg3, doublereal *abd, doublereal *p1ip, doublereal *p2ip, 
	integer *ld4, integer *ldnk, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, p1ip_dim1, p1ip_offset, p2ip_dim1, 
	    p2ip_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static integer i__, j;
    static doublereal b0, b1, b2, b3, df, xv, eps, rss, work[16], sumw;
    extern /* Subroutine */ int dpbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *);
    static integer mflag, ileft;
    extern /* Subroutine */ int dpbsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *);
    static doublereal vnikx[4]	/* was [4][1] */;
    extern doublereal bvalue_(doublereal *, doublereal *, integer *, integer *
	    , doublereal *, integer *);
    static integer lenkno;
    extern /* Subroutine */ int bsplvd_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    sinerp_(doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *);
    extern integer interv_(doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);

/* Purpose : */
/*       Compute smoothing spline for smoothing parameter lambda */
/*       and compute one of three `criteria' (OCV , GCV , "df match"). */
/* See comments in ./sbart.c from which this is called */
/* local variables */

/* in ../../../appl/interv.c */
    /* Parameter adjustments */
    --lev;
    --sz;
    --w;
    --y;
    --x;
    --sg3;
    --sg2;
    --sg1;
    --sg0;
    --hs3;
    --hs2;
    --hs1;
    --hs0;
    --xwy;
    --coef;
    --knot;
    p1ip_dim1 = *ld4;
    p1ip_offset = 1 + p1ip_dim1;
    p1ip -= p1ip_offset;
    abd_dim1 = *ld4;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    p2ip_dim1 = *ldnk;
    p2ip_offset = 1 + p2ip_dim1;
    p2ip -= p2ip_offset;

    /* Function Body */
    lenkno = *nk + 4;
    ileft = 1;
    eps = 1e-11;
/* compute the coefficients coef() of estimated smooth */
    i__1 = *nk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	coef[i__] = xwy[i__];
	abd[i__ * abd_dim1 + 4] = hs0[i__] + *lambda * sg0[i__];
    }
    i__1 = *nk - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	abd[(i__ + 1) * abd_dim1 + 3] = hs1[i__] + *lambda * sg1[i__];
    }
    i__1 = *nk - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	abd[(i__ + 2) * abd_dim1 + 2] = hs2[i__] + *lambda * sg2[i__];
    }
    i__1 = *nk - 3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	abd[(i__ + 3) * abd_dim1 + 1] = hs3[i__] + *lambda * sg3[i__];
    }
/*     factorize banded matrix abd (into upper triangular): */
    dpbfa_(&abd[abd_offset], ld4, nk, &c__3, info);
    if (*info != 0) {
/*        matrix could not be factorized -> ier := info */
	return 0;
    }
/*     solve linear system (from factorized abd): */
    dpbsl_(&abd[abd_offset], ld4, nk, &c__3, &coef[1]);
/*     Value of smooth at the data points */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xv = x[i__];
	sz[i__] = bvalue_(&knot[1], &coef[1], nk, &c__4, &xv, &c__0);
    }
/*     Compute the criterion function if requested (icrit > 0) : */
    if (*icrit >= 1) {
/* --- Ordinary or Generalized CV or "df match" --- */
/*     Get Leverages First */
	sinerp_(&abd[abd_offset], ld4, nk, &p1ip[p1ip_offset], &p2ip[
		p2ip_offset], ldnk, &c__0);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xv = x[i__];
	    i__2 = *nk + 1;
	    ileft = interv_(&knot[1], &i__2, &xv, &c__0, &c__0, &ileft, &
		    mflag);
	    if (mflag == -1) {
		ileft = 4;
		xv = knot[4] + eps;
	    } else if (mflag == 1) {
		ileft = *nk;
		xv = knot[*nk + 1] - eps;
	    }
	    j = ileft - 3;
/*           call bspvd(knot,4,1,xv,ileft,4,vnikx,work) */
	    bsplvd_(&knot[1], &lenkno, &c__4, &xv, &ileft, work, vnikx, &c__1)
		    ;
	    b0 = vnikx[0];
	    b1 = vnikx[1];
	    b2 = vnikx[2];
	    b3 = vnikx[3];
/* Computing 2nd power */
	    d__1 = b0;
/* Computing 2nd power */
	    d__2 = b1;
/* Computing 2nd power */
	    d__3 = b2;
/* Computing 2nd power */
	    d__4 = b3;
/* Computing 2nd power */
	    d__5 = w[i__];
	    lev[i__] = (p1ip[j * p1ip_dim1 + 4] * (d__1 * d__1) + p1ip[j * 
		    p1ip_dim1 + 3] * 2. * b0 * b1 + p1ip[j * p1ip_dim1 + 2] * 
		    2. * b0 * b2 + p1ip[j * p1ip_dim1 + 1] * 2. * b0 * b3 + 
		    p1ip[(j + 1) * p1ip_dim1 + 4] * (d__2 * d__2) + p1ip[(j + 
		    1) * p1ip_dim1 + 3] * 2. * b1 * b2 + p1ip[(j + 1) * 
		    p1ip_dim1 + 2] * 2. * b1 * b3 + p1ip[(j + 2) * p1ip_dim1 
		    + 4] * (d__3 * d__3) + p1ip[(j + 2) * p1ip_dim1 + 3] * 2. 
		    * b2 * b3 + p1ip[(j + 3) * p1ip_dim1 + 4] * (d__4 * d__4))
		     * (d__5 * d__5);
	}
/*     Evaluate Criterion */
	df = 0.;
	if (*icrit == 1) {
/* Generalized CV -------------------- */
	    rss = *ssw;
	    sumw = 0.;
/*       w(i) are sqrt( wt[i] ) weights scaled in ../R/smspline.R such */
/*       that sumw =  number of observations with w(i) > 0 */
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = (y[i__] - sz[i__]) * w[i__];
		rss += d__1 * d__1;
		df += lev[i__];
/* Computing 2nd power */
		d__1 = w[i__];
		sumw += d__1 * d__1;
	    }
/* Computing 2nd power */
	    d__1 = 1. - (*dofoff + *penalt * df) / sumw;
	    *crit = rss / sumw / (d__1 * d__1);
/*            call dblepr("spar", 4, spar, 1) */
/*            call dblepr("crit", 4, crit, 1) */
	} else if (*icrit == 2) {
/* Ordinary CV ------------------ */
	    *crit = 0.;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = (y[i__] - sz[i__]) * w[i__] / (1 - lev[i__]);
		*crit += d__1 * d__1;
	    }
	    *crit /= *n;
/*            call dblepr("spar", 4, spar, 1) */
/*            call dblepr("crit", 4, crit, 1) */
	} else {
/* df := sum( lev[i] ) */
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		df += lev[i__];
	    }
	    if (*icrit == 3) {
/* df matching -------------------- */
/* Computing 2nd power */
		d__1 = *dofoff - df;
		*crit = d__1 * d__1 + 3;
	    } else {
/* if(icrit .eq. 4) then df - dofoff (=> zero finding) */
		*crit = df - *dofoff;
	    }
	}
    }
/*     Criterion evaluation */
    return 0;
} /* sslvrg_ */

