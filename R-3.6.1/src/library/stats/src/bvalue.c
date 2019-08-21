/* bvalue.f -- translated by f2c (version 20190311).
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

doublereal bvalue_(doublereal *t, doublereal *bcoef, integer *n, integer *k, 
	doublereal *x, integer *jderiv)
{
    /* Initialized data */

    static integer i__ = 1;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static integer j;
    static doublereal aj[20];
    static integer jc;
    static doublereal dm[20], dp[20];
    static integer jj, km1, imk, kmj, ilo, nmi;
    static doublereal fkmj;
    static integer mflag, jcmin, jcmax;
    //    extern /* Subroutine */ int rwarn_(char *, ftnlen);
    static integer jdrvp1;
    extern integer interv_(doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *);

/* Calculates value at  x  of  jderiv-th derivative of spline from B-repr. */
/* The spline is taken to be continuous from the right. */

/* calls  interv()  (from ../../../appl/interv.c ) */

/* ******  i n p u t ****** */
/*  t, bcoef, n, k......forms the b-representation of the spline  f  to */
/*        be evaluated. specifically, */
/*  t.....knot sequence, of length  n+k, assumed nondecreasing. */
/*  bcoef.....b-coefficient sequence, of length  n . */
/*  n.....length of  bcoef  and dimension of s(k,t), */
/*        a s s u m e d  positive . */
/*  k.....order of the spline . */

/*  w a r n i n g . . .   the restriction  k <= kmax (=20)  is imposed */
/*        arbitrarily by the dimension statement for  aj, dm, dm  below, */
/*        but is  n o w h e r e  c h e c k e d  for. */
/*  however in R, this is only called from bvalus() with k=4 anyway! */

/*  x.....the point at which to evaluate . */
/*  jderiv.....integer giving the order of the derivative to be evaluated */
/*        a s s u m e d  to be zero or positive. */

/* ******  o u t p u t  ****** */
/*  bvalue.....the value of the (jderiv)-th derivative of  f  at  x . */

/* ******  m e t h o d  ****** */
/*     the nontrivial knot interval  (t(i),t(i+1))  containing  x  is lo- */
/*  cated with the aid of  interv(). the  k  b-coeffs of  f  relevant for */
/*  this interval are then obtained from  bcoef (or taken to be zero if */
/*  not explicitly available) and are then differenced  jderiv  times to */
/*  obtain the b-coeffs of  (d^jderiv)f  relevant for that interval. */
/*  precisely, with  j = jderiv, we have from x.(12) of the text that */

/*     (d^j)f  =  sum ( bcoef(.,j)*b(.,k-j,t) ) */

/*  where */
/*                   / bcoef(.),                     ,  j .eq. 0 */
/*                   / */
/*    bcoef(.,j)  =  / bcoef(.,j-1) - bcoef(.-1,j-1) */
/*                   / ----------------------------- ,  j > 0 */
/*                   /    (t(.+k-j) - t(.))/(k-j) */

/*     then, we use repeatedly the fact that */

/*    sum ( a(.)*b(.,m,t)(x) )  =  sum ( a(.,x)*b(.,m-1,t)(x) ) */
/*  with */
/*                 (x - t(.))*a(.) + (t(.+m-1) - x)*a(.-1) */
/*    a(.,x)  =    --------------------------------------- */
/*                 (x - t(.))      + (t(.+m-1) - x) */

/*  to write  (d^j)f(x)  eventually as a linear combination of b-splines */
/*  of order  1 , and the coefficient for  b(i,1,t)(x)  must then */
/*  be the desired number  (d^j)f(x). (see x.(17)-(19) of text). */

/* Arguments */
/*     dimension t(n+k) */
/*  current fortran standard makes it impossible to specify the length of */
/*  t  precisely without the introduction of otherwise superfluous */
/*  additional arguments. */
/* Local Variables */

/*     initialize */
    /* Parameter adjustments */
    --t;
    --bcoef;

    /* Function Body */
    ret_val = 0.;
    if (*jderiv >= *k) {
	goto L99;
    }

/*  *** find  i  s.t.  1 <= i < n+k  and  t(i) < t(i+1) and */
/*      t(i) <= x < t(i+1) . if no such i can be found,  x  lies */
/*      outside the support of the spline  f and bvalue = 0. */
/*  {this case is handled in the calling R code} */
/*      (the asymmetry in this choice of  i  makes  f  rightcontinuous) */
    if (*x != t[*n + 1] || t[*n + 1] != t[*n + *k]) {
	i__1 = *n + *k;
	i__ = interv_(&t[1], &i__1, x, &c__0, &c__0, &i__, &mflag);
	if (mflag != 0) {
          //rwarn_("bvalue()  mflag != 0: should never happen!", (ftnlen)42);
	    goto L99;
	}
    } else {
	i__ = *n;
    }
/*  *** if k = 1 (and jderiv = 0), bvalue = bcoef(i). */
    km1 = *k - 1;
    if (km1 <= 0) {
	ret_val = bcoef[i__];
	goto L99;
    }

/*  *** store the k b-spline coefficients relevant for the knot interval */
/*     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dm(j) = x - t(i+1-j), */
/*     dp(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable */
/*     from input to zero. set any t.s not obtainable equal to t(1) or */
/*     to t(n+k) appropriately. */
    jcmin = 1;
    imk = i__ - *k;
    if (imk >= 0) {
	i__1 = km1;
	for (j = 1; j <= i__1; ++j) {
	    dm[j - 1] = *x - t[i__ + 1 - j];
/* L9: */
	}
    } else {
	jcmin = 1 - imk;
	i__1 = i__;
	for (j = 1; j <= i__1; ++j) {
	    dm[j - 1] = *x - t[i__ + 1 - j];
/* L5: */
	}
	i__1 = km1;
	for (j = i__; j <= i__1; ++j) {
	    aj[*k - j - 1] = 0.;
	    dm[j - 1] = dm[i__ - 1];
/* L6: */
	}
    }

    jcmax = *k;
    nmi = *n - i__;
    if (nmi >= 0) {
	i__1 = km1;
	for (j = 1; j <= i__1; ++j) {
/*     the following if() happens; e.g. in   pp <- predict(cars.spl, xx) */
/*     -       if( (i+j) .gt. lent) write(6,9911) i+j,lent */
/*     -  9911         format(' i+j, lent ',2(i6,1x)) */
	    dp[j - 1] = t[i__ + j] - *x;
/* L19: */
	}
    } else {
	jcmax = *k + nmi;
	i__1 = jcmax;
	for (j = 1; j <= i__1; ++j) {
	    dp[j - 1] = t[i__ + j] - *x;
/* L15: */
	}
	i__1 = km1;
	for (j = jcmax; j <= i__1; ++j) {
	    aj[j] = 0.;
	    dp[j - 1] = dp[jcmax - 1];
/* L16: */
	}
    }

    i__1 = jcmax;
    for (jc = jcmin; jc <= i__1; ++jc) {
	aj[jc - 1] = bcoef[imk + jc];
/* L21: */
    }

/*               *** difference the coefficients  jderiv  times. */
    if (*jderiv >= 1) {
	i__1 = *jderiv;
	for (j = 1; j <= i__1; ++j) {
	    kmj = *k - j;
	    fkmj = (doublereal) kmj;
	    ilo = kmj;
	    i__2 = kmj;
	    for (jj = 1; jj <= i__2; ++jj) {
		aj[jj - 1] = (aj[jj] - aj[jj - 1]) / (dm[ilo - 1] + dp[jj - 1]
			) * fkmj;
		--ilo;
/* L24: */
	    }
/* L23: */
	}
    }

/*  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative, */
/*     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv). */
    if (*jderiv != km1) {
	jdrvp1 = *jderiv + 1;
	i__1 = km1;
	for (j = jdrvp1; j <= i__1; ++j) {
	    kmj = *k - j;
	    ilo = kmj;
	    i__2 = kmj;
	    for (jj = 1; jj <= i__2; ++jj) {
		aj[jj - 1] = (aj[jj] * dm[ilo - 1] + aj[jj - 1] * dp[jj - 1]) 
			/ (dm[ilo - 1] + dp[jj - 1]);
		--ilo;
/* L34: */
	    }
/* L33: */
	}
    }
    ret_val = aj[0];

L99:
    return ret_val;
} /* bvalue_ */

