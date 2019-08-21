/* bvalus.f -- translated by f2c (version 20190311).
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

static integer c__4 = 4;

/* Subroutine */ int bvalus_(integer *n, doublereal *knot, doublereal *coef, 
	integer *nk, doublereal *x, doublereal *s, integer *order)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    extern doublereal bvalue_(doublereal *, doublereal *, integer *, integer *
	    , doublereal *, integer *);

/* Args */
/* Local */
    /* Parameter adjustments */
    --s;
    --x;
    --coef;
    --knot;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s[i__] = bvalue_(&knot[1], &coef[1], nk, &c__4, &x[i__], order);
/* L10: */
    }
    return 0;
} /* bvalus_ */

