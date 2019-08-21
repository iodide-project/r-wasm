/* dqrutl.f -- translated by f2c (version 20190311).
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

static integer c__1000 = 1000;
static integer c__10000 = 10000;
static integer c__100 = 100;
static integer c__10 = 10;
static integer c__1 = 1;

/* dqr Utilities:  Interface to the different "switches" of  dqrsl(). */

/* Subroutine */ int dqrqty_(doublereal *x, integer *n, integer *k, 
	doublereal *qraux, doublereal *y, integer *ny, doublereal *qty)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, qty_dim1, qty_offset, i__1;

    /* Local variables */
    static integer j, info;
    extern /* Subroutine */ int dqrsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal dummy[1];

    /* Parameter adjustments */
    --qraux;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    qty_dim1 = *n;
    qty_offset = 1 + qty_dim1;
    qty -= qty_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
	dqrsl_(&x[x_offset], n, n, k, &qraux[1], &y[j * y_dim1 + 1], dummy, &
		qty[j * qty_dim1 + 1], dummy, dummy, dummy, &c__1000, &info);
/* L10: */
    }
    return 0;
} /* dqrqty_ */


/* Subroutine */ int dqrqy_(doublereal *x, integer *n, integer *k, doublereal 
	*qraux, doublereal *y, integer *ny, doublereal *qy)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, qy_dim1, qy_offset, i__1;

    /* Local variables */
    static integer j, info;
    extern /* Subroutine */ int dqrsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal dummy[1];

    /* Parameter adjustments */
    --qraux;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    qy_dim1 = *n;
    qy_offset = 1 + qy_dim1;
    qy -= qy_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
	dqrsl_(&x[x_offset], n, n, k, &qraux[1], &y[j * y_dim1 + 1], &qy[j * 
		qy_dim1 + 1], dummy, dummy, dummy, dummy, &c__10000, &info);
/* L10: */
    }
    return 0;
} /* dqrqy_ */


/* Subroutine */ int dqrcf_(doublereal *x, integer *n, integer *k, doublereal 
	*qraux, doublereal *y, integer *ny, doublereal *b, integer *info)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int dqrsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal dummy[1];

    /* Parameter adjustments */
    --qraux;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    b_dim1 = *k;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
	dqrsl_(&x[x_offset], n, n, k, &qraux[1], &y[j * y_dim1 + 1], dummy, &
		y[j * y_dim1 + 1], &b[j * b_dim1 + 1], dummy, dummy, &c__100, 
		info);
/* L10: */
    }
    return 0;
} /* dqrcf_ */


/* Subroutine */ int dqrrsd_(doublereal *x, integer *n, integer *k, 
	doublereal *qraux, doublereal *y, integer *ny, doublereal *rsd)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, rsd_dim1, rsd_offset, i__1;

    /* Local variables */
    static integer j, info;
    extern /* Subroutine */ int dqrsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal dummy[1];

    /* Parameter adjustments */
    --qraux;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    rsd_dim1 = *n;
    rsd_offset = 1 + rsd_dim1;
    rsd -= rsd_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
	dqrsl_(&x[x_offset], n, n, k, &qraux[1], &y[j * y_dim1 + 1], dummy, &
		y[j * y_dim1 + 1], dummy, &rsd[j * rsd_dim1 + 1], dummy, &
		c__10, &info);
/* L10: */
    }
    return 0;
} /* dqrrsd_ */


/* Subroutine */ int dqrxb_(doublereal *x, integer *n, integer *k, doublereal 
	*qraux, doublereal *y, integer *ny, doublereal *xb)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, xb_dim1, xb_offset, i__1;

    /* Local variables */
    static integer j, info;
    extern /* Subroutine */ int dqrsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal dummy[1];

    /* Parameter adjustments */
    --qraux;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    xb_dim1 = *n;
    xb_offset = 1 + xb_dim1;
    xb -= xb_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
	dqrsl_(&x[x_offset], n, n, k, &qraux[1], &y[j * y_dim1 + 1], dummy, &
		y[j * y_dim1 + 1], dummy, dummy, &xb[j * xb_dim1 + 1], &c__1, 
		&info);
/* L10: */
    }
    return 0;
} /* dqrxb_ */

