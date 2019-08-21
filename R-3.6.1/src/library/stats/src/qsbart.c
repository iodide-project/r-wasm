/* qsbart.f -- translated by f2c (version 20190311).
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

/* ----------------------------------------------------------------------- */

/*  R : A Computer Language for Statistical Data Analysis */
/*  Copyright (C) 1998-2016 The R Core Team */

/*  This program is free software; you can redistribute it and/or modify */
/*  it under the terms of the GNU General Public License as published by */
/*  the Free Software Foundation; either version 2 of the License, or */
/*  (at your option) any later version. */

/*  This program is distributed in the hope that it will be useful, */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*  GNU General Public License for more details. */

/*  You should have received a copy of the GNU General Public License */
/*  along with this program; if not, a copy is available at */
/*  https://www.R-project.org/Licenses/ */

/* ----------------------------------------------------------------------- */
/* Called from R's smooth.spline in ../R/smspline.R  as .Fortran(C, ..) */
/*    and from C's */
/* An interface to sbart() --- fewer arguments BUT unspecified scrtch() dimension */

/* NB: this routine alters ws [and isetup]. */
/* renamed for safety */

/* Subroutine */ int rbart_(doublereal *penalt, doublereal *dofoff, 
	doublereal *xs, doublereal *ys, doublereal *ws, doublereal *ssw, 
	integer *n, doublereal *knot, integer *nk, doublereal *coef, 
	doublereal *sz, doublereal *lev, doublereal *crit, integer *iparms, 
	doublereal *spar, doublereal *parms, doublereal *scrtch, integer *ld4,
	 integer *ldnk, integer *ier)
{
    extern /* Subroutine */ int sbart_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *);
    static integer isetup;

/* Args: */
/*          ^^^^^^^^ dimension (9+2*ld4+ldnk)*nk = (17 + 1)*nk [last nk never accessed] */
/* Vars: */
    /* Parameter adjustments */
    --lev;
    --sz;
    --ws;
    --ys;
    --xs;
    --coef;
    --knot;
    --iparms;
    --parms;
    --scrtch;

    /* Function Body */
    if (iparms[4] == 1) {
/* spar is lambda */
	isetup = 2;
    } else {
	isetup = 0;
    }
    sbart_(penalt, dofoff, &xs[1], &ys[1], &ws[1], ssw, n, &knot[1], nk, &
	    coef[1], &sz[1], &lev[1], crit, &iparms[1], spar, &iparms[2], &
	    iparms[3], &parms[1], &parms[2], &parms[3], &parms[4], &parms[5], 
	    &isetup, &scrtch[1], &scrtch[*nk + 1], &scrtch[(*nk << 1) + 1], &
	    scrtch[*nk * 3 + 1], &scrtch[(*nk << 2) + 1], &scrtch[*nk * 5 + 1]
	    , &scrtch[*nk * 6 + 1], &scrtch[*nk * 7 + 1], &scrtch[(*nk << 3) 
	    + 1], &scrtch[*nk * 9 + 1], &scrtch[*nk * 9 + *ld4 * *nk + 1], &
	    scrtch[*nk * 9 + (*ld4 << 1) * *nk + 1], ld4, ldnk, ier);
/*          = icrit   spar   ispar    iter */
/*          = lspar   uspar    tol      eps      ratio */
/*          = 0|2    xwy  == X'W y */
/*          =   hs0	      hs1	     hs2	    hs3		==> X'W X */
/*          =   sg0	      sg1	     sg2	    sg3		==> SIGMA */
/*          =   abd [ld4 x nk]						==> R */
/*          =   p1ip[ld4 x nk]          p2ip [ldnk x nk] */
    return 0;
} /* rbart_ */

