/* kmns.f -- translated by f2c (version 20190311).
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

/* Code in this file based on Applied Statistics algorithms */
/* (C) Royal Statistical Society 1979 */

/* a minimal modification of AS136 to use double precision */
/* all variables are now declared. */
/* B.D. Ripley 1998/06/17 */
/* Declaration re-ordering to satisfy "f77 -ansi",  M.Maechler 2001/04/12 */

/* ~= R's  kmeans(x=A, centers=C, iter.max=ITER, algorithm = "Hartigan-Wong") */

/* Subroutine */ int kmns_(doublereal *a, integer *m, integer *n, doublereal *
	c__, integer *k, integer *ic1, integer *ic2, integer *nc, doublereal *
	an1, doublereal *an2, integer *ncp, doublereal *d__, integer *itran, 
	integer *live, integer *iter, doublereal *wss, integer *ifault)
{
    /* Initialized data */

    static doublereal big = 1e30;
    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l;
    static doublereal aa, da, db, dc;
    static integer ii, ij, il;
    static doublereal dt[2];
    static integer indx;
    static doublereal temp;
    extern /* Subroutine */ int kmns1_(integer *, integer *, integer *), 
	    optra_(doublereal *, integer *, integer *, doublereal *, integer *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *), qtran_(
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *);
    static integer itrace, imaxqtr;


/*     ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1 */

/*     Divide M points in N-dimensional space into K clusters so that */
/*     the within cluster sum of squares is minimized. */

/*                      ------        ------ */
/*                     data x[,]      centers[,] */

/*     Define BIG to be a very large positive number */

    /* Parameter adjustments */
    --d__;
    --ic2;
    --ic1;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --wss;
    --live;
    --itran;
    --ncp;
    --an2;
    --an1;
    --nc;
    c_dim1 = *k;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */

    itrace = *ifault;
    imaxqtr = itran[1];
    *ifault = 3;
    if (*k <= 1 || *k >= *m) {
	return 0;
    }
    *ifault = 0;

/*     For each point I, find its two closest centres, IC1(I) and */
/*     IC2(I).     Assign it to IC1(I). */

    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ic1[i__] = 1;
	ic2[i__] = 2;
	for (il = 1; il <= 2; ++il) {
	    dt[il - 1] = zero;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		da = a[i__ + j * a_dim1] - c__[il + j * c_dim1];
		dt[il - 1] += da * da;
	    }
	}
/* IL */
	if (dt[0] > dt[1]) {
	    ic1[i__] = 2;
	    ic2[i__] = 1;
	    temp = dt[0];
	    dt[0] = dt[1];
	    dt[1] = temp;
	}
	i__2 = *k;
	for (l = 3; l <= i__2; ++l) {
	    db = zero;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		dc = a[i__ + j * a_dim1] - c__[l + j * c_dim1];
		db += dc * dc;
		if (db >= dt[1]) {
		    goto L50;
		}
	    }
	    if (db >= dt[0]) {
		dt[1] = db;
		ic2[i__] = l;
	    } else {
		dt[1] = dt[0];
		ic2[i__] = ic1[i__];
		dt[0] = db;
		ic1[i__] = l;
	    }
L50:
	    ;
	}
/* L60: */
    }

/*     Update cluster centres to be the average of points contained */
/*     within them. */
/*     NC(L) := #{units in cluster L},  L = 1..K */
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	nc[l] = 0;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    c__[l + j * c_dim1] = zero;
	}
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = ic1[i__];
	++nc[l];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    c__[l + j * c_dim1] += a[i__ + j * a_dim1];
	}
    }

/*     Check to see if there is any empty cluster at this stage */

    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	if (nc[l] == 0) {
	    *ifault = 1;
	    return 0;
	}
	aa = (doublereal) nc[l];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    c__[l + j * c_dim1] /= aa;
	}

/*     Initialize AN1, AN2, ITRAN & NCP */
/*     AN1(L) = NC(L) / (NC(L) - 1) */
/*     AN2(L) = NC(L) / (NC(L) + 1) */
/*     ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage, */
/*              = 0 otherwise */
/*     In the optimal-transfer stage, NCP(L) stores the step at which */
/*     cluster L is last updated. */
/*     In the quick-transfer stage, NCP(L) stores the step at which */
/*     cluster L is last updated plus M. */

	an2[l] = aa / (aa + one);
	an1[l] = big;
	if (aa > one) {
	    an1[l] = aa / (aa - one);
	}
	itran[l] = 1;
	ncp[l] = -1;
    }
    indx = 0;
    i__1 = *iter;
    for (ij = 1; ij <= i__1; ++ij) {

/* OPtimal-TRAnsfer stage: there is only one pass through the data. */
/*     Each point is re-allocated, if necessary, to the cluster that will */
/*     induce the maximum reduction in within-cluster sum of squares. */

	optra_(&a[a_offset], m, n, &c__[c_offset], k, &ic1[1], &ic2[1], &nc[1]
		, &an1[1], &an2[1], &ncp[1], &d__[1], &itran[1], &live[1], &
		indx);
	if (itrace > 0) {
	    kmns1_(k, &ij, &indx);
	}

/*     Stop if no transfer took place in the last M optimal transfer steps. */
	if (indx == *m) {
	    goto L150;
	}

/* Quick-TRANSfer stage: Each point is tested in turn to see if it should */
/*     be re-allocated to the cluster to which it is most likely to be */
/*     transferred, IC2(I), from its present cluster, IC1(I). */
/*     Loop through the data until no further change is to take place. */

	qtran_(&a[a_offset], m, n, &c__[c_offset], k, &ic1[1], &ic2[1], &nc[1]
		, &an1[1], &an2[1], &ncp[1], &d__[1], &itran[1], &indx, &
		itrace, &imaxqtr);

	if (imaxqtr < 0) {
	    *ifault = 4;
	    goto L150;
	}

/*     If there are only two clusters, there is no need to re-enter the */
/*     optimal transfer stage. */

	if (*k == 2) {
	    goto L150;
	}

/*     NCP has to be set to 0 before entering OPTRA. */

	i__2 = *k;
	for (l = 1; l <= i__2; ++l) {
	    ncp[l] = 0;
	}
    }

/*     Since the specified number of iterations has been exceeded, set */
/*     IFAULT = 2.   This may indicate unforeseen looping. */

/* iter -------------------------------------- */
    *ifault = 2;
L150:
    *iter = ij;

/*     Compute within-cluster sum of squares for each cluster. */

    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	wss[l] = zero;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    c__[l + j * c_dim1] = zero;
	}
    }
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = ic1[i__];
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    c__[ii + j * c_dim1] += a[i__ + j * a_dim1];
	}
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *k;
	for (l = 1; l <= i__2; ++l) {
	    c__[l + j * c_dim1] /= (doublereal) nc[l];
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ii = ic1[i__];
	    da = a[i__ + j * a_dim1] - c__[ii + j * c_dim1];
	    wss[ii] += da * da;
	}
    }

    return 0;
} /* kmns_ */



/* Subroutine */ int optra_(doublereal *a, integer *m, integer *n, doublereal 
	*c__, integer *k, integer *ic1, integer *ic2, integer *nc, doublereal 
	*an1, doublereal *an2, integer *ncp, doublereal *d__, integer *itran, 
	integer *live, integer *indx)
{
    /* Initialized data */

    static doublereal big = 1e30;
    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, l1, l2;
    static doublereal r2, da, db, dc, dd, de, df;
    static integer ll;
    static doublereal rr, al1, al2, alt, alw;


/*     ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1 */

/*     This is the OPtimal TRAnsfer stage. */
/*                 ---------------------- */
/*     Each point is re-allocated, if necessary, to the cluster that */
/*     will induce a maximum reduction in the within-cluster sum of */
/*     squares. */


    /* Parameter adjustments */
    --d__;
    --ic2;
    --ic1;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --live;
    --itran;
    --ncp;
    --an2;
    --an1;
    --nc;
    c_dim1 = *k;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */

/*     If cluster L is updated in the last quick-transfer stage, it */
/*     belongs to the live set throughout this stage.   Otherwise, at */
/*     each step, it is not in the live set if it has not been updated */
/*     in the last M optimal transfer steps. */

/* BIG := a very large positive number */
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	if (itran[l] == 1) {
	    live[l] = *m + 1;
	}
    }
/*     ---------------------- Loop over each point ------------------- */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	++(*indx);
	l1 = ic1[i__];
	l2 = ic2[i__];
	ll = l2;

/*     If point I is the only member of cluster L1, no transfer. */

	if (nc[l1] == 1) {
	    goto L90;
	}

/*     If L1 has not yet been updated in this stage, no need to */
/*     re-compute D(I). */

	if (ncp[l1] != 0) {
	    de = zero;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		df = a[i__ + j * a_dim1] - c__[l1 + j * c_dim1];
		de += df * df;
	    }
	    d__[i__] = de * an1[l1];
	}

/*     Find the cluster with minimum R2. */

	da = zero;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    db = a[i__ + j * a_dim1] - c__[l2 + j * c_dim1];
	    da += db * db;
	}
	r2 = da * an2[l2];
	i__2 = *k;
	for (l = 1; l <= i__2; ++l) {

/*     If I >= LIVE(L1), then L1 is not in the live set.   If this is */
/*     true, we only need to consider clusters that are in the live set */
/*     for possible transfer of point I.   Otherwise, we need to consider */
/*     all possible clusters. */

	    if (i__ >= live[l1] && i__ >= live[l] || l == l1 || l == ll) {
		goto L60;
	    }
	    rr = r2 / an2[l];
	    dc = zero;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		dd = a[i__ + j * a_dim1] - c__[l + j * c_dim1];
		dc += dd * dd;
		if (dc >= rr) {
		    goto L60;
		}
	    }
	    r2 = dc * an2[l];
	    l2 = l;
L60:
	    ;
	}
	if (r2 >= d__[i__]) {
/*     If no transfer is necessary, L2 is the new IC2(I). */
	    ic2[i__] = l2;
	} else {

/*     Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and */
/*     L2, and update IC1(I) & IC2(I). */
	    *indx = 0;
	    live[l1] = *m + i__;
	    live[l2] = *m + i__;
	    ncp[l1] = i__;
	    ncp[l2] = i__;
	    al1 = (doublereal) nc[l1];
	    alw = al1 - one;
	    al2 = (doublereal) nc[l2];
	    alt = al2 + one;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		c__[l1 + j * c_dim1] = (c__[l1 + j * c_dim1] * al1 - a[i__ + 
			j * a_dim1]) / alw;
		c__[l2 + j * c_dim1] = (c__[l2 + j * c_dim1] * al2 + a[i__ + 
			j * a_dim1]) / alt;
	    }
	    --nc[l1];
	    ++nc[l2];
	    an2[l1] = alw / al1;
	    an1[l1] = big;
	    if (alw > one) {
		an1[l1] = alw / (alw - one);
	    }
	    an1[l2] = alt / al2;
	    an2[l2] = alt / (alt + one);
	    ic1[i__] = l2;
	    ic2[i__] = l1;
	}
L90:
	if (*indx == *m) {
	    return 0;
	}
    }
/*     ---------------------- each point ------------------- */
    i__1 = *k;
    for (l = 1; l <= i__1; ++l) {
	itran[l] = 0;
/*       LIVE(L) has to be decreased by M before re-entering OPTRA */
/* before entering QTRAN. */
	live[l] -= *m;
    }

    return 0;
} /* optra_ */



/* Subroutine */ int qtran_(doublereal *a, integer *m, integer *n, doublereal 
	*c__, integer *k, integer *ic1, integer *ic2, integer *nc, doublereal 
	*an1, doublereal *an2, integer *ncp, doublereal *d__, integer *itran, 
	integer *indx, integer *itrace, integer *imaxqtr)
{
    /* Initialized data */

    static doublereal big = 1e30;
    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, l1, l2;
    static doublereal r2, da, db, dd, de, al1, al2, alt, alw;
    static integer icoun, istep;
    extern /* Subroutine */ int rchkusr_(void), kmnsqpr_(integer *, integer *,
	     integer *, integer *, integer *);


/*     ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1 */

/*     This is the Quick TRANsfer stage. */
/*                 -------------------- */
/*     IC1(I) is the cluster which point I belongs to. */
/*     IC2(I) is the cluster which point I is most likely to be */
/*         transferred to. */
/*     For each point I, IC1(I) & IC2(I) are switched, if necessary, to */
/*     reduce within-cluster sum of squares.  The cluster centres are */
/*     updated after each step. */


/*     Define BIG to be a very large positive number */

    /* Parameter adjustments */
    --d__;
    --ic2;
    --ic1;
    a_dim1 = *m;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --itran;
    --ncp;
    --an2;
    --an1;
    --nc;
    c_dim1 = *k;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */

/*     In the quick transfer stage, NCP(L) */
/*     is equal to the step at which cluster L is last updated plus M. */

    icoun = 0;
    istep = 0;
/*   Repeat { */
L10:
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itrace > 0 && istep >= 1 && i__ == 1) {
	    kmnsqpr_(&istep, &icoun, &ncp[1], k, itrace);
	}
/* only from */
	++icoun;
	++istep;
	if (istep >= *imaxqtr) {
	    *imaxqtr = -1;
	    return 0;
	}
	l1 = ic1[i__];
	l2 = ic2[i__];

/*     If point I is the only member of cluster L1, no transfer. */

	if (nc[l1] == 1) {
	    goto L60;
	}

/*     If ISTEP > NCP(L1), no need to re-compute distance from point I to */
/*     cluster L1.   Note that if cluster L1 is last updated exactly M */
/*     steps ago, we still need to compute the distance from point I to */
/*     cluster L1. */

	if (istep <= ncp[l1]) {
	    da = zero;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		db = a[i__ + j * a_dim1] - c__[l1 + j * c_dim1];
		da += db * db;
	    }
	    d__[i__] = da * an1[l1];
	}

/*     If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer of */
/*     point I at this step. */

	if (istep < ncp[l1] || istep < ncp[l2]) {
	    r2 = d__[i__] / an2[l2];
	    dd = zero;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		de = a[i__ + j * a_dim1] - c__[l2 + j * c_dim1];
		dd += de * de;
		if (dd >= r2) {
		    goto L60;
		}
	    }

/*     Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters */
/*     L1 & L2.   Also update IC1(I) & IC2(I).   Note that if any */
/*     updating occurs in this stage, INDX is set back to 0. */

	    icoun = 0;
	    *indx = 0;
	    itran[l1] = 1;
	    itran[l2] = 1;
	    ncp[l1] = istep + *m;
	    ncp[l2] = istep + *m;
	    al1 = (doublereal) nc[l1];
	    alw = al1 - one;
	    al2 = (doublereal) nc[l2];
	    alt = al2 + one;
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		c__[l1 + j * c_dim1] = (c__[l1 + j * c_dim1] * al1 - a[i__ + 
			j * a_dim1]) / alw;
		c__[l2 + j * c_dim1] = (c__[l2 + j * c_dim1] * al2 + a[i__ + 
			j * a_dim1]) / alt;
	    }
	    --nc[l1];
	    ++nc[l2];
	    an2[l1] = alw / al1;
	    an1[l1] = big;
	    if (alw > one) {
		an1[l1] = alw / (alw - one);
	    }
	    an1[l2] = alt / al2;
	    an2[l2] = alt / (alt + one);
	    ic1[i__] = l2;
	    ic2[i__] = l1;
	}

/*     If no re-allocation took place in the last M steps, return. */

L60:
	if (icoun == *m) {
	    return 0;
	}
    }
    rchkusr_();
/* allow user interrupt */
    goto L10;
/*     -------- */
} /* qtran_ */

