/*
 *  tree/src/grow.c by B. D. Ripley  Copyright (C) 1994-2005
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 or 3 of the License
 *  (at your option).
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 */
#include <R.h>

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("tree", String)
#else
#define _(String) (String)
#endif

static double XLOGX(double x) 
{
    return (x > 0)?x*log(x):0;
}

#define True 1
#define False 0
#define EPS 1e-4
#define NALEVEL -99999
#ifndef max
# define max(a, b) ((a>b)?a:b)
#endif

#define DEBUG False
#define Printf if (DEBUG) printf

static void scat(char *s, char c)
{
    int i;
    i = strlen(s);
    s[i++] = c;
    s[i] = '\0';
}

static double *X, *y, *w, *dev, *yval, *yprob, mindev,  devtarget,
    *tvar, *cprob, *scprob, *tyc, *w1;
static int  nobs, nvar, minsize, mincut, nnode, nmax,*twhere, *ttw, *ty, Gini;
static Sint *levels, *node, *var, *where, *ordered;

static char **cutleft, **cutright;
static int nc, *indl, *indr, *ind, exists, offset, maxnl;
static double *yp;
static double *tab, *cnt, *n, *ys;


static void fillin_node(int inode)
{
    int     j, k, nl, yparent;
    double  yl, sum, t, n1;
    char   *labl, *labr;

    labl = (char *) S_alloc(maxnl, sizeof(char));
    labr = (char *) S_alloc(maxnl, sizeof(char));
    *labl = *labr = '\0';
    cutleft[inode] = labl;
    cutright[inode] = labr;
    var[inode] = 0;
    if (nc) {
	n1 = 0;
	for (k = 0; k < nc; k++) yprob[nc * inode + k] = 0.0;
	for (j = 0; j < nobs; j++)
	    if (where[j] == inode) {
		n1 += w[j];
		yprob[nc * inode + (int) y[j]-1] += w[j];
	    }
	n[inode] = n1;
	yparent = -1;
	if (inode > 0) {
	    for(j = 0; j < inode; j++) 
		if(node[j] == node[inode]/2) yparent = y[j] - 1;
	}
	nl = 0;
	yl = -1.0;
	for (k = 0; k < nc; k++) {
/*	 if (yprob[nc*inode + k] > yl) {
	 nl = k;
	 yl = yprob[nc * inode + k];
	 } */
	    if (yprob[nc*inode + k] >= yl) {
		if (yprob[nc*inode + k] == yl) {
		    if(k == yparent) nl = k;
		} else {
		    nl = k;
		    yl = yprob[nc * inode + k];
		}
	    }
	    if (n1 > 0) yprob[nc * inode + k] /= n1;
	    else yprob[nc * inode + k] = 1.0/nc;
	}
/*for(k = 0; k < nc; k++) Printf(" %g", yprob[nc * inode + k]); Printf("\n");*/
	nl++;
	if(inode >= exists + offset) yval[inode] = nl;
	sum = 0.0;
	for (j = 0; j < nobs; j++)
	    if (where[j] == inode)
		sum += w[j] * log(yprob[nc * inode + (int) y[j] - 1]);
	dev[inode] = -2 * sum;
    }
    else {
	n1 = 0;
	sum = 0.0;
	for (j = 0; j < nobs; j++)
	    if (where[j] == inode) {
		n1 += w[j];
		sum += w[j] * y[j];
	    }
	n[inode] = n1;
	t = sum / n1;
	yval[inode] = t;
	sum = 0.0;
	for (j = 0; j < nobs; j++)
	    if (where[j] == inode) sum +=  w[j] * (y[j] - t) * (y[j] - t);
	dev[inode] = sum;
    }
}


/* Check split. If best so far, write splits, var, left partition */

/* corrected from R. Sedgewick 'Algorithms in C' */
static void shellsort(double *a, int *b, double *w, int N)
{
    int     i, j, h, at;
    double  v, wt;
    for (h = 1; h <= N / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
	for (i = h; i < N; i++) {
	    v = a[i]; at = b[i]; wt = w[i];
	    j = i;
	    while (j >= h && a[j - h] > v) {
		a[j] = a[j - h]; b[j] = b[j - h]; w[j] = w[j - h];
		j -= h;
	    }
	    a[j] = v; b[j] = at; w[j] = wt;
	}
}
static void shelldsort(double *a, double *b, double *w, int N)
{
    int     i, j, h;
    double  v, at, wt;
    for (h = 1; h <= N / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
	for (i = h; i < N; i++) {
	    v = a[i]; at = b[i]; wt = w[i];
	    j = i;
	    while (j >= h && a[j - h] > v) {
		a[j] = a[j - h]; b[j] = b[j - h]; w[j] = w[j - h];
		j -= h;
	    }
	    a[j] = v; b[j] = at; w[j] = wt;
	}
}


static void split_cont(int inode, int iv, double *bval)
{
    int     i,j, js, k, ns, lo, hi;
    double  ldev, bdev, sdev, tmp, split, bsplit, cntl, totw, ysum = 0.0,
	    ytot = 0.0, y2 = 0.0, psum;

    Printf("..trying split on var %d ", iv);
    ns = 0;
    sdev = 0.0;
    totw = 0.0;
    for (j = 0; j < nobs; j++)
	if (where[j] == inode) {
	    tmp = X[j + nobs * iv]; 
	    if (!ISNA(tmp)) {
		if (nc) ty[ns] = y[j] - 1;
		else tyc[ns] = y[j];
		w1[ns] = w[j];
		tvar[ns++] = tmp;
		totw += w[j];
	    } else {
		if (nc) sdev -= 2*w[j]*log(yprob[nc * inode + (int) y[j] - 1]);
		else {
		    tmp = y[j] - yval[inode];
		    sdev += w[j]*tmp*tmp;
		}
	    }
	}
    if ( Gini && sdev > 0) 
	error(_("cannot use 'Gini' with missing values"));
    Printf(" count %d", ns);
    if ( ns < 2 || totw < EPS ) { Printf("\n"); return;}
    cntl = 0;
    if (nc) {
	shellsort(tvar, ty, w1, ns);
	for (k = 0; k < 2 * nc; k++)
	    tab[k] = 0;		/* left then right cnt */
    } else {
	shelldsort(tvar, tyc, w1, ns);
	ysum = ytot = y2 = 0.0;
	for (j = 0; j < ns; j++) {
	    ytot += w1[j]*tyc[j];
	    y2 += w1[j]*tyc[j]*tyc[j];
	}
    }
    psum = 0.0;
    lo = hi = -1;
    for(i = 0; i < ns; i++) {
	psum += w1[i];
	if(lo < 0 && psum >= mincut) lo = i;
	if(hi < 0 && ns - psum <= mincut) hi = i;
    }
    lo = mincut - 1;
    hi = ns - mincut;
    js = lo;
    tmp = tvar[js];
    if (tvar[ns - 1] == tmp)
    {Printf("\n"); return;}			/* no valid split */
    while (tvar[js + 1] == tmp)
	js++;
    if (js >= hi) {Printf("\n"); return;}
    split = 0.5 * (tmp + tvar[js + 1]);
    for (j = 0; j < ns; j++)
	if (tvar[j] < split) {
	    cntl += w1[j];
	    if (!nc) ysum += w1[j]*tyc[j];
	    else tab[ty[j]] += w1[j];
	} else  if (nc) tab[ty[j] + nc] += w1[j];
    if (nc) {
	if (Gini) {
	    ysum = 0.0;
	    for (k = 0; k < nc; k++) {
		tmp = tab[k] / cntl;
		ysum += tmp*tmp;
	    }
	    ldev = totw - cntl*ysum;
	    ysum = 0.0;
	    for (k = 0; k < nc; k++) {
		tmp = tab[k + nc] / (totw-cntl);
		ysum += tmp*tmp;
	    }
	    ldev = (totw - cntl)*ysum;
	} else {
	    ldev = XLOGX(cntl) + XLOGX((totw - cntl));
	    for (k = 0; k < nc; k++) {
		ldev -= XLOGX(tab[k]) + XLOGX(tab[k + nc]);
	    }
	}
	ldev *= 2;
    } else {
	ldev = y2 - ysum*ysum/cntl - (ytot-ysum)*(ytot-ysum)/(totw-cntl);
    }

/*   Printf("split %g dev %g counts %g %g\n", split, ldev, cntl, totw-cntl);*/
    bdev = ldev;
    bsplit = split;
    while (js < hi - 1) {
	js++;
	tmp = tvar[js];
	if (tvar[ns - 1] == tmp) break;
	cntl += w[js];
	if (nc) {
	    tab[ty[js]] += w1[js];
	    tab[ty[js] + nc] -= w1[js];
	} else ysum += w1[js]*tyc[js];
	while (tvar[js + 1] == tmp) {
	    js++;
	    cntl += w[js];
	    if (nc) {
		tab[ty[js]] += w1[js];
		tab[ty[js] + nc] -= w1[js];
	    } else ysum += w1[js]*tyc[js];
	}
	if (js >= hi) break;
	split = 0.5 * (tmp + tvar[js + 1]);
	if (nc) {
	    if (Gini) {
		ysum = 0.0;
		for (k = 0; k < nc; k++) {
		    tmp = tab[k] / cntl;
		    ysum += tmp*tmp;
		}
		ldev = totw - cntl*ysum;
		ysum = 0.0;
		for (k = 0; k < nc; k++) {
		    tmp = tab[k + nc] / (totw-cntl);
		    ysum += tmp*tmp;
		}
		ldev -= (totw - cntl) * ysum;
	    } else {
		ldev = XLOGX(cntl) + XLOGX((totw - cntl));
		for (k = 0; k < nc; k++)
		    ldev -= XLOGX(tab[k]) + XLOGX(tab[k + nc]);
	    }
	    ldev *= 2;
	} else {
	    ldev = y2 - ysum*ysum/cntl - (ytot-ysum)*(ytot-ysum)/(totw-cntl);
	}
/*  Printf("split %g dev %g counts %g %g\n", split, ldev, cntl, totw-cntl);*/
	if (ldev < bdev) {
	    bdev = ldev;
	    bsplit = split;
	}
    }
    bdev = bdev + sdev;
    Printf(" val %f, split %g\n", bdev, bsplit);
    if (bdev >= *bval) return;
    if (bdev >= devtarget) return;
    *bval = bdev;
    var[inode] = iv + 1;
    sprintf(cutleft[inode], "<%g", bsplit);
    sprintf(cutright[inode], ">%g", bsplit);
    for (j = 0; j < nobs; j++)
	if (where[j] == inode) {
	    tmp = X[j + nobs * iv];
	    if (ISNA(tmp)) ttw[j] = NALEVEL;
	    else ttw[j] =  tmp > bsplit;
	} else ttw[j] = -1;
}


static void split_disc(int inode, int iv, double *bval)
{
    int     i, ii, iis, j, k, l, mi, nl = levels[iv], nll;
    double  bdev, ldev, sdev, val, fence, bfence, cntl, 
	cntr, cntl1, cntr1, tmp, ysum, ytot, y2;
    char    *labl, *labr;

    Printf("..trying split on var %d ", iv);

    for (l = 0; l < nl; l++) ind[l] = False;
    sdev = 0.0;
    for (j = 0; j < nobs; j++)
	if (where[j] == inode) {
	    if (ISNA(X[j + nobs * iv])) {
		twhere[j] = NALEVEL;
		if (nc) sdev -= 2*w[j]*log(yprob[nc * inode + (int) y[j] - 1]);
		else {
		    tmp = y[j] - yval[inode];
		    sdev += w[j]*tmp*tmp;
		}
	    } else {
		twhere[j] = (int) X[j + nobs * iv] - 1;
		if (w[j] > 0) ind[twhere[j]] = True;
	    }
	} else twhere[j] = -1;
    if ( Gini && sdev > 0) 
	error(_("cannot use 'Gini' with missing values"));
    ytot = y2 = 0.0;
    nll = 0;
    for (l = 0; l < nl; l++) {
	nll += ind[l];
	if (ind[l]) ind[nll - 1] = l;
	cnt[l] = 0;
	if (nc) for (k = 0; k < nc; k++) tab[k + nc * l] = 0;
	else ys[l] = 0;
    }
    if (nll < 2) {
	Printf(" no split\n");
	return;
    }

    for (j = 0; j < nobs; j++) {
	l = twhere[j];
	if (l >= 0) {
	    cnt[l] += w[j];
	    if (nc) tab[(int) y[j] - 1 + nc * l] += w[j];
	    else {
		ys[l] += w[j] * y[j];
		y2 += w[j] * y[j] * y[j];
		ytot += w[j] * y[j];
	    }
	}
    }
    /* remove empty levels */
    for (l = 0; l < nll; l++) {
	cnt[l] = cnt[ind[l]];
	if (nc) for (k = 0; k < nc; k++) 
	    tab[k + nc * l] = tab[k + nc * ind[l]];
	else ys[l] = ys[ind[l]];
    }
    if (nll == 2) {  /* Only 2 levels */
	Printf(" counts %g %g", cnt[0], cnt[1]);
	for (l = 0; l < nll; l++)
	    if (cnt[l] < mincut) {
		Printf("\n");
		return;
	    }
	ldev = 0.0;
	for (l = 0; l < nll; l++) {
	    if (nc) {
		if (Gini) {
		    ysum = 0.0;
		    for (k = 0; k < nc; k++) {
			tmp = tab[k + nc * l] / cnt[l];
			ysum += tmp*tmp;
		    }
		    ldev += cnt[l]*(1 - ysum);
		} else {
		    for (k = 0; k < nc; k++) ldev -= XLOGX(tab[k + nc * l]);
		    ldev += XLOGX(cnt[l]);
		}
	    } else {
		ldev += ys[l]*ys[l]/cnt[l];
	    }
	}
	if (!nc) ldev = y2 - ldev;
	else ldev *= 2;
	val = ldev + sdev;
	Printf(" val %f\n", val);
	if (val >= devtarget || val >= *bval) return;
	*bval = val;
	var[inode] = iv + 1;
	labl = cutleft[inode];
	labr = cutright[inode];
	strcpy(labl, ":");
	strcpy(labr, ":");
	scat(labl, 'a' + ind[0]);
	scat(labr, 'a' + ind[1]);
	for (j = 0; j < nobs; j++)
	    if (twhere[j] < 0) ttw[j] = twhere[j];
	    else ttw[j] = (twhere[j] != ind[0]);

    } else {

	/* Treat 2 levels  and regression separately, also ordered */

	if (nc <= 2 || ordered[iv]) {
	    if(ordered[iv]) {
		for (l = 0; l < nll; l++) scprob[l] = cprob[l] = l;
	    } else {
		if (nc) {
		    for (l = 0; l < nll; l++) {
			cprob[l] = (double) tab[1+nc*l] / cnt[l];
			scprob[l] = cprob[l];
		    }
		    shellsort(scprob, indl, w1, nll);
		} else {
		    for (l = 0; l < nll; l++) {
			cprob[l] = ys[l]/cnt[l];
			scprob[l] = cprob[l];
		    }
		    shellsort(scprob, indl, w1, nll);
		}
	    }
	    bdev = devtarget;
	    bfence = -1;

	    Printf(" cnts "); for(l = 0; l < nll; l++) Printf(" %g", cnt[l]);

	    for (i = 1; i < nll; i++) {
		fence = scprob[i];
		if (scprob[i-1] == fence) continue;
		cntl = cntr = 0;
		for (l = 0; l < nll; l++)
		    if (cprob[l] < fence) cntl += cnt[l];
		    else cntr += cnt[l];
		if (cntl < mincut || cntr < mincut) continue;
		if (nc) {
		    if (Gini) {
			ldev = n[inode];
			for (k = 0; k < nc; k++) {
			    cntl1 = cntr1 = 0;
			    for (l = 0; l < nll; l++)
				if (cprob[l] < fence) cntl1 += tab[k + nc * l];
				else cntr1 += tab[k + nc * l];
			    ldev -= cntl1*cntl1/cntl + cntr1*cntr1/cntr;
			}
		    } else {
			ldev = XLOGX(cntl) + XLOGX(cntr);
			for (k = 0; k < nc; k++) {
			    cntl1 = cntr1 = 0;
			    for (l = 0; l < nll; l++)
				if (cprob[l] < fence) cntl1 += tab[k + nc * l];
				else cntr1 += tab[k + nc * l];
			    ldev -= XLOGX(cntl1) + XLOGX(cntr1);
			}
		    }
		    ldev *= 2;
		} else {
		    ysum = 0.0;
		    for (l = 0; l < nll; l++)
			if (cprob[l] < fence) ysum += ys[l];
		    ldev = y2 - ysum*ysum/cntl - (ytot-ysum)*(ytot-ysum)/cntr;
		}
/*  Printf("fence %f dev %f #l %g #r %g\n", fence, ldev, cntl, cntr); */ 
		if (ldev < bdev) {
		    bdev = ldev;
		    bfence = fence;
		}
	    }
	    val = bdev + sdev;
	    Printf(" val %f fence %f\n", val, bfence);
	    if (val >= devtarget || val >= *bval) return;
	    *bval = val;
	    var[inode] = iv + 1;
	    labl = cutleft[inode];
	    labr = cutright[inode];
	    strcpy(labl, ":");
	    strcpy(labr, ":");
	    for (l = 0; l < nll; l++)
		if (cprob[l] < bfence) scat(labl, 'a' + ind[l]);
		else scat(labr, 'a' + ind[l]);
/* Printf("%s %s\n", labl, labr); */
	    for (l = 0; l < nl; l++) indl[l] = False;
	    for (l = 0; l < nll; l++)
		if (cprob[l] < bfence) indl[ind[l]] = True;
/*	 for(l = 0; l < nl; l++) Printf(" %d", indl[l]); Printf("\n");*/
	    for (j = 0; j < nobs; j++)
		if (twhere[j] < 0) ttw[j] = twhere[j];
		else ttw[j] = !indl[twhere[j]];

	} else {

	    Printf(" cnts "); for(l = 0; l < nll; l++) Printf(" %g", cnt[l]);
	    Printf("\n"); 
	    indl[0] = True;
	    for (l = 1; l < nll; l++) indl[l] = False;
	    bdev = devtarget;
	    mi = 1;
	    iis = -1;
	    for(l = 1; l < nll; l++) mi *= 2;
	    for(i = 0; i < mi-1; i++) { /* go though all splits */
		ii = i;
		for(l = 1; l < nll; l++) {
		    indl[l] = (ii%2);
		    ii /= 2;
		}
		cntl = cntr = 0;
		for (l = 0; l < nll; l++)
		    if (indl[l])
			cntl += cnt[l];
		    else
			cntr += cnt[l];
		if (cntl < mincut || cntr < mincut) continue;
		if (Gini) {
		    ldev = n[inode];
		    for (k = 0; k < nc; k++) {
			cntl1 = cntr1 = 0;
			for (l = 0; l < nll; l++)
			    if (indl[l]) cntl1 += tab[k + nc * l];
			    else cntr1 += tab[k + nc * l];
			ldev -= cntl1*cntl1/cntl + cntr1*cntr1/cntr;
		    }
		} else {
		    ldev = XLOGX(cntl) + XLOGX(cntr);
		    for (k = 0; k < nc; k++) {
			cntl1 = cntr1 = 0;
			for (l = 0; l < nll; l++)
			    if (indl[l]) cntl1 += tab[k + nc * l];
			    else cntr1 += tab[k + nc * l];
			ldev -= XLOGX(cntl1) + XLOGX(cntr1);
		    }
		}
		ldev *= 2;
/*	     for(l = 0; l < nll; l++) Printf("%d ", indl[l]);
	     Printf(": %d %f\n", i, ldev); */
		if (ldev < bdev) {
		    bdev = ldev;
		    iis = i;
		}
	    }
	    val = bdev + sdev;
	    Printf(" val %f at bin val %d\n", val, iis);
	    if (val >= *bval || val >= devtarget) return;
	    *bval = val;
	    for(l = 1; l < nll; l++) {
		indl[l] = (iis%2);
		iis /= 2;
	    }
	    var[inode] = iv + 1;
	    labl = cutleft[inode];
	    labr = cutright[inode];
	    strcpy(labl, ":");
	    strcpy(labr, ":");
	    for (l = 0; l < nll; l++)
		if (indl[l]) scat(labl, 'a' + ind[l]);
		else scat(labr, 'a' + ind[l]);
	    for (l = 0; l < nl; l++) indr[l] = False;
	    for (l = 0; l < nll; l++) indr[ind[l]] = indl[l];
	    for (j = 0; j < nobs; j++)
		if (twhere[j] < 0) ttw[j] = twhere[j];
		else ttw[j] = !indr[twhere[j]];
	}
    }
}

static void shift_up_node(int i, int N)
{
    int j, k;
/*    Printf("shifting %d to %d\n", i, i+N); */
    var[i+N] = var[i];
    cutleft[i+N] = cutleft[i];
    cutright[i+N] = cutright[i];
/*    Printf("(%d) %d to %d %s %s %p\n", node[i], i, i+N, cutleft[i+N], 
      cutright[i+N], *(cutleft+i+N));*/
    n[i+N] = n[i];
    dev[i+N] = dev[i];
    yval[i+N] = yval[i];
    node[i+N] = node[i];
    for (k = 0; k < nc; k++) yprob[(i+N)*nc+k] = yprob[i*nc+k];
    for (j = 0; j < nobs; j++) if (where[j] == i) where[j] +=N;
}

static void shift_down_node(int i, int N)
{
    int j, k;
/*    Printf("shifting %d to %d %p\n", i+N, i); */
    var[i] = var[i+N];
    cutleft[i] = cutleft[i+N];
    cutright[i] = cutright[i+N];
    n[i] = n[i+N];
    dev[i] = dev[i+N];
    yval[i] = yval[i+N];
    node[i] = node[i+N];
/*    Printf("(%d) %d to %d %s %s %p\n", node[i], i+N, i, cutleft[i], 
      cutright[i], *(cutleft+i)); */
    for (k = 0; k < nc; k++) yprob[i*nc+k] = yprob[(i+N)*nc+k];
    for (j = 0; j < nobs; j++) if (where[j] == i+N) where[j] -=N; 
}

static void divide_node(int inode)
{
    int     i, iv, j, k, shift, shifted = False;
    double  bval, tmp;

    if (inode >= nmax) error(_("tree is too big"));

    fillin_node(inode);
    if ( n[inode] < minsize ) return;

    if (Gini) {
	bval = 0.0;
	for (k = 0; k < nc; k++) {
	    tmp = yprob[inode*nc + k];
	    bval += tmp *tmp;
	}
	bval = n[inode] * (1 - bval);
	bval *= 2.0;
	Printf("gini = %g\n", bval);
	devtarget = bval;
    } else {
	bval = dev[inode];
	devtarget = dev[inode] - mindev*dev[0];
    }
    if(devtarget <= (1e-6)*dev[0]) return;
    Printf("\n--evaluating node %d(%d) size %g\n", inode, 
	   (int)node[inode], n[inode]);

    for (iv = 0; iv < nvar; iv++)
	if (levels[iv])
	    split_disc(inode, iv, &bval);
	else
	    split_cont(inode, iv, &bval);

    Printf("..best value is %g\n", bval);
   
    if (bval < devtarget) {
        Printf("..splitting\n");
	if ( node[inode] >=  1073741824 ) {
	    error(_("maximum depth reached\n"));
	    return;
	}
   
	if (inode < nnode-1) {
	    shifted = nnode;
	    for (i= nnode-1; i > inode; i--) shift_up_node(i, nmax-nnode);
	    nnode = inode + 1;
/*Printf("..shifted up\n");*/
	} else shifted = False;
	/* write left as nnode */
	for (j = 0; j < nobs; j++) {
	    if (ttw[j] == 0) where[j] = nnode;
	    if (ttw[j] == NALEVEL) where[j] += NALEVEL;
	}
	node[nnode++] = 2 * node[inode];
	divide_node(nnode-1);
	Printf("..done left at %d\n", inode);
	/* write right as nnode */
	for (j = 0; j < nobs; j++)
	    if (where[j] == inode) where[j] = nnode;
	node[nnode++] = 2 * node[inode] + 1;
	divide_node(nnode-1);
	Printf("..done right at %d\n", inode);
	if (shifted) {
	    shift = nnode - inode -1;
	    for (i = inode+1; i < shifted; i++) 
		shift_down_node(i+shift, nmax-shifted-shift);
	    offset += shift;
	    nnode = shifted + shift;
/*Printf("..shifted down\n");*/
	}
    }
}

void 
BDRgrow1(double *pX, double *pY, double *pw, Sint *plevels, Sint *junk1, 
	 Sint *pnobs, Sint *pncol, Sint *pnode, Sint *pvar, char **pcutleft, 
	 char **pcutright, double *pn, double *pdev, double *pyval, 
	 double *pyprob, Sint *pminsize, Sint *pmincut, double *pmindev, 
	 Sint *pnnode, Sint *pwhere, Sint *pnmax, Sint *stype, Sint *pordered)
{
    int i, nl;

    X = pX; y = pY; w = pw; dev = pdev; yval = pyval; yprob = pyprob;
    nobs = *pnobs; nvar = *pncol;
    levels = plevels; node = pnode; var = pvar; n = pn; mindev = *pmindev;
    minsize = *pminsize; mincut = *pmincut; nmax = *pnmax; nnode = *pnnode;
    where = pwhere; cutleft = pcutleft; cutright = pcutright; 
    ordered= pordered; Gini = *stype;
    nc = levels[nvar];
    Printf("nnode: %d\n", nnode);
    Printf("nvar: %d\n", nvar);
    for(i = 0; i <= nvar; i++) Printf("%d ", (int)levels[i]);
    Printf("\n");
    /* allocate scratch storage */
    nl = 0;
    for(i = 0; i <= nvar; i++)
	if (levels[i] > nl) nl = levels[i];
    maxnl = max(nl, 10);
    twhere = (int *) S_alloc(nobs, sizeof(int));
    ttw = (int *) S_alloc(nobs, sizeof(int));
    tvar = (double *) S_alloc(nobs, sizeof(double));
    ind = (int *) S_alloc(nl, sizeof(int));
    w1 = (double *) S_alloc(nobs, sizeof(double));
    cnt = (double *) S_alloc(nl, sizeof(double));
    cprob = (double*) S_alloc(nl, sizeof(double));
    scprob = (double*) S_alloc(nl, sizeof(double));
    indl = (int*) S_alloc(nl, sizeof(int));
    if (nc > 0) {
	yp = (double *) S_alloc(nc, sizeof(double));
	tab = (double*) S_alloc(nl*(1+nc), sizeof(double));
	indr = (int*) S_alloc(nl, sizeof(int));
	ty = (int *) S_alloc(nobs, sizeof(int));
    } else {
	tyc = (double *) S_alloc(nobs, sizeof(double));
	ys = (double *) S_alloc(nl, sizeof(double));
    }
    exists = nnode;
    offset = 0;
    if (exists <= 1) {
	for(i = 0; i < nobs; i++) where[i] = 0;
	nnode = 1;
	node[0] = 1;
	divide_node(0);
    } else {
	/* Adjust from S indexing */
	for(i = 0; i < nobs; i++) where[i]--;
	for(i = 0; i < exists; i++)
	    if (!var[i+offset]) {
/* Printf("trying node %d at offset %d, nnode %d\n", i, offset, nnode);*/
		divide_node(i + offset);
	    }
    }
    /* Adjust to S indexing */

    for(i = 0; i < nobs; i++) {
	if(where[i] < 0) where[i] -= NALEVEL;
	where[i]++;
    }
    *pnnode = nnode;
    Printf("Finished!\n");
}

