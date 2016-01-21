/*
 *  treefix/treefix.c by B. D. Ripley  Copyright (C) 1994-2009
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

#define True 1
#define False 0
#define EPS 1e-4
#include <math.h>
#include <R.h>
#include <stdio.h>
#include <string.h> /* for strchr */

# include <R_ext/Utils.h>

#ifndef Salloc
#  define Salloc(n, t) (t *)S_alloc(n, sizeof(t))
#endif

/* workaround for log(0) */
#define safe_log(result, x) result = log(x + 1e-200)

static int mysum(Sint *vec, int n)
{
    int     i, s = 0;
    for (i = 0; i < n; i++) s += vec[i];
    return (s);
}

void VR_dev1(Sint *nnode, Sint *nodes, Sint *parent, 
	     double *dev, double *sdev,
	     Sint *y, Sint *ny, Sint *yf, Sint *where, double *wt,
	     Sint *nc, double *loss)
{
    int i, j, wh, nn, nr = *nnode, nclass = *nc;
	
    parent[0] = -1;
    for (i = 1; i < nr; i++)
    {
	nn = nodes[i]/2;
	for (j = 0; j < i; j++)
	    if(nodes[j] == nn) {parent[i] = j; break;}
    }
    for (i = 0; i < nr; i++) dev[i] = sdev[i] = 0.0;
    for (j = 0; j < *ny; j++)
    {
	wh = where[j] - 1; /* C indexed */
	sdev[wh] += wt[j] * loss[y[j] - 1 + nclass*(yf[wh]-1)];
	dev[wh] += wt[j] * loss[y[j] - 1 + nclass*(yf[wh]-1)];
	while((wh = parent[wh]) >= 0) 
	    dev[wh] += wt[j] * loss[y[j] - 1 + nclass*(yf[wh]-1)];
    }
}

void VR_dev2(Sint *nnode, Sint *nodes, Sint *parent, 
	     double *dev, double *sdev,
	     Sint *y, Sint *ny, double *yprob, Sint* where, double *wt)
{
    int i, j, wh, nn, nr = *nnode;
    double tmp;
	
    parent[0] = -1;
    for (i = 1; i < nr; i++)
    {
	nn = nodes[i]/2;
	for (j = 0; j < i; j++)
	    if(nodes[j] == nn) {parent[i] = j; break;}
    }
    for (i = 0; i < nr; i++) dev[i] = sdev[i] = 0.0;
    for (j = 0; j < *ny; j++)
    {
	wh = where[j] - 1; /* C indexed */
	safe_log(tmp, yprob[wh + nr*(y[j]-1)]);
	sdev[wh] += wt[j] * tmp;
	dev[wh] += wt[j] * tmp;
	while((wh = parent[wh]) >= 0) {	  
	    safe_log(tmp, yprob[wh + nr*(y[j]-1)]);
	    dev[wh] += wt[j] * tmp;
	}
		
    }
}

# define sqr(x) (x) * (x)
void VR_dev3(Sint *nnode, Sint *nodes, Sint *parent, 
	     double *dev, double *sdev,
	     double *y, Sint *ny, double *yf, Sint* where, double *wt)
{
    int i, j, wh, nn, nr = *nnode;
	
    parent[0] = -1;
    for (i = 1; i < nr; i++)
    {
	nn = nodes[i]/2;
	for (j = 0; j < i; j++)
	    if(nodes[j] == nn) {parent[i] = j; break;}
    }
    for (i = 0; i < nr; i++) dev[i] = sdev[i] = 0.0;
    for (j = 0; j < *ny; j++)
    {
	wh = where[j] - 1; /* C indexed */
	sdev[wh] += wt[j] * sqr(y[j] - yf[wh]);
	dev[wh] += wt[j] * sqr(y[j] - yf[wh]);
	while((wh = parent[wh]) >= 0) {
	    dev[wh] += wt[j] * sqr(y[j] - yf[wh]);
	}
    }
}

void    
VR_prune2(Sint *nnode, Sint *nodes, Sint *leaf, double *dev, double *sdev,
	  double *ndev, double *nsdev, Sint *keep, Sint *ord, double *g,
	  Sint *size, double *cdev, double *alph, Sint *inodes, Sint *tsize,
	  double *tdev, double *ntdev)
{
    int     i, in, ir, j, k, cur, nr = *nnode, sz, First, na = 0;
    double  alpha = 0.0, rt, sum;

    for (i = 0; i < nr; i++) keep[i] = True;
    /* start with full tree */
    inodes[na] = 0;
    alph[na] = -1.0e+200;
    tsize[na] = mysum(leaf, nr);
    sum = 0;
    for (j = 0; j < nr; j++)
	if (leaf[j]) sum += ndev[j];
	else if (keep[j]) sum += nsdev[j];	/* presumably all kept */
    ntdev[na] = sum;
    sum = 0;
    for (j = 0; j < nr; j++)
	if (leaf[j]) sum += dev[j];
	else if (keep[j]) sum += sdev[j];
    tdev[na++] = sum;

    do {
	First = True;
	for (i = 0; i < nr; i++) {
	    cdev[i] = dev[i];
	    size[i] = 1;
	}
	for (k = nr - 1; k >= 0; k--) {
	    cur = ord[k] - 1;
	    if (keep[cur] && !leaf[cur]) {
		/* sum over descendants */
		rt = sdev[cur];
		sz = 0;
		ir = nodes[cur];
		for (j = 0; j < nr; j++)
		    if (nodes[j] / 2 == ir) {
			rt += cdev[j];
			sz += size[j];
		    }
		size[cur] = sz;
		g[cur] = (dev[cur] - rt) / (sz - 1);
		cdev[cur] = rt;
		if (First) {
		    First = False;
		    alpha = g[cur];
		}
		else if (g[cur] < alpha) alpha = g[cur];
	    }
	}
	for (k = 0; k < nr; k++) {
	    cur = ord[k] - 1;
	    if (keep[cur] && !leaf[cur] &&
		fabs(g[cur] - alpha) < EPS * (1 + fabs(alpha))) {
		/* prune at cur */
		leaf[cur] = True;
		alph[na] = alpha;
		ir = nodes[cur];
		inodes[na] = ir;
		for (j = 0; j < nr; j++)
		    if (keep[j]) {
			in = nodes[j];
			if (in >= 2 * ir) {
			    while (in >= 2 * ir) in = in / 2;
			    if (in == ir) {
				keep[j] = False;
				leaf[j] = False;
			    }
			}
		    }
		tsize[na] = mysum(leaf, nr);
		sum = 0;
		for (j = 0; j < nr; j++)
		    if (leaf[j]) sum += ndev[j];
		    else if (keep[j]) sum += nsdev[j];
		ntdev[na] = sum;
		sum = 0;
		for (j = 0; j < nr; j++)
		    if (leaf[j])
			sum += dev[j];
		    else if (keep[j])
			sum += sdev[j];
		tdev[na++] = sum;
	    }
	}
    } while (mysum(keep, nr) > 1);
    *nnode = na;
}

/* Take each case and drop it down the tree as needed */

void    
VR_pred1(double *x, Sint *vars, char **lsplit, char **rsplit,
	 Sint *nlevels, Sint *nodes, Sint *fn, Sint *nnode,
	 Sint *nr, Sint *nc, Sint *where)
{
    int     nobs = *nr, cur, i, k, goleft, ival, cnode, var, *left, *right;
    double  val, sp;

    left = Salloc(*nnode, int);
    right = Salloc(*nnode, int);
    for (i = 0; i < *nnode; i++)
	if (vars[i] > 0) {
	    cnode = nodes[i];
            left[i] = *nnode; /* protect against corrupted tree */
            right[i] = *nnode;
	    for (k = i + 1; k < *nnode; k++) {
		if (nodes[k] == 2 * cnode) left[i] = k;
		if (nodes[k] == 2 * cnode + 1) right[i] = k;
	    }
	}

    for (i = 0; i < nobs; i++) {
	cur = 0;			/* current node, C indexed */
	while (True) {
	    if (cur >= *nnode)
		PROBLEM "corrupt tree" RECOVER(NULL_ENTRY);
	    if (vars[cur] == 0) {	/* at a leaf */
		where[i] = cur + 1;
		break;
	    }
	    var = vars[cur] - 1;	/* C indexing */
	    val = x[i + nobs * var];
#ifdef USING_R
	    if (ISNA(val)) {
#else
	    if (is_na(&val, DOUBLE)) {
#endif
		where[i] = cur + 1;
		break;
	    }
	    if (nlevels[var] == 0) {
		sp = R_atof(lsplit[cur] + 1);
		goleft = (val < sp);
	    }
	    else {
		ival = 'a' + (int) val - 1;
		if (strchr(lsplit[cur], ival) != NULL)
		    goleft = True;
		else if (strchr(rsplit[cur], ival) != NULL)
		    goleft = False;
		else {		/* unforeseen level */
		    where[i] = cur + 1;
		    break;
		}
	    }
	    if (goleft) cur = left[cur];
	    else cur = right[cur];
	}
    }
}

static int nobs, nnode, *left, *right;
static double *lprob, *where;
static double *x;
static char **lsplit, **rsplit;
static Sint *vars, *nlevels, *nodes;

static void 
downtree(int i, int cur, double prob)
    {
	int     k, ival, cnode, var;
	double   goleft;
	double  val, sp;
	if (cur >= nnode) PROBLEM "corrupt tree" RECOVER(NULL_ENTRY);
	where[cur + nnode * i] += prob;
	if (vars[cur] == 0)		/* at a leaf */
	    return;
	var = vars[cur] - 1;		/* C indexing */
	val = x[i + nobs * var];
#ifdef USING_R
	if (ISNA(val)) {
#else
	if (is_na(&val, DOUBLE)) {
#endif
	    goleft = lprob[cur];
	}
	else if (nlevels[var] == 0) {
	    sp = R_atof(lsplit[cur] + 1);
	    goleft = (val < sp);
	}
	else {
	    ival = 'a' + (int) val - 1;
	    if (strchr(lsplit[cur], ival) != NULL)
		goleft = True;
	    else if (strchr(rsplit[cur], ival) != NULL)
		goleft = False;
	    else
		goleft = lprob[cur];
	}
	cnode = nodes[cur];
	if (goleft > 0) {
	    for (k = cur + 1; k < nnode; k++)
		if (nodes[k] == 2 * cnode)
		    break;
	    downtree(i, k, prob * goleft);
	}
	if (goleft < 1) {
	    for (k = cur + 1; k < nnode; k++)
		if (nodes[k] == 2 * cnode + 1)
		    break;
	    downtree(i, k, prob * (1 - goleft));
	}
    }

void    
VR_pred2(double *px, Sint *pvars, char **plsplit, char **prsplit,
	 Sint *pnlevels, Sint *pnodes, Sint *fn, Sint *pnnode,
	 Sint *nr, double *pwhere)
    {
	int     cnode, i, k;

	nobs = *nr;
	nnode = *pnnode;
	x = px;
	vars = pvars;
	lsplit = plsplit;
	rsplit = prsplit;
	nlevels = pnlevels;
	nodes = pnodes;
	where = pwhere;
	lprob = Salloc((long)nnode, double);
	left = Salloc((long)nnode, int);
	right = Salloc((long)nnode, int);
	for (i = 0; i < nnode; i++)
	    if (vars[i] > 0) {
		cnode = nodes[i];
		for (k = i + 1; k < nnode; k++) {
		    if (nodes[k] == 2 * cnode) left[i] = k;
		    if (nodes[k] == 2 * cnode + 1) right[i] = k;
		}
		lprob[i] = (double)fn[left[i]] / (fn[left[i]] + fn[right[i]]);
	    }

	for (i = 0; i < nobs; i++) {
	    for (k = 0; k < nnode; k++) where[k + nnode * i] = 0.0;
	    downtree(i, 0, 1.0);
	}
    }
