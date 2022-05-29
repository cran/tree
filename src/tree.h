/*
 * tree/src/tree.h  Copyright (C) 2005-2022 B. D. Ripley
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
void 
BDRgrow1(double *pX, double *pY, double *pw, int *plevels, int *junk1, 
	 int *pnobs, int *pncol, int *pnode, int *pvar, char **pcutleft, 
	 char **pcutright, double *pn, double *pdev, double *pyval, 
	 double *pyprob, int *pminsize, int *pmincut, double *pmindev, 
	 int *pnnode, int *pwhere, int *pnmax, int *stype, int *pordered);


void VR_dev1(int *nnode, int *nodes, int *parent, 
	     double *dev, double *sdev,
	     int *y, int *ny, int *yf, int *where, double *wt,
	     int *nc, double *loss);

void VR_dev2(int *nnode, int *nodes, int *parent, 
	     double *dev, double *sdev,
	     int *y, int *ny, double *yprob, int* where, double *wt);

void VR_dev3(int *nnode, int *nodes, int *parent, 
	     double *dev, double *sdev,
	     double *y, int *ny, double *yf, int* where, double *wt);

void    
VR_prune2(int *nnode, int *nodes, int *leaf, double *dev, double *sdev,
	  double *ndev, double *nsdev, int *keep, int *ord, double *g,
	  int *size, double *cdev, double *alph, int *inodes, int *tsize,
	  double *tdev, double *ntdev);

void    
VR_pred1(double *x, int *vars, char **lsplit, char **rsplit,
	 int *nlevels, int *nodes, int *fn, int *nnode,
	 int *nr, int *nc, int *where);

void    
VR_pred2(double *px, int *pvars, char **plsplit, char **prsplit,
	 int *pnlevels, int *pnodes, int *fn, int *pnnode,
	 int *nr, double *pwhere);







