#include <R.h>
void 
BDRgrow1(double *pX, double *pY, double *pw, Sint *plevels, Sint *junk1, 
	 Sint *pnobs, Sint *pncol, Sint *pnode, Sint *pvar, char **pcutleft, 
	 char **pcutright, double *pn, double *pdev, double *pyval, 
	 double *pyprob, Sint *pminsize, Sint *pmincut, double *pmindev, 
	 Sint *pnnode, Sint *pwhere, Sint *pnmax, Sint *stype, Sint *pordered);


void VR_dev1(Sint *nnode, Sint *nodes, Sint *parent, 
	     double *dev, double *sdev,
	     Sint *y, Sint *ny, Sint *yf, Sint *where, double *wt,
	     Sint *nc, double *loss);

void VR_dev2(Sint *nnode, Sint *nodes, Sint *parent, 
	     double *dev, double *sdev,
	     Sint *y, Sint *ny, double *yprob, Sint* where, double *wt);

void VR_dev3(Sint *nnode, Sint *nodes, Sint *parent, 
	     double *dev, double *sdev,
	     double *y, Sint *ny, double *yf, Sint* where, double *wt);

void    
VR_prune2(Sint *nnode, Sint *nodes, Sint *leaf, double *dev, double *sdev,
	  double *ndev, double *nsdev, Sint *keep, Sint *ord, double *g,
	  Sint *size, double *cdev, double *alph, Sint *inodes, Sint *tsize,
	  double *tdev, double *ntdev);

void    
VR_pred1(double *x, Sint *vars, char **lsplit, char **rsplit,
	 Sint *nlevels, Sint *nodes, Sint *fn, Sint *nnode,
	 Sint *nr, Sint *nc, Sint *where);

void    
VR_pred2(double *px, Sint *pvars, char **plsplit, char **prsplit,
	 Sint *pnlevels, Sint *pnodes, Sint *fn, Sint *pnnode,
	 Sint *nr, double *pwhere);







