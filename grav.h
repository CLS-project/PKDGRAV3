#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED

#include "pkd.h"
#include "moments.h"

static inline double softmassweight(double m1,double h12,double m2,double h22){
    return((m1+m2)*(h12*h22)/(h22*m1+h12*m2));
}

void PPInteractSIMD( int nPart, ILP *ilp, const FLOAT *r, const FLOAT *a,
		     FLOAT fMass, FLOAT fSoft,
		     momFloat *ax, momFloat *ay, momFloat *az,
		     momFloat *fPot, momFloat *rhosum, momFloat *maisum );

int pkdGravInteract(PKD pkd,KDN *pBucket,LOCR *l,ILP *ilp,int nPart,ILC *ilc,int nCell,ILPB *ilpb,int nPartBucket,double dirLsum,double normLsum,int bEwaldKicking,double *pdFlop);
#ifdef HERMITE
double pkdRho(PKD pkd, double rhopmaxlocal,double summ, double sumr, double *dir2, 
	      double *dir, double x, double y, double z, double vx, double vy, double vz, 
	      double rv, double v2, double a3, int iOrder,int jOrder);
#endif
#endif



