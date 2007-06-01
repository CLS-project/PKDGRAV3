#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED

#include "pkd.h"
#include "moments.h"

inline double softmassweight(double m1,double h12,double m2,double h22);

void PPInteractSIMD( int nPart, ILP *ilp, const FLOAT *r, const FLOAT *a,
		     FLOAT fMass, FLOAT fSoft,
		     float *ax, float *ay, float *az, float *fPot,
		     float *rhosum, float *maisum );

int pkdGravInteract(PKD pkd,KDN *pBucket,LOCR *l,ILP *ilp,int nPart,ILC *ilc,int nCell,ILPB *ilpb,int nPartBucket,double *pdFlop);

#endif



