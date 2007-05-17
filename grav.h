#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED

#include "pkd.h"
#include "moments.h"

void PPInteractSIMD( int nPart, ILP *ilp, const FLOAT *r, const FLOAT *a,
		     FLOAT fMass, FLOAT fSoft,
		     float *ax, float *ay, float *az, float *fPot );

int pkdGravInteract(PKD pkd,KDN *pBucket,LOCR *l,ILP *ilp,int nPart,ILC *ilc,int nCell,ILPB *ilpb,int nPartBucket,double *pdFlop);

#endif



