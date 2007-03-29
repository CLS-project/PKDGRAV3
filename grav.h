#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED

#include "pkd.h"
#include "moments.h"

int pkdGravInteract(PKD pkd,KDN *pBucket,LOCR *l,ILP *ilp,int nPart,ILC *ilc,int nCell,ILPB *ilpb,int nPartBucket,double *pdFlop);

#endif



