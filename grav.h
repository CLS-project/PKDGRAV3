#ifndef GRAV_HINCLUDED
#define GRAV_HINCLUDED

#include "pkd.h"
#include "moments.h"

static inline double softmassweight(double m1,double h12,double m2,double h22) {
    double tmp = h12*h22;
    if (m1 == 0.0) return(h22);
    if (m2 == 0.0) return(h12);
    if (tmp > 0.0) return((m1+m2)*tmp/(h22*m1+h12*m2));
    else return(0.0);
    }

void PPInteractSIMD( int nPart, ILP ilp, const FLOAT *r, const FLOAT *a,
		     FLOAT fMass, FLOAT fSoft,
		     momFloat *ax, momFloat *ay, momFloat *az,
		     momFloat *fPot, momFloat *rhosum, momFloat *maisum );

#ifdef LOCAL_EXPANSION
int pkdGravInteract(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,KDN *pBucket,LOCR *pLoc,ILP ilp,ILC ilc,
		    double dirLsum,double normLsum,int bEwald,double *pdFlop,double *pdEwFlop,double dRhoFac);
#else
int pkdGravInteract(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,KDN *pBucket,LOCR *pLoc,ILP *ilp,int nPart,ILC *ilc,int nCell,
		    double dirLsum,double normLsum,int bEwald,double *pdFlop,double *pdEwFlop,double dRhoFac);
#endif

double pkdRho1(double rhopmaxlocal, double summ, double dir, double x, double y, double z, double vx, double vy, double vz);

#ifdef HERMITE
double pkdRho3(PKD pkd, double rhopmaxlocal,double summ, double sumr, double *dir2,
	      double *dir, double x, double y, double z, double vx, double vy, double vz,
	      double rv, double v2, double a3, int iOrder,int jOrder);
#endif
#endif



