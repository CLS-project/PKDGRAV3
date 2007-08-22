#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#ifdef USE_SIMD
#include "simd.h"
#endif

#ifdef HAVE_CONFIG_H
typedef MOMFLOAT momFloat;
#else
#ifdef MOMQUAD
typedef long double momFloat;
#define sqrt(x)	sqrtl(x)
#else
typedef double momFloat;
#endif
#endif

/*
 ** moment tensor components for reduced multipoles.
 */
/* IMPORTANT: The order of the fields MUST match with VMOMR */
typedef struct momReduced {
	momFloat m;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	} MOMR;

#ifdef USE_SIMD
/*
 ** vector moment tensor components for reduced multipoles.
 */

typedef union {
    v4sf p;
    float f[4];
} momPacked;

/* IMPORTANT: The order of the fields MUST match with MOMR */
typedef struct momVectorReduced {
    momPacked m;
    momPacked xx,yy,xy,xz,yz;
    momPacked xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    momPacked xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
} VMOMR;

/* IMPORTANT: The order of the fields MUST match with GLAM  */
/* IMPORTANT: The size of this structure must aligned by 16 */
typedef struct momGenLocrAddMomrSIMDArray {
    VMOMR q;
    momPacked dir;
    momPacked g0;
    momPacked t1;
    momPacked t2;
    momPacked t3r;
    momPacked t4r;
    momPacked x;
    momPacked y;
    momPacked z;
    momPacked zero;
} VGLAM;
#endif

/* IMPORTANT: The order of the fields MUST match with VGLAM */
/* IMPORTANT: The size of this structure must aligned by 16 */
typedef struct momGenLocrAddMomrArray {
    MOMR q;
    momFloat dir;
    momFloat g0;
    momFloat t1;
    momFloat t2;
    momFloat t3r;
    momFloat t4r;
    momFloat x;
    momFloat y;
    momFloat z;
    momFloat zero; /* 32 floats here */
} GLAM;

/*
** components required for evaluating a multipole interaction.
*/

#ifdef USE_SIMD_MOMR
typedef struct ilCell {   
    double x[4], y[4], z[4];
    momPacked xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
    momPacked xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    momPacked xx,xy,xz,yy,yz;
    momPacked m;
} ILC;
#else
typedef struct ilCell {
    double x,y,z; 
#ifdef HERMITE
    double vx,vy,vz;
#endif
    MOMR mom;
    } ILC;
#endif

/*
 ** moment tensor components for complete multipoles.
 */
typedef struct momComplete {
	momFloat m;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	momFloat zz;
	momFloat xzz,yzz,zzz;
	momFloat xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} MOMC;


/*
 ** moment tensor components for reduced local expansion.
 ** note that we have the 5th-order terms here now!
 */
typedef struct locReduced {
	momFloat m;
	momFloat x,y,z;
	momFloat xx,xy,yy,xz,yz;
	momFloat xxx,xxy,xyy,yyy,xxz,xyz,yyz;
	momFloat xxxx,xxxy,xxyy,xyyy,yyyy,xxxz,xxyz,xyyz,yyyz;
        momFloat xxxxx,xxxxy,xxxyy,xxyyy,xyyyy,yyyyy,xxxxz,xxxyz,xxyyz,xyyyz,yyyyz;        
	} LOCR;

/*
 ** moment tensor components for complete local expansion.
 */
typedef struct locComplete {
	momFloat m;
	momFloat x,y,z;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	momFloat zz;
	momFloat xzz,yzz,zzz;
	momFloat xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} LOCC;


void momAddMomc(MOMC *,MOMC *);
void momAddMomr(MOMR *,MOMR *);
void momMulAddMomc(MOMC *,momFloat,MOMC *);
void momMulAddMomr(MOMR *,momFloat,MOMR *);
void momSubMomc(MOMC *,MOMC *);
void momSubMomr(MOMR *,MOMR *);
void momMakeMomc(MOMC *,momFloat,momFloat,momFloat,momFloat);
momFloat momMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momOldMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momShiftMomc(MOMC *,momFloat,momFloat,momFloat);
void momShiftMomr(MOMR *,momFloat,momFloat,momFloat);
double momShiftLocr(LOCR *,momFloat,momFloat,momFloat);
void momReduceMomc(MOMC *,MOMR *);
/*
** All the variants of EvalMomr...
*/
void momEvalMomr(MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z,
		 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az,
		 momFloat *magai);
void momGenEvalMomr(MOMR *m,momFloat g0,momFloat g1,momFloat g2,momFloat g3,momFloat g4,momFloat g5,
		    momFloat x,momFloat y,momFloat z,
		    momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az,momFloat *magai);
#ifdef USE_SIMD
double momGenEvalVMomr(int n,GLAM *p,momFloat ax,momFloat ay,momFloat az,
		     momFloat *fPot,momFloat *aix,momFloat *aiy,momFloat *aiz,
		     momFloat *rhosum,momFloat *maisum);
void momPadSIMDMomr( int *nCell, ILC *ilc );
void momEvalSIMDMomr( int nCell, const ILC *ilc, const double *r, const FLOAT *a,
		      float *ax, float *ay, float *az, float *fPot,
		      float *rhosum, float *maisum );
#endif
void momClearLocc(LOCC *);
void momClearLocr(LOCR *);
void momClearMomr(MOMR *);
void momMomr2Momc(MOMR *,MOMC *);
/*
** All the variants of LocrAddMomr...
*/
void momSymLocrAddMomr(LOCR *l1,LOCR *l2,MOMR *q1,MOMR *q2,momFloat dir,momFloat x,momFloat y,momFloat z);
double momLocrAddMomr(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat);
void momGenLocrAddMomr(LOCR *l,MOMR *q,momFloat dir,
		       momFloat g0,momFloat t1,momFloat t2,momFloat t3r,momFloat t4r,
		       momFloat x,momFloat y,momFloat z);
#ifdef USE_SIMD
double momGenLocrAddVMomr(LOCR *l,int n,GLAM *p,momFloat ax,momFloat ay,momFloat az,momFloat *rhosum,momFloat *maisum);
double momGenLocrAddSIMDMomr(LOCR *l,int n,GLAM *p,momFloat ax,momFloat ay,momFloat az,momFloat *rhosum,momFloat *maisum);
#endif
void momEwaldLocrAddMomr(LOCR *l,MOMR *m,momFloat r2,int bInHole,momFloat x,momFloat y,momFloat z);
void momNooptLocrAddMomr(LOCR *l,MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z);
void momLoccAddMomrAccurate(LOCC *l,MOMC *m,momFloat g0,momFloat x,momFloat y,momFloat z);
void momLocrAddMomrAccurate(LOCR *l,MOMR *m,momFloat g0,momFloat x,momFloat y,momFloat z);
/*
** These are the prefered versions that should be used in pkdgrav2.
*/
double momLocrAddMomr5(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat,double *,double *,double *);
double momLocrAddMomr5Noopt(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat);
void momEvalLocr(LOCR *,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *); 
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);
void momPrintLocc(LOCC *);
void momPrintLocr(LOCR *);


#endif
