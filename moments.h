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

/*
** V4SUM takes a momPacked and adds horizontally all the elements of the packed vector.
*/
#define V4SUM(A) ((A).f[0] + (A).f[1] + (A).f[2] + (A).f[3])

typedef struct momVectorReduced {
	momPacked m;
	momPacked xx,yy,xy,xz,yz;
	momPacked xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momPacked xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	} VMOMR;

typedef struct momGenLocrAddMomrArray {
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
} GLAM;
#endif

/*
** components required for evaluating a multipole interaction.
*/

#ifdef USE_SIMD
typedef struct ilCell {   
    double x[4], y[4], z[4];
    momPacked xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
    momPacked xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    momPacked xx,xy,xz,yy,yz;
    momPacked m;
} ILC;
#else
typedef struct ilCell {
    double x,y,z; /*,vx,vy,vz;*/
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
 */
typedef struct locReduced {
	momFloat m;
	momFloat x,y,z;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
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
void momShiftLocr(LOCR *,momFloat,momFloat,momFloat);
void momShiftLocr2(LOCR *,momFloat,momFloat,momFloat);
void momReduceMomc(MOMC *,MOMR *);
void momReduceLocc(LOCC *,LOCR *);
/*
** All the variants of EvalMomr...
*/
void momEvalMomr(MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z,
		 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az,momFloat *magai);
void momGenEvalMomr(MOMR *m,momFloat dir,momFloat g0,momFloat t1,momFloat t2,
		    momFloat t3r,momFloat t4r,momFloat x,momFloat y,momFloat z,
		    momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az,
		    momFloat *magai);
#ifdef USE_SIMD
double momGenEvalVMomr(int n,GLAM *p,momFloat ax,momFloat ay,momFloat az,
		     momFloat *fPot,momFloat *aix,momFloat *aiy,momFloat *aiz,
		     momFloat *rhosum,momFloat *maisum);
void momPadSIMDMomr( int *nCell, ILC *ilc );
void momEvalSIMDMomr( int nCell, ILC *ilc, const double *r,
		      float *ax, float *ay, float *az,
		      float *fPot );
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


void momEvalLocc(LOCC *,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *); 
void momEvalLocr(LOCR *,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *); 
void momSubLocc(LOCC *,LOCC *);
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);
void momPrintLocc(LOCC *);
void momPrintLocr(LOCR *);


#endif
