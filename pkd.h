#ifndef PKD_HINCLUDED
#define PKD_HINCLUDED

#include <stdint.h>
#include <string.h>

#include "mdl.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#ifdef COOLING
#include "cooling.h" /* before parameters.h */
#endif
#include "parameters.h"
#include "ilp.h"
#include "ilc.h"
#include "cl.h"
#include "moments.h"
#include "cosmo.h"
#include "fio.h"
#ifdef USE_HDF5
#include "iohdf5.h"
#endif
#ifdef USE_GRAFIC
#include "grafic.h"
#endif
#ifdef PLANETS
#include "ssio.h"
#endif

#ifdef __cplusplus
#define CAST(T,V) reinterpret_cast<T>(V)
#else
#define CAST(T,V) (V)
#endif

typedef uint_fast32_t local_t; /* Count of particles locally (per processor) */
typedef uint_fast64_t total_t; /* Count of particles globally (total number) */

static inline int d2i(double d)  {
    return (int)d;
}

static inline int64_t d2u64(double d) {
    return (uint64_t)d;
}

/*
** Handy type punning macro.
*/
#define UNION_CAST(x, sourceType, destType) \
	(((union {sourceType a; destType b;})x).b)
/*
** The following sort of definition should really be in a global
** configuration header file -- someday...
*/

#define CID_PARTICLE	0
#define CID_CELL	1
#define CID_GROUP	2
#define CID_RM		3
#define CID_BIN		4
#define CID_SHAPES	5
#define CID_PK          2
#define CID_PNG         2

/*
** These macros implement the index manipulation tricks we use for moving
** around in the top-tree. Note: these do NOT apply to the local-tree!
*/
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define SIBLING(i) 	(i^1)
#define PARENT(i)	(i>>1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

/*
** This is useful for debugging the very-active force calculation.
*/
#define A_VERY_ACTIVE  1


#define MAX_TIMERS		10

/*
** Memory models.  Each is a bit mask that indicates that additional fields should be
** added to the particle structure.
*/
#define PKD_MODEL_VELOCITY     (1<<0)  /* Velocity Required */
#define PKD_MODEL_ACCELERATION (1<<1)  /* Acceleration Required */
#define PKD_MODEL_POTENTIAL    (1<<2)  /* Potential Required */
#define PKD_MODEL_GROUPS       (1<<3)  /* Group profiling */
#define PKD_MODEL_HERMITE      (1<<4)  /* Hermite integrator */
#define PKD_MODEL_RELAXATION   (1<<5)  /* Trace relaxation */
#define PKD_MODEL_MASS         (1<<6)  /* Mass for each particle */
#define PKD_MODEL_SOFTENING    (1<<7)  /* Softening for each particle */
#define PKD_MODEL_VELSMOOTH    (1<<8)  /* Velocity Smoothing */
#define PKD_MODEL_SPH          (1<<9)  /* Sph Fields */
#define PKD_MODEL_STAR         (1<<10) /* Star Fields */
#define PKD_MODEL_PSMETRIC     (1<<11) /* Phase-space metric*/
#define PKD_MODEL_RUNGDEST     (1<<12) /* New Domain Decomposition */

#define PKD_MODEL_NODE_MOMENT  (1<<24) /* Include moment in the tree */
#define PKD_MODEL_NODE_ACCEL   (1<<25) /* mean accel on cell (for grav step) */
#define PKD_MODEL_NODE_VEL     (1<<26) /* center of mass velocity for cells */
#define PKD_MODEL_NODE_SPHBNDS (1<<27) /* Include 3 extra bounds in tree */

#define PKD_MODEL_NODE_BND     (1<<28) /* Include normal bounds in tree */
#define PKD_MODEL_NODE_BND6    (1<<29) /* Include phase-space bounds in tree */


#ifdef USE_PSD
typedef struct pLite {
    FLOAT r[6];
    float fMass;
    int i;
    uint8_t uRung;
    } PLITE;
#else
typedef struct pLite {
    FLOAT r[3];
    int i;
    uint8_t uRung;
    } PLITE;
#endif


typedef struct pIO {
    total_t iOrder;
    double r[3];
    double v[3];
    float fMass;
    float fSoft;
    float  fDensity;
    float  fPot;
    } PIO;

#define PKD_MAX_CLASSES 256
#define MAX_RUNG     63

/*
** Here we define some special reserved nodes. Node-0 is a sentinel or null node, node-1
** is here defined as the ROOT of the local tree (or top tree), node-2 is unused and
** node-3 is the root node for the very active tree.
*/
#define VAROOT          3
#define ROOT		1
#define NRESERVED_NODES MAX_RUNG+1

typedef struct partclass {
    float       fMass;    /* Particle mass */
    float       fSoft;    /* Current softening */
    FIO_SPECIES eSpecies; /* Species: dark, star, etc. */
    } PARTCLASS;

typedef struct hermitefields {
    double ad[3];
    double r0[3];
    double v0[3];
    double a0[3];
    double ad0[3];
    double rp[3];
    double vp[3];
    double app[3];
    double adpp[3];
    double dTime0;
    } HERMITEFIELDS;

typedef struct velsmooth {
    float vmean[3];
    float divv;
    float veldisp2;
    } VELSMOOTH;

typedef struct sphfields {
    char *pNeighborList; /* pointer to nearest neighbor list - compressed */
    double vPred[3];
    float u;	        /* thermal energy */ 
    float uPred;	/* predicted thermal energy */
    float uDot;
    float c;		/* sound speed */
    float divv;		
    float BalsaraSwitch;    /* Balsara viscosity reduction */
    float fMetals;	    /* mass fraction in metals, a.k.a, Z - tipsy output variable */

    /* diffusion */
    float diff; 
    float fMetalsPred;
    float fMetalsDot;

    } SPHFIELDS;

typedef struct starfields {
    float fTimer;  /* For gas -- cooling shutoff, for stars -- when formed */
    double totaltime; /* diagnostic -- get rid of it */
    } STARFIELDS;   

#define IORDERBITS 42    
#define IORDERMAX ((((uint64_t) 1)<<IORDERBITS)-1)

typedef struct particle {
    /*-----Base-Particle-Data----*/
    uint64_t iOrder     :  IORDERBITS;
    uint8_t  uNewRung   :  6;
    uint8_t  uRung      :  6;
    uint8_t  bSrcActive :  1;
    uint8_t  bDstActive :  1;
    uint8_t  iClass     :  8;
    double r[3];
    /*-----Used-for-Smooth-------*/
    float fBall;
    float fDensity;
#if USE_PSD
    float fDensityGrad[6];
#endif
    /* a, fPot, v, pGroup, pBin moved to memory models */

#ifdef PLANETS
    /* (collision stuff) */
    int iOrgIdx;		/* for tracking of mergers, aggregates etc. */
    FLOAT w[3];			/* spin vector */
    int iColor;			/* handy color tag */
    int iColflag;	        /* handy collision tag 1 for c1, 2 for c2*/
    uint64_t iOrderCol;              /* iOrder of colliding oponent.*/
    FLOAT dtCol;
    /* end (collision stuff) */
#ifdef SYMBA
    FLOAT rb[3]; /* position before drift */
    FLOAT vb[3]; /* velocity before drift */
    FLOAT drmin; /* minimum distance from neighbors normalized by Hill*/
    FLOAT drmin2; /* min. dis. during drift */
    uint64_t iOrder_VA[5]; /* iOrder's of particles within 3 hill radius*/
    int   i_VA[5];    /* pointers of particles */
    int   n_VA;       /* number of particles */
    double  hill_VA[5]; /* mutual hill radius calculated in grav.c */
    double a_VA[3];          /* accralation due to close encounters */
    /* int   iKickRung; */
#endif
#endif/* PLANETS */
    } PARTICLE;

#ifdef USE_PSD
typedef struct {
    double fCenter[6];
    double fMax[6];
    double size;
    } BND;

typedef struct {
    double scale[6];
    } PSMETRIC;
#define BND_COMBINE(b,b1,b2)\
{\
	int BND_COMBINE_j;\
	for (BND_COMBINE_j=0;BND_COMBINE_j<6;++BND_COMBINE_j) {\
		FLOAT BND_COMBINE_t1,BND_COMBINE_t2,BND_COMBINE_max,BND_COMBINE_min;\
		BND_COMBINE_t1 = (b1).fCenter[BND_COMBINE_j] + (b1).fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2).fCenter[BND_COMBINE_j] + (b2).fMax[BND_COMBINE_j];\
		BND_COMBINE_max = (BND_COMBINE_t1 > BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		BND_COMBINE_t1 = (b1).fCenter[BND_COMBINE_j] - (b1).fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2).fCenter[BND_COMBINE_j] - (b2).fMax[BND_COMBINE_j];\
		BND_COMBINE_min = (BND_COMBINE_t1 < BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		(b).fCenter[BND_COMBINE_j] = 0.5*(BND_COMBINE_max + BND_COMBINE_min);\
		(b).fMax[BND_COMBINE_j] = 0.5*(BND_COMBINE_max - BND_COMBINE_min);\
		}\
	}

#else
typedef struct bndBound {
    double fCenter[3];
    double fMax[3];
    double size;
    } BND;

#define BND_COMBINE(b,b1,b2)\
{\
	int BND_COMBINE_j;\
	for (BND_COMBINE_j=0;BND_COMBINE_j<3;++BND_COMBINE_j) {\
		FLOAT BND_COMBINE_t1,BND_COMBINE_t2,BND_COMBINE_max,BND_COMBINE_min;\
		BND_COMBINE_t1 = (b1).fCenter[BND_COMBINE_j] + (b1).fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2).fCenter[BND_COMBINE_j] + (b2).fMax[BND_COMBINE_j];\
		BND_COMBINE_max = (BND_COMBINE_t1 > BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		BND_COMBINE_t1 = (b1).fCenter[BND_COMBINE_j] - (b1).fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2).fCenter[BND_COMBINE_j] - (b2).fMax[BND_COMBINE_j];\
		BND_COMBINE_min = (BND_COMBINE_t1 < BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		(b).fCenter[BND_COMBINE_j] = 0.5*(BND_COMBINE_max + BND_COMBINE_min);\
		(b).fMax[BND_COMBINE_j] = 0.5*(BND_COMBINE_max - BND_COMBINE_min);\
		}\
	}

#endif

typedef struct {
    double *fCenter;
    double *fMax;
    double *size;
    } pBND;

#define MINDIST(bnd,pos,min2) {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;					\
    for (BND_j=0;BND_j<3;++BND_j) {\
	BND_dMin = fabs((bnd).fCenter[BND_j] - (pos)[BND_j]) - (bnd).fMax[BND_j]; \
	if (BND_dMin > 0) (min2) += BND_dMin*BND_dMin;			\
	}\
    }

#if USE_PSD
#define MINDIST6(bnd,pos,min2, scale) {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;					\
    for (BND_j=0;BND_j<6;++BND_j) {\
	BND_dMin = fabs((bnd).fCenter[BND_j] - (pos)[BND_j]) - (bnd).fMax[BND_j]; \
	if (BND_dMin > 0) (min2) += BND_dMin*BND_dMin*scale[BND_j]*scale[BND_j];			\
	}\
    }
#endif

static inline int IN_BND(const FLOAT *R,const pBND *b) {
    int i;
    for( i=0; i<3; i++ )
	if ( R[i]<b->fCenter[i]-b->fMax[i] || R[i]>=b->fCenter[i]+b->fMax[i] )
	    return 0;
    return 1;
    }


/*
** General partition macro
** LT,LE: Compare less-than/less-than or equal
** ii,dj: Increment i and decrement j
** SWAP: Swap the i'th and j'th element
** LOWER,UPPER: comparison predicates
** e.g.,
** PARTICLE *pi = pkdParticle(pkd,i);
** PARTICLE *pj = pkdParticle(pkd,j);
**    PARTITION(pi<pj,pi<=pj,
**              pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
**              pkdSwapParticle(pkd,pi,pj),
**	        pi->r[d] >= fSplit,pj->r[d] < fSplit);
** When finished, the 'I' variable points to the first element of
** the upper partition (or one past the end).
** NOTE: Because this function supports tiled data structures,
**       the LT, followed by "if (LE)" needs to remain this way.
*/
#define PARTITION(LT,LE,INCI,DECJ,SWAP,LOWER,UPPER)	\
    {							\
    while ((LT) && (LOWER)) { INCI; }			\
    if ((LE) && (LOWER)) { INCI; }			\
    else {						\
	while ((LT) && (UPPER)) { DECJ; }		\
	while (LT) {					\
		    { SWAP; }				\
		    do { DECJ; } while (UPPER);	\
		    do { INCI; } while (LOWER);		\
	    }						\
	}						\
    }

#define NEW_STACK(S,inc) \
    int S ## __s=0, S ## __ns=inc, S ## __inc=inc; \
    do { S = malloc((S ## __ns)*sizeof(*S)); assert(S != NULL); } while (0)
#define FREE_STACK(S) do { free(S); } while (0)

#define PUSH(S,v) do { S[(S ## __s)++] = (v); } while (0)
#define POP(S) (S[--(S ## __s)])
#define STACK_EMPTY(S) (S ## __s == 0)
#define EXTEND_STACK(S) do { \
    if ( (S ## __s)+1 >= (S ## __ns) ) { \
        assert( (S ## __s)+1 == (S ## __ns) ); \
        (S ## __ns) += (S ## __inc); \
        S = realloc(S,(S ## __ns)*sizeof(*S)); \
        } \
} while (0)
#define CLEAR_STACK(S) do { S ## __s=0; } while (0)


typedef struct kdNode {
    double r[3];
    int iLower;
    int iParent;
    int pLower;		/* also serves as thread id for the LTT */
    int pUpper;		/* pUpper < 0 indicates no particles in tree! */
    float bMax;
    float fSoft2;
    uint32_t nActive; /* local active count used for walk2 */
    uint8_t uMinRung;
    uint8_t uMaxRung;
    uint8_t bSrcActive;
    uint8_t bDstActive;
    } KDN;


#if 0
/*
** This structure should allow for both k-D trees as well as spatial binary 
** binary trees and finally Peano-Hilbert SFC type binary trees. For SFC type
** binary trees we don't really need the double precision split if we have the 
** key of a cell (which is not included here).
*/
typedef struct kdNode {
    union uCellFields {	
	struct bucketFields { /* 32 bytes for M=8 */
	    uint32_t iPart[M];  /* loop until iPart[i] == 0xffffffff */
	    } bf;
	struct cellFields {  /* 32 bytes */
	    /*
	    ** iDummy will always be == 0xffffffff for a cell.
	    */
	    uint32_t iDummy;
	    /*
	    ** Since cells will be allocated from high memory downward the 
	    ** indecies here will all be negative. Cells will be allocated
	    ** until they collide with the particles which are allocated in 
	    ** a heap fashion from low to high memory.
	    */
	    int iLower;
	    int iUpper;
	    int iParent;
	    double dSplit;
	    uint16_t iDim;
	    uint16_t pidLower;
	    uint16_t pidUpper;
	    uint16_t pidParent;
	    } cf;
	} u;
    } KDN;
#endif

typedef struct sphBounds {
    struct minmaxBound {
	double min[3];
	double max[3];
    } A,B,BI;
} SPHBNDS;


#define NMAX_OPENCALC	1000

#define MAXSIDE(fMax,b) {\
    if ((fMax)[0] > (fMax)[1]) {\
	if ((fMax)[0] > (fMax)[2]) b = 2.0*(fMax)[0];\
	else b = 2.0*(fMax)[2];\
	}\
    else {\
	if ((fMax)[1] > (fMax)[2]) b = 2.0*(fMax)[1];\
	else b = 2.0*(fMax)[2];\
	}\
    }

#define MINSIDE(fMax,b) {\
    if ((fMax)[0] < (fMax)[1]) {\
	if ((fMax)[0] < (fMax)[2]) b = 2.0*(fMax)[0];\
	else b = 2.0*(fMax)[2];\
	}\
    else {\
	if ((fMax)[1] < (fMax)[2]) b = 2.0*(fMax)[1];\
	else b = 2.0*(fMax)[2];\
	}\
    }

#define CALCAXR(fMax,axr) {					\
    if ((fMax)[0] < (fMax)[1]) {				\
	if ((fMax)[1] < (fMax)[2]) {				\
	    if ((fMax)[0] > 0) axr = (fMax)[2]/(fMax)[0];	\
	    else axr = 1e6;					\
	}							\
	else if ((fMax)[0] < (fMax)[2]) {			\
	    if ((fMax)[0] > 0) axr = (fMax)[1]/(fMax)[0];	\
	    else axr = 1e6;					\
	}							\
	else if ((fMax)[2] > 0) axr = (fMax)[1]/(fMax)[2];	\
	else axr = 1e6;						\
    }								\
    else if ((fMax)[0] < (fMax)[2]) {				\
	if ((fMax)[1] > 0) axr = (fMax)[2]/(fMax)[1];		\
	else axr = 1e6;						\
    }								\
    else if ((fMax)[1] < (fMax)[2]) {				\
	if ((fMax)[1] > 0) axr = (fMax)[0]/(fMax)[1];		\
	else axr = 1e6;						\
    }								\
    else if ((fMax)[2] > 0) axr = (fMax)[0]/(fMax)[2];		\
    else axr = 1e6;						\
}


#define CALCOPEN(pkdn,minside) {					\
        FLOAT CALCOPEN_d2 = 0;						\
	FLOAT CALCOPEN_b;						\
        int CALCOPEN_j;							\
	pBND CALCOPEN_bnd;						\
        pkdNodeBnd(pkd, pkdn, &CALCOPEN_bnd);				\
        for (CALCOPEN_j=0;CALCOPEN_j<3;++CALCOPEN_j) {                  \
            FLOAT CALCOPEN_d = fabs(CALCOPEN_bnd.fCenter[CALCOPEN_j] - (pkdn)->r[CALCOPEN_j]) + \
                CALCOPEN_bnd.fMax[CALCOPEN_j];                          \
            CALCOPEN_d2 += CALCOPEN_d*CALCOPEN_d;                       \
            }								\
	MAXSIDE(CALCOPEN_bnd.fMax,CALCOPEN_b);				\
	if (CALCOPEN_b < minside) CALCOPEN_b = minside;			\
	if (CALCOPEN_b*CALCOPEN_b < CALCOPEN_d2) CALCOPEN_b = sqrt(CALCOPEN_d2); \
	(pkdn)->bMax = CALCOPEN_b;					\
	}

#if (0)
#define CALCOPEN(pkdn) {					\
        FLOAT CALCOPEN_d2 = 0;						\
        int CALCOPEN_j;							\
	pBND CALCOPEN_bnd; \
	pkdNodeBnd(pkd, pkdn, &CALCOPEN_bnd);			\
        for (CALCOPEN_j=0;CALCOPEN_j<3;++CALCOPEN_j) {                  \
            FLOAT CALCOPEN_d = fabs(CALCOPEN_bnd.fCenter[CALCOPEN_j] - (pkdn)->r[CALCOPEN_j]) + \
                CALCOPEN_bnd.fMax[CALCOPEN_j];                          \
            CALCOPEN_d2 += CALCOPEN_d*CALCOPEN_d;                       \
            }\
        CALCOPEN_d2 = sqrt(CALCOPEN_d2);	  \
        if (CALCOPEN_d2 < (pkdn)->bMax) (pkdn)->bMax = CALCOPEN_d2;	  \
	}
#endif

/*
** Components required for tree walking.
*/
#ifndef LOCAL_EXPANSION
typedef struct CheckElt {
    int iCell;
    int id;
    double cOpen;
    FLOAT rOffset[3];
    } CELT;
#endif

typedef struct CheckStack {
#ifdef LOCAL_EXPANSION
    ILPCHECKPT PartChkPt;
    ILCCHECKPT CellChkPt;
    CL cl;
#else
    int nPart;
    int nCell;
    int nCheck;
    CELT *Check;
#endif
    FLOCR L;
    float dirLsum;
    float normLsum;
    float fWeight;
    } CSTACK;

/*
** components required for time-step calculation (only grav.c)
*/

typedef struct RhoLocalArray {
    double d2;
    double m;
    } RHOLOCAL;

typedef struct ewaldTable {
    double hx,hy,hz;
    double hCfac,hSfac;
    } EWT;

struct EwaldVariables {
    double fEwCut2,fInner2,alpha,alpha2,k1,ka;
    double Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;
    EWT *ewt;
    int nMaxEwhLoop;
    int nEwhLoop;
    int nReps,nEwReps;
    };

/*
** components required for groupfinder:  --J.D.--
*/
typedef struct remoteMember {
    int iPid;
    int iIndex;
    } FOFRM;

typedef struct groupData {
    int iLocalId;
    int iGlobalId;
    FLOAT fMass;
    FLOAT fRMSRadius;
    FLOAT r[3];
    FLOAT rcom[3];
    FLOAT potordenmax;
    FLOAT v[3];
    int nLocal;
    int nTotal;
    int bMyGroup;
    int nLSubGroups;
    void *lSubGroups;
    int nSubGroups;
    void *subGroups;
    int nRemoteMembers;
    int iFirstRm;
} FOFGD;

typedef struct groupBin {
  FLOAT fRadius;
  int nMembers;
  FLOAT fMassInBin;
} FOFBIN;

typedef struct profileBin {
    double dRadius;
    double dMassInBin;
    double dVolume;
    double L[3];
    double vel_radial;
    double vel_radial_sigma;
    double vel_tang_sigma;
    uint64_t nParticles;
  } PROFILEBIN;

typedef struct shapesBin {
    double com[3];
    double dMassEnclosed;
    double dInertia[3][3];
    double ell_matrix[3][3];
    double ell_center[3];
    } SHAPESBIN;

#define PKD_GROUP_SIZE 64

typedef struct {
    float r[3];
    float v[3];
    float a[3];
    float fMass;
    float fSoft;
    float fSmooth2;
    } PINFOIN;

typedef struct {
    float a[3];
    float p;
    } PINFOOUT;

typedef struct pkdContext {
    MDL mdl;
    int idSelf;
    int nThreads;
    int nStore;
    int nRejects;
    int nLocal;
    int nVeryActive;
    int nActive;
    int nTreeBitsLo;
    int nTreeBitsHi;
    int iTreeMask;
    int nTreeTiles;
    int nMaxNodes;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
#ifdef USE_PSD
    FLOAT fPeriod[6];
#else
    FLOAT fPeriod[3];
#endif
    char *kdTopPRIVATE; /* Because this is a variable size, we use a char pointer, not a KDN pointer! */
    char **kdNodeListPRIVATE; /* BEWARE: also char instead of KDN */
    int iTopRoot;
    int nNodes;
    int nNodesFull;     /* number of nodes in the full tree (including very active particles) */
    int nNonVANodes;    /* number of nodes *not* in Very Active Tree, or index to the start of the VA nodes (except VAROOT) */
    BND bnd;
    size_t iTreeNodeSize;
    size_t iParticleSize;
    PARTICLE *pStorePRIVATE, *pStorePRIVATE2;
    PARTICLE *pTempPRIVATE;
    PARTCLASS *pClass;
    float fSoftFix;
    float fSoftFac;
    float fSoftMax;
    int nClasses;
    int nMaxBucketActive;
    PARTICLE **piActive;
    PARTICLE **piInactive;
    PLITE *pLite;
#ifdef COOLING
    COOL *Cool; /* Cooling Context */
#endif

#ifdef MPI_VERSION
    MDL_Datatype typeParticle;
#endif

    int nDomainRungs;
    int iFirstDomainRung;

    /*
    ** Advanced memory models
    */
    int oAcceleration; /* Three doubles */
    int oVelocity; /* Three doubles */
    int oPotential; /* One float */
    int oGroup; /* One int32 */
    int oMass; /* One float */
    int oSoft; /* One float */
    int oSph; /* Sph structure */
    int oStar; /* Star structure */
    int oHermite; /* Hermite structure */
    int oRelaxation;
    int oVelSmooth;
    int oPsMetric; /* Phase-space metric for density/group finding */
    int oRungDest; /* Destination processor for each rung */

    /*
    ** Advanced memory models - Tree Nodes
    */
    int oNodePosition; /* Three doubles */
    int oNodeVelocity; /* Three doubles */
    int oNodeAcceleration; /* Three doubles */
    int oNodeBmax;
    int oNodeSoft;
    int oNodeMom; /* an FMOMR */
    int oNodeBnd;
    int oNodeSphBounds; /* Three Bounds */
    int oNodeBnd6;

    /*
    ** Tree walk variables.
    */
    int nMaxStack;
    CSTACK *S;
#ifdef LOCAL_EXPANSION
    ILP ilp;
    ILC ilc;
    CLTILE clFreeList;
    CL cl;
    CL clNew;
#else
    ILP *ilp;
    ILC *ilc;
    CELT *Check;
    int nMaxPart, nMaxCell;
    int nMaxCheck;
#endif

    /*
    ** Opening angle table for mass weighting.
    */
#ifdef USE_DEHNEN_THETA
    float *fCritTheta;
    float *fCritMass;
    int nCritBins;
    float dCritLogDelta;
    float dCritThetaMin;
    float dCritThetaMax;
#else
    float fiCritTheta;
#endif
    /*
    ** New activation methods
    */
    uint8_t uMinRungActive;
    uint8_t uMaxRungActive;
    uint8_t uRungVeryActive;    /* NOTE: The first very active particle is at iRungVeryActive + 1 */

    /*
    ** Ewald summation setup.
    */
    MOMC momRoot;
    struct EwaldVariables ew;

    /*
    ** Timers stuff.
    */
    struct timer {
	double sec;
	double stamp;
	double system_sec;
	double system_stamp;
	double wallclock_sec;
	double wallclock_stamp;
	int iActive;
	} ti[MAX_TIMERS];
    int nGroups;
    FOFGD *groupData;
    int nRm;
    int nMaxRm;
    FOFRM *remoteMember;
    int nBins;

    FOFBIN *groupBin;

    PROFILEBIN *profileBins;

    /*
    ** Oh heck, just put all the parameters in here!
    ** This is set in pkdInitStep.
    */
    struct parameters param;
#ifdef PLANETS
    double dDeltaEcoll;
    double dSunMass;
    int    iCollisionflag; /*call pkddocollisionveryactive if iCollisionflag=1*/
#endif

#ifdef USE_CUDA
    void *cudaCtx;
#endif

    MDLGRID grid;
    float *gridData;
    } * PKD;

static inline void pkdMinMax( double *dVal, double *dMin, double *dMax ) {
    dMin[0] = dVal[0] < dMin[0] ? dVal[0] : dMin[0];
    dMin[1] = dVal[1] < dMin[1] ? dVal[1] : dMin[1];
    dMin[2] = dVal[2] < dMin[2] ? dVal[2] : dMin[2];
    dMax[0] = dVal[0] > dMax[0] ? dVal[0] : dMax[0];
    dMax[1] = dVal[1] > dMax[1] ? dVal[1] : dMax[1];
    dMax[2] = dVal[2] > dMax[2] ? dVal[2] : dMax[2];
    }

static inline void pkdMinMax6( double *dVal0, double *dVal1, double *dMin, double *dMax ) {
    dMin[0] = dVal0[0] < dMin[0] ? dVal0[0] : dMin[0];
    dMin[1] = dVal0[1] < dMin[1] ? dVal0[1] : dMin[1];
    dMin[2] = dVal0[2] < dMin[2] ? dVal0[2] : dMin[2];
    dMin[3] = dVal1[0] < dMin[3] ? dVal1[0] : dMin[3];
    dMin[4] = dVal1[1] < dMin[4] ? dVal1[1] : dMin[4];
    dMin[5] = dVal1[2] < dMin[5] ? dVal1[2] : dMin[5];

    dMax[0] = dVal0[0] > dMax[0] ? dVal0[0] : dMax[0];
    dMax[1] = dVal0[1] > dMax[1] ? dVal0[1] : dMax[1];
    dMax[2] = dVal0[2] > dMax[2] ? dVal0[2] : dMax[2];
    dMax[3] = dVal1[0] > dMax[3] ? dVal1[0] : dMax[3];
    dMax[4] = dVal1[1] > dMax[4] ? dVal1[1] : dMax[4];
    dMax[5] = dVal1[2] > dMax[5] ? dVal1[2] : dMax[5];
    }

/* New, rung based ACTIVE/INACTIVE routines */
static inline int pkdIsDstActive(PARTICLE *p,uint8_t uRungLo,uint8_t uRungHi) {
    return((p->uRung >= uRungLo)&&(p->uRung <= uRungHi)&&p->bDstActive);
    }

static inline int pkdIsSrcActive(PARTICLE *p,uint8_t uRungLo,uint8_t uRungHi) {
    return((p->uRung >= uRungLo)&&(p->uRung <= uRungHi)&&p->bSrcActive);
    }

static inline int pkdIsRungRange(PARTICLE *p,uint8_t uRungLo,uint8_t uRungHi) {
    return((p->uRung >= uRungLo)&&(p->uRung <= uRungHi));
    }

static inline int pkdRungVeryActive(PKD pkd) {
    return pkd->uRungVeryActive;
    }
static inline int pkdIsVeryActive(PKD pkd, PARTICLE *p) {
    return p->uRung > pkd->uRungVeryActive;
    }

static inline int pkdIsRungActive(PKD pkd, uint8_t uRung ) {
    return uRung >= pkd->uMinRungActive && uRung <= pkd->uMaxRungActive;
    }
static inline int pkdIsActive(PKD pkd, PARTICLE *p ) {
    return pkdIsRungActive(pkd,p->uRung);
    }

/*
** A tree node is of variable size.  The following routines are used to
** access individual fields.
*/
static inline KDN *pkdTreeBase( PKD pkd ) {
    return (KDN *)pkd->kdNodeListPRIVATE;
    }
static inline size_t pkdNodeSize( PKD pkd ) {
    return pkd->iTreeNodeSize;
    }
static inline size_t pkdMaxNodeSize() {
#ifdef USE_PSD
    return sizeof(KDN) + 13*sizeof(double) + sizeof(FMOMR) + 6*sizeof(double) + sizeof(SPHBNDS);
#else
    return sizeof(KDN) +  7*sizeof(double) + sizeof(FMOMR) + 6*sizeof(double) + sizeof(SPHBNDS);
#endif
    }
static inline void pkdCopyNode(PKD pkd, KDN *a, KDN *b) {
    memcpy(a,b,pkdNodeSize(pkd));
    }
static inline void *pkdNodeField( KDN *n, int iOffset ) {
    char *v = (char *)n;
    /*assert(iOffset);*/ /* Remove this for better performance */
    return (void *)(v + iOffset);
    }

static inline FMOMR *pkdNodeMom(PKD pkd,KDN *n) {
    return CAST(FMOMR *,pkdNodeField(n,pkd->oNodeMom));
    }
static inline double *pkdNodeVel( PKD pkd, KDN *n ) {
    return CAST(double *,pkdNodeField(n,pkd->oNodeVelocity));
    }
static inline double *pkdNodeAccel( PKD pkd, KDN *n ) {
    return CAST(double *,pkdNodeField(n,pkd->oNodeAcceleration));
    }
static inline SPHBNDS *pkdNodeSphBounds( PKD pkd, KDN *n ) {
    return CAST(SPHBNDS *,pkdNodeField(n,pkd->oNodeSphBounds));
    }

static inline void pkdNodeBnd( PKD pkd, KDN *n, pBND *bnd ) {
    const int o = pkd->oNodeBnd;
    const int e = 3*sizeof(double)*(1+(pkd->oNodeBnd6!=0));
    bnd->fCenter = CAST(double *,pkdNodeField(n,o));
    bnd->fMax = CAST(double *,pkdNodeField(n,o+e));
    bnd->size = CAST(double *,pkdNodeField(n,o+e+e));
    }

static inline KDN *pkdNode(PKD pkd,KDN *pBase,int iNode) {
    return (KDN *)&((char *)pBase)[pkd->iTreeNodeSize*iNode];
    }
int pkdNodes(PKD pkd);
void pkdExtendTree(PKD pkd);
static inline KDN *pkdTreeNode(PKD pkd,int iNode) {
    return (KDN *)&pkd->kdNodeListPRIVATE[(iNode>>pkd->nTreeBitsLo)][pkd->iTreeNodeSize*(iNode&pkd->iTreeMask)];
    }
void *pkdTreeNodeGetElement(void *vData,int i,int iDataSize);
static inline KDN *pkdTopNode(PKD pkd,int iNode) {
    return (KDN *)&pkd->kdTopPRIVATE[pkd->iTreeNodeSize*iNode];
    }
void pkdAllocateTopTree(PKD pkd,int nCell);

/*
** The size of a particle is variable based on the memory model.
** The following three routines must be used instead of accessing pStore
** directly.  pkdParticle will return a pointer to the i'th particle.
** The Size and Base functions are intended for cache routines; no other
** code should care about sizes of the particle structure.
*/
static inline PARTICLE *pkdParticleBase( PKD pkd ) {
    return pkd->pStorePRIVATE;
    }
static inline size_t pkdParticleSize( PKD pkd ) {
    return pkd->iParticleSize;
    }
static inline size_t pkdParticleMemory(PKD pkd) {
    return (pkd->iParticleSize + sizeof(PLITE)) * (pkd->nStore+1);
    }
static inline PARTICLE *pkdParticle( PKD pkd, int i ) {
    char *v = (char *)pkd->pStorePRIVATE;
    PARTICLE *p = (PARTICLE *)(v + ((uint64_t)i)*pkd->iParticleSize);
    return p;
    }
static inline PARTICLE *pkdParticle2( PKD pkd, int i ) {
    char *v = (char *)pkd->pStorePRIVATE2;
    PARTICLE *p = (PARTICLE *)(v + ((uint64_t)i)*pkd->iParticleSize);
    return p;
    }
static inline void pkdSaveParticle(PKD pkd, PARTICLE *a) {
    memcpy(pkd->pTempPRIVATE,a,pkdParticleSize(pkd));
    }
static inline void pkdLoadParticle(PKD pkd, PARTICLE *a) {
    memcpy(a,pkd->pTempPRIVATE,pkdParticleSize(pkd));
    }
static inline void pkdCopyParticle(PKD pkd, PARTICLE *a, PARTICLE *b) {
    memcpy(a,b,pkdParticleSize(pkd));
    }
static inline void pkdSwapParticle(PKD pkd, PARTICLE *a, PARTICLE *b) {
    pkdSaveParticle(pkd,a);
    pkdCopyParticle(pkd,a,b);
    pkdLoadParticle(pkd,b);
    }

static inline void *pkdField( PARTICLE *p, int iOffset ) {
    char *v = (char *)p;
    /*assert(iOffset);*/ /* Remove this for better performance */
    return (void *)(v + iOffset);
    }

static inline int32_t *pkdInt32( PARTICLE *p, int iOffset ) {
    char *v = (char *)p;
    return (int32_t *)(v + iOffset);
    }

/* Here is the new way of getting mass and softening */
static inline float pkdMass( PKD pkd, PARTICLE *p ) {
    if ( pkd->oMass ) {
	float *pMass = CAST(float *,pkdField(p,pkd->oMass));
	return *pMass;
	}
    return pkd->pClass[p->iClass].fMass;
    }
static inline float pkdSoft0( PKD pkd, PARTICLE *p ) {
    if ( pkd->oSoft ) {
	float *pSoft = CAST(float *,pkdField(p,pkd->oSoft));
	return *pSoft;
	}
    return pkd->pClass[p->iClass].fSoft;
    }
static inline float pkdSoft( PKD pkd, PARTICLE *p ) {
    float fSoft;

    if ( pkd->fSoftFix >= 0.0 ) fSoft = pkd->fSoftFix;
    else fSoft = pkdSoft0(pkd,p);
    fSoft *= pkd->fSoftFac;
    if ( fSoft > pkd->fSoftMax ) fSoft = pkd->fSoftMax;
    return fSoft;
    }
static inline FIO_SPECIES pkdSpecies( PKD pkd, PARTICLE *p ) {
    return pkd->pClass[p->iClass].eSpecies;
    }
static inline double *pkdVel( PKD pkd, PARTICLE *p ) {
    return CAST(double *,pkdField(p,pkd->oVelocity));
    }
static inline float *pkdAccel( PKD pkd, PARTICLE *p ) {
    return CAST(float *,pkdField(p,pkd->oAcceleration));
    }
static inline float *pkdPot( PKD pkd, PARTICLE *p ) {
    return CAST(float *,pkdField(p,pkd->oPotential));
    }
static inline uint16_t *pkdRungDest( PKD pkd, PARTICLE *p ) {
    return CAST(uint16_t *,pkdField(p,pkd->oRungDest));
    }
/* Sph variables */
static inline SPHFIELDS *pkdSph( PKD pkd, PARTICLE *p ) {
    return ((SPHFIELDS *) pkdField(p,pkd->oSph));
    }
static inline STARFIELDS *pkdStar( PKD pkd, PARTICLE *p ) {
    return ((STARFIELDS *) pkdField(p,pkd->oStar));
    }
static inline double *pkd_vPred( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->vPred[0]);
    }
static inline float *pkd_u( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->u);
    }
static inline float *pkd_uPred( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->uPred);
    }
static inline float *pkd_uDot( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->uDot);
    }
static inline float *pkd_c( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->c);
    }
static inline float *pkd_divv( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->divv);
    }
static inline float *pkd_fMetals( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->fMetals);
    }
static inline float *pkd_diff( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->diff);
    }
static inline float *pkd_fMetalsDot( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->fMetalsDot);
    }
static inline float *pkd_fMetalsPred( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->fMetalsPred);
    }
static inline char **pkd_pNeighborList( PKD pkd, PARTICLE *p ) {
    return &(((SPHFIELDS *) pkdField(p,pkd->oSph))->pNeighborList);
    }

#ifdef USE_PSD
/* Phase-space variables */
static inline PSMETRIC *pkdPsMetric( PKD pkd, PARTICLE *p ) {
    return ((PSMETRIC *) pkdField(p,pkd->oPsMetric));
    }
#endif

static inline float *pkd_Timer( PKD pkd, PARTICLE *p ) {
    return &(((STARFIELDS *) pkdField(p,pkd->oStar))->fTimer);
    }

static inline int pkdIsDeleted(PKD pkd,PARTICLE *p) {
    return (pkdSpecies(pkd,p) == FIO_SPECIES_LAST);
    }

static inline int pkdIsNew(PKD pkd,PARTICLE *p) {
    return (p->iOrder == IORDERMAX);
    }


typedef struct CacheStatistics {
    double dpNumAccess;
    double dpMissRatio;
    double dpCollRatio;
    double dcNumAccess;
    double dcMissRatio;
    double dcCollRatio;
    } CASTAT;

/*
** From tree.c:
*/
void pkdVATreeBuild(PKD pkd,int nBucket);
void pkdTreeBuild(PKD pkd,int nBucket,KDN *pkdn,int bExcludeVeryActive);
void pkdCombineCells1(PKD,KDN *pkdn,KDN *p1,KDN *p2);
void pkdCombineCells2(PKD,KDN *pkdn,KDN *p1,KDN *p2);
void pkdDistribCells(PKD,int,KDN *);
void pkdCalcRoot(PKD,MOMC *);
void pkdDistribRoot(PKD,MOMC *);
void pkdTreeNumSrcActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi);
void pkdBoundWalk(PKD pkd,BND *pbnd,uint8_t uRungLo,uint8_t uRungHi,uint32_t *pnActive,uint32_t *pnContained);

#include "parameters.h"
/*
** From pkd.c:
*/
double pkdGetTimer(PKD,int);
double pkdGetSystemTimer(PKD,int);
double pkdGetWallClockTimer(PKD,int);
void pkdClearTimer(PKD,int);
void pkdStartTimer(PKD,int);
void pkdStopTimer(PKD,int);
void pkdInitialize(
    PKD *ppkd,MDL mdl,int nStore,int nBucket,int nTreeBitsLo, int nTreeBitsHi,
    int iCacheSize,FLOAT *fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,
    uint64_t mMemoryModel, int nDomainRungs);
void pkdFinish(PKD);
size_t pkdClCount(PKD pkd);
size_t pkdClMemory(PKD pkd);
size_t pkdIlcMemory(PKD pkd);
size_t pkdIlpMemory(PKD pkd);
size_t pkdTreeMemory(PKD pkd);
void pkdReadFIO(PKD pkd,FIO fio,uint64_t iFirst,int nLocal,double dvFac, double dTuFac);
#ifdef USE_MDL_IO
void pkdIOInitialize( PKD pkd, int nLocal);
#endif

void pkdSetSoft(PKD pkd,double dSoft);
void pkdSetCrit(PKD pkd,double dCrit);
void pkdCalcBound(PKD,BND *);
void pkdEnforcePeriodic(PKD,BND *);
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul);

void pkdBucketWeight(PKD pkd,int iBucket,FLOAT fWeight);
int pkdWeight(PKD,int,FLOAT,int,int,int,int *,int *,FLOAT *,FLOAT *);
void pkdCountVA(PKD,int,FLOAT,int *,int *);
double pkdTotalMass(PKD pkd);
int pkdLowerPart(PKD,int,FLOAT,int,int);
int pkdUpperPart(PKD,int,FLOAT,int,int);
int pkdWeightWrap(PKD,int,FLOAT,FLOAT,int,int,int,int,int *,int *);
int pkdLowerPartWrap(PKD,int,FLOAT,FLOAT,int,int,int);
int pkdUpperPartWrap(PKD,int,FLOAT,FLOAT,int,int,int);
int pkdLowerOrdPart(PKD,uint64_t,int,int);
int pkdUpperOrdPart(PKD,uint64_t,int,int);
int pkdActiveOrder(PKD);

/*#define PEANO_HILBERT_KEY_MAX 0x3ffffffffffull*/ /* 2d */
#define PEANO_HILBERT_KEY_MAX 0x7fffffffffffffffull /* 3d */
void pkdPeanoHilbertDecomp(PKD pkd, uint64_t nTotal, int nRungs, int iMethod);
void pkdRungOrder(PKD pkd, int iRung, total_t *nMoved);
int pkdColRejects(PKD,int);
int pkdColRejects_Old(PKD,int,FLOAT,FLOAT,int);

int pkdSwapRejects(PKD,int);
int pkdSwapSpace(PKD);
int pkdFreeStore(PKD);
int pkdLocal(PKD);
int pkdActive(PKD);
int pkdInactive(PKD);

int pkdNumSrcActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi);
int pkdNumDstActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi);
int pkdColOrdRejects(PKD,uint64_t,int);
void pkdLocalOrder(PKD);
uint32_t pkdWriteFIO(PKD pkd,FIO fio,double dvFac);
uint32_t pkdWriteTipsy(PKD,char *,uint64_t,int,double,int);
void
pkdGravAll(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bPeriodic,
	   int iOrder,int bEwald,double fEwCut,double fEwhCut,double dThetaMin,double dThetaMax,
	   int *nActive,double *pdPartSum, double *pdCellSum,CASTAT *pcs, double *pdFlop);
void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double *L,double *F,double *W);
void pkdDrift(PKD pkd,double dDelta,double,double,uint8_t uRungLo,uint8_t uRungHi);
void pkdScaleVel(PKD pkd,double dvFac);
void pkdStepVeryActiveKDK(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dStep, double dTime, double dDelta,
			  int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double diCrit2,
			  int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass);
#ifdef HERMITE
void
pkdStepVeryActiveHermite(PKD pkd, double dStep, double dTime, double dDelta,
			 int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double diCrit2,
			 int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass);
void pkdCopy0(PKD pkd,double dTime);
void pkdPredictor(PKD pkd,double dTime);
void pkdCorrector(PKD pkd,double dTime);
void pkdSunCorrector(PKD pkd,double dTime,double dSunMass);
void pkdPredictorInactive(PKD pkd,double dTime);
void pkdAarsethStep(PKD pkd, double dEta);
void pkdFirstDt(PKD pkd);
#endif /* Hermite */
void pkdKickKDKOpen(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
void pkdKickKDKClose(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
void pkdKick(PKD pkd,double dTime,double dDelta,double,double,double,uint8_t uRungLo,uint8_t uRungHi);
void pkdSwapAll(PKD pkd, int idSwap);
void pkdInitStep(PKD pkd,struct parameters *p,CSM csm);
void pkdSetRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung);
void pkdZeroNewRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung);
void pkdActiveRung(PKD pkd, int iRung, int bGreater);
int pkdCurrRung(PKD pkd,uint8_t uRung);
void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
		  double dEta,double dVelFac,double dAccFac,
		  int bDoGravity,int bEpsAcc,int bSqrtPhi,double dhMinOverSoft);
void pkdSphStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,double dAccFac);
void pkdStarForm(PKD pkd, double dRateCoeff, double dTMax, double dDenMin,
		 double dDelta, double dTime,
		 double dInitStarMass, double dESNPerStarMass, double dtCoolingShutoff,
		 double dtFeedbackDelay,  double dMassLossPerStarMass,    
		 double dZMassPerStarMass, double dMinGasMass,
		 int bdivv, int *nFormed, double *dMassFormed,
		 int *nDeleted);
void pkdCooling(PKD pkd,double,double,int,int,int,int);
#define CORRECTENERGY_IN 1
#define CORRECTENERGY_OUT 2
#define CORRECTENERGY_SPECIAL 3
void pkdCorrectEnergy(PKD pkd, double dTuFac, double z, double dTime, int iType );
void pkdDensityStep(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, double dEta, double dRhoFac);
int pkdUpdateRung(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
		  uint8_t uRung,int iMaxRung,int *nRungCount);
uint8_t pkdDtToRung(double dT, double dDelta, uint8_t uMaxRung);
int pkdOrdWeight(PKD pkd,uint64_t iOrdSplit,int iSplitSide,int iFrom,int iTo,
		 int *pnLow,int *pnHigh);
void pkdDeleteParticle(PKD pkd, PARTICLE *p);
void pkdNewParticle(PKD pkd, PARTICLE *p);
int pkdResetTouchRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
void pkdSetRungVeryActive(PKD pkd, int iRung );
int pkdIsGas(PKD,PARTICLE *);
int pkdIsDark(PKD,PARTICLE *);
int pkdIsStar(PKD,PARTICLE *);
void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
		  int *nDeltaStar);
void pkdNewOrder(PKD pkd, int nStart);

struct outGetNParts { 
    local_t n;
    local_t nGas;
    local_t nDark;
    local_t nStar;
    total_t nMaxOrder;
    };

void pkdGetNParts(PKD pkd, struct outGetNParts *out );
void pkdSetNParts(PKD pkd, int nGas, int nDark, int nStar);
void pkdInitRelaxation(PKD pkd);

int pkdPackIO(PKD pkd,
	      PIO *io, int nMax,
	      local_t *iIndex,
	      total_t iMinOrder, total_t iMaxOrder,
	      double dvFac);

int pkdUnpackIO(PKD pkd,
		PIO *io, int nMax,
		local_t *iIndex,
		total_t iMinOrder, total_t iMaxOrder,
		double dvFac);

#ifdef PLANETS
void pkdSunIndirect(PKD pkd,double aSun[],double adSun[],int iFlag);
void pkdGravSun(PKD pkd,double aSun[],double adSun[],double dSunMass);
void pkdHandSunMass(PKD pkd,double dSunMass);
void pkdReadSS(PKD pkd,char *pszFileName,int nStart,int nLocal);
void pkdWriteSS(PKD pkd,char *pszFileName,int nStart);

#ifdef SYMBA
#define symfac 1.925925925925926 /* 2.08/1.08 */
#define rsym2  1.442307692307692 /* 3.0/2.08 */

void
pkdStepVeryActiveSymba(PKD pkd, double dStep, double dTime, double dDelta,
		       int iRung, int iKickRung, int iRungVeryActive,
		       int iAdjust, double diCrit2,
		       int *pnMaxRung, double dSunMass, int);
int pkdSortVA(PKD pkd, int iRung);
void pkdGravVA(PKD pkd, int iRung);
int pkdCheckDrminVA(PKD pkd, int iRung, int multiflag, int nMaxRung);
void pkdKickVA(PKD pkd, double dt);
int pkdDrminToRung(PKD pkd, int iRung, int iMaxRung, int *nRungCount);
int pkdGetPointa(PKD pkd);
int pkdDrminToRungVA(PKD pkd, int iRung, int iMaxRung, int multiflag);
void pkdMomSun(PKD pkd,double momSun[]);
void pkdDriftSun(PKD pkd,double vSun[],double dt);
void pkdKeplerDrift(PKD pkd,double dt,double mu,int tag_VA);



#endif /* SYMBA */
#endif /* PLANETS*/

#ifdef USE_GRAFIC
void pkdGenerateIC(PKD pkd, GRAFICCTX gctx, int iDim,
		   double fSoft, double fMass, int bCannonical);
#endif
int pkdGetClasses( PKD pkd, int nMax, PARTCLASS *pClass );
void pkdSetClasses( PKD pkd, int n, PARTCLASS *pClass, int bUpdate );

int pkdSelSrcAll(PKD pkd);
int pkdSelDstAll(PKD pkd);
int pkdSelSrcGas(PKD pkd);
int pkdSelDstGas(PKD pkd);
int pkdSelSrcStar(PKD pkd);
int pkdSelDstStar(PKD pkd, int, double);
int pkdSelSrcDeleted(PKD pkd);
int pkdSelDstDeleted(PKD pkd);

int pkdSelSrcMass(PKD pkd,double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse );
int pkdSelDstMass(PKD pkd,double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse );
int pkdSelSrcById(PKD pkd,uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse );
int pkdSelDstById(PKD pkd,uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse );
int pkdSelSrcPhaseDensity(PKD pkd,double dMinDensity, double dMaxDensity, int setIfTrue, int clearIfFalse );
int pkdSelDstPhaseDensity(PKD pkd,double dMinDensity, double dMaxDensity, int setIfTrue, int clearIfFalse );
int pkdSelSrcBox(PKD pkd,double *dCenter, double *dSize, int setIfTrue, int clearIfFalse );
int pkdSelDstBox(PKD pkd,double *dCenter, double *dSize, int setIfTrue, int clearIfFalse );
int pkdSelSrcSphere(PKD pkd,double *r, double dRadius, int setIfTrue, int clearIfFalse );
int pkdSelDstSphere(PKD pkd,double *r, double dRadius, int setIfTrue, int clearIfFalse );
int pkdSelSrcCylinder(PKD pkd,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse );
int pkdSelDstCylinder(PKD pkd,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse );
int pkdDeepestPot(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
    double *r, float *fPot);
void pkdProfile(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
		const double *dCenter, const double *dRadii, int nBins,
		const double *com, const double *vcm, const double *L);
int pkdFindProcessor(const PKD pkd, const FLOAT *R);
void pkdCalcDistance(PKD pkd, double *dCenter);
uint_fast32_t pkdCountDistance(PKD pkd, double r2i, double r2o );
void pkdCalcCOM(PKD pkd, double *dCenter, double dRadius,
		double *com, double *vcm, double *L,
		double *M, uint64_t *N);
void pkdGridInitialize(PKD pkd, int n1, int n2, int n3, int a1, int s, int n);
void pkdGridProject(PKD pkd);
#ifdef MDL_FFTW
void pkdMeasurePk(PKD pkd, double dCenter[3], double dRadius,
		  int nGrid, float *fPower, int *nPower);
#endif

#ifdef USE_CUDA
#ifdef __cplusplus
extern "C" {
#endif
    void *pkdGravCudaPPAllocate(void);
    void pkdGravCudaPPFree(void *cudaCtx);
    ILP_BLK * pkdGravCudaPPAllocateBlk();
    void pkdGravCudaPPFreeBlk(ILP_BLK *blk);
    extern void pkdGravCudaPP(PKD pkd, void *cudaCtx, int nP, PARTICLE **pPart, PINFOIN *pInfoIn, ILP ilp );
#ifdef __cplusplus
}
#endif
#endif

static inline void vec_sub(double *r,const double *a,const double *b ) {
    int i;
    for (i=0; i<3; i++) r[i] = a[i] - b[i];
}

static inline void vec_add_const_mult(double *r,const double *a,double c,const double *b) {
    int i;
    for (i=0; i<3; i++) r[i] = a[i] + c * b[i];
}

static inline void matrix_vector_mult(double *b,double mat[3][3], const double *a) {
    int i,j ;
    for (i=0; i<3; i++) {
        b[i] = 0.0;
        for (j=0; j<3; j++) b[i] += mat[i][j] * a[j];
    }
}

static inline double dot_product(const double *a,const double *b) {
    int i;
    double r = 0.0;
    for(i=0; i<3; i++) r += a[i]*b[i];
    return r;
    }

static inline void cross_product(double *r,const double *a,const double *b) {
    r[0] = a[1] * b[2] - a[2] * b[1] ;
    r[1] = a[2] * b[0] - a[0] * b[2] ;
    r[2] = a[0] * b[1] - a[1] * b[0] ;
}

static inline void mat_transpose(double mat[3][3], double trans_mat[3][3]) {
    int i,j ;
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            trans_mat[i][j] = mat[j][i];
	    }
	}
    }




#endif
