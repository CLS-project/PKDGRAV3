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
#ifdef USE_GRAFIC
#include "grafic.h"
#endif
#include "basetype.h"

#if defined(HAVE_LIBAIO_H)
#include <libaio.h>
#elif defined(HAVE_AIO_H)
#include <aio.h>
#endif

#ifdef __cplusplus
#define CAST(T,V) reinterpret_cast<T>(V)
#else
#define CAST(T,V) ((T)(V))
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
#define CID_PARTICLE2	8
#define CID_CELL2	9
#define CID_HEALPIX     7
#define CID_GROUP	2
#define CID_RM		3
#define CID_BIN		4
#define CID_SHAPES	5
#define CID_PK          2
#define CID_PNG         2
#define CID_SADDLE_BUF  3
#define CID_TREE_ROOT   3

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
#define PKD_MODEL_RELAXATION   (1<<5)  /* Trace relaxation */
#define PKD_MODEL_MASS         (1<<6)  /* Mass for each particle */
#define PKD_MODEL_DENSITY      (1<<7)  /* Density for each particle */
#define PKD_MODEL_BALL         (1<<8)  /* Ball for each particle */
#define PKD_MODEL_SOFTENING    (1<<9)  /* Softening for each particle */
#define PKD_MODEL_VELSMOOTH    (1<<10) /* Velocity Smoothing */
#define PKD_MODEL_SPH          (1<<11) /* Sph Fields */
#define PKD_MODEL_STAR         (1<<12) /* Star Fields */
#define PKD_MODEL_PARTICLE_ID  (1<<13) /* Particles have a unique ID */

#define PKD_MODEL_NODE_MOMENT  (1<<24) /* Include moment in the tree */
#define PKD_MODEL_NODE_ACCEL   (1<<25) /* mean accel on cell (for grav step) */
#define PKD_MODEL_NODE_VEL     (1<<26) /* center of mass velocity for cells */
#define PKD_MODEL_NODE_SPHBNDS (1<<27) /* Include 3 extra bounds in tree */

#define PKD_MODEL_NODE_BND     (1<<28) /* Include normal bounds in tree */
#define PKD_MODEL_NODE_VBND    (1<<29) /* Include velocity bounds in tree for phase-space density*/

#define EPHEMERAL_BYTES 8

typedef struct {
    FLOAT rscale[3];
    FLOAT vscale[3];
    } PSMETRIC;

#define PKD_MAX_CLASSES 256
#define MAX_RUNG     63

/*
** Here we define some special reserved nodes. Node-0 is a sentinel or null node, node-1
** is here defined as the ROOT of the local tree (or top tree), node-2 is unused and
** node-3 is the root node "fixed" tree.
*/
#define FIXROOT         3
#define VAROOT          3
#define ROOT		1
#define NRESERVED_NODES MAX_RUNG+1

typedef struct partclass {
    float       fMass;    /* Particle mass */
    float       fSoft;    /* Current softening */
    FIO_SPECIES eSpecies; /* Species: dark, star, etc. */
    } PARTCLASS;

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


typedef struct partLightCone {
    float pos[3];
    float vel[3];
    } LIGHTCONEP;


typedef struct bndBound {
    double fCenter[3];
    double fMax[3];
    } BND;

#define BND_COMBINE(b,b1,b2)\
{\
	int BND_COMBINE_j;\
	for (BND_COMBINE_j=0;BND_COMBINE_j<3;++BND_COMBINE_j) {\
		FLOAT BND_COMBINE_t1,BND_COMBINE_t2,BND_COMBINE_max,BND_COMBINE_min;\
		BND_COMBINE_t1 = (b1)->fCenter[BND_COMBINE_j] + (b1)->fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2)->fCenter[BND_COMBINE_j] + (b2)->fMax[BND_COMBINE_j];\
		BND_COMBINE_max = (BND_COMBINE_t1 > BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		BND_COMBINE_t1 = (b1)->fCenter[BND_COMBINE_j] - (b1)->fMax[BND_COMBINE_j];\
		BND_COMBINE_t2 = (b2)->fCenter[BND_COMBINE_j] - (b2)->fMax[BND_COMBINE_j];\
		BND_COMBINE_min = (BND_COMBINE_t1 < BND_COMBINE_t2)?BND_COMBINE_t1:BND_COMBINE_t2;\
		(b)->fCenter[BND_COMBINE_j] = 0.5*(BND_COMBINE_max + BND_COMBINE_min);\
		(b)->fMax[BND_COMBINE_j] = 0.5*(BND_COMBINE_max - BND_COMBINE_min);\
		}\
	}

static inline int IN_BND(const FLOAT *R,const BND *b) {
    int i;
    for( i=0; i<3; i++ )
	if ( R[i]<b->fCenter[i]-b->fMax[i] || R[i]>=b->fCenter[i]+b->fMax[i] )
	    return 0;
    return 1;
    }


#if defined(USE_SIMD) && defined(__SSE2__)
static inline FLOAT mindist(const BND *bnd,const FLOAT *pos) {
#ifdef __AVX__
    typedef union {
	uint64_t i[4];
	__m256d p;
	} vint64;
    static const vint64 isignmask = {{0x8000000000000000,0x8000000000000000,0x8000000000000000,0x8000000000000000}};
    static const vint64 zero = {{0,0,0,0}};
    __m256i justthree = {-1,-1,-1,0};
    __m256d fCenter = _mm256_maskload_pd(bnd->fCenter,justthree);
    __m256d fMax = _mm256_maskload_pd(bnd->fMax,justthree);
    __m256d ppos = _mm256_maskload_pd(pos,justthree);
    __m128d d;
    __m256d m = _mm256_max_pd(zero.p,_mm256_sub_pd(_mm256_andnot_pd(isignmask.p,_mm256_sub_pd(fCenter,ppos)),fMax));
    m = _mm256_mul_pd(m,m);
    d = _mm_hadd_pd(_mm256_extractf128_pd(m,1),_mm256_castpd256_pd128(m));
    return _mm_cvtsd_f64(_mm_hadd_pd(d,d));
#else
    typedef union {
	uint64_t i[2];
	__m128d p;
	} vint64;
    static const vint64 isignmask = {{0x8000000000000000,0x8000000000000000}};
    static const vint64 zero = {{0,0}};
    __m128d m2,m1;
    m1 = _mm_max_sd(zero.p,_mm_sub_sd(_mm_andnot_pd(isignmask.p,_mm_sub_sd(_mm_load_sd(bnd->fCenter+0),_mm_load_sd(pos+0))),_mm_load_sd(bnd->fMax+0)));
    m2 = _mm_mul_sd(m1,m1);
    m1 = _mm_max_sd(zero.p,_mm_sub_sd(_mm_andnot_pd(isignmask.p,_mm_sub_sd(_mm_load_sd(bnd->fCenter+1),_mm_load_sd(pos+1))),_mm_load_sd(bnd->fMax+1)));
    m2 = _mm_add_sd(m2,_mm_mul_sd(m1,m1));
    m1 = _mm_max_sd(zero.p,_mm_sub_sd(_mm_andnot_pd(isignmask.p,_mm_sub_sd(_mm_load_sd(bnd->fCenter+2),_mm_load_sd(pos+2))),_mm_load_sd(bnd->fMax+2)));
    m2 = _mm_add_sd(m2,_mm_mul_sd(m1,m1));
    return _mm_cvtsd_f64(m2);
#endif
    }
#define MINDIST(bnd,pos,min2) ((min2) = mindist(bnd,pos))
#else
#define MINDIST(bnd,pos,min2) {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;					\
    for (BND_j=0;BND_j<3;++BND_j) {\
	BND_dMin = fabs((bnd)->fCenter[BND_j] - (pos)[BND_j]) - (bnd)->fMax[BND_j]; \
	if (BND_dMin > 0) (min2) += BND_dMin*BND_dMin;			\
	}\
    }
#endif

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
#define STACK_SIZE(S) (S ## __s)

/* This macro give the id (processor) and index of the two child cells */
#define pkdGetChildCells(c,id,idLower,idxLower,idUpper,idxUpper)	\
    do {								\
        idxLower = c->iLower;   /* This is always true */		\
	idLower = idUpper = id; /* Default */				\
	if (c->bTopTree) {						\
	    /* We may point off node, but then nodes are adjacent */	\
	    if (c->bRemote) {						\
		idLower = c->pLower;					\
		idUpper = idLower;					\
		idxUpper = idxLower+1;					\
		}							\
	    else {							\
		idxUpper = c->pUpper;					\
		}							\
	    }								\
	else { idxUpper = idxLower+1; }					\
	} while(0)							\

typedef struct kdNode {
    double r[3];
    int pLower;		     /* also serves as thread id for the LTT */
    int pUpper;		     /* pUpper < 0 indicates no particles in tree! */
    uint32_t iLower;         /* Local lower node (or remote processor w/bRemote=1) */
    uint16_t iDepth;
    uint16_t uMinRung   : 6;
    uint16_t uMaxRung   : 6;
    uint16_t bSrcActive : 1;
    uint16_t bDstActive : 1;
    uint16_t bTopTree   : 1; /* This is a top tree node: pLower,pUpper are node indexes */
    uint16_t bRemote    : 1; /* children are remote */
    float bMax;
    float fSoft2;
    } KDN;

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
	BND *CALCOPEN_bnd = pkdNodeBnd(pkd, pkdn);			\
        for (CALCOPEN_j=0;CALCOPEN_j<3;++CALCOPEN_j) {                  \
            FLOAT CALCOPEN_d = fabs(CALCOPEN_bnd->fCenter[CALCOPEN_j] - (pkdn)->r[CALCOPEN_j]) + \
                CALCOPEN_bnd->fMax[CALCOPEN_j];                          \
            CALCOPEN_d2 += CALCOPEN_d*CALCOPEN_d;                       \
            }								\
	MAXSIDE(CALCOPEN_bnd->fMax,CALCOPEN_b);				\
	if (CALCOPEN_b < minside) CALCOPEN_b = minside;			\
	if (CALCOPEN_b*CALCOPEN_b < CALCOPEN_d2) CALCOPEN_b = sqrt(CALCOPEN_d2); \
	(pkdn)->bMax = CALCOPEN_b;					\
	}

#if (0)
#define CALCOPEN(pkdn) {						\
        FLOAT CALCOPEN_d2 = 0;						\
        int CALCOPEN_j;							\
	BND *CALCOPEN_bnd = pkdNodeBnd(pkd, pkdn);			\
        for (CALCOPEN_j=0;CALCOPEN_j<3;++CALCOPEN_j) {                  \
            FLOAT CALCOPEN_d = fabs(CALCOPEN_bnd->fCenter[CALCOPEN_j] - (pkdn)->r[CALCOPEN_j]) + \
                CALCOPEN_bnd->fMax[CALCOPEN_j];                          \
            CALCOPEN_d2 += CALCOPEN_d*CALCOPEN_d;                       \
            }\
        CALCOPEN_d2 = sqrt(CALCOPEN_d2);	  \
        if (CALCOPEN_d2 < (pkdn)->bMax) (pkdn)->bMax = CALCOPEN_d2;	  \
	}
#endif

/*
** Components required for tree walking.
*/
typedef struct CheckStack {
    ILPCHECKPT PartChkPt;
    ILCCHECKPT CellChkPt;
    CL cl;
    LOCR L;
    float dirLsum;
    float normLsum;
    int iNodeIndex;
    } CSTACK;

/*
** components required for time-step calculation (only grav.c)
*/

typedef struct RhoLocalArray {
    double d2;
    double m;
    } RHOLOCAL;

typedef union {
    double *d;
#if defined(USE_SIMD) && !defined(__CUDACC__)
    v_df *p;
#endif
    uint64_t *i;
    } ewaldDouble;

#if defined(USE_SIMD) && !defined(__CUDACC__)
typedef struct {
    struct {
	vdouble m;
	vdouble xx,yy,xy,xz,yz;
	vdouble xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	vdouble xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	vdouble zz;
	vdouble xzz,yzz,zzz;
	vdouble xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} ewm;
    struct {
	vdouble fEwCut2,fInner2,alpha,alpha2,ialpha,k1,ka;
	vdouble Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z,Q2;
	} ewp;
    } ewaldSIMD;
#endif

/*
** This is the temporary group table used when Grasshopping.
** We eventually contruct a proper table.
*/
typedef remoteID GHtmpGroupTable;

typedef struct {
    remoteID key;
    remoteID name;
    uint32_t iLink;
    } FOFRemote;

typedef struct {
    uint64_t nTotal;      /* Total particles in this group */
    remoteID id;          /* Owner (or myself) */
    remoteID rmt;
    uint32_t iGlobalId;   /* Global unique group id */
    uint32_t nLocal;      /* Local to this processor */
    uint32_t iTreeRoot;   /* Our local tree root */
    uint32_t iAllRoots;
    uint16_t nRemote;     /* Number of remote partners */
    uint8_t  bNeedGrav : 1;
    uint8_t  bComplete : 1;
    float fMass;
    float fRMSRadius;
    double dEnergy;
    double rref[3];
    double ravg[3];
    double rmin[3];
    double rmax[3];
    double rcom[3];
    double vcom[3];
    } HopGroupTable;

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
    float fMass;
    float fRMSRadius;
    double r[3];
    double rcom[3];
    float potordenmax;
    float v[3];
    int nLocal;
    int nTotal;
    int bMyGroup;
    int nRemoteMembers;
    int iFirstRm;
} FOFGD;

struct remote_root_id
{
    int iPid;
    int iLocalGroupId;
    int iLocalRootId;
};

struct saddle_point_group
{
    int iGlobalId;
    int iLocalId;
    int iPid;
    FLOAT fDensity;
};

struct saddle_point
{
    /* Information about the particle that is the saddle point */
    FLOAT fDensity;
    int iLocalId;
    int iPid;

    /* The group that owns the saddle point */
    struct saddle_point_group owner;
    /* The other group joined by the saddle point */
    struct saddle_point_group nbr;
    struct saddle_point_group parent;
};

struct saddle_point_buffer
{
    /* Destination Group Id. I.e., who should get the following saddle points. */
    int iLocalId;
    int iGlobalId;

    int nSaddlePoints;
    struct saddle_point sp[32];
};

struct saddle_point_list
{
    /* Number of saddle points in the list */
    int n;
    /* Size of allocated array */
    int size;

    struct saddle_point *sp;
    struct saddle_point_buffer *buf;
};

struct tree_root_cache_msg
{
    struct remote_root_id tr;
    int iPid;
    int iLocalId;
};

struct psGroup {
    int iPid;
    int iLocalId;
    int iGlobalId;
    /*-----Unique-Local-Data-----*/
    int bridge;
    int dup;
    int nLocal;
    int nTreeRoots;
    struct remote_root_id *treeRoots;
    int nSaddlePoints;
    int *sp;
    /*-----Shared-Data-----------*/
    uint64_t nTotal;
    FLOAT fDensity;
    FLOAT fMass;
    FLOAT fRMSRadius;
    FLOAT r[3], rcom[3];
    FLOAT v[3], vcom[3];
    FLOAT fMass_com;
};

struct psGroupTable
{
    int nGroups;
    struct psGroup *pGroup;
};

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

/*#define PKD_GROUP_SIZE 256*/

typedef struct {
    uint64_t nActiveBelow;
    uint64_t nActiveAbove;
    uint64_t nTotalBelow;
    uint64_t nTotalAbove;
    } ORBCOUNT;

struct psContext;

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
    int nTreeTilesReserved;
    int nMaxNodes;
    uint64_t nRung[IRUNGMAX+1];
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    FLOAT fPeriod[3];
    char **kdNodeListPRIVATE; /* BEWARE: also char instead of KDN */
    int iTopTree[NRESERVED_NODES];
    int nNodes;
    int nNodesFull;     /* number of nodes in the full tree (including very active particles) */
    int nNonVANodes;    /* number of nodes *not* in Very Active Tree, or index to the start of the VA nodes (except VAROOT) */
    BND bnd;
    BND vbnd;
    BND fixbnd;
    size_t iTreeNodeSize;
    size_t iParticleSize;
    PARTICLE *pStorePRIVATE;
    PARTICLE *pTempPRIVATE;
    double dTimeRedshift0;

#define NUMLCBUFS 2
#if defined(HAVE_LIBAIO_H)
    struct iocb cbLightCone[NUMLCBUFS];
    struct io_event eventsLightCone[NUMLCBUFS];
    io_context_t ctxLightCone;
#elif defined(HAVE_AIO_H)
    struct aiocb cbLightCone[NUMLCBUFS];
    struct aiocb const * pcbLightCone[NUMLCBUFS];
#endif
    LIGHTCONEP *pLightCone[NUMLCBUFS];
    off_t iFilePositionLightCone;
    int fdLightCone;
    int iLightConeBuffer;
    int nLightCone, nLightConeMax;
    int64_t nHealpixPerDomain;
    int64_t nSideHealpix;
    uint32_t *pHealpixCounts;
    PARTCLASS *pClass;
    float fSoftFix;
    float fSoftFac;
    float fSoftMax;
    int nClasses;
    void *pLite;
#ifdef COOLING
    COOL *Cool; /* Cooling Context */
#endif
    /*
    ** Advanced memory models
    */
    int oPosition;
    int oAcceleration; /* Three float */
    int oVelocity; /* Three vel_t */
    int oPotential; /* One float */
    int oGroup; /* One int32 */
    int oMass; /* One float */
    int oSoft; /* One float */
    int oDensity; /* One float */
    int oBall; /* One float */
    int oSph; /* Sph structure */
    int oStar; /* Star structure */
    int oRelaxation;
    int oVelSmooth;
    int oRungDest; /* Destination processor for each rung */
    int oParticleID;

    /*
    ** Advanced memory models - Tree Nodes
    */
    int oNodePosition; /* Three vel_t */
    int oNodeVelocity; /* Three vel_t */
    int oNodeAcceleration; /* Three doubles */
    int oNodeBmax;
    int oNodeSoft;
    int oNodeMom; /* an FMOMR */
    int oNodeBnd;
    int oNodeSphBounds; /* Three Bounds */
    int oNodeVBnd; /* Velocity bounds */

    /*
    ** Tree walk variables.
    */
    int nMaxStack;
    CSTACK *S;
    ILP ilp;
    ILC ilc;
    LSTFREELIST clFreeList;
    CL cl;
    CL clNew;
    double dFlop;

    /*
    ** Opening angle table for mass weighting.
    */
    float fiCritTheta;
    /* Potential Energy for when potential is not in the particle */
    double dEnergyU;
    /* also put kinetic energy here to calculate it on the fly */
    double dEnergyT;
    double dEnergyW;
    double dEnergyF[3];
    double dEnergyL[3];

    /*
    ** New activation methods
    */
    uint8_t uMinRungActive;
    uint8_t uMaxRungActive;
    uint8_t uRungVeryActive;    /* NOTE: The first very active particle is at iRungVeryActive + 1 */

    /*
    ** Ewald summation setup.
    */
    struct EwaldVariables ew;
    EwaldTable ewt;
    ewaldSIMD es;
    workEwald *ewWork;

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
    struct psGroupTable psGroupTable;

    FOFGD *groupData;



    uint64_t iStartGID;
    int nGroups, nLocalGroups;
    /*
    ** Some variables needed for pkdNewFof().
    */
    uint32_t iHead;
    uint32_t iTail;
    uint32_t  *Fifo;
    int bCurrGroupContained;
    uint32_t nCurrFofParticles;
    FLOAT fMinFofContained[3];
    FLOAT fMaxFofContained[3];    
    struct smGroupArray *ga;
    uint32_t iRemoteGroup,nMaxRemoteGroups;
    FOFRemote *tmpFofRemote;
    

    GHtmpGroupTable *tmpHopGroups;
    HopGroupTable *hopGroups;
    int *hopRootIndex;
    int hopSavedRoots;
    remoteID *hopRoots;

    struct saddle_point_list saddle_points;
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

#ifdef USE_CUDA
    void *cudaCtx;
#endif

    MDLGRID grid;
    float *gridData;
    struct psContext *psx;

    } * PKD;


#ifdef __SSE2__
#define pkdMinMax(dVal,dMin,dMax) do {					\
    (dMin)[0] = _mm_cvtsd_f64(_mm_min_sd(_mm_set_sd((dMin)[0]),_mm_set_sd((dVal)[0]))); \
    (dMin)[1] = _mm_cvtsd_f64(_mm_min_sd(_mm_set_sd((dMin)[1]),_mm_set_sd((dVal)[1])));	\
    (dMin)[2] = _mm_cvtsd_f64(_mm_min_sd(_mm_set_sd((dMin)[2]),_mm_set_sd((dVal)[2])));	\
    (dMax)[0] = _mm_cvtsd_f64(_mm_max_sd(_mm_set_sd((dMax)[0]),_mm_set_sd((dVal)[0])));	\
    (dMax)[1] = _mm_cvtsd_f64(_mm_max_sd(_mm_set_sd((dMax)[1]),_mm_set_sd((dVal)[1])));	\
    (dMax)[2] = _mm_cvtsd_f64(_mm_max_sd(_mm_set_sd((dMax)[2]),_mm_set_sd((dVal)[2]))); \
    } while(0)
#else
#define pkdMinMax(dVal,dMin,dMax) {\
    (dMin)[0] = (dVal)[0] < (dMin)[0] ? (dVal)[0] : (dMin)[0];	\
    (dMin)[1] = (dVal)[1] < (dMin)[1] ? (dVal)[1] : (dMin)[1];		\
    (dMin)[2] = (dVal)[2] < (dMin)[2] ? (dVal)[2] : (dMin)[2];		\
    (dMax)[0] = (dVal)[0] > (dMax)[0] ? (dVal)[0] : (dMax)[0];		\
    (dMax)[1] = (dVal)[1] > (dMax)[1] ? (dVal)[1] : (dMax)[1];		\
    (dMax)[2] = (dVal)[2] > (dMax)[2] ? (dVal)[2] : (dMax)[2];		\
    }
#endif
#define pkdMinMax6(dVal0,dVal1,dMin,dMax) {\
    (dMin)[0] = (dVal0)[0] < (dMin)[0] ? (dVal0)[0] : (dMin)[0];	\
    (dMin)[1] = (dVal0)[1] < (dMin)[1] ? (dVal0)[1] : (dMin)[1];	\
    (dMin)[2] = (dVal0)[2] < (dMin)[2] ? (dVal0)[2] : (dMin)[2];	\
    (dMin)[3] = (dVal1)[0] < (dMin)[3] ? (dVal1)[0] : (dMin)[3];	\
    (dMin)[4] = (dVal1)[1] < (dMin)[4] ? (dVal1)[1] : (dMin)[4];	\
    (dMin)[5] = (dVal1)[2] < (dMin)[5] ? (dVal1)[2] : (dMin)[5];	\
    (dMax)[0] = (dVal0)[0] > (dMax)[0] ? (dVal0)[0] : (dMax)[0];	\
    (dMax)[1] = (dVal0)[1] > (dMax)[1] ? (dVal0)[1] : (dMax)[1];	\
    (dMax)[2] = (dVal0)[2] > (dMax)[2] ? (dVal0)[2] : (dMax)[2];	\
    (dMax)[3] = (dVal1)[0] > (dMax)[3] ? (dVal1)[0] : (dMax)[3];	\
    (dMax)[4] = (dVal1)[1] > (dMax)[4] ? (dVal1)[1] : (dMax)[4];	\
    (dMax)[5] = (dVal1)[2] > (dMax)[5] ? (dVal1)[2] : (dMax)[5];	\
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
    return sizeof(KDN) + 2*sizeof(BND) + sizeof(FMOMR) + 6*sizeof(double) + sizeof(SPHBNDS);
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
static inline vel_t *pkdNodeVel( PKD pkd, KDN *n ) {
    return CAST(vel_t *,pkdNodeField(n,pkd->oNodeVelocity));
    }
static inline float *pkdNodeAccel( PKD pkd, KDN *n ) {
    return CAST(float *,pkdNodeField(n,pkd->oNodeAcceleration));
    }
static inline SPHBNDS *pkdNodeSphBounds( PKD pkd, KDN *n ) {
    return CAST(SPHBNDS *,pkdNodeField(n,pkd->oNodeSphBounds));
    }

static inline BND *pkdNodeBnd( PKD pkd, KDN *n ) {
    return CAST(BND *,pkdNodeField(n,pkd->oNodeBnd));
    }

static inline BND *pkdNodeVBnd( PKD pkd, KDN *n ) {
    return CAST(BND *,pkdNodeField(n,pkd->oNodeVBnd));
    }

static inline KDN *pkdNode(PKD pkd,KDN *pBase,int iNode) {
    return (KDN *)&((char *)pBase)[pkd->iTreeNodeSize*iNode];
    }
int pkdNodes(PKD pkd);
void pkdExtendTree(PKD pkd);
static inline KDN *pkdTreeNode(PKD pkd,int iNode) {
    char *kdn = &pkd->kdNodeListPRIVATE[(iNode>>pkd->nTreeBitsLo)][pkd->iTreeNodeSize*(iNode&pkd->iTreeMask)];
    return (KDN *)kdn;
    }
static inline int pkdTreeAlignNode(PKD pkd) {
    if (pkd->nNodes&1) ++pkd->nNodes;
    return pkd->nNodes;
    }

static inline int pkdTreeAllocNodes(PKD pkd, int nNodes) {
    int iNode = pkd->nNodes;
    pkd->nNodes += nNodes;
    while(pkd->nNodes > pkd->nMaxNodes) pkdExtendTree(pkd);
    return iNode;
    }
static inline int pkdTreeAllocNode(PKD pkd) {
    return pkdTreeAllocNodes(pkd,1);
    }
static inline void pkdTreeAllocNodePair(PKD pkd,int *iLeft, int *iRight) {
    *iLeft = pkdTreeAllocNode(pkd);
    *iRight = pkdTreeAllocNode(pkd);
    }
static inline void pkdTreeAllocRootNode(PKD pkd,int *iRoot) {
    pkdTreeAllocNodePair(pkd,NULL,iRoot);
    }

void *pkdTreeNodeGetElement(void *vData,int i,int iDataSize);
static inline KDN *pkdTopNode(PKD pkd,int iNode) {
    assert(0); // no top tree now
    }
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
    return (pkd->iParticleSize + EPHEMERAL_BYTES) * (pkd->nStore+1);
    }
static inline PARTICLE *pkdParticleGet( PKD pkd, void *pBase, int i ) {
    char *v = (char *)pBase;
    PARTICLE *p = (PARTICLE *)(v + ((uint64_t)i)*pkd->iParticleSize);
    return p;
    }
static inline PARTICLE *pkdParticle( PKD pkd, int i ) {
    return pkdParticleGet(pkd,pkd->pStorePRIVATE,i);
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

static inline int32_t *pkdGroup( PKD pkd, PARTICLE *p ) {
    assert(pkd->oGroup);
    return CAST(int32_t *, pkdField(p,pkd->oGroup));
    }

static inline float pkdDensity( PKD pkd, PARTICLE *p ) {
    assert(pkd->oDensity);
    return * CAST(float *, pkdField(p,pkd->oDensity));
    }
static inline void pkdSetDensity( PKD pkd, PARTICLE *p, float fDensity ) {
    if (pkd->oDensity) *CAST(float *, pkdField(p,pkd->oDensity)) = fDensity;
    }

static inline float pkdBall( PKD pkd, PARTICLE *p ) {
    assert(pkd->oBall);
    return *CAST(float *, pkdField(p,pkd->oBall));
    }
static inline void pkdSetBall(PKD pkd, PARTICLE *p, float fBall) {
    if (pkd->oBall) *CAST(float *, pkdField(p,pkd->oBall)) = fBall;
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


/*
** Integerized coordinates: signed integer -0x7fffffff to +0x7fffffff
** We assume a periodic box of width 1 so a simple multiple will convert.
** The situation is more complicated with non-periodic boxes, or for boxes
** with a different period so this is not currently supported.
*/
#define pkdPosRaw(pkd,p,d) (CAST(pos_t *,pkdField(p,pkd->oPosition))[d])
#define pkdSetPosRaw(pkd,p,d,v) (CAST(pos_t *,pkdField(p,pkd->oPosition))[d]) = (v);
#ifdef INTEGER_POSITION
#define pkdDblToPos(pkd,d) (pos_t)((d)*0x80000000u)
#define pkdPos(pkd,p,d) ((CAST(pos_t *,pkdField(p,pkd->oPosition))[d]) * (1.0/0x80000000u))
#define pkdSetPos(pkd,p,d,v) (void)((CAST(pos_t *,pkdField(p,pkd->oPosition))[d]) = (v)*0x80000000u)
#ifdef __AVX__
#define pkdGetPos3(pkd,p,d1,d2,d3) do {					\
	union { __m256d p; double d[4]; } r_pkdGetPos3;			\
	r_pkdGetPos3.p = _mm256_mul_pd(_mm256_cvtepi32_pd(*(__m128i *)(CAST(pos_t *,pkdField(p,pkd->oPosition)))),_mm256_set1_pd(1.0/0x80000000u) ); \
	d1 = r_pkdGetPos3.d[0];						\
	d2 = r_pkdGetPos3.d[1];						\
	d3 = r_pkdGetPos3.d[2];						\
	} while(0)
#else
#define pkdGetPos3(pkd,p,d1,d2,d3) do { d1=pkdPos(pkd,p,0); d2=pkdPos(pkd,p,1); d3=pkdPos(pkd,p,2); } while(0)
#endif
#else
#define pkdDblToPos(pkd,d) (d)
#define pkdPos(pkd,p,d) (CAST(pos_t *,pkdField(p,pkd->oPosition))[d])
#define pkdSetPos(pkd,p,d,v) (void)((CAST(pos_t *,pkdField(p,pkd->oPosition))[d]) = (v))
#define pkdGetPos3(pkd,p,d1,d2,d3) ((d1)=pkdPos(pkd,p,0),(d2)=pkdPos(pkd,p,1),(d3)=pkdPos(pkd,p,2))
#endif
#define pkdGetPos1(pkd,p,d) pkdGetPos3(pkd,p,(d)[0],(d)[1],(d)[2])

static inline vel_t *pkdVel( PKD pkd, PARTICLE *p ) {
    return CAST(vel_t *,pkdField(p,pkd->oVelocity));
    }
static inline float *pkdAccel( PKD pkd, PARTICLE *p ) {
    return CAST(float *,pkdField(p,pkd->oAcceleration));
    }
static inline float *pkdPot( PKD pkd, PARTICLE *p ) {
    return pkd->oPotential ? CAST(float *,pkdField(p,pkd->oPotential)) : NULL;
    }
static inline uint16_t *pkdRungDest( PKD pkd, PARTICLE *p ) {
    return CAST(uint16_t *,pkdField(p,pkd->oRungDest));
    }
static inline uint64_t *pkdParticleID( PKD pkd, PARTICLE *p ) {
    return CAST(uint64_t *,pkdField(p,pkd->oParticleID));
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

static inline float *pkd_Timer( PKD pkd, PARTICLE *p ) {
    return &(((STARFIELDS *) pkdField(p,pkd->oStar))->fTimer);
    }

static inline int pkdIsDeleted(PKD pkd,PARTICLE *p) {
    return (pkdSpecies(pkd,p) == FIO_SPECIES_LAST);
    }

static inline int pkdIsNew(PKD pkd,PARTICLE *p) {
    return (p->iOrder == IORDERMAX);
    }

/*
** From tree.c:
*/
void pkdVATreeBuild(PKD pkd,int nBucket);
void pkdTreeBuild(PKD pkd,int nBucket,uint32_t uRoot,uint32_t uTemp);
uint32_t pkdDistribTopTree(PKD pkd, uint32_t uRoot, uint32_t nTop, KDN *pTop);
void pkdOpenCloseCaches(PKD pkd,int bOpen,int bFixed);
void pkdTreeInitMarked(PKD pkd);
void pkdDumpTrees(PKD pkd,int bOnlyVA,uint8_t uRungDD);
void pkdCombineCells1(PKD,KDN *pkdn,KDN *p1,KDN *p2);
void pkdCombineCells2(PKD,KDN *pkdn,KDN *p1,KDN *p2);
void pkdCalcRoot(PKD,uint32_t,double *,MOMC *);
void pkdDistribRoot(PKD,double *,MOMC *);
void pkdTreeBuildByGroup(PKD pkd, int nBucket);

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
    PKD *ppkd,MDL mdl,int nStore,uint64_t nMinTotalStore,uint64_t nMinEphemeral,
    int nBucket,int nGroup,int nTreeBitsLo, int nTreeBitsHi,
    int iCacheSize,int iWorkQueueSize,int iCUDAQueueSize,FLOAT *fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,
    uint64_t mMemoryModel, int bLightCone, int bLightConeParticles);
void pkdFinish(PKD);
size_t pkdClCount(PKD pkd);
size_t pkdClMemory(PKD pkd);
size_t pkdIlcMemory(PKD pkd);
size_t pkdIlpMemory(PKD pkd);
size_t pkdTreeMemory(PKD pkd);
void pkdReadFIO(PKD pkd,FIO fio,uint64_t iFirst,int nLocal,double dvFac, double dTuFac);
void pkdSetSoft(PKD pkd,double dSoft);
void pkdSetCrit(PKD pkd,double dCrit);
void pkdCalcBound(PKD,BND *);
void pkdCalcVBound(PKD,BND *);
void pkdEnforcePeriodic(PKD,BND *);
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul);
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

void pkdOrbBegin(PKD pkd, int nRungs);
int pkdOrbSelectRung(PKD pkd, int iRung);
void pkdOrbUpdateRung(PKD pkd);
void pkdOrbFinish(PKD pkd);
void pkdOrbSplit(PKD pkd,int iDomain);
int pkdOrbRootFind(
    PKD pkd,double dFraction,uint64_t nLowerMax, uint64_t nUpperMax,
    double dReserveFraction, BND *bnd, double *dSplitOut, int *iDim);
void pkdOrbUpdateRung(PKD pkd);
/*#define PEANO_HILBERT_KEY_MAX 0x3ffffffffffull*/ /* 2d */
#define PEANO_HILBERT_KEY_MAX 0x7fffffffffffffffull /* 3d */
void pkdPeanoHilbertDecomp(PKD pkd, int nRungs, int iMethod);
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
void pkdLocalOrder(PKD,uint64_t iMinOrder,uint64_t iMaxOrder);
void pkdCheckpoint(PKD pkd,const char *fname);
void pkdRestore(PKD pkd,const char *fname);
uint32_t pkdWriteFIO(PKD pkd,FIO fio,double dvFac,BND *bnd);
void pkdWriteFromNode(PKD pkd,int iNode, FIO fio,double dvFac,BND *bnd);
void pkdWriteViaNode(PKD pkd, int iNode);
void pkdGravAll(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
    int bKickClose,int bKickOpen,double *dtClose,double *dtOpen,
    double dAccFac,double dTime,int nReps,int bPeriodic,
    int iOrder,int bEwald,int nGroup,int iRoot1, int iRoot2,
    double fEwCut,double fEwhCut,double dThetaMin,
    uint64_t *pnActive,
    double *pdPart,double *pdPartNumAccess,double *pdPartMissRatio,
    double *pdCell,double *pdCellNumAccess,double *pdCellMissRatio,
    double *pdFlop,uint64_t *pnRung);
void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double *L,double *F,double *W);
void pkdDrift(PKD pkd,int iRoot,double dTime,double dDelta,double,double);
void pkdScaleVel(PKD pkd,double dvFac);
void pkdStepVeryActiveKDK(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dStep, double dTime, double dDelta,
			  int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double diCrit2,
			  int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass);
void pkdKickKDKOpen(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
void pkdKickKDKClose(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
void pkdKick(PKD pkd,double dTime,double dDelta,double,double,double,uint8_t uRungLo,uint8_t uRungHi);
void pkdKickTree(PKD pkd,double dTime,double dDelta,double,double,double,int iRoot);
void pkdSwapAll(PKD pkd, int idSwap);
void pkdInitStep(PKD pkd,struct parameters *p,CSM csm);
void pkdSetRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung);
void pkdZeroNewRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung);
void pkdActiveRung(PKD pkd, int iRung, int bGreater);
void pkdCountRungs(PKD pkd,uint64_t *nRungs);
void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
		  double dEta,double dVelFac,double dAccFac,
		  int bDoGravity,int bEpsAcc,double dhMinOverSoft);
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
		  uint8_t uRung,int iMaxRung,uint64_t *nRungCount);
void pkdUpdateRungByTree(PKD pkd,int iRoot,uint8_t uMinRung,int iMaxRung,uint64_t *nRungCount);
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
    total_t n;
    total_t nGas;
    total_t nDark;
    total_t nStar;
    total_t nMaxOrder;
    };

void pkdGetNParts(PKD pkd, struct outGetNParts *out );
void pkdSetNParts(PKD pkd, int nGas, int nDark, int nStar);
void pkdInitRelaxation(PKD pkd);

#ifdef USE_GRAFIC
void pkdGenerateIC(PKD pkd, GRAFICCTX gctx, int iDim,
		   double fSoft, double fMass, int bCannonical);
#endif
int pkdGetClasses( PKD pkd, int nMax, PARTCLASS *pClass );
void pkdSetClasses( PKD pkd, int n, PARTCLASS *pClass, int bUpdate );
void pkdSetClass( PKD pkd, float fMass, float fSoft, FIO_SPECIES eSpecies, PARTICLE *p );

int pkdSelSrcAll(PKD pkd);
int pkdSelDstAll(PKD pkd);
int pkdSelSrcGas(PKD pkd);
int pkdSelDstGas(PKD pkd);
int pkdSelSrcStar(PKD pkd);
int pkdSelDstStar(PKD pkd, int, double);
int pkdSelSrcDeleted(PKD pkd);
int pkdSelDstDeleted(PKD pkd);
int pkdSelSrcGroup(PKD pkd, int iGroup);
int pkdSelDstGroup(PKD pkd, int iGroup);

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
void pkdCalcDistance(PKD pkd, double *dCenter);
uint_fast32_t pkdCountDistance(PKD pkd, double r2i, double r2o );
void pkdCalcCOM(PKD pkd, double *dCenter, double dRadius,
		double *com, double *vcm, double *L,
		double *M, uint64_t *N);
void pkdGridInitialize(PKD pkd, int n1, int n2, int n3, int a1, int s, int n);
void pkdGridProject(PKD pkd);
#ifdef MDL_FFTW
void pkdMeasurePk(PKD pkd, double dCenter[3], double dRadius, double dTotalMass,
    int nGrid, int nBins, float *fK, float *fPower, int *nPower);
#endif
void pkdOutPsGroup(PKD pkd,char *pszFileName,int iType);

void pkdLightConeOpen(PKD pkd, const char *fname,int nSideHealpix);
void pkdLightConeClose(PKD pkd, const char *healpixname);


#ifdef USE_CUDA
#ifdef __cplusplus
extern "C" {
#endif
    extern int CUDAinitWorkPP( void *vpp, void *vwork );
    extern int CUDAcheckWorkPP( void *vpp, void *vwork );
    extern int CUDAinitWorkPC( void *vpp, void *vwork );
    extern int CUDAcheckWorkPC( void *vpp, void *vwork );
    extern int CUDAinitWorkEwald( void *vpp, void *vwork );
    extern int CUDAcheckWorkEwald( void *vpp, void *vwork );
    extern void cudaEwaldInit(struct EwaldVariables *ewIn, EwaldTable *ewt );
#ifdef __cplusplus
}
#endif
#endif

#define vec_sub(r,a,b) do {\
    int i;\
    for (i=0; i<3; i++) (r)[i] = (a)[i] - (b)[i];	\
} while(0)

#define vec_add_const_mult(r,a,c,b) do {\
    int i;\
    for (i=0; i<3; i++) (r)[i] = (a)[i] + (c) * (b)[i];	\
} while(0)

#define matrix_vector_mult(b,mat,a) do {\
    int i;\
    for (i=0; i<3; i++) {\
        int j;\
	(b)[i] = 0.0;					\
        for (j=0; j<3; j++) (b)[i] += (mat)[i][j] * (a)[j];	\
    }\
} while(0)

static inline double dot_product(const double *a,const double *b) {
    double r = 0.0;
    int i;
    for(i=0; i<3; i++) r += a[i]*b[i];
    return r;
    }

#define cross_product(r,a,b) do {\
    (r)[0] = (a)[1] * (b)[2] - (a)[2] * (b)[1] ;	\
    (r)[1] = (a)[2] * (b)[0] - (a)[0] * (b)[2] ;	\
    (r)[2] = (a)[0] * (b)[1] - (a)[1] * (b)[0] ;	\
} while(0)

#define mat_transpose(mat,trans_mat) do {\
    int i;				 \
    for (i=0; i<3; i++) {			\
	int j;					\
        for (j=0; j<3; j++) {			\
            (trans_mat)[i][j] = (mat)[j][i];	\
	    }					\
	}					\
} while(0)


#endif
