#ifndef ILP_H
#define ILP_H
#ifdef LOCAL_EXPANSION
#include <stdint.h>

#ifndef ILP_PART_PER_TILE
#define ILP_PART_PER_TILE 4096 /* 4096*24 ~ 100k */
#endif

#define ILP_ALIGN_BITS 2
#define ILP_ALIGN_SIZE (1<<ILP_ALIGN_BITS)
#define ILP_ALIGN_MASK (ILP_ALIGN_SIZE-1)

#include "simd.h"

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILP_PART_PER_TILE];
#ifdef USE_SIMD_PP
    v4sf p[ILP_PART_PER_TILE/4];
#endif
    } ilpFloat;

typedef union {
    uint64_t i[ILP_PART_PER_TILE];
    } ilpInt64;


typedef struct ilpTile {
    struct ilpTile *next;
    struct ilpTile *prev;
    uint32_t nMaxPart;          /* Maximum number of particles in this tile */
    uint32_t nPart;             /* Current number of particles */

    struct {
	ilpFloat dx, dy, dz;        /* Offset from ilp->cx, cy, cz */
	ilpFloat d2;                /* Distance squared: calculated */
	ilpFloat m;                 /* Mass */
	ilpFloat fourh2;            /* Softening: calculated */
#ifdef HERMITE
	ilpFloat vx, vy, vz;
#endif
#if defined(SYMBA) || defined(PLANETS)
	ilpInt64 iOrder;
#endif
	} d;
    /* Everything in this structure is sorted */
    struct {
	ilpFloat d2;                /* Distance squared: calculated */
	ilpFloat m;                 /* Mass */
	} s;

    } *ILPTILE;

typedef struct ilpContext {
    ILPTILE first;              /* first tile in the chain */
    ILPTILE tile;               /* Current tile in the chain */
    double cx, cy, cz;          /* Center coordinates */
    double d2cmax;
    uint32_t nPrevious;         /* Particles in tiles prior to "tile" */
    } *ILP;

typedef struct {
    ILPTILE  tile;
    uint32_t nPart;
    uint32_t nPrevious;
    } ILPCHECKPT;

ILPTILE ilpExtend(ILP ilp);    /* Add tile and return new tile */
ILPTILE ilpClear(ILP ilp);     /* Go back to, and return first tile (empty) */
void ilpInitialize(ILP *ilp);
void ilpFinish(ILP ilp);

static inline void ilpCheckPt(ILP ilp,ILPCHECKPT *cp) {
    cp->tile = ilp->tile;
    cp->nPart = ilp->tile->nPart;
    cp->nPrevious = ilp->nPrevious;
    }

static inline void ilpRestore(ILP ilp,ILPCHECKPT *cp) {
    ilp->tile = cp->tile;
    ilp->nPrevious = cp->nPrevious;
    ilp->tile->nPart = cp->nPart;
    }

static inline uint32_t ilpCount(ILP ilp) {
    return ilp->nPrevious + ilp->tile->nPart;
    }

float ilpSelect(ILP ilp,uint32_t n,float *rMax);

#if defined(SYMBA) || defined(PLANETS)
#define ilpAppend_1(ilp,I) tile->d.iOrder.i[i] = (I);
#else
#define ilpAppend_1(ilp,I)
#endif

#if defined(HERMITE)
#define ilpAppend_2(ilp,VX,VY,VZ)					\
    tile->d.vx.f[i] = (VX);					\
    tile->d.vy.f[i] = (VY);					\
    tile->d.vz.f[i] = (VZ);
#else
#define ilpAppend_2(ilp,VX,VY,VZ)
#endif


#define ilpAppend(ilp,X,Y,Z,M,S,I,VX,VY,VZ)				\
    {									\
	ILPTILE tile = (ilp)->tile;					\
	uint_fast32_t i;						\
	if ( tile->nPart == tile->nMaxPart ) tile = ilpExtend((ilp));	\
	i = tile->nPart;						\
	tile->d.dx.f[i] = (ilp)->cx - (X);				\
	tile->d.dy.f[i] = (ilp)->cy - (Y);				\
	tile->d.dz.f[i] = (ilp)->cz - (Z);				\
	assert( (M) > 0.0 );						\
	tile->d.m.f[i] = (M);						\
	tile->d.fourh2.f[i] = (S);					\
	ilpAppend_1((ilp),I);						\
	ilpAppend_2((ilp),VX,VY,VZ);					\
	++tile->nPart;							\
    }

#define ILP_LOOP(ilp,ptile) for( ptile=(ilp)->first; ptile!=(ilp)->tile->next; ptile=ptile->next )



static inline void ilpCompute(ILP ilp, float fx, float fy, float fz ) {
    ILPTILE tile;
    uint32_t j;

#if defined(USE_SIMD_PP)
    v4sf px, py, pz, t1, t2, t3;
    px = SIMD_SPLAT(fx);
    py = SIMD_SPLAT(fy);
    pz = SIMD_SPLAT(fz);

    ILP_LOOP(ilp,tile) {
	uint32_t n = tile->nPart >> ILP_ALIGN_BITS; /* # Packed floats */
	uint32_t r = tile->nPart &  ILP_ALIGN_MASK; /* Remaining */

	if ( r != 0 ) {
	    for ( j=r; j<ILP_ALIGN_SIZE; j++ ) {
		int o = (n<<ILP_ALIGN_BITS) + j;
		tile->d.dx.f[o] = tile->d.dy.f[o] = tile->d.dz.f[o] = 0;
		tile->d.m.f[o] = 0;
		tile->d.fourh2.f[o] = tile->d.fourh2.f[0];
		}
	    n++;
	    }
	for ( j=0; j<n; j++ ) {
	    tile->d.dx.p[j] = t1 = SIMD_ADD(tile->d.dx.p[j],px);
	    tile->d.dy.p[j] = t2 = SIMD_ADD(tile->d.dy.p[j],py);
	    tile->d.dz.p[j] = t3 = SIMD_ADD(tile->d.dz.p[j],pz);
	    tile->d.d2.p[j] = SIMD_MADD(t3,t3,SIMD_MADD(t2,t2,SIMD_MUL(t1,t1)));
	    }
	}
#else
    ILP_LOOP(ilp,tile) {
	for (j=0;j<tile->nPart;++j) {
	    tile->d.dx.f[j] += fx;
	    tile->d.dy.f[j] += fy;
	    tile->d.dz.f[j] += fz;
	    tile->d.d2.f[j] = tile->d.dx.f[j]*tile->d.dx.f[j]
			      + tile->d.dy.f[j]*tile->d.dy.f[j] + tile->d.dz.f[j]*tile->d.dz.f[j];
	    }
	}
#endif
    }
#else /* LOCAL_EXPANSION */
typedef struct ilPart {
    double m,x,y,z;
#if defined(SOFTLINEAR)
    double h;
#elif defined(SOFTSQUARE)
    double twoh2;
#else
    double fourh2;
#endif  
#if defined(SYMBA) || defined(PLANETS) || !defined(LOCAL_EXPANSION)
    uint64_t iOrder;
#endif
#if defined(HERMITE) || !defined(LOCAL_EXPANSION)
    double vx,vy,vz;
#endif
    } ILP;
#endif /* LOCAL_EXPANSION */

#endif
