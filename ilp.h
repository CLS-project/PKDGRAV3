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
    int64_t i[ILP_PART_PER_TILE]; /* Signed because negative marks softened cells */
    } ilpInt64;


typedef struct ilpTile {
    struct ilpTile *next;
    struct ilpTile *prev;
    uint32_t nMaxPart;          /* Maximum number of particles in this tile */
    uint32_t nPart;             /* Current number of particles */

    ilpFloat dx, dy, dz;        /* Offset from ilp->cx, cy, cz */
    ilpFloat d2;                /* Distance squared: calculated */
    ilpFloat m;                 /* Mass */
    ilpFloat fourh2;            /* Softening: calculated */
/* #ifdef HERMITE */
    ilpFloat vx, vy, vz;
/* #endif */
/* #if defined(SYMBA) || defined(PLANETS) */
    ilpInt64 iOrder;
/* #endif */
    } *ILPTILE;

typedef struct ilpContext {
    ILPTILE first;              /* first tile in the chain */
    ILPTILE tile;               /* Current tile in the chain */
    double cx, cy, cz;          /* Center coordinates */
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

/* #if defined(SYMBA) || defined(PLANETS) */
#define ilpAppend_1(ilp,I) tile->iOrder.i[ILP_APPEND_i] = (I);
/* #else */
/* #define ilpAppend_1(ilp,I) */
/* #endif */

/* #if defined(HERMITE) */
#define ilpAppend_2(ilp,VX,VY,VZ)					\
    tile->vx.f[ILP_APPEND_i] = (VX);					\
    tile->vy.f[ILP_APPEND_i] = (VY);					\
    tile->vz.f[ILP_APPEND_i] = (VZ);
/* #else */
/* #define ilpAppend_2(ilp,VX,VY,VZ) */
/* #endif */


#define ilpAppend(ilp,X,Y,Z,M,S,I,VX,VY,VZ)				\
    {									\
    ILPTILE tile = (ilp)->tile;						\
    uint_fast32_t ILP_APPEND_i;						\
    if ( tile->nPart == tile->nMaxPart ) tile = ilpExtend((ilp));	\
    ILP_APPEND_i = tile->nPart;						\
    tile->dx.f[ILP_APPEND_i] = (ilp)->cx - (X);			\
    tile->dy.f[ILP_APPEND_i] = (ilp)->cy - (Y);			\
    tile->dz.f[ILP_APPEND_i] = (ilp)->cz - (Z);			\
    assert( (M) > 0.0 );						\
    tile->m.f[ILP_APPEND_i] = (M);					\
    tile->fourh2.f[ILP_APPEND_i] = (S);				\
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

	/*
	** This is a little trick to speed up the calculation. By setting
	** unused entries in the list to have a zero mass, the resulting
	** forces are zero. Setting the distance to a large value avoids
	** softening the non-existent forces which is slightly faster.
	*/
	if ( r != 0 ) {
	    for ( j=r; j<ILP_ALIGN_SIZE; j++ ) {
		int o = (n<<ILP_ALIGN_BITS) + j;
		tile->dx.f[o] = tile->dy.f[o] = tile->dz.f[o] = 1e18f;
		tile->m.f[o] = 0.0f;
		tile->fourh2.f[o] = tile->fourh2.f[0];
		}
	    n++;
	    }
	for ( j=0; j<n; j++ ) {
	    tile->dx.p[j] = t1 = SIMD_ADD(tile->dx.p[j],px);
	    tile->dy.p[j] = t2 = SIMD_ADD(tile->dy.p[j],py);
	    tile->dz.p[j] = t3 = SIMD_ADD(tile->dz.p[j],pz);
	    tile->d2.p[j] = SIMD_MADD(t3,t3,SIMD_MADD(t2,t2,SIMD_MUL(t1,t1)));
	    }
	}
#else
    ILP_LOOP(ilp,tile) {
	for (j=0;j<tile->nPart;++j) {
	    tile->dx.f[j] += fx;
	    tile->dy.f[j] += fy;
	    tile->dz.f[j] += fz;
	    tile->d2.f[j] = tile->dx.f[j]*tile->dx.f[j]
			      + tile->dy.f[j]*tile->dy.f[j] + tile->dz.f[j]*tile->dz.f[j];
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
