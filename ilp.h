#ifndef ILP_H
#define ILP_H

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
	ilpFloat m;                 /* Mass */
	ilpFloat fourh2;            /* Softening: calculated */
	ilpFloat d2;                /* Distance squared: calculated */
#ifdef HERMITE
	ilpFloat vx, vy, vz;
#endif
#if defined(SYMBA) || defined(PLANETS)
	ilpInt64 iOrder;
#endif
    } d;
} *ILPTILE;

typedef struct ilpContext
{
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

static inline void ilpCheckPt(ILP ilp,ILPCHECKPT *cp)
{
    cp->tile = ilp->tile;
    cp->nPart = ilp->tile->nPart;
    cp->nPrevious = ilp->nPrevious;
}

static inline void ilpRestore(ILP ilp,ILPCHECKPT *cp)
{
    ilp->tile = cp->tile;
    ilp->nPrevious = cp->nPrevious;
    ilp->tile->nPart = cp->nPart;
}

static inline uint32_t ilpCount(ILP ilp) {
    return ilp->nPrevious + ilp->tile->nPart;
}


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
	tile->d.m.f[i] = (M);						\
	tile->d.fourh2.f[i] = (S);					\
	ilpAppend_1((ilp),I);						\
	ilpAppend_2((ilp),VX,VY,VZ);					\
	++tile->nPart;							\
    }

#define ILP_LOOP(ilp,ptile) for( ptile=(ilp)->first; ptile!=(ilp)->tile->next; ptile=ptile->next )

#endif
