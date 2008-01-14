#ifndef ILP_H
#define ILP_H

#ifndef ILP_PART_PER_TILE
#define ILP_PART_PER_TILE 2048
#endif

#define ILP_ALIGN_BITS 2
#define ILP_ALIGN_SIZE (1<<ILP_ALIGN_BITS)

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
    uint64_t i;
} ilpInt64;


typedef struct ilpTile {
    struct ilpTile *next;
    struct ilpTile *prev;
    uint32_t nMaxPart;          /* Maximum number of particles in this tile */
    uint32_t nPart;             /* Current number of particles */

    ilpFloat dx, dy, dz;        /* Offset from ilp->cx, cy, cz */
    ilpFloat m;                 /* Mass */
    ilpFloat fourh2;            /* Softening: calculated */
    ilpFloat d2;                /* Distance squared: calculated */
#ifdef HERMITE
    ilpFloat vx, vy, vz;
#endif
#if defined(SYMBA) || defined(PLANETS)
    ilpInt64 iOrder[ILP_PART_PER_TILE];
#endif
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

#endif
