#ifndef ILP_H
#define ILP_H
#define ILP_PART_PER_TILE 8192

#if defined(__SSE__)
#include <xmmintrin.h>
#elif defined(__ALTIVEC__)
#include <altivec.h>
#endif

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILP_PART_PER_TILE];
#ifdef USE_SIMD_PP
#if defined(__SSE__)
    __m128 p[ILP_PART_PER_TILE/4];
#elif defined(__ALTIVEC__)
    vector float p[ILP_PART_PER_TILE/4];
#endif
#endif
} ilpFloat;

typedef union {
    uint64_t i;
} ilpInt64;


typedef struct ilpTile {
    struct ilpTile *next;
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

ILPTILE ilpExtend(ILP ilp);    /* Add tile and return new tile */
ILPTILE ilpClear(ILP ilp);     /* Go back to, and return first tile (empty) */
void ilpInitialize(ILP *ilp);
void ilpFinish(ILP ilp);
void ilpSetCount(ILP ilp,uint32_t count);

static inline uint32_t ilpCount(ILP ilp) {
    return ilp->nPrevious + ilp->tile->nPart;
}

#endif
