#ifndef ILC_H
#define ILC_H
#ifdef LOCAL_EXPANSION
#ifndef ILC_PART_PER_TILE
#define ILC_PART_PER_TILE 1024 /* 1024*100 ~ 100k */
#endif

#define ILC_ALIGN_BITS 2
#define ILC_ALIGN_SIZE (1<<ILC_ALIGN_BITS)
#define ILC_ALIGN_MASK (ILC_ALIGN_SIZE-1)

#include "simd.h"

/*
** We use a union here so that the compiler can properly align the values.
*/

#ifdef USE_SIMD_MOMR
#define ILC_BLOCKS (ILC_PART_PER_TILE/4)
typedef union {
    MOMFLOAT f[4];
    v4sf p;
    } ilcFloat;
#else
#define ILC_BLOCKS ILC_PART_PER_TILE
typedef union {
    MOMFLOAT f;
    } ilcFloat;
#endif

typedef union {
    uint64_t i;
    } ilcInt64;


typedef struct ilcTile {
    struct ilcTile *next;
    struct ilcTile *prev;
    uint32_t nMaxCell;          /* Maximum number of cells in this tile */
    uint32_t nCell;             /* Current number of cells */

    struct {
	ilcFloat x,y,z;
#ifdef HERMITE
	ilcFloat vx,vy,vz;
#endif
	ilcFloat xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
	ilcFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	ilcFloat xx,xy,xz,yy,yz;
	ilcFloat m;
	} d[ILC_BLOCKS];
    } *ILCTILE;

typedef struct ilcContext {
    ILCTILE first;              /* first tile in the chain */
    ILCTILE tile;               /* Current tile in the chain */
    double cx, cy, cz;          /* Center coordinates */
    double d2cmax;
    uint32_t nPrevious;         /* Particles in tiles prior to "tile" */
    } *ILC;

typedef struct {
    ILCTILE  tile;
    uint32_t nCell;
    uint32_t nPrevious;
    } ILCCHECKPT;

ILCTILE ilcExtend(ILC ilc);    /* Add tile and return new tile */
ILCTILE ilcClear(ILC ilc);     /* Go back to, and return first tile (empty) */
void ilcInitialize(ILC *ilc);
void ilcFinish(ILC ilc);

static inline void ilcCheckPt(ILC ilc,ILCCHECKPT *cp) {
    cp->tile = ilc->tile;
    cp->nCell = ilc->tile->nCell;
    cp->nPrevious = ilc->nPrevious;
    }

static inline void ilcRestore(ILC ilc,ILCCHECKPT *cp) {
    ilc->tile = cp->tile;
    ilc->nPrevious = cp->nPrevious;
    ilc->tile->nCell = cp->nCell;
    }

static inline uint32_t ilcCount(ILC ilc) {
    return ilc->nPrevious + ilc->tile->nCell;
    }

#define ILC_LOOP(ilc,ctile) for( ctile=(ilc)->first; ctile!=(ilc)->tile->next; ctile=ctile->next )
#else /* LOCAL_EXPANSION */
#include "moments.h"
/*
** components required for evaluating a multipole interaction.
*/

typedef struct ilCell {
    double x,y,z;
    MOMR mom;
    } ILC;
#endif /* LOCAL_EXPANSION */
#endif
