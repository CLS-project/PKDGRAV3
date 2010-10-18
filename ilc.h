#ifndef ILC_H
#define ILC_H
#ifdef LOCAL_EXPANSION
#include <stdint.h>

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
typedef union {
    float f[ILC_PART_PER_TILE];
#ifdef USE_SIMD_PC
    v4sf p[ILC_PART_PER_TILE/4];
#endif
    } ilcFloat;

typedef union {
    uint64_t i[ILC_PART_PER_TILE];
    } ilcInt64;


typedef struct ilcTile {
    struct ilcTile *next;
    struct ilcTile *prev;
    uint32_t nMaxCell;          /* Maximum number of cells in this tile */
    uint32_t nCell;             /* Current number of cells */

    struct {
	ilcFloat dx,dy,dz;
	ilcFloat d2;
#ifdef HERMITE
	ilcFloat vx,vy,vz;
#endif
	ilcFloat xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
	ilcFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	ilcFloat xx,xy,xz,yy,yz;
	ilcFloat m,u;
	} d;
/*
** ILP has some extra data here that gets partitioned to find the closest
** particles. ILC doesn't need this (at least not yet).
*/
    } *ILCTILE;

typedef struct ilcContext {
    ILCTILE first;              /* first tile in the chain */
    ILCTILE tile;               /* Current tile in the chain */
    double cx, cy, cz;          /* Center coordinates */
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

#if defined(HERMITE)
#define ilcAppend_2(ilc,VX,VY,VZ)					\
    tile->d.vx.f[ILC_APPEND_i] = (VX);					\
    tile->d.vy.f[ILC_APPEND_i] = (VY);					\
    tile->d.vz.f[ILC_APPEND_i] = (VZ);
#else
#define ilcAppend_2(ilc,VX,VY,VZ)
#endif


#define ilcAppend(ilc,X,Y,Z,M,U,VX,VY,VZ)				\
    {									\
    ILCTILE tile = (ilc)->tile;						\
    uint_fast32_t ILC_APPEND_i;						\
    if ( tile->nCell == tile->nMaxCell ) tile = ilcExtend((ilc));	\
    ILC_APPEND_i = tile->nCell;						\
    tile->d.dx.f[ILC_APPEND_i] = (ilc)->cx - (X);			\
    tile->d.dy.f[ILC_APPEND_i] = (ilc)->cy - (Y);			\
    tile->d.dz.f[ILC_APPEND_i] = (ilc)->cz - (Z);			\
    assert( (M)->m > 0.0 );						\
    tile->d.m.f[ILC_APPEND_i] = (M)->m;					\
    tile->d.u.f[ILC_APPEND_i] = (U);					\
    tile->d.xx.f[ILC_APPEND_i] = (M)->xx;      				\
    tile->d.xy.f[ILC_APPEND_i] = (M)->xy;	       			\
    tile->d.xz.f[ILC_APPEND_i] = (M)->xz;	       			\
    tile->d.yy.f[ILC_APPEND_i] = (M)->yy;      				\
    tile->d.yz.f[ILC_APPEND_i] = (M)->yz;				\
    tile->d.xxx.f[ILC_APPEND_i] = (M)->xxx;            			\
    tile->d.xxy.f[ILC_APPEND_i] = (M)->xxy;            			\
    tile->d.xxz.f[ILC_APPEND_i] = (M)->xxz;            			\
    tile->d.xyy.f[ILC_APPEND_i] = (M)->xyy;            			\
    tile->d.xyz.f[ILC_APPEND_i] = (M)->xyz;            			\
    tile->d.yyy.f[ILC_APPEND_i] = (M)->yyy;            			\
    tile->d.yyz.f[ILC_APPEND_i] = (M)->yyz;            			\
    tile->d.xxxx.f[ILC_APPEND_i] = (M)->xxxx;                  		\
    tile->d.xxxy.f[ILC_APPEND_i] = (M)->xxxy;                  		\
    tile->d.xxxz.f[ILC_APPEND_i] = (M)->xxxz;                  		\
    tile->d.xxyy.f[ILC_APPEND_i] = (M)->xxyy;                  		\
    tile->d.xxyz.f[ILC_APPEND_i] = (M)->xxyz;                  		\
    tile->d.xyyy.f[ILC_APPEND_i] = (M)->xyyy;                  		\
    tile->d.xyyz.f[ILC_APPEND_i] = (M)->xyyz;                  		\
    tile->d.yyyy.f[ILC_APPEND_i] = (M)->yyyy;                  		\
    tile->d.yyyz.f[ILC_APPEND_i] = (M)->yyyz;                  		\
    ilcAppend_2((ilc),VX,VY,VZ);					\
    ++tile->nCell;							\
    }

#define ILC_LOOP(ilc,ctile) for( ctile=(ilc)->first; ctile!=(ilc)->tile->next; ctile=ctile->next )

static inline void ilcCompute(ILC ilc, float fx, float fy, float fz ) {
    ILCTILE tile;
    uint32_t j;

#if defined(USE_SIMD_PC)
    v4sf px, py, pz, t1, t2, t3;
    px = SIMD_SPLAT(fx);
    py = SIMD_SPLAT(fy);
    pz = SIMD_SPLAT(fz);

    ILC_LOOP(ilp,tile) {
	uint32_t n = tile->nCell >> ILC_ALIGN_BITS; /* # Packed floats */
	uint32_t r = tile->nCell &  ILC_ALIGN_MASK; /* Remaining */

	if ( r != 0 ) {
	    for ( j=r; j<ILC_ALIGN_SIZE; j++ ) {
		int o = (n<<ILC_ALIGN_BITS) + j;
		tile->d.dx.f[o] = tile->d.dy.f[o] = tile->d.dz.f[o] = 1e18;
		tile->d.m.f[o] = 0;
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
    ILC_LOOP(ilc,tile) {
	for (j=0;j<tile->nCell;++j) {
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
