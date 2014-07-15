#ifndef ILC_H
#define ILC_H
#include <stdint.h>
#include <assert.h>
#include <stdint.h>

#include "lst.h"

#define ILC_TILE_SIZE (100*1024) /* 100k */
#ifndef ILC_PART_PER_BLK
#define ILC_PART_PER_BLK (32) /* Don't mess with this: see CUDA */
#endif

#if !defined(__CUDACC__)
#include "simd.h"
#endif

#include "moments.h"

/*
** We use a union here so that the compiler can properly align the values.
*/
typedef union {
    float f[ILC_PART_PER_BLK];
#if defined(USE_SIMD_PC) && !defined(__CUDACC__)
    v_sf p[ILC_PART_PER_BLK/SIMD_WIDTH];
#endif
    } ilcFloat;

typedef struct {
    ilcFloat dx,dy,dz;
    ilcFloat xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
    ilcFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    ilcFloat xx,xy,xz,yy,yz;
    ilcFloat x,y,z;
    ilcFloat m,u;
    } ILC_BLK;

typedef struct {
    char dummy[0];
    } ILC_XTR;

typedef struct ilcTile {
    LSTTILE lstTile;
    ILC_BLK *blk;
    ILC_XTR *xtr;
    } *ILCTILE;

typedef struct ilcContext {
    LST lst;
    double cx, cy, cz;          /* Center coordinates */
    } *ILC;

#define ILCCHECKPT LSTCHECKPT

#define ilcClear(ilc) lstClear(&(ilc)->lst)
#define ilcExtend(ilc) lstExtend(&(ilc)->lst)
#define ilcMemory(ilc) lstMemory(&(ilc)->lst)
#define ilcCheckPt(ilc,cp) lstCheckPt(&(ilc)->lst,(cp))
#define ilcRestore(ilc,cp) lstRestore(&(ilc)->lst,(cp))
#define ilcCount(ilc) lstCount(&(ilc)->lst)

void ilcInitialize(ILC *ilc);
void ilcFinish(ILC ilc);

/*
** The X, Y and Z coordinates must already be relative to cx, cy and cz!
*/
static inline void ilcAppendFloat(ILC ilc,float X,float Y,float Z,FMOMR *M,float U,float VX,float VY,float VZ) {
    ILCTILE tile = (ILCTILE)lstReposition(&ilc->lst);
    uint_fast32_t blk = tile->lstTile.nBlocks;
    uint_fast32_t prt = tile->lstTile.nInLast;

    tile->blk[blk].dx.f[prt] = (X);
    tile->blk[blk].dy.f[prt] = (Y);
    tile->blk[blk].dz.f[prt] = (Z);
    assert( (M)->m > 0.0 );
    tile->blk[blk].m.f[prt] = (M)->m;
    tile->blk[blk].u.f[prt] = (U);
    tile->blk[blk].x.f[prt] = (M)->x;
    tile->blk[blk].y.f[prt] = (M)->y;
    tile->blk[blk].z.f[prt] = (M)->z;
    tile->blk[blk].xx.f[prt] = (M)->xx;
    tile->blk[blk].xy.f[prt] = (M)->xy;
    tile->blk[blk].xz.f[prt] = (M)->xz;
    tile->blk[blk].yy.f[prt] = (M)->yy;
    tile->blk[blk].yz.f[prt] = (M)->yz;
    tile->blk[blk].xxx.f[prt] = (M)->xxx;
    tile->blk[blk].xxy.f[prt] = (M)->xxy;
    tile->blk[blk].xxz.f[prt] = (M)->xxz;
    tile->blk[blk].xyy.f[prt] = (M)->xyy;
    tile->blk[blk].xyz.f[prt] = (M)->xyz;
    tile->blk[blk].yyy.f[prt] = (M)->yyy;
    tile->blk[blk].yyz.f[prt] = (M)->yyz;
    tile->blk[blk].xxxx.f[prt] = (M)->xxxx;
    tile->blk[blk].xxxy.f[prt] = (M)->xxxy;
    tile->blk[blk].xxxz.f[prt] = (M)->xxxz;
    tile->blk[blk].xxyy.f[prt] = (M)->xxyy;
    tile->blk[blk].xxyz.f[prt] = (M)->xxyz;
    tile->blk[blk].xyyy.f[prt] = (M)->xyyy;
    tile->blk[blk].xyyz.f[prt] = (M)->xyyz;
    tile->blk[blk].yyyy.f[prt] = (M)->yyyy;
    tile->blk[blk].yyyz.f[prt] = (M)->yyyz;

    ++tile->lstTile.nInLast;
    }

static inline void ilcAppend(ILC ilc,double X,double Y,double Z,FMOMR *M,float U,float VX,float VY,float VZ) {
    ilcAppendFloat(ilc,(ilc)->cx-(X),(ilc)->cy-(Y),(ilc)->cz-(Z),M,U,VX,VY,VZ);
    }

#define ILC_LOOP(ilc,ctile) for( ctile=(ILCTILE)((ilc)->lst.list); ctile!=NULL; ctile=(ILCTILE)(ctile->lstTile.next) )

#endif
