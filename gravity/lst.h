/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
** Generic List Management
**
*/

#ifndef LST_H
#define LST_H
#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <stdlib.h>
#include <stdint.h>

typedef struct lstTile {
    struct lstTile *next;
    uint16_t nBlocks;
    uint16_t nInLast;
    uint32_t nRefs;
} LSTTILE;
/* N.B. The Area pointers immediately follow */

struct lstContext;
typedef void *(*lstAreaAllocate)(size_t);
typedef void  (*lstAreaFree)(void *);

typedef struct {
    lstAreaAllocate fnAllocate; /* Function to allocate memory */
    lstAreaFree fnFree;         /* Function to free memory */
    size_t nAreaSize;           /* Size of each block: a tile will have many */
} LSTAREAINFO;

typedef struct {
    LSTTILE *list;     /* List of available tiles */
    size_t nTiles;     /* Current number of allocated tiles */
    int nRefs;         /* Number of lists using this free list */
} LSTFREELIST;

typedef struct lstContext {
    LSTAREAINFO *info;      /* Information on areas in each tile */
    LSTFREELIST *freeList;  /* Tiles that are not currently being used */
    LSTFREELIST defFreeList;
    LSTTILE *list;          /* The first tile in our list */
    LSTTILE *tile;          /* The last (current) tile in the list */
    size_t nTileSize;       /* How much memory a single tile uses */
    int nAreas;             /* Number of areas in info[] */
    int nBlocksPerTile;     /* The number of blocks in each tile */
    int nPerBlock;          /* Number of items in each block */
    int nPrevious;
} LST;

typedef struct {
    uint16_t nBlocks;
    uint16_t nInLast;
    uint32_t nPrevious;
} LSTCHECKPT;

#ifdef __cplusplus
extern "C" {
#endif
void *lstSIMDAllocate(size_t nBytes);
void lstSIMDFree(void *data);
#ifdef USE_CUDA
void *lstCUDAAllocate(size_t nBytes);
void lstCUDAFree(void *data);
#endif
void lstInitialize(LST *lst, LSTFREELIST *freeList, int nBlocksPerTile, int nPerBlock, int nAreas, ...);
void lstFree(LST *lst);
void lstFreeTile(LST *lst,LSTTILE *tile);
size_t lstMemory(LST *lst);
void lstCheckPt(LST *lst,LSTCHECKPT *cp);
void lstRestore(LST *lst,LSTCHECKPT *cp);
void lstClone(LST *dst,LST *src);
LSTTILE *lstSplit(LST *lst);

void lstClear(LST *lst);
void *lstExtend(LST *lst);
#ifdef __cplusplus
}
#endif
static inline uint32_t lstCount(LST *lst) {
    return lst->nPrevious + lst->tile->nBlocks*lst->nPerBlock  + lst->tile->nInLast;
}

static inline void *lstReposition(LST *lst) {
    LSTTILE *tile = lst->tile;
    if (tile->nRefs > 1) tile=lstSplit(lst);
    if (tile->nInLast == lst->nPerBlock ) {
        if ( ++tile->nBlocks == lst->nBlocksPerTile ) {
            --tile->nBlocks;
            tile = (LSTTILE *)lstExtend(lst);
        }
        else tile->nInLast = 0;
    }
    return tile;
}
#endif
