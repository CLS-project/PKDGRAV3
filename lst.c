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

#include "lst.h"
#include "core/simd.h"

#include <stdarg.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*
** The default allocate and free for areas is malloc and free
*/
static void *defaultAllocate(size_t nBytes) {
    return malloc(nBytes);
}

static void defaultFree(void *data) {
    free(data);
}

static void cloneTile(LST *lst, LSTTILE *dst, LSTTILE *src) {
    void **sb, **db;
    int i;
    int nBlocks;

    dst->nBlocks = src->nBlocks;
    dst->nInLast = src->nInLast;

    sb = (void **)(src+1);
    db = (void **)(dst+1);
    nBlocks = src->nBlocks + (src->nInLast ? 1 : 0);

    for ( i=0; i<lst->nAreas; ++i) {
        memcpy(db[i],sb[i],lst->info[i].nAreaSize * nBlocks);
    }
}

LSTTILE *lstNewTile(LST *lst) {
    LSTTILE *tile;
    void **blks;
    int i;

    if ( lst->freeList->list != NULL ) {
        tile = lst->freeList->list;
        lst->freeList->list = tile->next;
    }
    else {
        tile = malloc(sizeof(struct lstTile) + lst->nAreas*sizeof(void *));
        assert(tile!=NULL);
        blks = (void **)(tile+1);
        for ( i=0; i<lst->nAreas; ++i) {
            blks[i] = (*lst->info[i].fnAllocate)(lst->info[i].nAreaSize * lst->nBlocksPerTile);
            assert(blks[i] != NULL);
            /*
            ** This memset is important. For SIMD versions, the lists must contain valid
            ** floating point data so we set everything to zero.
            */
            memset(blks[i],0,lst->info[i].nAreaSize * lst->nBlocksPerTile);
        }
        lst->freeList->nTiles++;
    }

    tile->next = NULL;
    tile->nBlocks = tile->nInLast = 0;
    tile->nRefs = 1;

    return tile;
}

void lstFreeTile(LST *lst,LSTTILE *tile) {
    /* If this is also owned by someone else then they need to give it back */
    if ( --tile->nRefs == 0 ) {
        tile->next = lst->freeList->list;
        lst->freeList->list = tile;
    }
}

size_t lstMemory(LST *lst) {
    return sizeof(struct lstContext) +
           lst->nAreas * sizeof(LSTAREAINFO) +
           lst->freeList->nTiles * lst->nTileSize;
}

static void moveToFreeList(LST *lst,LSTTILE *tile) {
    LSTTILE *next;
    for (; tile!=NULL; tile=next) {
        next = tile->next;
        lstFreeTile(lst,tile);
    }
}

void lstClear(LST *lst) {
    assert( lst != NULL );
    moveToFreeList(lst,lst->list);
    lst->tile = lst->list = lstNewTile(lst);
    lst->nPrevious = 0;
}

void lstCheckPt(LST *lst,LSTCHECKPT *cp) {
    cp->nBlocks = lst->tile->nBlocks;
    cp->nInLast = lst->tile->nInLast;
    cp->nPrevious = lst->nPrevious;
}

/*
** This is called if we need to split off the last block
** because we want to add to it but it is in use.
*/
LSTTILE *lstSplit(LST *lst) {
    LSTTILE *tile, *prev;
    assert(lst->tile->nRefs > 1);
    tile = lstNewTile(lst);
    cloneTile(lst, tile, lst->tile);
    if (lst->tile == lst->list) lst->list = tile;
    else {
        for (prev=lst->list; prev->next != lst->tile; prev=prev->next) {}
        prev->next = tile;
    }
    lstFreeTile(lst,lst->tile);
    lst->tile = tile;
    return tile;
}

void lstRestore(LST *lst,LSTCHECKPT *cp) {
    LSTTILE *tile = lst->list;
    uint32_t nPrevious = 0;

    while ( nPrevious < cp->nPrevious ) {
        assert(tile!=NULL);
        nPrevious += tile->nBlocks * lst->nPerBlock + tile->nInLast;
        tile = tile->next;
    }
    assert(nPrevious == cp->nPrevious);

    /* Set the latest tile */
    lst->tile = tile;
    lst->nPrevious = cp->nPrevious;
    lst->tile->nBlocks = cp->nBlocks;
    lst->tile->nInLast = cp->nInLast;

    /* Dump any tiles that are not needed */
    moveToFreeList(lst,lst->tile->next);
    lst->tile->next = NULL;
}

/* Add a tile to the list */
void *lstExtend(LST *lst) {
    assert( lst != NULL );
    assert( lst->tile->nInLast==lst->nPerBlock ); /* It would be silly to extend a perfectly good list */
    lst->nPrevious += lst->tile->nBlocks*lst->nPerBlock  + lst->tile->nInLast;
    lst->tile = lst->tile->next = lstNewTile(lst);
    return lst->tile;
}

void lstClone(LST *dst,LST *src) {
    LSTTILE *stile, *dtile;

    moveToFreeList(dst,dst->list);
    dst->tile = dst->list = lstNewTile(dst);

    stile = src->list;
    dtile = dst->tile;
    while (stile != NULL) {
        cloneTile(src, dtile, stile);
        stile = stile->next;
        if (stile) {
            dtile->next = lstNewTile(dst);
            dtile = dtile->next;
        }
    }
    dst->tile = dtile;
    dst->nPrevious = src->nPrevious;
}

void lstInitialize(LST *lst, LSTFREELIST *freeList, int nBlocksPerTile, int nPerBlock, int nAreas, ...) {
    va_list args;
    int i;
    assert(lst!=NULL);

    if (nBlocksPerTile==0) nBlocksPerTile=1;
    lst->nAreas = nAreas;
    lst->nBlocksPerTile = nBlocksPerTile;
    lst->nPerBlock = nPerBlock;
    lst->nTileSize = sizeof(struct lstTile);
    lst->defFreeList.list = NULL;
    lst->defFreeList.nRefs = 0;
    lst->defFreeList.nTiles = 0;
    if (freeList) lst->freeList = freeList;
    else lst->freeList = &lst->defFreeList;
    lst->freeList->nRefs++;
    lst->info = malloc( lst->nAreas * sizeof(LSTAREAINFO) );
    va_start(args,nAreas);
    for ( i=0; i<lst->nAreas; ++i) {
        lst->info[i].nAreaSize = va_arg(args,int);
        lst->info[i].fnAllocate = va_arg(args,lstAreaAllocate);
        lst->info[i].fnFree     = va_arg(args,lstAreaFree );
        if (lst->info[i].fnAllocate==NULL && lst->info[i].fnFree==NULL) {
            lst->info[i].fnAllocate = defaultAllocate;
            lst->info[i].fnFree = defaultFree;
        }
        assert(lst->info[i].fnAllocate!=NULL);
        assert(lst->info[i].fnFree!=NULL);
        lst->nTileSize += lst->info[i].nAreaSize * nBlocksPerTile;
    }
    va_end(args);
    lst->list = lst->tile = lstNewTile(lst);
}

void lstFree(LST *lst) {
    LSTTILE *tile, *next;
    int i;
    void **blks;

    moveToFreeList(lst,lst->list);
    if (--lst->freeList->nRefs == 0) {
        for ( tile= lst->freeList->list; tile!=NULL; tile=next) {
            next = tile->next;
            blks = (void **)(tile+1);
            for ( i=0; i<lst->nAreas; ++i) {
                (*lst->info[i].fnFree)(blks[i]);
            }
            lst->freeList->nTiles--;
            free(tile);
        }
        assert(lst->freeList->nTiles == 0);
        lst->freeList->list = NULL;
    }
    free(lst->info);
    lst->info = NULL;
}
