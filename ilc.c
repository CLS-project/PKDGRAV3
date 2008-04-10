#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include "ilc.h"

/*
** Private: Create a new tile
*/
static ILCTILE newTile(ILCTILE prev) {
    ILCTILE tile = SIMD_malloc(sizeof(struct ilcTile));
    assert( tile != NULL );
    assert(ILC_PART_PER_TILE%4 == 0 );

    tile->next = NULL;
    tile->prev = prev;
    tile->nMaxCell = ILC_PART_PER_TILE;
    tile->nCell = 0;
    return tile;
    }

/*
** If the current tile is full (nCell == nMaxCell), then
** this function is called to get a new, empty tile.
*/
ILCTILE ilcExtend(ILC ilc) {
    assert( ilc != NULL );
    assert( ilc->tile != NULL );
    assert( ilc->first != NULL );
    assert( ilc->tile->nCell == ilc->tile->nMaxCell );

    ilc->nPrevious += ilc->tile->nCell;

    /* Use the next tile if it exists, or create a new one */
    if ( ilc->tile->next != NULL ) {
	ilc->tile = ilc->tile->next;
	ilc->tile->nCell = 0;
	}
    else {
	ilc->tile = ilc->tile->next = newTile(ilc->tile);
	}

    return ilc->tile;
    }

/*
** Empty the list of particles (go back to the first tile)
*/
ILCTILE ilcClear(ILC ilc) {
    assert( ilc != NULL );
    ilc->tile = ilc->first;
    ilc->nPrevious = 0;
    assert( ilc->tile != NULL );
    ilc->tile->nCell = 0;
    return ilc->tile;
    }

void ilcInitialize(ILC *ilc) {
    *ilc = malloc(sizeof(struct ilcContext));
    assert( *ilc != NULL );
    (*ilc)->first = (*ilc)->tile = newTile(NULL);
    (*ilc)->nPrevious = 0;
    }

void ilcFinish(ILC ilc) {
    ILCTILE tile, next;

    assert( ilc != NULL );

    /* Free all allocated tiles first */
    for ( tile=ilc->first; tile!=NULL; tile=next ) {
	next = tile->next;
	SIMD_free(tile);
	}

    free(ilc);
    }
