#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "ilp.h"
#include "pkd.h" /* for PARTITION macro */

/*
** Private: Create a new tile
*/
static ILPTILE newTile(ILPTILE prev) {
    ILPTILE tile = SIMD_malloc(sizeof(struct ilpTile));
    assert( tile != NULL );
    assert(ILP_PART_PER_TILE%4 == 0 );

    tile->next = NULL;
    tile->prev = prev;
    tile->nMaxPart = ILP_PART_PER_TILE;
    tile->nPart = 0;
    return tile;
    }

/*
** If the current tile is full (nPart == nMaxPart), then
** this function is called to get a new, empty tile.
*/
ILPTILE ilpExtend(ILP ilp) {
    assert( ilp != NULL );
    assert( ilp->tile != NULL );
    assert( ilp->first != NULL );
    assert( ilp->tile->nPart == ilp->tile->nMaxPart );

    ilp->nPrevious += ilp->tile->nPart;

    /* Use the next tile if it exists, or create a new one */
    if ( ilp->tile->next != NULL ) {
	ilp->tile = ilp->tile->next;
	ilp->tile->nPart = 0;
	}
    else {
	ilp->tile = ilp->tile->next = newTile(ilp->tile);
	}

    return ilp->tile;
    }

/*
** Empty the list of particles (go back to the first tile)
*/
ILPTILE ilpClear(ILP ilp) {
    assert( ilp != NULL );
    ilp->tile = ilp->first;
    ilp->nPrevious = 0;
    assert( ilp->tile != NULL );
    ilp->tile->nPart = 0;
    return ilp->tile;
    }

void ilpInitialize(ILP *ilp) {
    *ilp = malloc(sizeof(struct ilpContext));
    assert( *ilp != NULL );
    (*ilp)->first = (*ilp)->tile = newTile(NULL);
    (*ilp)->nPrevious = 0;
    }

void ilpFinish(ILP ilp) {
    ILPTILE tile, next;

    assert( ilp != NULL );

    /* Free all allocated tiles first */
    for ( tile=ilp->first; tile!=NULL; tile=next ) {
	next = tile->next;
	SIMD_free(tile);
	}

    free(ilp);
    }
