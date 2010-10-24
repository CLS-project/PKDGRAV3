#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include "cl.h"

/*
** Private: Create a new tile
*/
static CLTILE newTile(CLTILE prev) {
    CLTILE tile = SIMD_malloc(sizeof(struct clTile));

    assert( tile != NULL );
    assert(CL_PART_PER_TILE%4 == 0 );

    tile->next = NULL;
    tile->prev = prev;
    tile->nMaxItems = CL_PART_PER_TILE;
    tile->nItems = 0;

    return tile;
    }

/*
** If the current tile is full (nItems == nMaxItems), then
** this function is called to get a new, empty tile.
*/
CLTILE clExtend(CL cl) {
    assert( cl != NULL );
    assert( cl->tile != NULL );
    assert( cl->first != NULL );
    assert( cl->tile->nItems == cl->tile->nMaxItems );

    cl->nPrevious += cl->tile->nItems;

    /* Use the next tile if it exists, or create a new one */
    if ( cl->tile->next != NULL ) {
	cl->tile = cl->tile->next;
	cl->tile->nItems = 0;
	}
    else {
	cl->tile = cl->tile->next = newTile(cl->tile);
	}

    return cl->tile;
    }

/*
** Empty the list of particles (go back to the first tile)
*/
CLTILE clClear(CL cl) {
    assert( cl != NULL );
    cl->tile = cl->first;
    cl->nPrevious = 0;
    assert( cl->tile != NULL );
    cl->tile->nItems = 0;
    return cl->tile;
    }

void clInitialize(CL *cl) {
    *cl = malloc(sizeof(struct clContext));
    assert( *cl != NULL );
    (*cl)->first = (*cl)->tile = newTile(NULL);
    (*cl)->nPrevious = 0;
    }

void clFinish(CL cl) {
    CLTILE tile, next;

    assert( cl != NULL );

    /* Free all allocated tiles first */
    for ( tile=cl->first; tile!=NULL; tile=next ) {
	next = tile->next;
	SIMD_free(tile);
	}

    free(cl);
    }

CLTILE clClone(CL cl,CL src) {
    CLTILE tile;
    float r[3], fOffset[3], fCenter[3], fMax[3];
    int j;

    clClear(cl);
    CL_LOOP(src,tile) {
	for ( j=0; j<tile->nItems; ++j ) {
	    r[0] = tile->x.f[j];
	    r[1] = tile->y.f[j];
	    r[2] = tile->z.f[j];
	    fOffset[0] = tile->xOffset.f[j];
	    fOffset[1] = tile->yOffset.f[j];
	    fOffset[2] = tile->zOffset.f[j];
	    fCenter[0] = tile->xCenter.f[j];
	    fCenter[1] = tile->yCenter.f[j];
	    fCenter[2] = tile->zCenter.f[j];
	    fMax[0] = tile->xMax.f[j];
	    fMax[1] = tile->yMax.f[j];
	    fMax[2] = tile->zMax.f[j];
	    clAppend(cl,tile->iCell.i[j],tile->id.i[j],tile->iLower.i[j],tile->nc.i[j],
		tile->cOpen.f[j],tile->m.f[j],tile->fourh2.f[j],
		r,fOffset,fCenter,fMax);
	}
    }
    return cl->tile;
}
