#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include "ilp.h"

void ilpInitialize(ILP *ilp) {
    *ilp = malloc(sizeof(struct ilpContext));
    assert( *ilp != NULL );
    lstInitialize(&(*ilp)->lst,NULL,ILP_BLK_PER_TILE, ILP_PART_PER_BLK, 2, /* two areas */
	sizeof(ILP_BLK),  lstSIMDAllocate,lstSIMDFree,
	sizeof(ILP_EXTRA),NULL,           NULL);
    (*ilp)->cx = (*ilp)->cy = (*ilp)->cz = 0.0;
    }

void ilpFinish(ILP ilp) {
    lstFree(&ilp->lst);
    free(ilp);
    }
