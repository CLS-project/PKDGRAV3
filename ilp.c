#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>

#include "ilp.h"
#include "cudautil.h"

void ilpInitialize(ILP *ilp) {
    *ilp = malloc(sizeof(struct ilpContext));
    assert( *ilp != NULL );
    lstInitialize(&(*ilp)->lst,NULL,
	ILP_TILE_SIZE / sizeof(ILP_BLK), ILP_PART_PER_BLK,
#ifdef TIMESTEP_CRITICAL
	2, /* two areas */
	sizeof(ILP_BLK),SIMD_malloc,SIMD_free);
	sizeof(ILP_EXTRA),NULL,           NULL);
#else
	1, /* one area */
	sizeof(ILP_BLK),SIMD_malloc,SIMD_free);
#endif
    (*ilp)->cx = (*ilp)->cy = (*ilp)->cz = 0.0;
    }

void ilpFinish(ILP ilp) {
    lstFree(&ilp->lst);
    free(ilp);
    }
