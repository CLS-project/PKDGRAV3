#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>

#include "ilp.h"
#include "cudautil.h"

void ilpInitialize(ILP *ilp) {
    *ilp = malloc(sizeof(struct ilpContext));
    assert( *ilp != NULL );
    lstInitialize(&(*ilp)->lst,NULL,ILP_BLK_PER_TILE, ILP_PART_PER_BLK, 2, /* two areas */
#ifdef USE_CUDA
	sizeof(ILP_BLK),CUDA_malloc,CUDA_free,
#else
	sizeof(ILP_BLK),SIMD_malloc,SIMD_free,
#endif
	sizeof(ILP_EXTRA),NULL,           NULL);
    (*ilp)->cx = (*ilp)->cy = (*ilp)->cz = 0.0;
    }

void ilpFinish(ILP ilp) {
    lstFree(&ilp->lst);
    free(ilp);
    }
