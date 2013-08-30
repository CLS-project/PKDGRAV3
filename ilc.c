#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include "ilc.h"
#include "cudautil.h"

void ilcInitialize(ILC *ilc) {
    *ilc = malloc(sizeof(struct ilcContext));
    assert( *ilc != NULL );
    lstInitialize(&(*ilc)->lst,NULL,
	ILC_TILE_SIZE / sizeof(ILC_BLK), ILC_PART_PER_BLK,
	2, /* two areas */
#ifdef USE_CUDA
	sizeof(ILC_BLK),CUDA_malloc,CUDA_free,
#else
	sizeof(ILC_BLK),SIMD_malloc,SIMD_free,
#endif
	sizeof(ILC_XTR),NULL,           NULL);
    (*ilc)->cx = (*ilc)->cy = (*ilc)->cz = 0.0;
    }

void ilcFinish(ILC ilc) {
    lstFree(&ilc->lst);
    free(ilc);
    }
