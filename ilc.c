#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include "ilc.h"

void ilcInitialize(ILC *ilc) {
    *ilc = malloc(sizeof(struct ilcContext));
    assert( *ilc != NULL );
    lstInitialize(&(*ilc)->lst,NULL,ILC_BLK_PER_TILE, ILC_PART_PER_BLK, 2, /* two areas */
	sizeof(ILC_BLK),lstSIMDAllocate,lstSIMDFree,
	sizeof(ILC_XTR),NULL,           NULL);
    (*ilc)->cx = (*ilc)->cy = (*ilc)->cz = 0.0;
    }

void ilcFinish(ILC ilc) {
    lstFree(&ilc->lst);
    free(ilc);
    }
