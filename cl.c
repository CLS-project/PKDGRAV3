#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>

#include "cl.h"

void clInitialize(CL *cl,LSTFREELIST *freeList) {
    *cl = malloc(sizeof(struct clContext));
    assert( *cl != NULL );
    lstInitialize(&(*cl)->lst,freeList,CL_BLK_PER_TILE, CL_PART_PER_BLK, 1,
	sizeof(CL_BLK),  lstSIMDAllocate,lstSIMDFree);
    }

void clFinish(CL cl) {
    lstFree(&cl->lst);
    free(cl);
    }

