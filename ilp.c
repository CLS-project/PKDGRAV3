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

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <assert.h>

#include "ilp.h"

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
