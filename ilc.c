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

#include "ilc.h"

void ilcInitialize(ILC *ilc) {
    *ilc = malloc(sizeof(struct ilcContext));
    assert( *ilc != NULL );
	lstInitialize(&(*ilc)->lst, NULL,
		ILC_TILE_SIZE / sizeof(ILC_BLK), ILC_PART_PER_BLK,
		1, /* one area */
		sizeof(ILC_BLK), SIMD_malloc, SIMD_free);
    (*ilc)->cx = (*ilc)->cy = (*ilc)->cz = 0.0;
    }

void ilcFinish(ILC ilc) {
    lstFree(&ilc->lst);
    free(ilc);
    }
