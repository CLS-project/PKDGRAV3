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
#include "freestore.h"

static_assert(std::is_void<ServiceFreeStore::input>()  || std::is_trivial<ServiceFreeStore::input>());
static_assert(std::is_void<ServiceFreeStore::output>() || std::is_trivial<ServiceFreeStore::output>());
int ServiceFreeStore::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    static_assert(std::is_void<input>());
    auto out = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nOut==sizeof(output));
    *out = pkd->FreeStore();
    return sizeof(output);
}
