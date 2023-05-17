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
#ifndef CORE_REMOTE_H
#define CORE_REMOTE_H
#include "pkd_config.h"
#include "treenode.h"
#include "particle.h"
#include "mdl.h"

template<class STORE>
class remoteStore {
protected:
    mdl::mdlClass *mdl;
    STORE *localStore;
    int cache_id;
    auto &local() { return *localStore; }
public:
    remoteStore(mdl::mdlClass *mdl,STORE &localStore,int cache_id)
        : mdl(mdl), localStore(&localStore), cache_id(cache_id) {}

    auto operator()(int iIndex) {
        return local()[iIndex];
    }
    auto operator()(int iIndex,int iID) {
        if (iID==mdl->Self()) return (*this)(iIndex);
        else return local()[reinterpret_cast<typename STORE::value_type *>(mdlFetch(mdl,cache_id,iIndex,iID))];
    }
};

using remoteTree = remoteStore<treeStore>;
using remoteParticles = remoteStore<particleStore>;

#endif
