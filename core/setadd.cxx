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
#include "setadd.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void<ServiceSetAdd::input>()  || std::is_trivial<ServiceSetAdd::input>());
static_assert(std::is_void<ServiceSetAdd::output>() || std::is_trivial<ServiceSetAdd::output>());

// This is a weird service. Normally we use TraversePST, but here
// we are actually setting up the PST structure.
int ServiceSetAdd::operator()(int nIn, void *pIn, void *pOut) {
    assert(nIn == sizeof(input));
    SetAdd(node_pst,static_cast<input *>(pIn));
    return 0;
}

void ServiceSetAdd::SetAdd(PST pst,input *in) {
    PST pstNew;
    int n, idMiddle,iProcLower,iProcUpper;
    mdlassert(pst->mdl,pst->nLeaves==1);
    mdlassert(pst->mdl,in->idLower==mdlSelf(pst->mdl));
    n = in->idUpper - in->idLower;
    idMiddle = (in->idUpper + in->idLower) / 2;
    if ( n > 1 ) {
        int rID;
        /* Make sure that the pst lands on core zero */
        iProcLower = mdlThreadToProc(pst->mdl,in->idLower);
        iProcUpper = mdlThreadToProc(pst->mdl,in->idUpper-1);
        if (iProcLower!=iProcUpper) {
            idMiddle = mdlProcToThread(pst->mdl,mdlThreadToProc(pst->mdl,idMiddle));
        }
        pst->nLeaves += n - 1;
        pst->nLower = idMiddle - in->idLower;
        pst->nUpper = in->idUpper - idMiddle;

        in->idLower = idMiddle;
        pst->idUpper = in->idLower;
        rID = pst->mdl->ReqService(pst->idUpper,getServiceID(),in,sizeof(input));
        in->idLower = mdlSelf(pst->mdl);
        in->idUpper = idMiddle;
        pstInitialize(&pstNew,pst->mdl,pst->plcl);
        pst->pstLower = pstNew;
        pstNew->iLvl = pst->iLvl + 1;
        SetAdd(pst->pstLower,in);
        pst->mdl->GetReply(rID);
    }
}