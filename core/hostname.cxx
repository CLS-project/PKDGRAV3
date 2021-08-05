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
#include "hostname.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void<ServiceHostname::input>()  || std::is_trivial<ServiceHostname::input>());
static_assert(std::is_void<ServiceHostname::output>() || std::is_trivial<ServiceHostname::output>());

int ServiceHostname::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto out = static_cast<output *>(vout);
    auto outUp = out + pst->idUpper-pst->idSelf;
    auto rID = mdlReqService(pst->mdl,pst->idUpper,getServiceID(),vin,nIn);
    Traverse(pst->pstLower,vin,nIn,out,nOut);
    mdlGetReply(pst->mdl,rID,outUp,NULL);
    return pst->nLeaves*sizeof(output);
    }

int ServiceHostname::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto out = static_cast<output *>(vout);
    out->iMpiID = mdlSelf(pst->mdl);
    strncpy(out->szHostname,mdlName(pst->mdl),sizeof(out->szHostname));
    out->szHostname[sizeof(out->szHostname)-1] = 0;
    auto p = strchr(out->szHostname,'.');
    if (p) *p = 0;
    return sizeof(output);
    }
