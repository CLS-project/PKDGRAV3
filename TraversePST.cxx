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
#include "TraversePST.h"

int TraversePST::operator()(int nIn, void *pIn, void *pOut) {
    return Traverse(node_pst,pIn,nIn,pOut,getMaxBytesOut());
    }

int TraversePST::Traverse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    if (pstAmCore(pst))       return Service(pst,vin,nIn,vout,nOut);
    else if (pstOffNode(pst)) return OffNode(pst,vin,nIn,vout,nOut);
    else if (pstAmNode(pst))  return  AtNode(pst,vin,nIn,vout,nOut);
    else                      return Recurse(pst,vin,nIn,vout,nOut);
    }

// This is the default function for when we are still traversing the PST
// (OffNode,AtNode,Recurse) and will just naively recurse. Once we reach
// the leaf (AmCore) then Service() is called and must be provided.
int TraversePST::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = reinterpret_cast<mdl::mdlClass *>(pst->mdl);
    auto rID = mdl->ReqService(pst->idUpper,getServiceID(),vin,nIn);
    Traverse(pst->pstLower,vin,nIn,vout,nOut);
    mdl->GetReply(rID,vout,&nOut);
    return nOut;
    }

// This class is used for services that pass the same input (or none),
// and need to combine output results of the same size and type.
int TraverseCombinePST::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = reinterpret_cast<mdl::mdlClass *>(pst->mdl);
    auto rID = mdl->ReqService(pst->idUpper,getServiceID(),vin,nIn);
    Traverse(pst->pstLower,vin,nIn,vout,nOut);
    void *vout2 = alloca(nOut);
    mdl->GetReply(rID,vout2,&nOut);
    Combine(vout,vout2);
    return nOut;
    }
