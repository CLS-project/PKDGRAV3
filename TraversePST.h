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
#ifndef TRAVERSEPST_H
#define TRAVERSEPST_H
#include "pst.h"
#include <string>
#include <type_traits>

// This class will traverse the PST calling Recurse() at each level.
// The Recurse method is responsible for making a request,
// calling and calling Traverse() recursively. When a leaf is reach
// the Service() function is called.
class TraversePST : public mdl::BasicService {
    PST node_pst;
public:
    explicit TraversePST(PST node_pst,int service_id,
                         int nInBytes, int nOutBytes, const char *service_name="")
        : BasicService(service_id, nInBytes, nOutBytes, service_name), node_pst(node_pst) {}
    explicit TraversePST(PST node_pst,int service_id,
                         int nInBytes, const char *service_name="")
        : BasicService(service_id, nInBytes, 0, service_name), node_pst(node_pst) {}
    explicit TraversePST(PST node_pst,int service_id,const char *service_name="")
        : BasicService(service_id, 0, 0, service_name), node_pst(node_pst) {}
    virtual ~TraversePST() = default;
    int ReqService(PST pst,void *vin,int nIn);
protected:
    virtual int operator()(int nIn, void *pIn, void *pOut) final;
    virtual int Traverse(PST pst,void *vin,int nIn,void *vout,int nOut);

    virtual int OffNode(PST pst,void *vin,int nIn,void *vout,int nOut) {return Recurse(pst,vin,nIn,vout,nOut);}
    virtual int  AtNode(PST pst,void *vin,int nIn,void *vout,int nOut) {return Recurse(pst,vin,nIn,vout,nOut);}
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut);
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) = 0;
protected:
    static  int Traverse(unsigned sid, PST pst,void *vin,int nIn,void *vout,int nOut);
};

// This class is for services where the input does not change as the PST is walked,
// but where a customized Combine for a fixed size output is needed.
class TraverseCombinePST : public TraversePST {
public:
    explicit TraverseCombinePST(PST node_pst,int service_id,
                                int nInBytes=0, int nOutBytes=0, const char *service_name="")
        : TraversePST(node_pst,service_id, nInBytes, nOutBytes, service_name) {}
    virtual ~TraverseCombinePST() = default;
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) final;
    virtual int Combine(void *vout,void *vout2) = 0;
};

// This class is for services where the input does not change as the PST is walked,
// but where a field is gathered from each thread.
template<typename OUTPUT>
class TraverseGatherPST : public TraversePST {
public:
    explicit TraverseGatherPST(PST node_pst,int service_id,
                               int nInBytes=0, int nOutBytes=0, const char *service_name="")
        : TraversePST(node_pst,service_id, nInBytes, nOutBytes, service_name) {}
    virtual ~TraverseGatherPST() = default;
    using output = OUTPUT;
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) final {
        auto mdl = reinterpret_cast<mdl::mdlClass *>(pst->mdl);
        auto out1 = static_cast<output *>(vout);
        auto out2 = out1 + pst->nLower;
        int nBytesLower = pst->nLower * sizeof(output);
        int nBytesUpper = pst->nUpper * sizeof(output);
        assert(nBytesLower + nBytesUpper <= nOut);
        nOut = nBytesLower + nBytesUpper;
        auto rID = mdl->ReqService(pst->idUpper,getServiceID(),vin,nIn);
        Traverse(pst->pstLower,vin,nIn,out1,nBytesLower);
        nBytesUpper = mdl->GetReply(rID,out2);
        assert(nBytesLower + nBytesUpper == nOut);
        return nOut;
    }
};

// This function more generally handles the common case where the return
// needs to simply add the lower and upper values. This works for any
// type that provides a "+=" function.
template<class TYPENAME>
class TraverseCount : public TraverseCombinePST {
public:
    typedef TYPENAME output;
    static_assert(std::is_standard_layout<output>());
    explicit TraverseCount(PST pst,int service_id,int nInBytes, const char *service_name="")
        : TraverseCombinePST(pst,service_id,nInBytes,sizeof(output),service_name) {}
    explicit TraverseCount(PST pst,int service_id,const char *service_name="")
        : TraverseCombinePST(pst,service_id,0,sizeof(output),service_name) {}
protected:
    int Combine(void *vout,void *vout2) final {
        auto out  = static_cast<output *>(vout);
        auto out2 = static_cast<output *>(vout2);
        *out += *out2;
        return sizeof(output);
    }
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override = 0;
};



#endif
