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
#include "pkd_config.h"
#include "TraversePST.h"

namespace OldDD {
class ServiceColRejects : public TraversePST {
public:
    typedef void input;
    struct output {
        int id;
        int nRejects;
        int nSpace;
        int nLocal;
    };
    explicit ServiceColRejects(PST pst)
        : TraversePST(pst,PST_COLREJECTS,0,mdlThreads(pst->mdl)*sizeof(output),"ColRejects") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceSwapRejects : public TraversePST {
public:
    typedef int input;
    typedef ServiceColRejects::output output;
    explicit ServiceSwapRejects(PST pst)
        : TraversePST(pst,PST_SWAPREJECTS,
                      mdlThreads(pst->mdl)*sizeof(input),
                      mdlThreads(pst->mdl)*sizeof(output),"SwapRejects") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceColOrdRejects : public TraversePST {
public:
    struct input {
        uint64_t iOrdSplit;
        int iSplitSide;
    };
    typedef ServiceColRejects::output output;
    explicit ServiceColOrdRejects(PST pst)
        : TraversePST(pst,PST_COLORDREJECTS,0,mdlThreads(pst->mdl)*sizeof(output),"ColOrdRejects") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceDomain : public TraversePST {
public:
    explicit ServiceDomain(PST node_pst,int service_id,
                           int nInBytes, const char *service_name="")
        : TraversePST(node_pst,service_id, nInBytes, service_name) {}
protected:
    int RejMatch(PST pst,int n1,ServiceColRejects::output *p1,int n2,ServiceColRejects::output *p2,int *pidSwap);
};

class ServiceDomainDecomp : public ServiceDomain {
public:
    struct input {
        Bound bnd;
        int nBndWrap[3];
        int bDoRootFind;
        int bDoSplitDimFind;
        uint64_t nActive;
        uint64_t nTotal;
    };
    typedef void output;
    explicit ServiceDomainDecomp(PST pst)
        : ServiceDomain(pst,PST_DOMAINDECOMP,sizeof(input),"DomainDecomp") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
private:
    void RootSplit(PST pst,int iSplitDim,int bDoRootFind,int bDoSplitDimFind);
};

#ifndef NEW_REORDER
class ServiceDomainOrder : public ServiceDomain {
public:
    struct input {
        uint64_t iMinOrder;
        uint64_t iMaxOrder;
        input() = default;
        input(uint64_t iMinOrder,uint64_t iMaxOrder) : iMinOrder(iMinOrder), iMaxOrder(iMaxOrder) {}
        input(uint64_t iMaxOrder) : iMinOrder(0), iMaxOrder(iMaxOrder) {}
    };
    typedef void output;
    explicit ServiceDomainOrder(PST pst)
        : ServiceDomain(pst,PST_DOMAINORDER,sizeof(input),"DomainOrder") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
private:
    uint64_t OrdSplit(PST pst,uint64_t iMinOrder,uint64_t iMaxOrder);
};

class ServiceLocalOrder : public TraversePST {
public:
    struct input {
        uint64_t iMinOrder;
        uint64_t iMaxOrder;
        input() = default;
        input(uint64_t iMinOrder,uint64_t iMaxOrder) : iMinOrder(iMinOrder), iMaxOrder(iMaxOrder) {}
        input(uint64_t iMaxOrder) : iMinOrder(0), iMaxOrder(iMaxOrder) {}
    };
    typedef void output;
    explicit ServiceLocalOrder(PST pst)
        : TraversePST(pst,PST_LOCALORDER,sizeof(input),"LocalOrder") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};
#endif

class ServiceWeight : public TraversePST {
public:
    struct input {
        int iSplitDim;
        double fSplit;
        int iSplitSide;
        int ittr;
        int pFlag;
    };
    struct output {
        uint64_t nLow;
        uint64_t nHigh;
        double fLow;
        double fHigh;
    };
    explicit ServiceWeight(PST pst)
        : TraversePST(pst,PST_WEIGHT,sizeof(input),sizeof(output),"Weight") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceWeightWrap : public TraversePST {
public:
    struct input {
        int iSplitDim;
        double fSplit;
        double fSplit2;
        int iSplitSide;
        int ittr;
    };
    struct output {
        uint64_t nLow;
        uint64_t nHigh;
    };
    explicit ServiceWeightWrap(PST pst)
        : TraversePST(pst,PST_WEIGHTWRAP,sizeof(input),sizeof(output),"WeightWrap") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

class ServiceOrdWeight : public TraversePST {
public:
    struct input {
        uint64_t iOrdSplit;
        int iSplitSide;
        int ittr;
    };
    struct output {
        uint64_t nLow;
        uint64_t nHigh;
    };
    explicit ServiceOrdWeight(PST pst)
        : TraversePST(pst,PST_ORDWEIGHT,sizeof(input),sizeof(output),"OrdWeight") {}
protected:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};
} // namespace OldDD