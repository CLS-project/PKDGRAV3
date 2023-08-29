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
#ifndef DOMAINS_REORDER_H
#define DOMAINS_REORDER_H

#include "TraversePST.h"

namespace NewDD {

class ServiceReorder : public TraversePST {
public:
    using dd_offset_type = mdl::mdlClass::dd_offset_type;
    struct input {
        uint64_t nPerProc;
        uint64_t nPerCore;
        uint64_t iMaxOrder;
        input() = default;
        input(uint64_t nPerProc,uint64_t iMaxOrder) : nPerProc(nPerProc), nPerCore(0), iMaxOrder(iMaxOrder) {}
    };
    typedef void output;
    explicit ServiceReorder(PST pst)
        : TraversePST(pst,PST_REORDER,sizeof(input),"Reorder") {}
protected:
    virtual int OffNode(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int  AtNode(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};

} // NewDD

#endif // DOMAINS_REORDER_H
