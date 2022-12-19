#ifndef BD2FB51E_EAF8_41BC_B8A5_6333435EB545
#define BD2FB51E_EAF8_41BC_B8A5_6333435EB545
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
#include "io/service.h"

class ServiceRsHaloLoadIds : public ServiceInput {
public:
    struct input {
        bool bAppend;
    };
    explicit ServiceRsHaloLoadIds(PST pst)
        : ServiceInput(pst,PST_RS_HALO_LOAD_IDS,sizeof(input)) {}
    virtual void Read(PST pst,uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd) override;
    virtual void start(PST pst,uint64_t nElements,void *vin,int nIn) override;
    virtual void finish(PST pst,uint64_t nElements,void *vin,int nIn) override;
};

#endif /* BD2FB51E_EAF8_41BC_B8A5_6333435EB545 */
