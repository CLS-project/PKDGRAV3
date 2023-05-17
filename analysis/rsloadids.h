#ifndef E64EF539_F17B_4CA8_91AA_C22D36E59AC9
#define E64EF539_F17B_4CA8_91AA_C22D36E59AC9
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

class ServiceRsLoadIds : public ServiceInput {
public:
    struct input {
        bool bAppend;
    };
    explicit ServiceRsLoadIds(PST pst)
        : ServiceInput(pst,PST_RS_LOAD_IDS,sizeof(input)) {}
    virtual void Read(PST pst,uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd) override;
    virtual void start(PST pst,uint64_t nElements,void *vin,int nIn) override;
    virtual void finish(PST pst,uint64_t nElements,void *vin,int nIn) override;
};

#endif /* E64EF539_F17B_4CA8_91AA_C22D36E59AC9 */
