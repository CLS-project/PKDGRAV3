#ifndef DF755BD8_0CD7_42F4_AE2E_201894E9ECAB
#define DF755BD8_0CD7_42F4_AE2E_201894E9ECAB
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

class ServiceRsSaveIds : public ServiceOutput {
public:
    using input = ServiceOutput::input;
    using output = void;

    explicit ServiceRsSaveIds(PST pst)
        : ServiceOutput(pst,PST_RS_SAVE_IDS,sizeof(input)) {}
protected:
    virtual void Write(PST pst,void *vin,int nIn,int iGroup,const std::string &filename,int iSegment,int nSegment) override;
};


#endif /* DF755BD8_0CD7_42F4_AE2E_201894E9ECAB */
