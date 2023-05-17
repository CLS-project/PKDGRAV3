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

#ifndef DB0087E1_11FE_4040_88EB_FA1695330E03
#define DB0087E1_11FE_4040_88EB_FA1695330E03

#include "TraversePST.h"
#include "io/service.h"

class ServiceRsExtract : public ServiceOutput {
public:
    using header = ServiceOutput::input;
    using input = uint64_t;
    using output = void;
    explicit ServiceRsExtract(PST pst)
        : ServiceOutput(pst,PST_RS_EXTRACT,sizeof(header)+sizeof(input)*(pst->mdl->Threads()+1)) {}
protected:
    virtual void Write(PST pst,void *vin,int nIn,int iGroup,const std::string &filename,int iSegment,int nSegment) override;
    // virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
};


#endif /* DB0087E1_11FE_4040_88EB_FA1695330E03 */
