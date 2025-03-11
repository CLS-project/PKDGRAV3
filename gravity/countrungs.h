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
#ifndef SERVICE_COUNTRUNGS_H
#define SERVICE_COUNTRUNGS_H
#include "TraversePST.h"

class ServiceCountRungs : public TraverseCount<blitz::TinyVector<uint64_t,MAX_RUNG+1>> {
public:
    typedef void input;
    explicit ServiceCountRungs(PST pst,int service_id=PST_COUNTRUNGS, int nInBytes=0, const char *service_name="CountRungs")
        : TraverseCount(pst,service_id,nInBytes,service_name) {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut);
};
#endif