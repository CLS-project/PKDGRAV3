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
#ifndef SERVICE_COUNTSPECIES_H
#define SERVICE_COUNTSPECIES_H
#include <cstdint>
#include "TraversePST.h"
#include "io/fio.h"

struct CountSpeciesCounts {
    blitz::TinyVector<std::uint_fast64_t,FIO_SPECIES_LAST> counts;
    std::uint_fast64_t nMaxOrder;
    CountSpeciesCounts() {}
    // In this context += is used to "combine" the results -- not necessarily to add them
    CountSpeciesCounts &operator+=(const CountSpeciesCounts &rhs) {
        counts += rhs.counts;
        nMaxOrder = std::max(nMaxOrder,rhs.nMaxOrder);
        return *this;
    }
};

class ServiceCountSpecies : public TraverseCount<CountSpeciesCounts> {
public:
    typedef void input;
    explicit ServiceCountSpecies(PST pst,int service_id=PST_COUNTSPECIES, int nInBytes=0, const char *service_name="CountSpecies")
        : TraverseCount(pst,service_id,nInBytes,service_name) {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut);
};
#endif