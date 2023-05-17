#ifndef FD041C3A_4183_478E_A8E7_CA209613BB18
#define FD041C3A_4183_478E_A8E7_CA209613BB18
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

class ServiceRsHaloCount : public ServiceFileSizes {
public:
    explicit ServiceRsHaloCount(PST pst)
        : ServiceFileSizes(pst,PST_RS_HALO_COUNT,"RsHaloCount") {}
    virtual uint64_t GetSize(const std::string &filename,uint64_t file_size) override;
};

#endif /* FD041C3A_4183_478E_A8E7_CA209613BB18 */
