#ifndef E5B60112_B4E2_46CB_AC8A_564BF342E5D3
#define E5B60112_B4E2_46CB_AC8A_564BF342E5D3

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

class ServiceRestore : public ServiceInput {
public:
    explicit ServiceRestore(PST pst)
        : ServiceInput(pst,PST_RESTORE) {}
    virtual void Read(PST pst,uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd) override;
};

#endif /* E5B60112_B4E2_46CB_AC8A_564BF342E5D3 */
