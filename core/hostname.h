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
#include "TraversePST.h"

class ServiceHostname : public TraversePST {
  public:
    static constexpr std::string_view name = "Hostname";
    typedef void input;
    struct output {
        int iMpiID;
        char szHostname[20];
    };
    explicit ServiceHostname(PST pst) : TraversePST(pst, name, 0, pst->mdl->Threads() * sizeof(output)) {}

  protected:
    virtual int Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) override;
    virtual int Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) override;
};
