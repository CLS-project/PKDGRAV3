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
#ifndef SERVICE_CACHESTATS_H
#define SERVICE_CACHESTATS_H
#include <cstdint>
#include "TraversePST.h"
#include "mdl.h"
#include "mdl2/WelfordStats.h"

class ServiceCacheStats : public TraverseCount<ServiceCacheStats> {
  public:
    static constexpr std::string_view name = "CacheStats";
    using input = void;
    struct output {
        WelfordStats<double> mpi_loop;
        WelfordStats<double> mpi_test;
        WelfordStats<double> gpu_test;
        WelfordStats<double> msg_test;
        // In this context += is used to "combine" the results -- not necessarily to add them
        auto operator+=(const output& rhs) -> output& {
            mpi_loop += rhs.mpi_loop;
            mpi_test += rhs.mpi_test;
            gpu_test += rhs.gpu_test;
            msg_test += rhs.msg_test;
            return *this;
        }
    };

    explicit ServiceCacheStats(PST pst, std::string_view service_name = name, int nInBytes = 0) : TraverseCount(pst, service_name, nInBytes) {}

  protected:
    // We are collecting MPI statistics at each node, so we use AtNode to do the work
    virtual int AtNode(PST pst, std::span<std::byte> in, std::span<std::byte> out) { return Service(pst, in, out); }
    virtual int Service(PST pst, std::span<std::byte> in, std::span<std::byte> out);
};
#endif