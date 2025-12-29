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

#ifndef FOF_H
#define FOF_H
#include "TraversePST.h"

class ServiceNewFof : public TraversePST {
  public:
    static constexpr std::string_view name = "NewFof";
    struct input {
        double dTau2;
        int nMinMembers;
        int bPeriodic;
        int nReplicas;
        int nBucket;
    };
    using output = void;
    explicit ServiceNewFof(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceFofPhases : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "FofPhases";
    using input = void;
    struct output {
        bool bMadeProgress;
    };
    explicit ServiceFofPhases(PST pst) : TraverseCombinePST(pst, name, 0, sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout, std::span<const std::byte> vin, PST) -> int override;
};

class ServiceFofFinishUp : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "FofFinishUp";
    struct input {
        int nMinGroupSize;
    };
    using output = uint64_t;
    explicit ServiceFofFinishUp(PST pst) : TraverseCombinePST(pst, name, sizeof(input), sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout, std::span<const std::byte> vin, PST) -> int override;
};

#endif
