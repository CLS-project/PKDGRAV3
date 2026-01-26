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

#ifndef HOP_H
#define HOP_H
#include "TraversePST.h"

class ServiceHopLink : public TraverseCount<ServiceHopLink> {
  public:
    static constexpr std::string_view name = "HopLink";
    struct input {
        int nSmooth;
        bool bPeriodic;
    };
    using output = uint64_t;
    explicit ServiceHopLink(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceHopJoin : public TraverseCount<ServiceHopJoin> {
  public:
    static constexpr std::string_view name = "HopJoin";
    struct input {
        double dHopTau;
        int nSmooth;
        bool bPeriodic;
    };
    struct output {
        uint64_t nGroups;
        bool bDone;
        auto operator+=(const output& rhs) -> output& {
            nGroups += rhs.nGroups;
            bDone = bDone && rhs.bDone;
            return *this;
        }
    };
    explicit ServiceHopJoin(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceHopFinishUp : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "HopFinishUp";
    struct input {
        blitz::TinyVector<double, 3> fPeriod;
        int nMinGroupSize;
        bool bPeriodic;
    };
    using output = uint64_t;
    explicit ServiceHopFinishUp(PST pst) : TraverseCombinePST(pst, name, sizeof(input), sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST) -> int override;
};

class ServiceHopTreeBuild : public TraversePST {
  public:
    static constexpr std::string_view name = "HopTreeBuild";
    struct input {
        int nBucket;
        int nGroup;
    };
    using output = void;
    explicit ServiceHopTreeBuild(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceHopGravity : public TraversePST {
  public:
    static constexpr std::string_view name = "HopGravity";
    struct input {
        double dTime;
        double dTheta;
    };
    using output = void;
    explicit ServiceHopGravity(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceHopUnbind : public TraverseCount<ServiceHopJoin> {
  public:
    static constexpr std::string_view name = "HopUnbind";
    struct input {
        double dTime;
        blitz::TinyVector<double, 3> fPeriod;
        int bPeriodic;
        int nMinGroupSize;
    };
    struct output {
        uint64_t nEvaporated;
        uint64_t nGroups;
        auto operator+=(const output& rhs) -> output& {
            nEvaporated += rhs.nEvaporated;
            nGroups += rhs.nGroups;
            return *this;
        }
    };
    explicit ServiceHopUnbind(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

#endif
