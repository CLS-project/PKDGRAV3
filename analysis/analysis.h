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
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "TraversePST.h"

class ServiceCountDistance : public TraverseCount<ServiceCountDistance> {
  public:
    static constexpr std::string_view name = "CountDistance";
    struct input {
        double dRadius2Inner;
        double dRadius2Outer;
    };
    using output = uint64_t;
    explicit ServiceCountDistance(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceCalcDistance : public TraversePST {
  public:
    static constexpr std::string_view name = "CalcDistance";
    struct input {
        blitz::TinyVector<double, 3> dCenter;
        double dRadius;
        int bPeriodic;
    };
    using output = void;
    explicit ServiceCalcDistance(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceCalcCOM : public TraversePST {
  public:
    static constexpr std::string_view name = "CalcCOM";
    struct input {
        blitz::TinyVector<double, 3> dCenter;
        double dRadius;
        int bPeriodic;
    };
    struct output {
        blitz::TinyVector<double, 3> com;
        blitz::TinyVector<double, 3> vcm;
        blitz::TinyVector<double, 3> L;
        double M;
        uint64_t N;
    };
    explicit ServiceCalcCOM(PST pst) : TraversePST(pst, name, sizeof(input), sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceTotalMass : public TraverseCount<ServiceTotalMass> {
  public:
    static constexpr std::string_view name = "TotalMass";
    using input = void;
    using output = double;
    explicit ServiceTotalMass(PST pst) : TraverseCount(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceProfile : public TraversePST {
  public:
    static constexpr std::string_view name = "Profile";
    static constexpr int max_profile_bins = 1000000;
    struct input {
        blitz::TinyVector<double, 3> dCenter;
        blitz::TinyVector<double, 3> com;
        blitz::TinyVector<double, 3> vcm;
        blitz::TinyVector<double, 3> L;
        uint32_t nBins;
        uint8_t uRungLo;
        uint8_t uRungHi;
        std::array<double, max_profile_bins> dRadii;
    };
    using output = void;
    explicit ServiceProfile(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

#endif