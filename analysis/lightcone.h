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
#ifndef LIGHTCONE_H
#define LIGHTCONE_H
#include "TraversePST.h"
#include <array>

class ServiceInitLightcone : public TraversePST {
  public:
    static constexpr std::string_view name = "InitLightcone";
    struct input {
        int bBowtie;
        int bLightConeParticles;
        int mOctants;
        double dBoxSize;
        double dRedshiftLCP;
        double alphaLCP;
        blitz::TinyVector<double, 3> hLCP;
    };
    using output = void;
    explicit ServiceInitLightcone(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceLightConeOpen : public TraversePST {
  public:
    static constexpr std::string_view name = "LightConeOpen";
    struct input {
        int nSideHealpix;
        char achOutFileLcp[PATH_MAX];
    };
    using output = void;
    explicit ServiceLightConeOpen(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceLightConeClose : public TraversePST {
  public:
    static constexpr std::string_view name = "LightConeClose";
    struct input {
        char achOutFileHpb[PATH_MAX];
    };
    using output = void;
    explicit ServiceLightConeClose(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceLightConeVel : public TraversePST {
  public:
    static constexpr std::string_view name = "LightConeVel";
    struct input {
        double dBoxSize;
    };
    using output = void;
    explicit ServiceLightConeVel(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

#endif // LIGHTCONE_H