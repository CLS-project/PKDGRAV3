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
#include "pkd_config.h"
#include "timestep.h"
#include <algorithm>
#include <numeric>

auto ServiceMaxRung::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->lcl()->pkd;
    auto& out = *reinterpret_cast<output*>(vout.data());
    out.uMaxRung = pkd->MaxRung();
    return sizeof(output);
}

auto pkdContext::MaxRung() -> uint8_t {
    if (particles.empty()) return 0;
    // return std::ranges::max(particles | std::views::transform([](const auto& p) { return p.new_rung(); }));
    return std::accumulate(particles.begin(), particles.end(), uint8_t{0}, [](uint8_t a, auto& p) { return std::max(a, p.new_rung()); });
}

auto ServiceSetGlobalRung::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->lcl()->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    pkd->SetGlobalRung(in->uMaxRung);
    return 0;
}

void pkdContext::SetGlobalRung(uint8_t minDt) {
    for (auto& p : particles) p.set_new_rung(minDt);
}

static inline auto DtToRung(double dT, double dDelta, uint8_t uMaxRung) -> uint8_t {
    union {
        double d;
        struct {
            uint64_t mantisa : 52;
            uint64_t exponent : 11;
            uint64_t sign : 1;
        } ieee;
    } T;
    T.d = fabs(dDelta) / dT;
    if (T.d <= 1.0) return 0;
    int iRung = T.ieee.exponent - 1023; /* log2(d) */
    if (iRung > uMaxRung) return uMaxRung;
    else return iRung;
}

auto ServiceAccelStep::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->lcl()->pkd;
    pkd->AccelStep(in->uRungLo, in->uRungHi, in->dDelta, in->iMaxRung, in->dEta, in->dVelFac, in->dAccFac, in->bDoGravity, in->bEpsAcc, in->dhMinOverSoft);
    return 0;
}

void pkdContext::AccelStep(uint8_t uRungLo, uint8_t uRungHi, double dDelta, int iMaxRung, double dEta, double dVelFac, double dAccFac, int bDoGravity,
                           int bEpsAcc, double dhMinOverSoft) {
    assert(particles.present(PKD_FIELD::oVelocity));
    assert(particles.present(PKD_FIELD::oAcceleration));
    assert(!bNoParticleOrder);

    for (auto& p : particles) {
        if (p.is_active()) {
            const auto& v = p.velocity();
            const auto& a = p.acceleration();
            double fSoft = p.soft();
            double vel = dot(v, v);
            double acc = dot(a, a);
            assert(vel >= 0);
            vel = sqrt(vel) * dVelFac;
            assert(acc >= 0);
            acc = sqrt(acc) * dAccFac;
            double dT = (acc > 0 && bEpsAcc) ? dEta * sqrt(fSoft / acc) : FLOAT_MAXVAL;
            uint8_t uNewRung = DtToRung(dT, dDelta, iMaxRung);
            if (uNewRung > p.new_rung()) p.set_new_rung(uNewRung);
        }
    }
}

auto ServiceDensityStep::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->lcl()->pkd;
    pkd->DensityStep(in->uRungLo, in->uRungHi, in->iMaxRung, in->dDelta, in->dEta, in->dRhoFac);
    return 0;
}

void pkdContext::DensityStep(uint8_t uRungLo, uint8_t uRungHi, int iMaxRung, double dDelta, double dEta, double dRhoFac) {
    assert(!bNoParticleOrder);
    for (auto& p : particles) {
        if (p.is_active()) {
            double dT = dEta / sqrt(p.density() * dRhoFac);
            p.set_new_rung(DtToRung(dT, dDelta, iMaxRung));
        }
    }
}
