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
using blitz::TinyVector;

auto ServiceMaxRung::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
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
    auto pkd = pst->pkd();
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
    auto pkd = pst->pkd();
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
    auto pkd = pst->pkd();
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

auto ServiceKick::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->pkd();
    pkd->Kick(in->dTime, in->dDelta, in->bDoGas, in->dDeltaVPred, in->dDeltaU, in->dDeltaUPred, in->uRungLo, in->uRungHi);
    return 0;
}

void pkdContext::Kick(double dTime, double dDelta, int bDoGas, double dDeltaVPred, double dDeltaU, double dDeltaUPred, uint8_t uRungLo, uint8_t uRungHi) {
    assert(particles.present(PKD_FIELD::oVelocity));
    assert(particles.present(PKD_FIELD::oAcceleration));

    if (bDoGas) { // pkd->param.bMeshlessHydro
        for (auto& p : particles) {
            if (!p.is_gas()) {
                if (p.is_rung_range(uRungLo, uRungHi)) {
                    auto& a = p.acceleration();
                    auto& v = p.velocity();
                    v += a * dDelta;
                }
            }
        }
    }
    else {
        for (auto& p : particles) {
            if (p.is_rung_range(uRungLo, uRungHi)) {
                auto& a = p.acceleration();
                auto& v = p.velocity();
                v += a * dDelta;
            }
        }
    }
}

auto ServiceDrift::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->pkd();
    pkd->Drift(in->iRoot, in->dTime, in->dDelta, in->dDeltaVPred, in->dDeltaUPred, in->bDoGas, in->bGasEvolveDensity);
    return 0;
}

/*
** Drift particles whose Rung falls between uRungLo (large step) and uRungHi (small step) inclusive,
** and those whose destination activity flag is set.
**
** Note that the drift funtion no longer wraps the particles around the periodic "unit" cell. This is
** now done by Domain Decomposition only.
*/
void pkdContext::Drift(int iRoot, double dTime, double dDelta, double dDeltaVPred, double dDeltaTime, int bDoGas, int bGasEvolveDensity) {
    int i;
    TinyVector<double, 3> rfinal, r0, dr;
    int pLower, pUpper;

    if (iRoot >= 0) {
        auto pRoot = tree[iRoot];
        pLower = pRoot->lower();
        pUpper = pRoot->upper();
    }
    else {
        pLower = 0;
        pUpper = Local();
    }

    assert(particles.present(PKD_FIELD::oVelocity));
    auto dMin = bnd.lower();
    auto dMax = bnd.upper();
    /*
    ** Update particle positions
    */
    if (bDoGas) {
        if (decorators.present(DECORATOR_FIELD::oMFM)) {
            assert(particles.present(PKD_FIELD::oAcceleration));
            for (i = pLower; i <= pUpper; ++i) {
                auto p = particles[i];
                const auto& v = p.velocity();

                if (p.is_gas()) {
                    // As for gas particles we use dx/dt = v/a, we must use the
                    // Kick factor provided by pkdgrav.
                    // See Stadel 2001, Appendix 3.8
                    dr = v * dDeltaVPred;
                }
                else {
                    dr = v * dDelta;
                }

#ifdef FORCE_1D
                dr[2] = 0.0;
                dr[1] = 0.0;
#endif
#ifdef FORCE_2D
                dr[2] = 0.0;
#endif
                r0 = p.position();
                p.set_position(rfinal = r0 + dr);
                assert(isfinite(rfinal[0]));
                assert(isfinite(rfinal[1]));
                assert(isfinite(rfinal[2]));
                dMin = min(dMin, rfinal);
                dMax = max(dMax, rfinal);
            }
        }
        else {
            assert(particles.present(PKD_FIELD::oNewSph));
            for (i = pLower; i <= pUpper; ++i) {
                auto p = particles[i];
                const auto& v = p.velocity();
                if (bGasEvolveDensity && p.is_gas()) {
                    auto& NewSph = p.newsph();
                    // float dfBalldt = 1.0f / 3.0f * p.ball() * p.density() * NewSph.divv;
                    // p.set_ball(p.ball() + dDelta * dfBalldt);
                    p.set_density(p.density() - dDelta * p.density() * NewSph.divv);
                }

                r0 = p.position();
                p.set_position(rfinal = r0 + dDelta * v);
                assert(isfinite(rfinal[0]));
                assert(isfinite(rfinal[1]));
                assert(isfinite(rfinal[2]));
                dMin = min(dMin, rfinal);
                dMax = max(dMax, rfinal);
            }
        }
    }
    else {
        for (i = pLower; i <= pUpper; ++i) {
            auto p = particles[i];
            const auto& v = p.velocity();
            r0 = p.position();
            p.set_position(rfinal = r0 + dDelta * v);
            dMin = min(dMin, rfinal);
            dMax = max(dMax, rfinal);
        }
    }
    bnd = Bound(dMin, dMax);
}
