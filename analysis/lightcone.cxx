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
#include <cmath>
#include "fmt/format.h"        // This will be part of c++20
using namespace fmt::literals; // Gives us ""_a and ""_format literals

#include "lightcone.h"
#include "pkd.h"
#include "core/simd.h"
extern "C" {
#include "core/healpix.h"
}
#include "core/abort_error.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceInitLightcone::input> || std::is_standard_layout_v<ServiceInitLightcone::input>);
static_assert(std::is_void_v<ServiceInitLightcone::output> || std::is_standard_layout_v<ServiceInitLightcone::output>);

int ServiceInitLightcone::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    abort_assert(1002, vin.size() == sizeof(input));
    pkdInitLightcone(pst->lcl()->pkd, in->bBowtie, in->bLightConeParticles, in->dBoxSize, in->dRedshiftLCP, in->alphaLCP, in->hLCP, in->mOctants);
    return 0;
}

int ServiceLightConeOpen::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    pst->lcl()->pkd->LightConeOpen(in->achOutFileLcp, in->nSideHealpix);
    return 0;
}

int ServiceLightConeClose::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    pst->lcl()->pkd->LightConeClose(in->achOutFileHpb);
    return 0;
}

int ServiceLightConeVel::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    pst->lcl()->pkd->LightConeVel(in->dBoxSize);
    return 0;
}

#define FILE_PROTECTION (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP)

void pkdContext::flushLightCone(int octant) {
    size_t count = nLightCone[octant] * sizeof(LIGHTCONEP);
    io_write(&afiLightCone[octant], pLightCone[octant], count);
    nLightCone[octant] = 0;
}

class HealpixHelper : public mdl::arrayCACHEhelper {
    virtual void init(void* dst) override final {
        auto m = static_cast<healpixData*>(dst);
        m->nGrouped = 0;
        m->nUngrouped = 0;
        m->fPotential = 0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto m1 = static_cast<healpixData*>(dst);
        auto m2 = static_cast<const healpixData*>(src);
        auto SumWithSaturate = [](uint32_t a, uint32_t b) -> uint32_t {
            uint64_t sum = static_cast<uint64_t>(a) + b;
            return (sum > 0xffffffffu) ? 0xffffffffu : static_cast<uint32_t>(sum);
        };
        m1->nGrouped = SumWithSaturate(m1->nGrouped, m2->nGrouped);
        m1->nUngrouped = SumWithSaturate(m1->nUngrouped, m2->nUngrouped);
        m1->fPotential += m2->fPotential;
    }

  public:
    explicit HealpixHelper(PKD pkd) : arrayCACHEhelper(pkd->pHealpixData, pkd->nHealpixPerDomain, true) {}
};

void pkdContext::LightConeClose(const char* healpixname_template) {
    for (auto octant = 0; octant < 8; ++octant) {
        if (afiLightCone[octant].fd > 0) {
            flushLightCone(octant);
            io_close(&afiLightCone[octant]);
        }
    }
    if (nSideHealpix) {
        assert(healpixname_template && healpixname_template[0]);
        auto healpixname = fmt::format(healpixname_template, fmt::arg("i", mdl->Self()));

        mdl->FinishCache(CID_HEALPIX);
        int fd = open(healpixname.c_str(), O_CREAT | O_WRONLY | O_TRUNC, FILE_PROTECTION);
        if (fd < 0) {
            perror(healpixname.c_str());
            abort();
        }
        for (auto i = 0; i < nHealpixPerDomain; ++i) {
            uint64_t sum = pHealpixData[i].nGrouped;
            sum += pHealpixData[i].nUngrouped;
            if (sum) pHealpixData[i].fPotential /= sum;
        }
        auto nWrite = nHealpixPerDomain * sizeof(*pHealpixData);
        if (write(fd, pHealpixData, nWrite) != nWrite) {
            perror("Wrong size writing Healpix");
        }
        close(fd);
    }
}

void pkdContext::LightConeOpen(const char* name_template, int nSideHealpix) {
    auto mask = (mLightConeOctants == 0) ? 1 : mLightConeOctants;
    for (auto octant = 0; octant < 8; ++octant) {
        if ((mask & (1 << octant)) == 0 || name_template[0] == '\0') {
            afiLightCone[octant].fd = -1;
            continue;
        }
        auto fname = fmt::format(name_template, "octant"_a = octant, "i"_a = mdl->Self());
        if (io_create(&afiLightCone[octant], fname.c_str()) < 0) {
            perror(fname.c_str());
            abort();
        }
    }

    this->nSideHealpix = nSideHealpix;
    if (nSideHealpix) {
        nHealpixPerDomain = (nside2npix64(nSideHealpix) + mdl->Threads() - 1) / mdl->Threads();
        nHealpixPerDomain = (nHealpixPerDomain + 15) & ~15;
        if (pHealpixData == NULL) {
            pHealpixData = static_cast<healpixData*>(malloc(nHealpixPerDomain * sizeof(*pHealpixData)));
            assert(pHealpixData != NULL);
        }
        for (auto i = 0; i < nHealpixPerDomain; ++i) {
            pHealpixData[i].nGrouped = 0;
            pHealpixData[i].nUngrouped = 0;
            pHealpixData[i].fPotential = 0;
        }
        mdl->CacheInitialize(CID_HEALPIX, std::make_shared<HealpixHelper>(this));
    }
}

void pkdContext::addToLightCone(double dvFac, double* r, float fPot, PARTICLE* p, int bParticleOutput) {
    auto P = particles[p];
    const auto& v = P.velocity();
    if (bParticleOutput) {
        auto octant = mLightConeOctants == 0 ? 0 : (((r[0] < 0) ? 4 : 0) + ((r[1] < 0) ? 2 : 0) + ((r[2] < 0) ? 1 : 0));
        if (afiLightCone[octant].fd > 0) {
            auto pLC = pLightCone[octant];
            auto n = nLightCone[octant];
            pLC[n].id = P.order();
            pLC[n].pos[0] = r[0];
            pLC[n].pos[1] = r[1];
            pLC[n].pos[2] = r[2];
            /*
            ** Convert momentum (a^2 * x_dot to physical peculiar velocity (a * x_dot in sim units)
            ** This is easy to do now, and somewhat more involved to do later.
            */
            pLC[n].vel[0] = v[0] * dvFac;
            pLC[n].vel[1] = v[1] * dvFac;
            pLC[n].vel[2] = v[2] * dvFac;
            pLC[n].pot = fPot;
            if (++nLightCone[octant] == nLightConeMax) flushLightCone(octant);
        }
    }
    if (nSideHealpix) {
        int64_t iPixel = vec2pix_ring64(nSideHealpix, r);
        assert(iPixel >= 0);
        int id = iPixel / nHealpixPerDomain;
        int idx = iPixel - id * nHealpixPerDomain;
        assert(id < mdl->Threads());
        assert(idx < nHealpixPerDomain);
        auto m = mdl->VirtualFetch<healpixData>(CID_HEALPIX, idx, id);
        if (P.group()) {
            if (m->nGrouped < 0xffffffffu) ++m->nGrouped; /* Increment with saturate */
        }
        else {
            if (m->nUngrouped < 0xffffffffu) ++m->nUngrouped; /* Increment with saturate */
        }
        m->fPotential += fPot;
    }
}

#define NBOX 184

void pkdContext::ProcessLightCone(PARTICLE* p, float fPot, double dLookbackFac, double dLookbackFacLCP, double dDriftDelta, double dKickDelta, double dBoxSize,
                                  int bLightConeParticles, blitz::TinyVector<double, 3> hlcp, double tanalpha_2) {
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    const double mrLCP = dLightSpeed * dLookbackFacLCP;
    const double depth = dLightSpeed * dLookbackFac;
    int bCone = 1;
    const int nLayerMax = ceil(mrLCP);
    double dxStart;
    auto P = particles[p];

    /*
    ** dxStart measures the fraction of the light surface that *just* enters
    ** a certain layer of replicas. If it lies well outside of this in the time step
    ** (i.e. it is greater than 1) then we don't need to check it (3 replicas was the
    ** max before). If it is negative then the light surface lies inside of 3
    ** replicas being checked during the whole time step interval, the xStart will
    ** always be 0 (we check the entire light surface progression during the entire
    ** time interval given by dKickDelta.
    */
    assert(depth >= 0);
    //    if (depth < 0) return; // timestep extends into the future, it seems.
    int l = floor(depth);
    if (l >= nLayerMax) l = nLayerMax - 1;
    int nBox = nBoxLC[l]; // check only this many boxes
    dxStart = (depth - nLayerMax) / (dKickDelta * dLightSpeed);
    if (dxStart > 1) return; // the timestep is still too deep!
    if (dxStart < 0) dxStart = 0;

    const auto& v = P.velocity();
    int j;

    auto r0 = P.position();
    for (j = 0; j < 3; ++j) {
        if (r0[j] < -0.5) r0[j] += 1.0;
        else if (r0[j] >= 0.5) r0[j] -= 1.0;
    }

    struct {
        double dt;
        double fOffset;
        int jPlane;
    } isect[4], temp;
    /*
    ** Figure out at what delta t the particle would cross each of the
    ** walls of the unit cell. Usually this delta t will be larger than
    ** the actual drift, but if not we need to treat each segment of the
    ** drift that is created by crossing a wall (there can be at most 4
    ** segments in 3-D).
    */
    for (j = 0; j < 3; ++j) {
        if (v[j] > 0) {
            /*
            ** Particle crosses the upper j-coordinate boundary of the unit cell at isect[j].dt.
            */
            isect[j].dt = (0.5 - r0[j]) / v[j];
            isect[j].fOffset = -1.0;
        }
        else if (v[j] < 0) {
            /*
            ** Particle crosses the lower j-coordinate boundary of the unit cell at isect[j].dt.
            */
            isect[j].dt = (-0.5 - r0[j]) / v[j];
            isect[j].fOffset = 1.0;
        }
        else {
            /*
            ** Particle is not moving in this dimension!
            */
            isect[j].dt = 2 * dDriftDelta; /* this will be ignored after the sort */
            isect[j].fOffset = 0.0;
        }
        assert(isect[j].dt >= 0);
        isect[j].jPlane = j;
    }
    /*
    ** Last option is that the particle does not cross any wall so we make an
    ** entry which contains just the drift of the particle.
    */
    isect[3].dt = dDriftDelta;
    isect[3].fOffset = 0.0;
    isect[3].jPlane = 3;
    /*
    ** Sort them by drift time, we will deal with each segment of the drift in this
    ** order.
    */
    if (isect[0].dt > isect[1].dt) {
        temp = isect[0];
        isect[0] = isect[1];
        isect[1] = temp;
    }
    if (isect[2].dt > isect[3].dt) {
        temp = isect[2];
        isect[2] = isect[3];
        isect[3] = temp;
    }
    temp = isect[1];
    isect[1] = isect[2];
    isect[2] = temp;
    if (isect[0].dt > isect[1].dt) {
        temp = isect[0];
        isect[0] = isect[1];
        isect[1] = temp;
    }
    if (isect[2].dt > isect[3].dt) {
        temp = isect[2];
        isect[2] = isect[3];
        isect[3] = temp;
    }
    if (isect[1].dt > isect[2].dt) {
        temp = isect[1];
        isect[1] = isect[2];
        isect[2] = temp;
    }

    dvec xStart = dxStart;
    dvec h[3];
    double ihm = 1.0 / sqrt(hlcp[0] * hlcp[0] + hlcp[1] * hlcp[1] + hlcp[2] * hlcp[2]); // make sure it is a unit vector!
    for (j = 0; j < 3; ++j) h[j] = hlcp[j] * ihm;                                       // normalize the hlcp unit vector
    /*
    ** If the input tan of alpha/2 is given as negative then we want all sky.
    */
    if (tanalpha_2 < 0) bCone = 0;
    else bCone = 1;
    dvec tan2alpha_2 = tanalpha_2 * tanalpha_2;
    nBox = (nBox + dvec::width() - 1) / dvec::width();

    int k;
    for (k = 0; k < 4; ++k) {
        double dtApprox, dt;
        dvec dlbt;
        assert(isect[k].dt <= dDriftDelta);
        if (k == 0) {
            /*
            ** Check lightcone from 0 <= dt < isect[k].dt
            */
            dt = isect[k].dt;
            dtApprox = dt / dDriftDelta * dKickDelta;
            dlbt = dLookbackFac;
        }
        else {
            /*
            ** Check lightcone from isect[k-1].dt <= dt < isect[k].dt
            */
            dt = isect[k].dt - isect[k - 1].dt;
            dtApprox = dt / dDriftDelta * dKickDelta;
            /*
            ** continue the light surface from the end of the previous segment
            */
            dlbt = dLookbackFac - isect[k - 1].dt / dDriftDelta * dKickDelta;
        }
        dvec t0 = dlbt * dlbt * dLightSpeed * dLightSpeed;
        dvec t1 = (dlbt - dtApprox) * (dlbt - dtApprox) * dLightSpeed * dLightSpeed;
        auto r1 = r0 + dt * v;
        for (int iOct = 0; iOct < nBox; ++iOct) {
            dvec off0, off1, off2;
            off0.load(&lcOffset0[iOct * dvec::width()]);
            off1.load(&lcOffset1[iOct * dvec::width()]);
            off2.load(&lcOffset2[iOct * dvec::width()]);
            dvec vrx0 = off0 + r0[0];
            dvec vry0 = off1 + r0[1];
            dvec vrz0 = off2 + r0[2];
            dvec vrx1 = off0 + r1[0];
            dvec vry1 = off1 + r1[1];
            dvec vrz1 = off2 + r1[2];
            dvec mr0 = vrx0 * vrx0 + vry0 * vry0 + vrz0 * vrz0;
            dvec mr1 = vrx1 * vrx1 + vry1 * vry1 + vrz1 * vrz1;
            dmask msk = (t1 <= max(mr0, mr1)) & (t0 >= min(mr0, mr1));
            if (!testz(msk)) {
                mr0 = sqrt(mr0);
                mr1 = sqrt(mr1);
                dvec vx = (dLightSpeed * dlbt - mr0) / (dLightSpeed * dtApprox - mr0 + mr1);
                msk = (vx >= xStart) & (vx < 1.0);
                if (!testz(msk)) {
                    dvec vr[3];
                    vr[0] = (1.0 - vx) * vrx0 + vx * vrx1;
                    vr[1] = (1.0 - vx) * vry0 + vx * vry1;
                    vr[2] = (1.0 - vx) * vrz0 + vx * vrz1;
                    if (bCone) {
                        /*
                        ** Now here we test for inclusion into a cone.
                        ** We need to have the unit vector of vr for this to calculate
                        ** the tangent of half the angle in the cone. For the full
                        ** sky lightcone we can skip this test.
                        ** (h is the direction unit vector of the cone)
                        */
                        dvec vrm = vr[0] * vr[0] + vr[1] * vr[1] + vr[2] * vr[2];
                        vrm = sqrt(vrm); // hopefully this is done as a vector operation
                        dvec a[3];       // difference between the 2 unit vectors
                        a[0] = vr[0] / vrm - h[0];
                        a[1] = vr[1] / vrm - h[1];
                        a[2] = vr[2] / vrm - h[2];
                        dvec am2 = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
                        dvec b[3]; // sum the 2 unit vectors
                        b[0] = vr[0] / vrm + h[0];
                        b[1] = vr[1] / vrm + h[1];
                        b[2] = vr[2] / vrm + h[2];
                        dvec bm2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
                        dvec tan2 = min(am2, bm2) / max(am2, bm2); // we have computed the tan-squared of the half angle
                        msk = msk & (tan2 < tan2alpha_2);          // we need to test if it falls in the thin shell, and in the correct angular region.
                    }
                    if (!testz(msk)) { // it tests twice the same thing for the full sky lightcone, sorry.
                        int m = movemask(msk);
                        for (int j = 0; j < dvec::width(); ++j) {
                            if (m & (1 << j)) {
                                double r[3];
                                r[0] = vr[0][j];
                                r[1] = vr[1][j];
                                r[2] = vr[2][j];
                                /*
                                ** Create a new light cone particle.
                                */
                                double mr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
                                /*
                                ** Lookup the expansion factor to convert the velocities into physical peculiar
                                ** in sim units. Input velocities are momenta p = a^2*x_dot and we want
                                ** v_pec = a*x_dot.
                                ** Use r -> 1/a spline table (but only if mr is in the interpolation domain!).
                                */
                                if (mr < mrLCP) {
                                    double dvFac = gsl_spline_eval(interp_scale, mr, interp_accel);
                                    addToLightCone(dvFac, r, fPot, p, bLightConeParticles);
                                }
                            }
                        }
                    } /* end of (!testz(msk)) */
                }
            }
        }
        if (isect[k].jPlane == 3) break;
        /*
        ** Now we need to reposition r0 to the new segment.
        */
        for (j = 0; j < 3; ++j) r0[j] = r1[j];
        r0[isect[k].jPlane] += isect[k].fOffset;
    }
}

void pkdContext::LightConeVel(double dBoxSize) {
    const int nTable = 1000;
    const double rMax = 3.0;
    gsl_spline* scale;
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    double dr, rt[nTable], at_inv[nTable];
    double dvFac;
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    int i;

    assert(particles.present(PKD_FIELD::oVelocity));
    /*
    ** Setup lookup table.
    */
    dr = rMax / (nTable - 1);
    for (i = 0; i < nTable; ++i) {
        rt[i] = i * dr;
        at_inv[i] = 1.0 / csmComoveLookbackTime2Exp(csm, rt[i] / dLightSpeed);
    }
    scale = gsl_spline_alloc(gsl_interp_cspline, nTable);
    gsl_spline_init(scale, rt, at_inv, nTable);
    /*
    ** Now loop over all particles using this table.
    */
    for (auto& p : particles) {
        auto& v = p.velocity();
        auto r = p.position();
        auto r2 = dot(r, r);
        /*
        ** Use r -> 1/a spline table.
        */
        dvFac = gsl_spline_eval(scale, sqrt(r2), acc);
        /* input velocities are momenta p = a^2*x_dot and we want v_pec = a*x_dot */
        v *= dvFac;
    }
    gsl_spline_free(scale);
    gsl_interp_accel_free(acc);
}
