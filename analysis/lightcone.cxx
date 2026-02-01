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

using blitz::TinyVector;
using LightConeInfo = pkdContext::LightConeInfo;

pkdContext::LightConeInfo::LightConeInfo() : pLightCone(static_cast<LIGHTCONEP*>(nullptr)) {
    for (auto& afi : afiLightCone) io_init(&afi, 0, 0, IO_AIO | IO_LIBAIO);
}
pkdContext::LightConeInfo::~LightConeInfo() {
    for (auto* lco : pLightCone) {
        if (lco) free(lco);
    }
    if (pHealpixData) free(pHealpixData);
    for (auto& afi : afiLightCone) io_free(&afi);
    if (interp_scale) gsl_spline_free(interp_scale);
    if (interp_accel) gsl_interp_accel_free(interp_accel);
}

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceInitLightcone::input> || std::is_standard_layout_v<ServiceInitLightcone::input>);
static_assert(std::is_void_v<ServiceInitLightcone::output> || std::is_standard_layout_v<ServiceInitLightcone::output>);

int ServiceInitLightcone::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    abort_assert(1002, vin.size() == sizeof(input));
    pst->pkd()->InitLightcone(in->bBowtie, in->bLightConeParticles, in->dBoxSize, in->dRedshiftLCP, in->alphaLCP, in->hLCP, in->mOctants);
    return 0;
}

int ServiceLightConeOpen::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    pst->pkd()->LightConeOpen(in->achOutFileLcp, in->nSideHealpix);
    return 0;
}

int ServiceLightConeClose::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    pst->pkd()->LightConeClose(in->achOutFileHpb);
    return 0;
}

int ServiceLightConeVel::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    pst->pkd()->LightConeVel(in->dBoxSize);
    return 0;
}

static void initLightBallOffsets(LightConeInfo& lcinfo, bool amThread0, double mrLCP) {
    Bound bnd(TinyVector<double, 3>(-0.5), TinyVector<double, 3>(0.5));
    double min2;
    int ix, iy, iz, l;

    lcinfo.nLayerMax = ceil(mrLCP);
    lcinfo.nBoxLC.clear();

    lcinfo.lcOffset0.clear();
    lcinfo.lcOffset1.clear();
    lcinfo.lcOffset2.clear();
    /*
    ** Set up the light cone offsets such that they proceed from the inner 8
    ** unit boxes outward layer by layer so that we can skip checks of the
    ** outer layers if we want.
    */
    for (l = 0; l < lcinfo.nLayerMax; ++l) {
        for (ix = -l; ix <= l + 1; ++ix) {
            for (iy = -l; iy <= l + 1; ++iy) {
                for (iz = -l; iz <= l + 1; ++iz) {
                    if ((ix > -l) && (ix < l + 1) && (iy > -l) && (iy < l + 1) && (iz > -l) && (iz < l + 1)) continue;
                    TinyVector<double, 3> r(ix - 0.5, iy - 0.5, iz - 0.5);
                    min2 = bnd.mindist(r);
                    if (min2 < mrLCP * mrLCP) {
                        if (amThread0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f> min2:%f\n", lcinfo.lcOffset0.size(), l, r[0], r[1], r[2], min2);
                        lcinfo.lcOffset0.push_back(r[0]);
                        lcinfo.lcOffset1.push_back(r[1]);
                        lcinfo.lcOffset2.push_back(r[2]);
                    }
                }
            }
        }
        /*
        ** Can save the nBox here for skipping the deeper parts of the lightcone checks!
        */
        lcinfo.nBoxLC.push_back(lcinfo.lcOffset0.size());
    }
}

static void initLightConeOffsets(LightConeInfo& lcinfo, bool amThread0, int bBowtie, blitz::TinyVector<double, 3> h, double alpha, double mrLCP) {
    int l, ix, iy, iz;
    double hm = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
    double min2;

    lcinfo.nLayerMax = ceil(mrLCP);
    for (int j = 0; j < 3; ++j) h[j] /= hm;
    if (amThread0) printf("mrLCP:%f nLayerMax:%d vector width:%d h:<%f %f %f>\n", mrLCP, lcinfo.nLayerMax, dvec::width(), h[0], h[1], h[2]);
    /*
    ** Check that all the components of h are positive and that h is not too close to the
    ** walls!
    */
    assert(h[0] > 0);
    assert(h[1] > 0);
    assert(h[2] > 0);

    lcinfo.nBoxLC.clear();

    // by default all projected cells need to be checked
    std::vector<char> xy(lcinfo.nLayerMax * lcinfo.nLayerMax, true);
    std::vector<char> xz(lcinfo.nLayerMax * lcinfo.nLayerMax, true);
    std::vector<char> yz(lcinfo.nLayerMax * lcinfo.nLayerMax, true);

    double alpha_xy = atan2(tan(alpha), sqrt(h[0] * h[0] + h[1] * h[1])); // angle projected onto the xy plane
    double alpha_xz = atan2(tan(alpha), sqrt(h[0] * h[0] + h[2] * h[2])); // angle projected onto the xz plane
    double alpha_yz = atan2(tan(alpha), sqrt(h[1] * h[1] + h[2] * h[2])); // angle projected onto the xy plane
    double beta_yx = atan2(h[0], h[1]);
    assert(beta_yx > alpha_xy); // make sure it is far enough from the walls
    double beta_zx = atan2(h[0], h[2]);
    assert(beta_zx > alpha_xz); // make sure it is far enough from the walls
    double beta_zy = atan2(h[1], h[2]);
    assert(beta_zy > alpha_yz); // make sure it is far enough from the walls

    double beta_xy = atan2(h[1], h[0]);
    assert(beta_xy > alpha_xy); // make sure it is far enough from the walls
    double beta_xz = atan2(h[2], h[0]);
    assert(beta_xz > alpha_xz); // make sure it is far enough from the walls
    double beta_yz = atan2(h[2], h[1]);
    assert(beta_yz > alpha_yz); // make sure it is far enough from the walls
    if (amThread0)
        printf("alpha:%f beta_xy:%f alpha_xy:%f beta_xz:%f alpha_xz:%f beta_yz:%f alpha_yz:%f\n", alpha, beta_xy, alpha_xy, beta_xz, alpha_xz, beta_yz,
               alpha_yz);
    for (int j = 0; j < lcinfo.nLayerMax; ++j) {
        for (int i = 0; i < lcinfo.nLayerMax; ++i) {
            double a1 = atan2(j, i + 1);
            double a2 = atan2(j + 1, i);
            if (beta_xy > a2 + alpha_xy || beta_xy < a1 - alpha_xy) xy[j * lcinfo.nLayerMax + i] = false; // doesn't intersect
            if (beta_xz > a2 + alpha_xz || beta_xz < a1 - alpha_xz) xz[j * lcinfo.nLayerMax + i] = false; // doesn't intersect
            if (beta_yz > a2 + alpha_yz || beta_yz < a1 - alpha_yz) yz[j * lcinfo.nLayerMax + i] = false; // doesn't intersect
        }
    }
    lcinfo.lcOffset0.clear();
    lcinfo.lcOffset1.clear();
    lcinfo.lcOffset2.clear();
    /*
    ** Set up the light cone offsets such that they proceed from the inner 8
    ** unit boxes outward layer by layer so that we can skip checks of the
    ** outer layers if we want.
    */
    for (l = 0; l < lcinfo.nLayerMax; ++l) {
        for (ix = 0; ix <= l; ++ix) {
            for (iy = 0; iy <= l; ++iy) {
                for (iz = 0; iz <= l; ++iz) {
                    if (ix < l && iy < l && iz < l) continue;
                    if (xy[iy * lcinfo.nLayerMax + ix] && xz[iz * lcinfo.nLayerMax + ix] && yz[iz * lcinfo.nLayerMax + iy]) {
                        double r[3] = {ix + 0.5, iy + 0.5, iz + 0.5};
                        Bound bnd(TinyVector<double, 3>(-0.5), TinyVector<double, 3>(0.5));
                        min2 = bnd.mindist(r);
                        if (min2 > mrLCP * mrLCP) continue;
                        if (amThread0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f>\n", lcinfo.lcOffset0.size(), l, r[0], r[1], r[2]);
                        lcinfo.lcOffset0.push_back(r[0]);
                        lcinfo.lcOffset1.push_back(r[1]);
                        lcinfo.lcOffset2.push_back(r[2]);
                        if (bBowtie) {
                            /*
                            ** Now include also the -,-,- octant for checks, as we will be producing a bowtie.
                            */
                            for (int j = 0; j < 3; ++j) r[j] = -r[j];
                            if (amThread0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f>\n", lcinfo.lcOffset0.size(), l, r[0], r[1], r[2]);
                            lcinfo.lcOffset0.push_back(r[0]);
                            lcinfo.lcOffset1.push_back(r[1]);
                            lcinfo.lcOffset2.push_back(r[2]);
                        }
                    }
                }
            }
        }
        /*
        ** Can save the nBox here for skipping the deeper parts of the lightcone checks!
        */
        lcinfo.nBoxLC.push_back(lcinfo.lcOffset0.size());
    }
    /*
    ** Add more boxes to get to a multiple of the vector width.
    ** Just make sure the added boxes far enough away to not get
    ** included.
    */
    l = lcinfo.nLayerMax - 1;
    double r[3] = {2 * l + 0.5, 2 * l + 0.5, 2 * l + 0.5};
    while (lcinfo.lcOffset0.size() % dvec::width()) {
        if (amThread0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f>\n", lcinfo.lcOffset0.size(), l, r[0], r[1], r[2]);
        lcinfo.lcOffset0.push_back(r[0]);
        lcinfo.lcOffset1.push_back(r[1]);
        lcinfo.lcOffset2.push_back(r[2]);
    }
    lcinfo.nBoxLC.push_back(lcinfo.lcOffset0.size());
}

static void SetupInterpScale(LightConeInfo& lcinfo, CSM csm, double dBoxSize, double mrMax) {
    const int nTable = 1000;
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    double dr, rt[nTable], at_inv[nTable];
    /*
    ** Setup lookup table.
    */
    dr = mrMax / (nTable - 1);
    for (int i = 0; i < nTable; ++i) {
        rt[i] = i * dr;
        at_inv[i] = 1.0 / csmComoveLookbackTime2Exp(csm, rt[i] / dLightSpeed);
    }
    lcinfo.interp_accel = gsl_interp_accel_alloc();
    lcinfo.interp_scale = gsl_spline_alloc(gsl_interp_cspline, nTable);
    gsl_spline_init(lcinfo.interp_scale, rt, at_inv, nTable);
}

/*
** Initialize Lightcone stuff.
*/
void pkdContext::InitLightcone(int bBowtie, int bLightConeParticles, double dBoxSize, double dRedshiftLCP, double alphaLCP, blitz::TinyVector<double, 3> hLCP,
                               int mOctants) {
#ifdef __linux__
    uint64_t nPageSize = sysconf(_SC_PAGESIZE);
#else
    uint64_t nPageSize = 512;
#endif
    lcinfo.mLightConeOctants = mOctants;
    /*
    ** Initialize Lookup table for converting lightcone velocities to physical (sim units)
    */
    double dTimeLCP = csmExp2Time(csm, 1.0 / (1.0 + dRedshiftLCP));
    double mrLCP = dLightSpeedSim(dBoxSize) * csmComoveKickFac(csm, dTimeLCP, (csmExp2Time(csm, 1.0) - dTimeLCP));
    SetupInterpScale(lcinfo, csm, dBoxSize, mrLCP);
    /*
    ** Initialize light cone offsets.
    */
    if (Self() == 0) printf("alphaLCP = %f\n", alphaLCP);
    if (alphaLCP < 0) {
        initLightBallOffsets(lcinfo, Self() == 0, mrLCP);
    }
    else {
        initLightConeOffsets(lcinfo, Self() == 0, bBowtie, hLCP, alphaLCP, mrLCP);
    }
    /*
    ** allocate enough space for light cone particle output
    */
    uint64_t nLightConeBytes = (1024 * 1024 * 2);
    lcinfo.nLightConeMax = nLightConeBytes / sizeof(LIGHTCONEP);
    lcinfo.nLightCone = 0; // array of 8 octants
    auto mask = lcinfo.mLightConeOctants == 0 ? 1 : lcinfo.mLightConeOctants;
    for (auto octant = 0; octant < 8; ++octant) {
        if (bLightConeParticles && (mask & (1 << octant))) {
            void* v;
            if (posix_memalign(&v, nPageSize, nLightConeBytes)) lcinfo.pLightCone[octant] = NULL;
            else lcinfo.pLightCone[octant] = static_cast<LIGHTCONEP*>(v);
            assert(lcinfo.pLightCone[octant] != NULL);
            io_init(&lcinfo.afiLightCone[octant], 8, 2 * 1024 * 1024, IO_AIO | IO_LIBAIO);
        }
        else {
            lcinfo.afiLightCone[octant].nBuffers = 0;
            lcinfo.pLightCone[octant] = NULL;
        }
        lcinfo.afiLightCone[octant].fd = -1;
    }
    lcinfo.pHealpixData = NULL;
}

#define FILE_PROTECTION (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP)

void flushLightCone(LightConeInfo& lcinfo, int octant) {
    size_t count = lcinfo.nLightCone[octant] * sizeof(LIGHTCONEP);
    io_write(&lcinfo.afiLightCone[octant], lcinfo.pLightCone[octant], count);
    lcinfo.nLightCone[octant] = 0;
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
    explicit HealpixHelper(LightConeInfo& lcinfo) : arrayCACHEhelper(lcinfo.pHealpixData, lcinfo.nHealpixPerDomain, true) {}
};

void pkdContext::LightConeClose(const char* healpixname_template) {
    for (auto octant = 0; octant < 8; ++octant) {
        if (lcinfo.afiLightCone[octant].fd > 0) {
            flushLightCone(lcinfo, octant);
            io_close(&lcinfo.afiLightCone[octant]);
        }
    }
    if (lcinfo.nSideHealpix) {
        assert(healpixname_template && healpixname_template[0]);
        auto healpixname = fmt::format(healpixname_template, fmt::arg("i", mdl->Self()));

        mdl->FinishCache(CID_HEALPIX);
        int fd = open(healpixname.c_str(), O_CREAT | O_WRONLY | O_TRUNC, FILE_PROTECTION);
        if (fd < 0) {
            perror(healpixname.c_str());
            abort();
        }
        for (auto i = 0; i < lcinfo.nHealpixPerDomain; ++i) {
            uint64_t sum = lcinfo.pHealpixData[i].nGrouped;
            sum += lcinfo.pHealpixData[i].nUngrouped;
            if (sum) lcinfo.pHealpixData[i].fPotential /= sum;
        }
        auto nWrite = lcinfo.nHealpixPerDomain * sizeof(*lcinfo.pHealpixData);
        if (write(fd, lcinfo.pHealpixData, nWrite) != nWrite) {
            perror("Wrong size writing Healpix");
        }
        close(fd);
    }
}

void pkdContext::LightConeOpen(const char* name_template, int nSideHealpix) {
    auto mask = (lcinfo.mLightConeOctants == 0) ? 1 : lcinfo.mLightConeOctants;
    for (auto octant = 0; octant < 8; ++octant) {
        if ((mask & (1 << octant)) == 0 || name_template[0] == '\0') {
            lcinfo.afiLightCone[octant].fd = -1;
            continue;
        }
        auto fname = fmt::format(name_template, "octant"_a = octant, "i"_a = mdl->Self());
        if (io_create(&lcinfo.afiLightCone[octant], fname.c_str()) < 0) {
            perror(fname.c_str());
            abort();
        }
    }

    this->lcinfo.nSideHealpix = nSideHealpix;
    if (nSideHealpix) {
        lcinfo.nHealpixPerDomain = (nside2npix64(nSideHealpix) + mdl->Threads() - 1) / mdl->Threads();
        lcinfo.nHealpixPerDomain = (lcinfo.nHealpixPerDomain + 15) & ~15;
        if (lcinfo.pHealpixData == NULL) {
            lcinfo.pHealpixData = static_cast<healpixData*>(malloc(lcinfo.nHealpixPerDomain * sizeof(*lcinfo.pHealpixData)));
            assert(lcinfo.pHealpixData != NULL);
        }
        for (auto i = 0; i < lcinfo.nHealpixPerDomain; ++i) {
            lcinfo.pHealpixData[i].nGrouped = 0;
            lcinfo.pHealpixData[i].nUngrouped = 0;
            lcinfo.pHealpixData[i].fPotential = 0;
        }
        mdl->CacheInitialize(CID_HEALPIX, std::make_shared<HealpixHelper>(lcinfo));
    }
}

static void addToLightCone(LightConeInfo& lcinfo, mdl::mdlClass* mdl, double dvFac, double* r, float fPot, particleStore::ParticleReference& P,
                           int bParticleOutput) {
    const auto& v = P.velocity();
    if (bParticleOutput) {
        auto octant = lcinfo.mLightConeOctants == 0 ? 0 : (((r[0] < 0) ? 4 : 0) + ((r[1] < 0) ? 2 : 0) + ((r[2] < 0) ? 1 : 0));
        if (lcinfo.afiLightCone[octant].fd > 0) {
            auto pLC = lcinfo.pLightCone[octant];
            auto n = lcinfo.nLightCone[octant];
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
            if (++lcinfo.nLightCone[octant] == lcinfo.nLightConeMax) flushLightCone(lcinfo, octant);
        }
    }
    if (lcinfo.nSideHealpix) {
        int64_t iPixel = vec2pix_ring64(lcinfo.nSideHealpix, r);
        assert(iPixel >= 0);
        int id = iPixel / lcinfo.nHealpixPerDomain;
        int idx = iPixel - id * lcinfo.nHealpixPerDomain;
        assert(id < mdl->Threads());
        assert(idx < lcinfo.nHealpixPerDomain);
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

void pkdContext::ProcessLightCone(particleStore::ParticleReference& P, float fPot, double dLookbackFac, double dLookbackFacLCP, double dDriftDelta,
                                  double dKickDelta, double dBoxSize, int bLightConeParticles, blitz::TinyVector<double, 3> hlcp, double tanalpha_2) {
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    const double mrLCP = dLightSpeed * dLookbackFacLCP;
    const double depth = dLightSpeed * dLookbackFac;
    int bCone = 1;
    const int nLayerMax = ceil(mrLCP);
    double dxStart;

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
    int nBox = lcinfo.nBoxLC[l]; // check only this many boxes
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
            off0.load(&lcinfo.lcOffset0[iOct * dvec::width()]);
            off1.load(&lcinfo.lcOffset1[iOct * dvec::width()]);
            off2.load(&lcinfo.lcOffset2[iOct * dvec::width()]);
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
                                    double dvFac = gsl_spline_eval(lcinfo.interp_scale, mr, lcinfo.interp_accel);
                                    addToLightCone(lcinfo, mdl, dvFac, r, fPot, P, bLightConeParticles);
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
