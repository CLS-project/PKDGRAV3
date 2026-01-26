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

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <cinttypes>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "smoothfcn.h"
using blitz::dot;
using blitz::TinyVector;

#ifdef M43D
    /* M43D Creates a 3D kernel by convolution of 3D tophats the way M4(1D) is made in 1D */
    #define BALL2(fBall) ((fBall) * (fBall))
    #define KERNEL(ak, ar2)                                                                                                                                    \
        {                                                                                                                                                      \
            ak = sqrt(ar2);                                                                                                                                    \
            if (ar2 < 1.0) ak = 6. * 0.25 / 350. / 3. * (1360 + ar2 * (-2880 + ar2 * (3528 + ak * (-1890 + ak * (-240 + ak * (270 - 6 * ar2))))));             \
            else if (ar2 < 4.0)                                                                                                                                \
                ak = 6. * 0.25 / 350. / 3. *                                                                                                                   \
                     (7040 - 1152 / ak + ak * (-10080 + ak * (2880 + ak * (4200 + ak * (-3528 + ak * (630 + ak * (240 + ak * (-90 + 2 * ar2))))))));           \
            else ak = 0.0;                                                                                                                                     \
        }
    #define DKERNEL(adk, ar2)                                                                                                                                  \
        {                                                                                                                                                      \
            adk = sqrt(ar2);                                                                                                                                   \
            if (ar2 < 1.0)                                                                                                                                     \
                adk = 6. * 0.25 / 350. / 3. * (-2880 * 2 + ar2 * (3528 * 4 + adk * (-1890 * 5 + adk * (-240 * 6 + adk * (270 * 7 - 6 * 9 * ar2)))));           \
            else if (ar2 < 4.0)                                                                                                                                \
                adk = 6. * 0.25 / 350. / 3. *                                                                                                                  \
                      ((1152 / ar2 - 10080) / adk +                                                                                                            \
                       (2880 * 2 + adk * (4200 * 3 + adk * (-3528 * 4 + adk * (630 * 5 + adk * (240 * 6 + adk * (-90 * 7 + 2 * 9 * ar2)))))));                 \
            else adk = 0.0;                                                                                                                                    \
        }

#else
    #ifdef HSHRINK
        /* HSHRINK M4 Kernel uses an effective h of (pi/6)^(1/3) times h for nSmooth neighbours */
        #define dSHRINKFACTOR 0.805995977
        #define BALL2(fBall) ((fBall) * (fBall) * (dSHRINKFACTOR * dSHRINKFACTOR))
        #define KERNEL(ak, ar2)                                                                                                                                \
            {                                                                                                                                                  \
                ak = 2.0 - sqrt(ar2);                                                                                                                          \
                if (ar2 < 1.0) ak = (1.0 - 0.75 * ak * ar2);                                                                                                   \
                else if (ar2 < 4.0) ak = 0.25 * ak * ak * ak;                                                                                                  \
                else ak = 0.0;                                                                                                                                 \
            }
        #define DKERNEL(adk, ar2)                                                                                                                              \
            {                                                                                                                                                  \
                adk = sqrt(ar2);                                                                                                                               \
                if (ar2 < 1.0) {                                                                                                                               \
                    adk = -3 + 2.25 * adk;                                                                                                                     \
                }                                                                                                                                              \
                else if (ar2 < 4.0) {                                                                                                                          \
                    adk = -0.75 * (2.0 - adk) * (2.0 - adk) / adk;                                                                                             \
                }                                                                                                                                              \
                else adk = 0.0;                                                                                                                                \
            }

    #else
        /* Standard M_4 Kernel */
        #define BALL2(fBall) ((fBall) * (fBall))
        #define KERNEL(ak, ar2)                                                                                                                                \
            {                                                                                                                                                  \
                ak = 2.0 - sqrt(ar2);                                                                                                                          \
                if (ar2 < 1.0) ak = (1.0 - 0.75 * ak * ar2);                                                                                                   \
                else if (ar2 < 4.0) ak = 0.25 * ak * ak * ak;                                                                                                  \
                else ak = 0.0;                                                                                                                                 \
            }
        #define DKERNEL(adk, ar2)                                                                                                                              \
            {                                                                                                                                                  \
                adk = sqrt(ar2);                                                                                                                               \
                if (ar2 < 1.0) {                                                                                                                               \
                    adk = -3 + 2.25 * adk;                                                                                                                     \
                }                                                                                                                                              \
                else if (ar2 < 4.0) {                                                                                                                          \
                    adk = -0.75 * (2.0 - adk) * (2.0 - adk) / adk;                                                                                             \
                }                                                                                                                                              \
                else adk = 0.0;                                                                                                                                \
            }
    #endif
#endif

void init_shared_gas_decorator(PKD pkd, void* dst) {
    auto p = pkd->particles[static_cast<PARTICLE*>(dst)];
    if (p.is_gas()) {
        auto src_decorator = p.decorator_field();
        p.set_adj_decorator();
        memcpy(p.decorator_field(), src_decorator, pkd->decorators.DecoratorSize());
    }
}

void init_shared_star_decorator(PKD pkd, void* dst) {
    auto p = pkd->particles[static_cast<PARTICLE*>(dst)];
    if (p.is_star()) {
        auto src_decorator = p.decorator_field();
        p.set_adj_decorator();
        memcpy(p.decorator_field(), src_decorator, pkd->decorators.DecoratorSize());
    }
}

void init_shared_bh_decorator(PKD pkd, void* dst) {
    auto p = pkd->particles[static_cast<PARTICLE*>(dst)];
    if (p.is_bh()) {
        auto src_decorator = p.decorator_field();
        p.set_adj_decorator();
        memcpy(p.decorator_field(), src_decorator, pkd->decorators.DecoratorSize());
    }
}

SmoothHelperBallSmooth::SmoothHelperBallSmooth(pkdContext* pkd) : SmoothHelper(pkd) {
    assert(pkd->particles.present(PKD_FIELD::oBall));
}
void SmoothHelperBallSmooth::init(void* pIn) {
    auto p = pkd->particles[static_cast<PARTICLE*>(pIn)];
    p.set_ball(0.0);
}
void SmoothHelperBallSmooth::initialize_local() {
    auto pRoot = pkd->tree[ROOT];
    for (auto& p : *pRoot) {
        if (p.is_active()) init(&p);
    }
}
void SmoothHelperBallSmooth::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    p.set_ball(p.is_gas() ? fBall : 0);
}

void SmoothHelperDensitySymmetric::init(void* pIn) {
    auto p = pkd->particles[static_cast<PARTICLE*>(pIn)];
    p.set_density(0.0);
}
void SmoothHelperDensitySymmetric::combine(void* dst, const void* src, const void* key) {
    auto P1 = pkd->particles[reinterpret_cast<PARTICLE*>(dst)];
    auto P2 = pkd->particles[reinterpret_cast<const PARTICLE*>(src)];
    P1.set_density(P1.density() + P2.density());
}
void SmoothHelperDensitySymmetric::initialize_local() {
    auto pRoot = pkd->tree[ROOT];
    for (auto& p : *pRoot) {
        if (p.is_active()) init(&p);
    }
}

void SmoothHelperDensityF1::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    auto ih2 = 1.0 / BALL2(fBall);
    auto fDensity = 0.0;
    for (auto& nn : nnList) {
        auto Q = pkd->particles[nn.pPart];
        auto fMass = Q.mass();
        auto r2 = nn.fDist2 * ih2;
        auto rs = 1 - r2;
        if (rs < 0) rs = 0.0;
        fDensity += rs * fMass;
    }
    fDensity *= 1.875f * M_1_PI * sqrtf(ih2) * ih2; /* F1 Kernel (15/8) */
    set_density(p, fDensity);
}

void SmoothHelperDensityM3::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    double ih2 = 1.0 / BALL2(fBall);
    double fDensity = 0.0;
    for (const auto& nn : nnList) {
        auto Q = pkd->particles[nn.pPart];
        auto fMass = Q.mass();
        auto r2 = nn.fDist2 * ih2;
        double rs;
        if (r2 < 1.0) {
            double r = sqrt(r2);
            rs = 1.0f - r;
            rs *= rs * rs; /* rs^3 */
            if (r < 0.5f) {
                double rs2 = 0.5f - r;
                rs2 *= rs2 * rs2; /* rs2^3 */
                rs -= 4.0f * rs2;
            }
        }
        else rs = 0.0;
        fDensity += rs * fMass;
    }
    p.set_density(16.0f * M_1_PI * sqrtf(ih2) * ih2 * fDensity);
}

void SmoothHelperGradientM3::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    double rs, idrho, r2min;
    int i, j;
    auto ih2 = 1.0 / BALL2(fBall);
    auto fNorm = 16.0f * M_1_PI * ih2 * ih2 * sqrtf(ih2);
    TinyVector<double, 3> frho(0.0);
    for (auto& nn : nnList) {
        auto q = pkd->particles[nn.pPart];
        double fMass = q.mass();
        auto r2 = nn.fDist2 * ih2;
        if (r2 < 1.0) {
            double r = sqrt(r2);
            if (r < 0.5f) {
                rs = -6.0f + 9.0f * r;
            }
            else {
                rs = 1.0f - r;
                rs *= rs; /* rs^2 */
                rs *= -3.0f;
                rs /= r;
            }
        }
        else rs = 0.0;
        rs *= fNorm * fMass;
        rs *= (q.density() - p.density()) / q.density();
        frho -= nnList[i].dr * rs;
    }
    idrho = 1.0 / sqrt(dot(frho, frho));
    for (j = 0; j < 3; ++j) frho[j] *= 0.5 * idrho * fBall;
    r2min = HUGE_VALF;
    if (nnList.size() == 0) p.set_group(-1);
    for (auto& nn : nnList) {
        TinyVector<double, 3> dr = nn.dr - frho;
        auto r2 = dot(dr, dr);
        if (r2 < r2min) {
            r2min = r2;
            hopParticleLink.iPid = nn.iPid;
            hopParticleLink.iIndex = nn.iIndex;
        }
    }
}

void SmoothHelperHopChains::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    auto mdl = pkd->mdl;
    int i, gid1, gid2;
    GHtmpGroupTable *g1, *g2, g;
    gid1 = p.group();
    g1 = &pkd->tmpHopGroups[gid1];
    bDone = true;
    for (auto& nn : nnList) {
        auto Q = pkd->particles[nn.pPart];
        gid2 = Q.group();
        if (nn.iPid == pkd->Self() && gid1 == gid2) continue;
        g.iPid = nn.iPid;
        g.iIndex = gid2;
        g2 = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, g.iIndex, g.iPid);

        /* Remote is authoratative. Update myself, but also what I currently link to. */
        if (g1->iPid > g2->iPid || (g1->iPid == g2->iPid && g1->iIndex > g2->iIndex)) {
            bDone = false;
            g = *g1;
            g1->iPid = g2->iPid;
            g1->iIndex = g2->iIndex;
            mdl->Release(CID_GROUP, g2);
            g2 = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, g.iIndex, g.iPid);
        }

        /* Update remote (or what we were pointing to) and what it points to if necessary. */
        while (g1->iPid < g2->iPid || (g1->iPid == g2->iPid && g1->iIndex < g2->iIndex)) {
            bDone = false;
            g = *g2;
            g2->iPid = g1->iPid;
            g2->iIndex = g1->iIndex;
            mdl->Release(CID_GROUP, g2);
            g2 = mdl->Acquire<GHtmpGroupTable>(CID_GROUP, g.iIndex, g.iPid);
        }
        mdl->Release(CID_GROUP, g2);
    }
}

void SmoothHelperDensity::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    auto ih2 = 4.0 / BALL2(fBall);
    auto fDensity = 0.0;
    for (auto& nn : nnList) {
        auto Q = pkd->particles[nn.pPart];
        auto r2 = nn.fDist2 * ih2;
        double rs;
        KERNEL(rs, r2);
        fDensity += rs * Q.mass();
    }
    p.set_density(M_1_PI * sqrt(ih2) * ih2 * fDensity);
}

void SmoothHelperDensitySymmetric::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    auto fMassP = p.mass();
    auto ih2 = 4.0 / (BALL2(fBall));
    auto fNorm = 0.5 * M_1_PI * sqrt(ih2) * ih2;
    for (auto& nn : nnList) {
        auto r2 = nn.fDist2 * ih2;
        double rs;
        KERNEL(rs, r2);
        rs *= fNorm;
        auto q = pkd->particles[nn.pPart];
        p.set_density(p.density() + rs * p.mass());
        q.set_density(q.density() + rs * fMassP);
    }
}

void SmoothHelperPrintNN::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    printf("%" PRIu64 ":", p.order());
    for (auto& nn : nnList) {
        if (pkdIsActive(pkd, nn.pPart)) printf("%" PRIu64 " ", (uint64_t)nn.pPart->iOrder);
        else printf("\033[7m%" PRIu64 "\033[0m ", (uint64_t)nn.pPart->iOrder);
    }
    printf("\n");
}