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
using blitz::TinyVector;
using blitz::dot;

#ifdef M43D
/* M43D Creates a 3D kernel by convolution of 3D tophats the way M4(1D) is made in 1D */
#define BALL2(fBall) ((fBall)*(fBall))
#define KERNEL(ak,ar2) { \
        ak = sqrt(ar2); \
        if (ar2 < 1.0) ak = 6.*0.25/350./3. *(1360+ar2*(-2880 \
             +ar2*(3528+ak*(-1890+ak*(-240+ak*(270-6*ar2)))))); \
        else if (ar2 < 4.0) ak = 6.*0.25/350./3. *(7040-1152/ak+ak*(-10080+ak*(2880+ak*(4200 \
             +ak*(-3528+ak*(630+ak*(240+ak*(-90+2*ar2)))))))); \
        else ak = 0.0;\
        }
#define DKERNEL(adk,ar2) { \
        adk = sqrt(ar2); \
        if (ar2 < 1.0) adk = 6.*0.25/350./3. * (-2880*2 \
             +ar2*(3528*4+ adk*(-1890*5 + adk*(-240*6+ adk*(270*7-6*9*ar2))))); \
        else if (ar2 < 4.0) adk = 6.*0.25/350./3. *((1152/ar2-10080)/adk+(2880*2+adk*(4200*3 \
             +adk*(-3528*4+adk*(630*5+adk*(240*6 +adk*(-90*7+2*9*ar2))))))); \
        else adk = 0.0;\
        }

#else
#ifdef HSHRINK
/* HSHRINK M4 Kernel uses an effective h of (pi/6)^(1/3) times h for nSmooth neighbours */
#define dSHRINKFACTOR        0.805995977
#define BALL2(fBall) ((fBall)*(fBall)*(dSHRINKFACTOR*dSHRINKFACTOR))
#define KERNEL(ak,ar2) { \
        ak = 2.0 - sqrt(ar2); \
        if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
        else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
        else ak = 0.0; \
        }
#define DKERNEL(adk,ar2) { \
        adk = sqrt(ar2); \
        if (ar2 < 1.0) { \
            adk = -3 + 2.25*adk; \
            } \
        else if (ar2 < 4.0) { \
            adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
            } \
        else adk = 0.0; \
        }

#else
/* Standard M_4 Kernel */
#define BALL2(fBall) ((fBall)*(fBall))
#define KERNEL(ak,ar2) { \
        ak = 2.0 - sqrt(ar2); \
        if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
        else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
        else ak = 0.0;\
        }
#define DKERNEL(adk,ar2) { \
        adk = sqrt(ar2); \
        if (ar2 < 1.0) { \
            adk = -3 + 2.25*adk; \
            } \
        else if (ar2 < 4.0) { \
            adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
            } \
        else adk = 0.0;\
        }
#endif
#endif

void NullSmooth(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
}

void initBall(void *vpkd, void *pIn) {
    auto pkd = static_cast<PKD>(vpkd);
    auto p = pkd->particles[static_cast<PARTICLE *>(pIn)];
    p.set_ball(0.0);
}

void BallSmooth(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    p.set_ball(fBall);
}

void initDensity(void *vpkd, void *p) {
    auto pkd = static_cast<PKD>(vpkd);
    auto P = pkd->particles[reinterpret_cast<PARTICLE *>(p)];
    P.set_density(0.0);
}

void combDensity(void *vpkd, void *p1,const void *p2) {
    auto pkd = static_cast<PKD>(vpkd);
    auto P1 = pkd->particles[reinterpret_cast<PARTICLE *>(p1)];
    auto P2 = pkd->particles[reinterpret_cast<const PARTICLE *>(p2)];
    P1.set_density(P1.density() + P2.density());
}

void DensityF1(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[p];
    double ih2,r2,rs,fDensity,fMass;
    int i;

    ih2 = 1.0/BALL2(fBall);
    fDensity = 0.0;
    for (i=0; i<nSmooth; ++i) {
        auto Q = pkd->particles[nnList[i].pPart];
        fMass = Q.mass();
        r2 = nnList[i].fDist2*ih2;
        rs = 1 - r2;
        if (rs < 0) rs = 0.0;
        fDensity += rs*fMass;
    }
    fDensity *= 1.875f*M_1_PI*sqrtf(ih2)*ih2; /* F1 Kernel (15/8) */
    if (smf->pfDensity) *smf->pfDensity = fDensity;
    else P.set_density(fDensity);
}

void DensityM3(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[p];
    double ih2,r2,rs,fDensity,fMass;
    int i;
    ih2 = 1.0f/BALL2(fBall);
    fDensity = 0.0;
    for (i=0; i<nSmooth; ++i) {
        auto Q = pkd->particles[nnList[i].pPart];
        fMass = Q.mass();
        r2 = nnList[i].fDist2*ih2;
        if (r2 < 1.0) {
            double r = sqrt(r2);
            rs = 1.0f - r;
            rs *= rs*rs; /* rs^3 */
            if (r < 0.5f) {
                double rs2 = 0.5f - r;
                rs2 *= rs2*rs2; /* rs2^3 */
                rs -= 4.0f*rs2;
            }
        }
        else rs = 0.0;
        fDensity += rs*fMass;
    }
    P.set_density(16.0f*M_1_PI*sqrtf(ih2)*ih2*fDensity);
}

void LinkGradientM3(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[p];
    double ih2,r2,rs,fMass,fNorm, idrho, r2min;
    int i, j;
    ih2 = 1.0/BALL2(fBall);
    fNorm = 16.0f*M_1_PI*ih2*ih2*sqrtf(ih2);
    TinyVector<double,3> frho(0.0);
    for (i=0; i<nSmooth; ++i) {
        auto Q = pkd->particles[nnList[i].pPart];
        fMass = Q.mass();
        r2 = nnList[i].fDist2*ih2;
        if (r2 < 1.0) {
            double r = sqrt(r2);
            if (r < 0.5f) {
                rs = -6.0f + 9.0f*r;
            }
            else {
                rs = 1.0f - r;
                rs *= rs; /* rs^2 */
                rs *= -3.0f;
                rs /= r;
            }
        }
        else rs = 0.0;
        rs *= fNorm*fMass;
        rs *= (Q.density() - P.density())/Q.density();
        frho -= nnList[i].dr*rs;
    }
    idrho = 1.0/sqrt(dot(frho,frho));
    for (j=0; j<3; ++j) frho[j] *= 0.5*idrho*fBall;
    r2min = HUGE_VALF;
    if (nSmooth==0) P.set_group(-1);
    for (i=0; i<nSmooth; ++i) {
        TinyVector<double,3> dr = nnList[i].dr - frho;
        r2 = dot(dr,dr);
        if (r2 < r2min) {
            r2min = r2;
            smf->hopParticleLink.iPid = nnList[i].iPid;
            smf->hopParticleLink.iIndex = nnList[i].iIndex;
        }
    }
}

void LinkHopChains(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    MDL mdl = pkd->mdl;
    auto P = pkd->particles[p];
    int i, gid1, gid2;
    GHtmpGroupTable *g1, *g2, g;
    gid1 = P.group();
    g1 = &pkd->tmpHopGroups[gid1];
    for (i=0; i<nSmooth; ++i) {
        auto Q = pkd->particles[nnList[i].pPart];
        gid2 = Q.group();
        if (nnList[i].iPid==pkd->Self() && gid1==gid2) continue;
        g.iPid = nnList[i].iPid;
        g.iIndex = gid2;
        g2 = static_cast<GHtmpGroupTable *>(mdlAcquire(mdl,CID_GROUP,g.iIndex,g.iPid));

        /* Remote is authoratative. Update myself, but also what I currently link to. */
        if (g1->iPid > g2->iPid || (g1->iPid == g2->iPid && g1->iIndex > g2->iIndex)) {
            smf->bDone = 0;
            g = *g1;
            g1->iPid = g2->iPid;
            g1->iIndex = g2->iIndex;
            mdlRelease(mdl,CID_GROUP,g2);
            g2 = static_cast<GHtmpGroupTable *>(mdlAcquire(mdl,CID_GROUP,g.iIndex,g.iPid));
        }

        /* Update remote (or what we were pointing to) and what it points to if necessary. */
        while (g1->iPid < g2->iPid || (g1->iPid == g2->iPid && g1->iIndex < g2->iIndex) ) {
            smf->bDone = 0;
            g = *g2;
            g2->iPid = g1->iPid;
            g2->iIndex = g1->iIndex;
            mdlRelease(mdl,CID_GROUP,g2);
            g2 = static_cast<GHtmpGroupTable *>(mdlAcquire(mdl,CID_GROUP,g.iIndex,g.iPid));
        }
        mdlRelease(mdl,CID_GROUP,g2);
    }
}

void Density(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[p];
    double ih2,r2,rs,fDensity;
    int i;

    ih2 = 4.0/BALL2(fBall);
    fDensity = 0.0;
    for (i=0; i<nSmooth; ++i) {
        auto Q = pkd->particles[nnList[i].pPart];
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        fDensity += rs*Q.mass();
    }
    P.set_density(M_1_PI*sqrt(ih2)*ih2*fDensity);
}

void DensitySym(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[p];
    double fNorm,ih2,r2,rs,fMassP;
    int i;
    fMassP = P.mass();
    ih2 = 4.0/(BALL2(fBall));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        rs *= fNorm;
        auto Q = pkd->particles[nnList[i].pPart];
        P.set_density(P.density() + rs*Q.mass());
        Q.set_density(Q.density() + rs*fMassP);
    }
}

void PrintNN(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    int i;

    printf("%" PRIu64 ":",(uint64_t)p->iOrder);
    for (i=0; i<nSmooth; ++i) {
        if (pkdIsActive(pkd,nnList[i].pPart))
            printf("%" PRIu64 " ",(uint64_t)nnList[i].pPart->iOrder);
        else
            printf("\033[7m%" PRIu64 "\033[0m ",(uint64_t)nnList[i].pPart->iOrder);
    }
    printf("\n");
}
