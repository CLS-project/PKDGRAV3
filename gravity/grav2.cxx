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

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stddef.h>
#include <assert.h>
#include <time.h>
#ifdef HAVE_SYS_TIME_H
    #include <sys/time.h>
#endif
#include "pkd.h"
#include "moments.h"
#include "qeval.h"
#include "ewald.h"
#include "grav.h"
#include "group/hop.h"

#include "cuda/cudautil.h"
#include "cuda/cudapppc.h"
#include "cuda/cudaewald.h"
#include "SPH/SPHOptions.h"
#include "SPH/SPHEOS.h"
#include "SPH/SPHpredict.h"

#include <algorithm>

#if 1
#if defined(USE_SIMD) && defined(__SSE2__)
/* Caution: This uses v/sqrt(v) so v cannot be zero! */
static inline float asqrtf(float v) {
    __m128 r2 = _mm_set_ss(v);
    __m128 r = _mm_rsqrt_ps(r2);
    r = _mm_mul_ss(r,_mm_sub_ss(_mm_set_ss(3.0/2.0),_mm_mul_ss(_mm_mul_ss(r,r),_mm_mul_ss(r2,_mm_set_ss(0.5)))));
    r = _mm_mul_ss(r,r2);
    v = _mm_cvtss_f32(r);
    return v;
}
static inline float rsqrtf(float v) {
    __m128 r2 = _mm_set_ss(v);
    __m128 r = _mm_rsqrt_ps(r2);
    r = _mm_mul_ss(r,_mm_sub_ss(_mm_set_ss(3.0/2.0),_mm_mul_ss(_mm_mul_ss(r,r),_mm_mul_ss(r2,_mm_set_ss(0.5)))));
    v =_mm_cvtss_f32(r);
    return v;
}
#else
static inline float asqrtf(float v) {
    return sqrtf(v);
}
static inline float rsqrtf(float v) {
    return 1.0f / sqrtf(v);
}
#endif
#endif

#define SQRT1(d2,dir)\
    {\
    dir = 1/sqrt(d2);\
    }

static void queueDensity( PKD pkd, workParticle *wp, ilpList &ilp, int bGravStep );
/*
** This is called after work has been done for this particle group.
** If everyone has finished, then the particle is updated.
*/
void pkdParticleWorkDone(workParticle *wp) {
    auto pkd = static_cast<PKD>(wp->ctx);
    int i,gid;
    vel_t v2;
    float fx, fy, fz;
    float maga, dT, dtGrav;
    unsigned char uNewRung;

    if ( --wp->nRefs == 0 ) {
        if (wp->SPHoptions->doDensity) {
            float maxkerneldeviation = 0.0f;
            // calculate maximum kernel mass deviation
            for (int i=0; i<wp->nP; i++) {
                float kerneldeviation = 0.0f;
                if (wp->SPHoptions->useNumDen) {
                    kerneldeviation = 4.0f/3.0f*M_PI*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoOut[i].nden - wp->SPHoptions->fKernelTarget;
                }
                else {
                    kerneldeviation = 4.0f/3.0f*M_PI*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoIn[i].fBall*wp->pInfoOut[i].rho - wp->SPHoptions->fKernelTarget;
                }
                if (wp->pInfoIn[i].isTooLarge) {
                    kerneldeviation = 0.0f;
                }
                kerneldeviation = (kerneldeviation > 0) ? kerneldeviation : -kerneldeviation;
                maxkerneldeviation = (kerneldeviation > maxkerneldeviation) ? kerneldeviation : maxkerneldeviation;
            }
            /*
            ** decide if loop has to continue
            ** if true, calculate new fBall for all particles
            ** else, exit loop
            */
            if (maxkerneldeviation/wp->SPHoptions->fKernelTarget > 1e-4f) {
                // do another loop
                for (int i=0; i<wp->nP; i++) {
                    float prefac = 4.0f/3.0f*M_PI;
                    float fBall = wp->pInfoIn[i].fBall;
                    float fx, dfdx;
                    if (wp->SPHoptions->useNumDen) {
                        fx = prefac * fBall * fBall * fBall * wp->pInfoOut[i].nden - wp->SPHoptions->fKernelTarget;
                        dfdx = prefac * 3.0f * fBall * fBall * wp->pInfoOut[i].nden + prefac * fBall * fBall * fBall * wp->pInfoOut[i].dndendfball;
                    }
                    else {
                        fx = prefac * fBall * fBall * fBall * wp->pInfoOut[i].rho - wp->SPHoptions->fKernelTarget;
                        dfdx = prefac * 3.0f * fBall * fBall * wp->pInfoOut[i].rho + prefac * fBall * fBall * fBall * wp->pInfoOut[i].drhodfball;
                    }
                    float newfBall = wp->pInfoIn[i].fBall - fx / dfdx;
                    if (fabsf(newfBall) >= wp->SPHoptions->ballSizeLimit || wp->pInfoIn[i].isTooLarge) {
                        wp->pInfoIn[i].fBall = wp->SPHoptions->ballSizeLimit;
                        wp->pInfoIn[i].isTooLarge = 1;
                    }
                    else if (newfBall < 0.5f * wp->pInfoIn[i].fBall) {
                        wp->pInfoIn[i].fBall = 0.5f * wp->pInfoIn[i].fBall;
                    }
                    else if (newfBall > 1.5f * wp->pInfoIn[i].fBall) {
                        wp->pInfoIn[i].fBall = 1.5f * wp->pInfoIn[i].fBall;
                    }
                    else {
                        wp->pInfoIn[i].fBall = newfBall;
                    }
                }
                wp->nRefs = 1;
                queueDensity(pkd,wp,*wp->ilp,wp->bGravStep);
                pkdParticleWorkDone(wp);
                return;
            }
            for ( int i=0; i<wp->nP; i++ ) {
                // save the new fBall for each particle
                wp->pInfoOut[i].fBall = wp->pInfoIn[i].fBall;
                // wp->pInfoOut[i].fBall = wp->pInfoOut[i].nSmooth;
            }
        }

        float fiDelta = 1.0/wp->ts->dDelta;
        float fEta = wp->ts->dEta;
        float fiAccFac = 1.0 / wp->ts->dAccFac;
        pkd->dFlop += wp->dFlop;
        pkd->dFlopSingleCPU += wp->dFlopSingleCPU;
        pkd->dFlopDoubleCPU += wp->dFlopDoubleCPU;
        pkd->dFlopSingleGPU += wp->dFlopSingleGPU;
        pkd->dFlopDoubleGPU += wp->dFlopDoubleGPU;
        auto p = pkd->particles[pkd->Local()];
        // char particle_buffer[pkd->particles.ParticleSize()];
        // auto p = reinterpret_cast<PARTICLE *>(particle_buffer);
        for ( i=0; i<wp->nP; i++ ) {
            //p = static_cast<PARTICLE *>(mdlAcquireWrite(pkd->mdl,CID_PARTICLE,wp->iPart[i]));
            p = pkd->particles[wp->pPart[i]];
            // pkd->CopyParticle(p,wp->pPart[i]);

            if (p.have_newsph()) {
                auto &NewSph = p.newsph();
                if (wp->SPHoptions->doDensity) {
                    p.set_density(wp->pInfoOut[i].rho);
                    p.set_ball(wp->pInfoOut[i].fBall);
                    if (wp->pInfoIn[i].fBall < wp->SPHoptions->ballSizeLimit) {
                        NewSph.Omega = 1.0f + wp->pInfoOut[i].fBall/(3.0f * wp->pInfoOut[i].rho)*wp->pInfoOut[i].drhodfball;
                    }
                    else {
                        NewSph.Omega = 1.0f;
                    }
                    if (wp->SPHoptions->doInterfaceCorrection) {
                        float imbalance = sqrtf(wp->pInfoOut[i].imbalanceX*wp->pInfoOut[i].imbalanceX + wp->pInfoOut[i].imbalanceY*wp->pInfoOut[i].imbalanceY + wp->pInfoOut[i].imbalanceZ*wp->pInfoOut[i].imbalanceZ) / (0.5f * wp->pInfoOut[i].fBall * wp->pInfoOut[i].rho);
                        imbalance *= calculateInterfaceCorrectionPrefactor(wp->SPHoptions->fKernelTarget,wp->SPHoptions->kernelType);
                        NewSph.expImb2 = expf(-imbalance*imbalance);
                    }
                    if (wp->SPHoptions->doUConversion && !wp->SPHoptions->doInterfaceCorrection) {
                        NewSph.u = SPHEOSUofRhoT(pkd,p.density(),NewSph.u,p.imaterial(),wp->SPHoptions);
                        NewSph.oldRho = p.density();
                    }
                    if (!wp->SPHoptions->doOnTheFlyPrediction) {
                        SPHpredictInDensity(pkd, p, wp->kick, wp->SPHoptions->nPredictRung, &NewSph.P, &NewSph.cs, &NewSph.T, wp->SPHoptions);
                    }
                }
                if (wp->SPHoptions->doDensityCorrection) {
                    float Tbar = wp->pInfoOut[i].corrT / wp->pInfoOut[i].corr;
                    float Pbar = wp->pInfoOut[i].corrP / wp->pInfoOut[i].corr;
                    float Ttilde = NewSph.expImb2 * NewSph.T + (1.0f - NewSph.expImb2) * Tbar;
                    float Ptilde = NewSph.expImb2 * NewSph.P + (1.0f - NewSph.expImb2) * Pbar;
                    float newRho = SPHEOSRhoofPT(pkd, Ptilde, Ttilde, p.imaterial(), wp->SPHoptions);
                    p.set_density(newRho);
                }
                if (wp->SPHoptions->doSPHForces) {
                    NewSph.divv = wp->pInfoOut[i].divv;
                    if (p.imaterial() == 0 && wp->SPHoptions->useIsentropic && wp->SPHoptions->useBuiltinIdeal) {
                        NewSph.uDot = (wp->SPHoptions->gamma - 1.0f) / pow(p.density(),wp->SPHoptions->gamma - 1.0f) * wp->pInfoOut[i].uDot;
                    }
                    else {
                        NewSph.uDot = wp->pInfoOut[i].uDot;
                    }
                }
            }

            if (wp->SPHoptions->doGravity) {
                auto r = p.position();
                auto m = p.mass();
#ifdef HERNQUIST_POTENTIAL
                // Hard coded just for the isolated galaxy test
                const double const_reduced_hubble_cgs = 3.2407789e-18;
                const double dMsolUnit = 1e10;
                const double dKpcUnit = 1.0;
                const double dKmPerSecUnit = sqrt(GCGS*dMsolUnit*MSOLG
                                                  /(dKpcUnit*KPCCM))/1e5;
                const double H0 = 70.4/ dKmPerSecUnit * ( dKpcUnit / 1e3);

                const double concentration = 9.0;
                const double M200 = 135.28423603962767;
                const double V200 = cbrt(10.*M200*H0);
                const double R200 = cbrt(M200/(100.*H0*H0));
                const double RS = R200 / concentration;

                const double al = RS * sqrt(2. * (log(1. + concentration) -
                                                  concentration / (1. + concentration)));

                const double mass = M200*(1.-0.041);
                const float sqrtgm_inv = 1.f / sqrtf(mass);
                const double epsilon =  0.2/dKpcUnit;
                const double epsilon2 = epsilon*epsilon;




                const float dx = pkdPos(pkd,p,0);
                const float dy = pkdPos(pkd,p,1);
                const float dz = pkdPos(pkd,p,2);

                /* Calculate the acceleration */
                const float rr = sqrtf(dx * dx + dy * dy + dz * dz + epsilon2);
                const float r_plus_a_inv = 1.f / (rr + al);
                const float r_plus_a_inv2 = r_plus_a_inv * r_plus_a_inv;
                const float term = -mass * r_plus_a_inv2 / rr;

                /* Calculate the circular orbital period */
                const float period = 2.f * M_PI * sqrtf(rr) * al *
                                     (1 + rr / al) * sqrtgm_inv;



                wp->pInfoOut[i].a[0] += term * dx;
                wp->pInfoOut[i].a[1] += term * dy;
                wp->pInfoOut[i].a[2] += term * dz;
#endif
                if (pkd->particles.present(PKD_FIELD::oAcceleration)) {
                    p.acceleration() = wp->pInfoOut[i].a;
                }
                if (pkd->particles.present(PKD_FIELD::oPotential)) {
                    p.potential() = wp->pInfoOut[i].fPot;
                }
                if (pkd->ga != NULL) {
                    gid = p.group();
                    if (gid && wp->pInfoOut[i].fPot < pkd->ga[gid].minPot) {
                        pkd->ga[gid].minPot = wp->pInfoOut[i].fPot;
                        pkd->ga[gid].iMinPart = wp->iPart[i];
                    }
                }
                pkd->dEnergyU += 0.5 * m * wp->pInfoOut[i].fPot;
                pkd->dEnergyW += m * blitz::dot(r,wp->pInfoOut[i].a);
                pkd->dEnergyF += m * wp->pInfoOut[i].a;

                if (wp->ts->bGravStep) {
                    float dirsum = wp->pInfoOut[i].dirsum;
                    float normsum = wp->pInfoOut[i].normsum;

                    /*
                    ** If this is the first time through, the accelerations will have
                    ** all been zero resulting in zero for normsum (and nan for dtGrav).
                    ** We repeat this process again, so dtGrav will be correct.
                    */
                    if (normsum > 0.0) {
                        /*
                        ** Use new acceleration here!
                        */
                        maga = blitz::dot(wp->pInfoOut[i].a,wp->pInfoOut[i].a);
                        if (maga>0.0f) maga = asqrtf(maga);
                        dtGrav = maga*dirsum/normsum;
                    }
                    else dtGrav = 0.0;
                    dtGrav += wp->ts->dPreFacRhoLoc*wp->pInfoIn[i].fDensity;
                    dtGrav = (wp->pInfoOut[i].rhopmax > dtGrav?wp->pInfoOut[i].rhopmax:dtGrav);
                    if (dtGrav > 0.0) {
                        dT = fEta * rsqrtf(dtGrav*wp->ts->dRhoFac);
                        if (pkd->particles.present(PKD_FIELD::oNewSph)) {
                            dT = std::min(dT,wp->pInfoOut[i].dtEst);
                        }
                        uNewRung = pkdDtToRungInverse(dT,fiDelta,wp->ts->uMaxRung-1);
                        if (pkd->particles.present(PKD_FIELD::oNewSph)) {
                            uNewRung = std::max(std::max((int)uNewRung,(int)round(wp->pInfoOut[i].maxRung) - wp->SPHoptions->nRungCorrection),0);
                        }
                    }
                    else {
                        uNewRung = 0; /* Assumes current uNewRung is outdated -- not ideal */
                        if (pkd->particles.present(PKD_FIELD::oNewSph)) {
                            uNewRung = pkdDtToRungInverse(wp->pInfoOut[i].dtEst,fiDelta,wp->ts->uMaxRung-1);
                            uNewRung = std::max(std::max((int)uNewRung,(int)round(wp->pInfoOut[i].maxRung) - wp->SPHoptions->nRungCorrection),0);
                        }
                    }
                } /* end of wp->bGravStep */
                else {
                    /*
                    ** We are doing eps/a timestepping.
                    */
                    fx = wp->pInfoOut[i].a[0];
                    fy = wp->pInfoOut[i].a[1];
                    fz = wp->pInfoOut[i].a[2];
                    maga = fx*fx + fy*fy + fz*fz;
                    if (maga > 0) {
                        float imaga = rsqrtf(maga) * fiAccFac;
                        dT = fEta*asqrtf(p.soft()*imaga);
                        if (p.have_newsph()) {
                            dT = std::min(dT,wp->pInfoOut[i].dtEst);
                        }
                        uNewRung = pkdDtToRungInverse(dT,fiDelta,wp->ts->uMaxRung-1);
                        if (p.have_newsph()) {
                            uNewRung = std::max(std::max((int)uNewRung,(int)round(wp->pInfoOut[i].maxRung) - wp->SPHoptions->nRungCorrection),0);
                        }
                    }
                    else {
                        uNewRung = 0;
                        if (p.have_newsph()) {
                            uNewRung = pkdDtToRungInverse(wp->pInfoOut[i].dtEst,fiDelta,wp->ts->uMaxRung-1);
                            uNewRung = std::max(std::max((int)uNewRung,(int)round(wp->pInfoOut[i].maxRung) - wp->SPHoptions->nRungCorrection),0);
                        }
                    }
                }
                /*
                ** Here we must make sure we do not try to take a larger opening
                ** timestep than the current active particles involved in the
                ** gravity calculation!
                */
                if (uNewRung < wp->ts->uRungLo) uNewRung = wp->ts->uRungLo;
                else if (uNewRung > wp->ts->uRungHi) uNewRung = wp->ts->uRungHi;
                if (!pkd->bNoParticleOrder) p.set_new_rung(uNewRung);
                /*
                ** Now we want to kick the particle velocities with a closing kick
                ** based on the old rung and an opening kick based on the new rung.
                ** However, we are not always allowed to decrease uNewRung by an
                ** arbitrary amount, as this depends on where we are in the
                ** substep hierarchy.
                ** This also requires the various timestep factors to have been
                ** pre-calculated for this dTime and all possible rungs at play.
                ** Note that the only persistent info needed is the old rung which
                ** gets updated at the end of this so that p.uRung indicates
                ** which particles are active in this time step.
                ** p.uRung = uNewRung could be done here, but we let the outside
                ** code handle this.
                */
                if (pkd->particles.present(PKD_FIELD::oVelocity)) {
                    auto &v = p.velocity();
                    if (wp->kick && wp->kick->bKickClose) {
                        v[0] += wp->kick->dtClose[p.rung()]*wp->pInfoOut[i].a[0];
                        v[1] += wp->kick->dtClose[p.rung()]*wp->pInfoOut[i].a[1];
                        v[2] += wp->kick->dtClose[p.rung()]*wp->pInfoOut[i].a[2];
                        if (wp->SPHoptions->doSPHForces) {
                            auto &NewSph = p.newsph();
                            if (wp->SPHoptions->useIsentropic) {
                                NewSph.u = SPHEOSIsentropic(pkd,NewSph.oldRho,NewSph.u,p.density(),p.imaterial(),wp->SPHoptions);
                                NewSph.oldRho = p.density();
                            }
                            NewSph.u += wp->kick->dtClose[p.rung()] * NewSph.uDot;
                            if (!wp->SPHoptions->doOnTheFlyPrediction && !wp->SPHoptions->doConsistentPrediction) {
                                NewSph.P = SPHEOSPCTofRhoU(pkd,p.density(),NewSph.u,&NewSph.cs,&NewSph.T,p.imaterial(),wp->SPHoptions);
                            }
                        }
                    }
                    v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
                    /*
                    ** Now calculate the kinetic energy term.
                    */
                    pkd->dEnergyT += 0.5*m*v2;
                    /* L is calculated with respect to the origin (0,0,0) */
                    pkd->dEnergyL[0] += m*(r[1]*v[2] - r[2]*v[1]);
                    pkd->dEnergyL[1] += m*(r[2]*v[0] - r[0]*v[2]);
                    pkd->dEnergyL[2] += m*(r[0]*v[1] - r[1]*v[0]);

                    if (wp->kick && wp->kick->bKickOpen) {
                        p.set_rung(uNewRung);
                        ++pkd->nRung[p.rung()];
                        if (wp->SPHoptions->VelocityDamper > 0.0f) {
                            v[0] *= exp(- wp->kick->dtOpen[p->uRung] * wp->SPHoptions->VelocityDamper);
                            v[1] *= exp(- wp->kick->dtOpen[p->uRung] * wp->SPHoptions->VelocityDamper);
                            v[2] *= exp(- wp->kick->dtOpen[p->uRung] * wp->SPHoptions->VelocityDamper);
                        }
                        v[0] += wp->kick->dtOpen[p.rung()]*wp->pInfoOut[i].a[0];
                        v[1] += wp->kick->dtOpen[p.rung()]*wp->pInfoOut[i].a[1];
                        v[2] += wp->kick->dtOpen[p.rung()]*wp->pInfoOut[i].a[2];
                        if (wp->SPHoptions->doSPHForces) {
                            auto &NewSph = p.newsph();
                            NewSph.u += wp->kick->dtOpen[p.rung()] * NewSph.uDot;
                        }
                        /*
                        ** On KickOpen we also always check for intersection with the lightcone
                        ** surface over the entire next timestep of the particle (not a half
                        ** timestep as is usual for kicking (we are drifting afterall).
                        */
                        if (wp->lc->dLookbackFac > 0) {
                            pkdProcessLightCone(pkd,&p,wp->pInfoOut[i].fPot,wp->lc->dLookbackFac,wp->lc->dLookbackFacLCP,
                                                wp->lc->dtLCDrift[p.rung()],wp->lc->dtLCKick[p.rung()],
                                                wp->lc->dBoxSize,wp->lc->bLightConeParticles,wp->lc->hLCP,wp->lc->tanalpha_2);
                        }
                    }
                    p.set_marked(true);
                }
            }
            else {
                ++pkd->nRung[p.rung()];
            }
            auto q = pkd->particles[static_cast<PARTICLE *>(mdlAcquireWrite(pkd->mdl,CID_PARTICLE,wp->iPart[i]))];
            q = p;
            mdlReleaseWrite(pkd->mdl,CID_PARTICLE,&q);
        }
        delete [] wp->pPart;
        delete [] wp->iPart;
        delete [] wp->pInfoIn;
        delete [] wp->pInfoOut;
        delete wp;
    }
}

static void queuePP( PKD pkd, workParticle *wp, ilpList &ilp, int bGravStep ) {
    for ( auto &tile : ilp ) {
#ifdef USE_CUDA
        if (pkd->cudaClient->queuePP(wp,tile,bGravStep)) continue;
#endif
#ifdef USE_METAL
        if (pkd->metalClient->queuePP(wp,tile,bGravStep)) continue;
#endif
        for (auto i=0; i<wp->nP; ++i) {
            pkdGravEvalPP(wp->pInfoIn[i],tile,wp->pInfoOut[i]);
            wp->dFlopSingleCPU += COST_FLOP_PP*tile.size();
        }
    }
}

static void queueDensity( PKD pkd, workParticle *wp, ilpList &ilp, int bGravStep ) {
    wp->ilp = &ilp;
    wp->bGravStep = bGravStep;
    for ( int i=0; i<wp->nP; i++ ) {
        wp->pInfoOut[i].rho = 0.0f;
        wp->pInfoOut[i].drhodfball = 0.0f;
        wp->pInfoOut[i].nden = 0.0f;
        wp->pInfoOut[i].dndendfball = 0.0f;
        wp->pInfoOut[i].nSmooth = 0.0f;
        wp->pInfoOut[i].imbalanceX = 0.0f;
        wp->pInfoOut[i].imbalanceY = 0.0f;
        wp->pInfoOut[i].imbalanceZ = 0.0f;
    }
    for ( auto &tile : ilp ) {
        for (auto i=0; i<wp->nP; ++i) {
            pkdDensityEval(wp->pInfoIn[i],tile,wp->pInfoOut[i],wp->SPHoptions);
        }
    }
}

static void queueDensityCorrection( PKD pkd, workParticle *wp, ilpList &ilp, int bGravStep ) {
    for ( auto &tile : ilp ) {
        for (auto i=0; i<wp->nP; ++i) {
            pkdDensityCorrectionEval(wp->pInfoIn[i],tile,wp->pInfoOut[i],wp->SPHoptions);
        }
    }
}

static void queueSPHForces( PKD pkd, workParticle *wp, ilpList &ilp, int bGravStep ) {
    for ( auto &tile : ilp ) {
        for (auto i=0; i<wp->nP; ++i) {
            pkdSPHForcesEval(wp->pInfoIn[i],tile,wp->pInfoOut[i],wp->SPHoptions);
        }
    }
}

static void addCentrifugalAcceleration(PKD pkd, workParticle *wp) {
    double dOmega, f;
    double r[3], rxy;

    if (wp->ts->dTime < wp->SPHoptions->CentrifugalT0) {
        return;
    }
    else if (wp->ts->dTime < wp->SPHoptions->CentrifugalT1) {
        dOmega = wp->SPHoptions->CentrifugalOmega0*(wp->ts->dTime - wp->SPHoptions->CentrifugalT0)/(wp->SPHoptions->CentrifugalT1 - wp->SPHoptions->CentrifugalT0);
    }
    else {
        dOmega = wp->SPHoptions->CentrifugalOmega0;
    }

    f = dOmega * dOmega;
    for (int i=0; i< wp->nP; i++) {
        pkdGetPos1(pkd,wp->pPart[i],r);
        rxy = sqrt(r[0] * r[0] + r[1] + r[1]);
        wp->pInfoOut[i].a[0] += f * r[0];
        wp->pInfoOut[i].a[1] += f * r[1];
    }
}

static void queuePC( PKD pkd,  workParticle *wp, ilcList &ilc, int bGravStep ) {
    for ( auto &tile : ilc ) {
#ifdef USE_CUDA
        if (pkd->cudaClient->queuePC(wp,tile,bGravStep)) continue;
#endif
#ifdef USE_METAL
        if (pkd->metalClient->queuePC(wp,tile,bGravStep)) continue;
#endif
        for (auto i=0; i<wp->nP; ++i) {
            pkdGravEvalPC(wp->pInfoIn[i],tile,wp->pInfoOut[i]);
            wp->dFlopSingleCPU += COST_FLOP_PC*tile.size();
        }
    }
}

static void queueEwald( PKD pkd, workParticle *wp ) {
    int i;
#ifdef USE_CUDA
    int nQueued = pkd->cudaClient->queueEwald(wp);
#else
    int nQueued = 0;
#endif
    ++wp->nRefs;
    for ( i=nQueued; i<wp->nP; ++i) {
        PINFOIN *in = &wp->pInfoIn[i];
        PINFOOUT *out = &wp->pInfoOut[i];
        double r[3];
        r[0] = wp->c[0] + in->r[0];
        r[1] = wp->c[1] + in->r[1];
        r[2] = wp->c[2] + in->r[2];
        wp->dFlop += pkdParticleEwald(pkd,r,out->a.data(),&out->fPot,&wp->dFlopSingleCPU,&wp->dFlopDoubleCPU);
    }
    pkdParticleWorkDone(wp);
}

/*
** This version of grav.c does all the operations inline, including
** v_sqrt's and such.
** Returns nActive.
*/

int pkdGravInteract(PKD pkd,
                    struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
                    treeStore::NodePointer pkdn,LOCR *pLoc,ilpList &ilp,ilcList &ilc,
                    float dirLsum,float normLsum,int bEwald,double *pdFlop,
                    SMX smx,SMF *smf,int iRoot1,int iRoot2,SPHOptions *SPHoptions) {
    float fBall;
    int i,nSoft,nActive;
    int nP;

    auto kdn_r = pkdn->position();

    /*
    ** Now process the two interaction lists for each active particle.
    */
    nActive = 0;
    nSoft = 0;

    /* Collect the bucket particle information */
    auto wp = new workParticle;
    assert(wp!=NULL);
    wp->nRefs = 1; /* I am using it currently */
    wp->ctx = pkd;
    wp->dFlop = 0.0;
    wp->dFlopSingleCPU = wp->dFlopSingleGPU = 0.0;
    wp->dFlopDoubleCPU = wp->dFlopDoubleGPU = 0.0;
    wp->nP = 0;

    nP = pkdn->count();
    wp->pPart = new PARTICLE *[nP];
    wp->iPart = new uint32_t[nP];
    wp->pInfoIn = new PINFOIN[nP];
    wp->pInfoOut = new PINFOOUT[nP];
    wp->c[0] = ilp.getReference(0); assert(wp->c[0] == ilc.getReference(0));
    wp->c[1] = ilp.getReference(1); assert(wp->c[1] == ilc.getReference(1));
    wp->c[2] = ilp.getReference(2); assert(wp->c[2] == ilc.getReference(2));

    wp->dirLsum = dirLsum; // not used here...
    wp->normLsum = normLsum; // not used here...

    wp->ts = ts;
    wp->lc = lc;
    wp->kick = kick;

    wp->SPHoptions = SPHoptions;

    for (auto &p : *pkdn) {
#ifdef NN_FLAG_IN_PARTICLE
        if (!p.is_rung_range(ts->uRungLo,ts->uRungHi) && !(SPHoptions->useDensityFlags && p.marked()) && !(SPHoptions->useNNflags && p.NN_flag())) continue;
#else
        if (!p.is_rung_range(ts->uRungLo,ts->uRungHi) && !(SPHoptions->useDensityFlags && p.marked())) continue;
#endif
        auto r = p.position();

        nP = wp->nP++;
        wp->pPart[nP] = &p;
        wp->iPart[nP] = &p - pkd->particles.begin();
        wp->pInfoIn[nP].r = pkd->ilc.get_dr(r);
        if (p.have_acceleration()) {
            const auto &ap = p.acceleration();
            wp->pInfoIn[nP].a[0]  = ap[0];
            wp->pInfoIn[nP].a[1]  = ap[1];
            wp->pInfoIn[nP].a[2]  = ap[2];
            //ap[0] = ap[1] = ap[2] = 0.0;
        }
        else {
            wp->pInfoIn[nP].a[0]  = 0;
            wp->pInfoIn[nP].a[1]  = 0;
            wp->pInfoIn[nP].a[2]  = 0;
        }

        if (p.have_newsph()) {
            const auto &NewSph = p.newsph();
            wp->pInfoIn[nP].fBall = p.ball();
            wp->pInfoIn[nP].isTooLarge = 0;
            wp->pInfoIn[nP].Omega = NewSph.Omega;
            wp->pInfoIn[nP].iMat = p.imaterial();
            SPHpredictOnTheFly(pkd, p, kick, ts->uRungLo, wp->pInfoIn[nP].v, &wp->pInfoIn[nP].P, &wp->pInfoIn[nP].cs, NULL, SPHoptions);
            wp->pInfoIn[nP].rho = p.density();
            wp->pInfoIn[nP].species = p.species();

            wp->pInfoOut[nP].rho = 0.0f;
            wp->pInfoOut[nP].drhodfball = 0.0f;
            wp->pInfoOut[nP].nden = 0.0f;
            wp->pInfoOut[nP].dndendfball = 0.0f;
            wp->pInfoOut[nP].fBall = 0.0f;
            wp->pInfoOut[nP].nSmooth = 0.0f;
            wp->pInfoOut[nP].imbalanceX = 0.0f;
            wp->pInfoOut[nP].imbalanceY = 0.0f;
            wp->pInfoOut[nP].imbalanceZ = 0.0f;
            wp->pInfoOut[nP].uDot = 0.0f;
            wp->pInfoOut[nP].divv = 0.0f;
            wp->pInfoOut[nP].dtEst = HUGE_VALF;
            wp->pInfoOut[nP].maxRung = 0.0f;
            wp->pInfoOut[nP].corrT = 0.0f;
            wp->pInfoOut[nP].corrP = 0.0f;
            wp->pInfoOut[nP].corr = 0.0f;
        }

        wp->pInfoOut[nP].a[0] = 0.0f;
        wp->pInfoOut[nP].a[1] = 0.0f;
        wp->pInfoOut[nP].a[2] = 0.0f;
        wp->pInfoOut[nP].fPot = 0.0f;
        wp->pInfoOut[nP].dirsum = dirLsum;
        wp->pInfoOut[nP].normsum = normLsum;
        wp->pInfoOut[nP].rhopmax = 0.0f;

        /*
        ** Calculate local density and kernel smoothing length for dynamical time-stepping
        */
        if (SPHoptions->doGravity) {
            if (ts->bGravStep) {
                /*
                ** Calculate local density using smooth; this is fast because the particles are
                ** likely to be cached already because they will be on the P-P list.
                */
                smf->pfDensity = &wp->pInfoIn[nP].fDensity;
                fBall = smSmoothSingle(smx,smf,p,iRoot1,iRoot2);
                wp->pInfoIn[nP].fSmooth2 = fBall * fBall;
            }
            else {
                /*
                ** We are not using GravStep!
                */
                wp->pInfoIn[nP].fSmooth2 = 0.0;
            }
        }
    }

    nActive += wp->nP;

    /*
    ** Evaluate the local expansion.
    */
    if (pLoc && SPHoptions->doGravity) {
        for ( i=0; i<wp->nP; i++ ) {
            momFloat ax,ay,az, dPot;
            blitz::TinyVector<double,3> r = wp->c + wp->pInfoIn[i].r;

            /*
            ** Evaluate local expansion.
            */
            blitz::TinyVector<double,3> dr = r - kdn_r;
            dPot = 0;
            ax = 0;
            ay = 0;
            az = 0;

            momEvalLocr(pLoc,dr[0],dr[1],dr[2],&dPot,&ax,&ay,&az);

            wp->pInfoOut[i].fPot += dPot;
            wp->pInfoOut[i].a[0] += ax;
            wp->pInfoOut[i].a[1] += ay;
            wp->pInfoOut[i].a[2] += az;
        }
    }

    if (SPHoptions->doGravity) {
        /*
        ** Evaluate the P-C interactions
        */
        queuePC( pkd,  wp, pkd->ilc, ts->bGravStep );

        /*
        ** Evaluate the P-P interactions
        */
        queuePP( pkd, wp, pkd->ilp, ts->bGravStep );
    }

    if (SPHoptions->doDensity) {
        /*
        ** Evaluate the Density on the P-P interactions
        */
        queueDensity( pkd, wp, pkd->ilp, ts->bGravStep );
    }

    if (SPHoptions->doDensityCorrection) {
        /*
        ** Evaluate the weighted averages of P and T
        */
        queueDensityCorrection( pkd, wp, pkd->ilp, ts->bGravStep );
    }

    if (SPHoptions->doSPHForces) {
        /*
        ** Evaluate the SPH forces on the P-P interactions
        */
        queueSPHForces( pkd, wp, pkd->ilp, ts->bGravStep );
        if (SPHoptions->doCentrifugal) {
            addCentrifugalAcceleration(pkd, wp);
        }
    }

    if (SPHoptions->doGravity) {
        /*
        ** Calculate the Ewald correction for this particle, if it is required.
        */
        if (bEwald) {
            queueEwald( pkd,  wp );
        }
    }

#ifdef TIMESTEP_CRITICAL
    if (SPHoptions->doGravity) {
        for ( i=0; i<wp->nP; i++ ) {
            double *c = wp->c;
            float *in = wp->pInfoIn[i].r;
            r[0] = c[0] + in[0];
            r[1] = c[1] + in[1];
            r[2] = c[2] + in[2];
            //p = wp->pPart[i];
            //pkdGetPos1(p.r,r);
            float fMass = p.mass();
            float fSoft = p.soft();
            const auto &v = p.velocity();
            /*
            ** Set value for time-step, note that we have the current ewald acceleration
            ** in this now as well!
            */
            if (ts->bGravStep && wp->ts->iTimeStepCrit == 1) {
                double vx,vy,vz;
                float d2,dir,dir2;
                float fx,fy,fz;
                float rhopmax,rhopmaxlocal;
                float summ;
                float fourh2;
                int n;
                ILPTILE tile;

                /*
                ** GravStep if iTimeStepCrit =
                ** 0: Mean field regime for dynamical time (normal/standard setting)
                ** 1: Gravitational scattering regime for dynamical time with eccentricity correction
                */
                rhopmax = 0.0;
                ILP_LOOP(ilp,tile) {
                    int blk,prt;
                    for ( blk=0; blk<=tile->lstTile.nBlocks; ++blk ) {
                        n = (blk==tile->lstTile.nBlocks ? tile->lstTile.nInLast : ilp->lst.nPerBlock);
                        for (prt=0; prt<n; ++prt) {
                            if (p.order() < wp->ts->nPartColl || tile->xtr[blk].iOrder.i[prt] < wp->ts->nPartColl) {
                                fx = r[0] - ilp->cx + tile->blk[blk].dx.f[prt];
                                fy = r[1] - ilp->cy + tile->blk[blk].dy.f[prt];
                                fz = r[2] - ilp->cz + tile->blk[blk].dz.f[prt];
                                d2 = fx*fx + fy*fy + fz*fz;
                                if (p.order() == tile->xtr[blk].iOrder.i[prt]) continue;
                                fourh2 = softmassweight(fMass,4*fSoft*fSoft,
                                                        tile->blk[blk].m.f[prt],tile->blk[blk].fourh2.f[prt]);
                                if (d2 > fourh2) {
                                    SQRT1(d2,dir);
                                    dir2 = dir*dir*dir;
                                }
                                else {
                                    /*
                                    ** This uses the Dehnen K1 kernel function now, it's fast!
                                    */
                                    SQRT1(fourh2,dir);
                                    dir2 = dir*dir;
                                    d2 *= dir2;
                                    dir2 *= dir;
                                    d2 = 1 - d2;
                                    dir *= 1.0f + d2*(0.5f + d2*(3.0f/8.0f + d2*(45.0f/32.0f)));
                                    dir2 *= 1.0f + d2*(1.5f + d2*(135.0f/16.0f));
                                }
                                summ = fMass+tile->blk[blk].m.f[prt];
                                rhopmaxlocal = summ*dir2;
                                vx = v[0] - tile->xtr[blk].vx.d[prt];
                                vy = v[1] - tile->xtr[blk].vy.d[prt];
                                vz = v[2] - tile->xtr[blk].vz.d[prt];
                                rhopmaxlocal = pkdRho1(rhopmaxlocal,summ,dir,
                                                       fx,fy,fz,vx,vy,vz,wp->ts->dEccFacMax);
                                rhopmax = (rhopmaxlocal > rhopmax)?rhopmaxlocal:rhopmax;
                            }
                        }
                    }
                }
                wp->pInfoOut[i].rhopmax = rhopmax;
            }
        } /* end of i-loop cells & particles */
    }
#endif

    pkdParticleWorkDone(wp);

    *pdFlop += nActive*(pkd->ilp.count()*COST_FLOP_PP + pkd->ilc.count()*COST_FLOP_PC) + nSoft*15;
    return (nActive);
}

#ifdef TIMESTEP_CRITICAL
/*
** Gravitational scattering regime (iTimeStepCrit=1)
*/
double pkdRho1(double rhopmaxlocal, double summ, double dir, double x, double y, double z, double vx, double vy, double vz, double EccFacMax) {

    double Etot, L2, ecc, eccfac, v2;
    /*
    ** Etot and L are normalized by the reduced mass
    */
    v2 = vx*vx + vy*vy + vz*vz;
    Etot = 0.5*v2 - summ*dir;
    L2 = (y*vz - z*vy)*(y*vz - z*vy) + (z*vx - x*vz)*(z*vx - x*vz) + (x*vy - y*vx)*(x*vy - y*vx);
    ecc = 1+2*Etot*L2/(summ*summ);
    ecc = (ecc <= 0)?0:sqrt(ecc);
    eccfac = (1 + 2*ecc)/fabs(1-ecc);
    eccfac = (eccfac > EccFacMax)?EccFacMax:eccfac;
    if (eccfac > 1.0) rhopmaxlocal *= eccfac;
    return rhopmaxlocal;
}
#endif
