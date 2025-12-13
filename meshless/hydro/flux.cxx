#include "eEOS/eEOS.h"
#include "hydro/hydro.h"
#include "hydro/limiters.h"
#include "interaction/mfm.h"
#include "master.h"
#include "walk/hydrowalk.h"

using blitz::dot;
using blitz::TinyVector;

void MSR::MeshlessFluxes(double dTime, double dDelta, uint8_t uRungLo,
                         uint8_t uRungHi) {
    double dsec;
    printf("Computing fluxes... ");
    TimerStart(TIMER_FLUXES);

    // /*
    // Mark those particles that are part of an interaction, either because
    // they are active or a neighbour of an active particle
    ReSmooth(dTime, 0.0f, SMX_MARK_ACTIVE_NN, 1);

    // Update the tree accordingly to have a correct "is_marked" for tree nodes
    SPHOptions SPHoptions;
    SPHoptions.dofBallFactor = false;
    SPHoptions.ballSizeLimit = INFINITY;
    TreeUpdateFlagBounds(parameters.get_bEwald(),ROOT,0,SPHoptions);

    HydroWalk<HydroInteraction>::OptionPack input;

    input.inter.uRungHi = uRungHi;
    input.inter.uRungLo = uRungLo;
    input.inter.bPeriodic = parameters.get_bPeriodic();
    input.inter.a = csmTime2Exp(csm, dTime);
    input.inter.H = csmTime2Hub(csm, dTime);
    input.inter.gamma = parameters.get_dConstGamma();
    input.inter.dTime = dTime;
    input.inter.dDelta = dDelta;
    input.inter.doExtensiveTest = false;

    input.walk.uRungHi = uRungHi;
    input.walk.uRungLo = uRungLo;
    input.walk.dTheta = 0.7; // Not really relevant
    input.walk.nReps = 1;
    input.walk.iRoot1 = ROOT;
    input.walk.iRoot2 = 0;
    input.walk.dofBallFactor = false;
    input.walk.fBallFactor = 1.0;
    input.walk.ballSizeLimit = INFINITY;

    HydroWalkTraits::output output;

    mdl->RunService("HydroWalkHydroInteraction", sizeof(input), &input,
                    &output);

    PrintStat(output.sLocal, "  particle  load:", 0);
    PrintStat(output.sActive, "  actives   load:", 0);
    PrintStat(output.sFlop, "  Gflop     load:", 1);
    PrintStat(output.sPart, "  P-P per active:", 2);
    PrintStat(output.sCell, "  P-C per active:", 2);
#ifdef INSTRUMENT
    PrintStat(output.sComputing, "     % computing:", 3);
    PrintStat(output.sWaiting, "     %   waiting:", 3);
    PrintStat(output.sSynchronizing, "     %   syncing:", 3);
#endif
#ifdef __linux__
    PrintStat(output.sFreeMemory, "free memory (GB):", 3);
    PrintStat(output.sRSS, "   resident size:", 3);
#endif
    print("  (cache access statistics are given per active particle)\n");
    PrintStat(output.sPartNumAccess, "  P-cache access:", 1);
    PrintStat(output.sCellNumAccess, "  C-cache access:", 1);
    PrintStat(output.sPartMissRatio, "  P-cache miss %:", 2);
    PrintStat(output.sCellMissRatio, "  C-cache miss %:", 2);
    print(
        "Total tiles processed: {:.5e}, on the GPU: {:.5e}, ratio: {:.2f} %\n",
        (double)output.nTilesTotal,
        (double)(output.nTilesTotal - output.nTilesCPU),
        100.0 -
        ((double)output.nTilesCPU) / ((double)output.nTilesTotal) * 100.0);
    // */

    /*
    #ifdef OPTIM_SMOOTH_NODE
    #ifdef OPTIM_AVOID_IS_ACTIVE
    SelActives();
    #endif
    #ifdef OPTIM_FLUX_VEC
    ReSmoothNode(dTime, dDelta, SMX_HYDRO_FLUX_VEC,1);
    #else
    ReSmoothNode(dTime, dDelta, SMX_HYDRO_FLUX,1);
    #endif
    #else // no OPTIM_SMOOTH_NODE
    ReSmooth(dTime, dDelta, SMX_HYDRO_FLUX,1);
    #endif
    */
    TimerStop(TIMER_FLUXES);
    dsec = TimerGet(TIMER_FLUXES);
    printf("took %.5f seconds\n", dsec);
}

void SmoothHelperHydroFluxes::pack(void *dst, const void *src) {
    auto p1 = static_cast<pack_type *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        const auto &sph = p2.sph();

        p1->position = p2.position();
        p1->velocity = p2.velocity();

        p1->B = sph.B;
        p1->gradRho = sph.gradRho;
        p1->gradVx = sph.gradVx;
        p1->gradVy = sph.gradVy;
        p1->gradVz = sph.gradVz;
        p1->gradP = sph.gradP;
        p1->lastUpdateTime = sph.lastUpdateTime;
        p1->lastAcc = sph.lastAcc;
        p1->omega = sph.omega;
        p1->spheom = sph.spheom;
        p1->P = sph.P;

        p1->fBall = p2.ball();
        p1->fDensity = p2.density();
        p1->uRung = p2.rung();
        p1->bMarked = p2.marked();
    }
}

void SmoothHelperHydroFluxes::unpack(void *dst, const void *src,
                                     const void *key) {
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const pack_type *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_adj_decorator();

        auto &sph = p1.sph();

        p1.set_position(p2->position);
        p1.velocity() = p2->velocity;

        sph.B = p2->B;
        sph.gradRho = p2->gradRho;
        sph.gradVx = p2->gradVx;
        sph.gradVy = p2->gradVy;
        sph.gradVz = p2->gradVz;
        sph.gradP = p2->gradP;
        sph.lastUpdateTime = p2->lastUpdateTime;
        sph.lastAcc = p2->lastAcc;
        sph.omega = p2->omega;
        sph.spheom = p2->spheom;
        sph.P = p2->P;

        p1.set_ball(p2->fBall);
        p1.set_density(p2->fDensity);
        p1.set_rung(p2->uRung);
        p1.set_marked(p2->bMarked);
    }
}

/* Zero all the conserved quantities in cached copies, which will be updated
 * during the hydro loop. They will be merged with the actual particle
 * information in combHydroFluxes
 */
void SmoothHelperHydroFluxes::init(void *dst) {
    auto p = pkd->particles[static_cast<PARTICLE *>(dst)];
    // For the init*Cached and comb functions we still have to explicitly
    // check if we are handling gas particles even ifdef OPTIM_REORDER_IN_NODES
    // because these operations are done in a cache-line basis, and a line may
    // contain other particles that are not of interest!
    if (p.is_gas()) {
        // assert((int)((char*)p.decorator_field()-(char*)dst) ==
        // pkd->particles.ParticleSize());
        auto &sph = p.sph();

        sph.Frho = 0.0;
        sph.Fmom = 0.0;
        sph.Fene = 0.0;
        sph.mom = 0.0;
        sph.E = 0.0;
        sph.Uint = 0.0;
        if (p.have_mfv()) {
            auto &mfv = p.mfv();
            mfv.drDotFrho = 0.0;
            p.set_mass(0.0);
        }
    }
}

void SmoothHelperHydroFluxes::flush(void *dst, const void *src) {
    auto p1 = static_cast<flush_type *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    if (p2.is_gas()) {
        const auto &sph = p2.sph();
        // assert((int)((char*)p2.decorator_field()-(char*)src) ==
        // pkd->particles.ParticleSize());

        p1->Frho = sph.Frho;
        p1->Fmom = sph.Fmom;
        p1->Fene = sph.Fene;
        p1->mom = sph.mom;
        p1->E = sph.E;
        p1->Uint = sph.Uint;
        if (p2.have_mfv()) {
            const auto &mfv = p2.mfv();
            p1->drDotFrho = mfv.drDotFrho;
            p1->fMass = p2.mass();
        }
    }
}

void SmoothHelperHydroFluxes::combine(void *dst, const void *src,
                                      const void *key) {
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const flush_type *>(src);

    if (p1.is_gas()) {
        auto &sph = p1.sph();

        sph.Frho += p2->Frho;
        sph.Fmom += p2->Fmom;
        sph.Fene += p2->Fene;
        sph.mom += p2->mom;
        sph.E += p2->E;
        sph.Uint += p2->Uint;
        if (p1.have_mfv()) {
            auto &mfv = p1.mfv();
            mfv.drDotFrho += p2->drDotFrho;
            p1.set_mass(p1.mass() + p2->fMass);
        }
    }
}


void hydroFluxFillBuffer(meshless::myreal *input_buffer, PARTICLE *qIn, int i,
                         int nBuff, double dr2, blitz::TinyVector<double, 3> dr,
                         SMF *smf) {
    PKD pkd = smf->pkd;
    auto Q = pkd->particles[qIn];
    double dDelta = smf->dDelta;
    auto &qsph = Q.sph();
    q(ball) = Q.ball();
    q(dx) = dr[0];
    q(dy) = dr[1];
    q(dz) = dr[2];
    q(dr) = sqrt(dr2);
    q(rung) = dDelta / (1 << Q.rung());
    q(rho) = Q.density();
    q(P) = qsph.P;
#ifdef ENTROPY_SWITCH
    q(S) = qsph.S;
#endif
    const auto &qv = Q.velocity();
    q(vx) = qv[0];
    q(vy) = qv[1];
    q(vz) = qv[2];

    q(gradRhoX) = qsph.gradRho[0];
    q(gradRhoY) = qsph.gradRho[1];
    q(gradRhoZ) = qsph.gradRho[2];

    q(gradPX) = qsph.gradP[0];
    q(gradPY) = qsph.gradP[1];
    q(gradPZ) = qsph.gradP[2];

    q(gradVxX) = qsph.gradVx[0];
    q(gradVxY) = qsph.gradVx[1];
    q(gradVxZ) = qsph.gradVx[2];

    q(gradVyX) = qsph.gradVy[0];
    q(gradVyY) = qsph.gradVy[1];
    q(gradVyZ) = qsph.gradVy[2];

    q(gradVzX) = qsph.gradVz[0];
    q(gradVzY) = qsph.gradVz[1];
    q(gradVzZ) = qsph.gradVz[2];

    q(lastUpdateTime) = qsph.lastUpdateTime;
    q(lastAccX) = qsph.lastAcc[0];
    q(lastAccY) = qsph.lastAcc[1];
    q(lastAccZ) = qsph.lastAcc[2];
    q(B_XX) = qsph.B[XX];
    q(B_YY) = qsph.B[YY];
    q(B_ZZ) = qsph.B[ZZ];
    q(B_XY) = qsph.B[XY];
    q(B_XZ) = qsph.B[XZ];
    q(B_YZ) = qsph.B[YZ];
    q(omega) = qsph.omega;
    q(spheom) = qsph.spheom;
}

void hydroFluxUpdateFromBuffer(meshless::myreal *output_buffer,
                               meshless::myreal *input_buffer, PARTICLE *pIn,
                               PARTICLE *qIn, int i, int nBuff, SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[pIn];
    auto Q = pkd->particles[qIn];
    const auto &pv = P.velocity();
    const auto &qv = Q.velocity();
    auto &psph = P.sph();
    auto &qsph = Q.sph();
    const auto &dDelta = smf->dDelta;
    const auto &aFac = smf->a;
    if (dDelta > 0) {
        if (P.have_mfv()) {
            P.set_mass(P.mass() - qout(minDt) * qout(Frho));
        }

        psph.mom[0] -= qout(minDt) * qout(FmomX);
        psph.mom[1] -= qout(minDt) * qout(FmomY);
        psph.mom[2] -= qout(minDt) * qout(FmomZ);

        psph.E -= qout(minDt) * qout(Fene);

        psph.Uint -= qout(minDt) *
                     (qout(Fene) - qout(FmomX) * pv[0] - qout(FmomY) * pv[1] -
                      qout(FmomZ) * pv[2] + 0.5 * dot(pv, pv) * qout(Frho));

#ifdef ENTROPY_SWITCH
        psph.S -= qout(minDt) * qout(FS);
#endif

        if (P.have_mfv()) {
            auto &mfv = P.mfv();
            mfv.drDotFrho[0] += qout(minDt) * qout(Frho) * q(dx) * aFac;
            mfv.drDotFrho[1] += qout(minDt) * qout(Frho) * q(dy) * aFac;
            mfv.drDotFrho[2] += qout(minDt) * qout(Frho) * q(dz) * aFac;
        }
        psph.Frho += qout(Frho);
        psph.Fene += qout(Fene);
        psph.Fmom[0] += qout(FmomX);
        psph.Fmom[1] += qout(FmomY);
        psph.Fmom[2] += qout(FmomZ);
    }
    else {
        psph.Frho += qout(Frho);
        psph.Fene += qout(Fene);
        psph.Fmom[0] += qout(FmomX);
        psph.Fmom[1] += qout(FmomY);
        psph.Fmom[2] += qout(FmomZ);
    }

#ifndef OPTIM_NO_REDUNDANT_FLUXES
#ifdef OPTIM_AVOID_IS_ACTIVE
    if (!Q.marked())
#else
    if (!Q.is_active())
#endif
#endif
    {

        // If this is not the case, something VERY odd must have happened
        assert(qsph.P == q(P));
        if (dDelta > 0) {
            if (Q.have_mfv()) {
                Q.set_mass(Q.mass() + qout(minDt) * qout(Frho));
            }

            qsph.mom[0] += qout(minDt) * qout(FmomX);
            qsph.mom[1] += qout(minDt) * qout(FmomY);
            qsph.mom[2] += qout(minDt) * qout(FmomZ);

            qsph.E += qout(minDt) * qout(Fene);

            qsph.Uint +=
                qout(minDt) *
                (qout(Fene) - qout(FmomX) * qv[0] - qout(FmomY) * qv[1] -
                 qout(FmomZ) * qv[2] + 0.5 * dot(qv, qv) * qout(Frho));
#ifdef ENTROPY_SWITCH
            qsph.S += qout(minDt) * qout(FS);
#endif

            if (Q.have_mfv()) {
                auto &mfv = P.mfv();
                mfv.drDotFrho[0] += qout(minDt) * qout(Frho) * q(dx) * aFac;
                mfv.drDotFrho[1] += qout(minDt) * qout(Frho) * q(dy) * aFac;
                mfv.drDotFrho[2] += qout(minDt) * qout(Frho) * q(dz) * aFac;
            }
            qsph.Frho -= qout(Frho);
            qsph.Fene -= qout(Fene);
            qsph.Fmom[0] -= qout(FmomX);
            qsph.Fmom[1] -= qout(FmomY);
            qsph.Fmom[2] -= qout(FmomZ);
        }
        else {
            qsph.Frho -= qout(Frho);
            qsph.Fene -= qout(Fene);
            qsph.Fmom[0] -= qout(FmomX);
            qsph.Fmom[1] -= qout(FmomY);
            qsph.Fmom[2] -= qout(FmomZ);
        }

    } // q marked/active
}

void hydroFluxGetBufferInfo(int *in, int *out) {
    *in = q_last;
    *out = out_last;
}

//#endif // OPTIM_FLUX_VEC

void hydroRiemann_wrapper(PARTICLE *p, float fBall, int nSmooth, int nBuff,
                          meshless::myreal *restrict input_buffer,
                          meshless::myreal *restrict output_buffer, SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[p];

#if defined(USE_SIMD_FLUX)
    MeshlessHydroSolver<dvec, dmask> solver(P.have_mfv());
#else
    MeshlessHydroSolver<vec<double, double>, mmask<bool>> solver(P.have_mfv());
#endif
    solver.hydroRiemann(P, fBall, nSmooth, nBuff, input_buffer, output_buffer,
                        pkd->csm->val.bComove, smf->dTime, smf->dDelta, smf->a,
                        smf->H, smf->dConstGamma, smf->eEOS);
}
