#include "hydro/hydro.h"
#include "master.h"
#include "eEOS/eEOS.h"
#ifdef OPTIM_FLUX_VEC
    #include "riemann_own.h"
#else
    #include "riemann.h"
#endif
using blitz::TinyVector;
using blitz::dot;

/* TODO: This function probably can be eliminated, and the fluxes reset when
 *   starting the ReSmooth
 */
void MSR::ResetFluxes(double dTime,double dDelta) {
    struct inDrift in;

    if (csm->val.bComove) {
        in.dDelta = csmComoveDriftFac(csm,dTime,dDelta);
        in.dDeltaVPred = csmComoveKickFac(csm,dTime,dDelta);
    }
    else {
        in.dDelta = dDelta;
        in.dDeltaVPred = dDelta;
    }
    in.dTime = dTime;
    in.dDeltaUPred = dDelta;
    pstResetFluxes(pst,&in,sizeof(in),NULL,0);
}

void pkdResetFluxes(PKD pkd, double dTime,double dDelta,double dDeltaVPred,double dDeltaTime) {
    /*
    ** Add the computed flux to the conserved variables for each gas particle
    */
    assert(pkd->particles.present(PKD_FIELD::oVelocity));
    assert(pkd->particles.present(PKD_FIELD::oSph));
    for (auto &P : pkd->particles) {
        if (P.is_gas() && P.is_active() ) {
            auto &sph = P.sph();
            sph.Frho = 0.0;
            sph.Fene = 0.0;
            sph.Fmom = 0.0;
        }
    }

}


void MSR::MeshlessFluxes(double dTime,double dDelta) {
    double dsec;
    printf("Computing fluxes... ");
    TimerStart(TIMER_FLUXES);
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

    TimerStop(TIMER_FLUXES);
    dsec = TimerGet(TIMER_FLUXES);
    printf("took %.5f seconds\n", dsec);
}


void packHydroFluxes(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<hydroFluxesPack *>(dst);
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
        p1->P = sph.P;

        p1->fBall = p2.ball();
        p1->fDensity = p2.density();
        p1->uRung = p2.rung();
        p1->bMarked = p2.marked();
    }
}

void unpackHydroFluxes(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const hydroFluxesPack *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
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
        sph.P = p2->P;

        p1.set_ball(p2->fBall);
        p1.set_density(p2->fDensity);
        p1.set_rung(p2->uRung);
        p1.set_marked(p2->bMarked);
    }
}

void initHydroFluxes(void *vpkd,void *dst) {
}

/* Zero all the conserved quantities in cached copies, which will be updated
 * during the hydro loop. They will be merged with the actual particle
 * information in combHydroFluxes
 */
void initHydroFluxesCached(void *vpkd,void *dst) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(dst)];
    assert(!pkd->bNoParticleOrder);
    // For the init*Cached and comb functions we still have to explicitly
    // check if we are handling gas particles even ifdef OPTIM_REORDER_IN_NODES
    // because these operations are done in a cache-line basis, and a line may
    // contain other particles that are not of interest!
    if (p.is_gas()) {
        auto &sph = p.sph();

        sph.Frho = 0.0;
        sph.Fmom = 0.0;
        sph.Fene = 0.0;
#ifndef USE_MFM
        sph.drDotFrho = 0.0;
#endif

        sph.mom = 0.0;
        sph.E = 0.0;
        sph.Uint = 0.0;

        p.set_mass(0.0);
    }
}

void flushHydroFluxes(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<hydroFluxesFlush *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    if (p2.is_gas()) {
        const auto &sph = p2.sph();

        p1->Frho = sph.Frho;
        p1->Fmom = sph.Fmom;
        p1->Fene = sph.Fene;
#ifndef USE_MFM
        p1->drDotFrho = sph.drDotFrho;
#endif

        p1->mom = sph.mom;
        p1->E = sph.E;
        p1->Uint = sph.Uint;

        p1->fMass = p2.mass();
    }
}

void combHydroFluxes(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const hydroFluxesFlush *>(src);

    assert(!pkd->bNoParticleOrder);
    if (p1.is_gas()) {
        auto &sph = p1.sph();

        sph.Frho += p2->Frho;
        sph.Fmom += p2->Fmom;
        sph.Fene += p2->Fene;
#ifndef USE_MFM
        sph.drDotFrho += p2->drDotFrho;
#endif

        sph.mom += p2->mom;
        sph.E += p2->E;
        sph.Uint += p2->Uint;

        p1.set_mass(p1.mass() + p2->fMass);
    }
}


/* This version is deprecated. It may be removed in future versions of the code
 * without notice.
 *
 * The maintained Riemann solver is the vectorized version of this function,
 * which is activated wit the OPTIM_SMOOTH_NODE and OPTIM_FLUX_VEC flags
 */
void hydroRiemann(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    //TODO Clean unused variables!
    PKD pkd = smf->pkd;
    auto P = pkd->particles[pIn];
    double minDt;
    double pvFrame[3], qvFrame[3], vFrame[3];
    double modApq, pDeltaHalf, qDeltaHalf;
    double pdivv, qdivv, psi;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3], dr[3];
    struct Input_vec_Riemann riemann_input;
    struct Riemann_outputs riemann_output;

    const auto &pv = P.velocity();
    auto &psph = P.sph();
    const double ph = 0.5*P.ball();
    const double pDensity = P.density();

    for (auto i = 0; i < nSmooth; ++i) {

        /* In the nnList there is a 'copy' of the own particle,
         * which we can omit as there are no fluxes to be computed here
         */
        if (nnList[i].pPart == pIn) continue;

        auto Q = pkd->particles[nnList[i].pPart];
        const auto &qv = Q.velocity();
        auto &qsph = Q.sph();
        const double qh = 0.5*Q.ball();

        const double dx = nnList[i].dr[0];
        const double dy = nnList[i].dr[1];
        const double dz = nnList[i].dr[2];
#ifdef FORCE_1D
        if (dz!=0) continue;
        if (dy!=0) continue;
#endif
#ifdef FORCE_2D
        if (dz!=0) continue;
#endif


        const double hpq = ph;
        const double rpq = sqrt(nnList[i].fDist2);
        // We only compute the fluxes if both particles are within the kernel
        // of each other
        if (2.*qh < rpq) continue;

        const double Wpq = cubicSplineKernel(rpq, hpq);
        if (Wpq==0.0) {
            continue;
        }

        /* We update the conservatives variables taking the minimum timestep
         * between the particles, as in AREPO */
        if (!Q.is_active()) {
            // If q is not active we now that p has the smallest dt
            minDt = smf->dDelta/(1<<P.rung()) ;
        }
        else {
            // Otherwise we need to explicitly check
            if (P.rung() > Q.rung()) {
                minDt = smf->dDelta/(1<<P.rung()) ;
            }
            else {
                minDt = smf->dDelta/(1<<Q.rung()) ;
            }
        }

        if (smf->dDelta > 0) {
            pDeltaHalf = smf->dTime - psph.lastUpdateTime + 0.5*smf->dDelta/(1<<P.rung());
            qDeltaHalf = smf->dTime - qsph.lastUpdateTime + 0.5*smf->dDelta/(1<<Q.rung());
        }
        else {
            /* For the initialization step we do not extrapolate because we
             * dont have a reliable dDelta
             */
            qDeltaHalf = 0.0;
            pDeltaHalf = 0.0;
        }
        if (pkd->csm->val.bComove) {
            qDeltaHalf /= smf->a;
            pDeltaHalf /= smf->a;
        }

        // DEBUG: Avoid temporal extrapolation
        //pDeltaHalf = 0.;
        //qDeltaHalf = 0.;



        // \tilde{\psi}_j (x_i)
        psi = -cubicSplineKernel(rpq, ph)/psph.omega;
        psiTilde_p[0] = (psph.B[XX]*dx + psph.B[XY]*dy + psph.B[XZ]*dz)*psi;
        psiTilde_p[1] = (psph.B[XY]*dx + psph.B[YY]*dy + psph.B[YZ]*dz)*psi;
        psiTilde_p[2] = (psph.B[XZ]*dx + psph.B[YZ]*dy + psph.B[ZZ]*dz)*psi;

        // \tilde{\psi}_i (x_j)
        psi = cubicSplineKernel(rpq, qh)/qsph.omega;
        psiTilde_q[0] = (qsph.B[XX]*dx + qsph.B[XY]*dy + qsph.B[XZ]*dz)*psi;
        psiTilde_q[1] = (qsph.B[XY]*dx + qsph.B[YY]*dy + qsph.B[YZ]*dz)*psi;
        psiTilde_q[2] = (qsph.B[XZ]*dx + qsph.B[YZ]*dy + qsph.B[ZZ]*dz)*psi;

        modApq = 0.0;
        for (auto j = 0; j < 3; ++j) {
            Apq[j] = psiTilde_p[j]/psph.omega - psiTilde_q[j]/qsph.omega;
            modApq += Apq[j]*Apq[j];
        }
        modApq = sqrt(modApq);

        /* DEBUG
        if (modApq<=0.0) {
           printf("dx %e \t dy %e \t dz %e \n", dx, dy, dz);
           printf("rpq %e hpq %e ratio %e Wpq %e \n", rpq, hpq, rpq/hpq, Wpq);
        }
        assert(modApq>0.0); // Area should be positive!
        */


        if (modApq > 0.) {
            for (auto j = 0; j < 3; ++j) face_unit[j] = Apq[j]/modApq;
        }
        else {
            for (auto j = 0; j < 3; ++j) face_unit[j] = 0.;
        }



        // Velocity of the quadrature mid-point
        for (auto j = 0; j < 3; ++j) {
            vFrame[j] = 0.5 * (pv[j] + qv[j]);

            // We boost to the reference of the p-q 'face'
            pvFrame[j] = pv[j] - vFrame[j];
            qvFrame[j] = qv[j] - vFrame[j];
        }

        // Mid-point rule
        dr[0] = -0.5*dx;
        dr[1] = -0.5*dy;
        dr[2] = -0.5*dz;

        // DEBUG: Avoid spatial extrapolation
        //dr[0] = 0.;
        //dr[1] = 0.;
        //dr[2] = 0.;

        // Divergence of the velocity field for the forward in time prediction
        pdivv = psph.gradVx[0] + psph.gradVy[1] + psph.gradVz[2];
        qdivv = qsph.gradVx[0] + qsph.gradVy[1] + qsph.gradVz[2];

        pdivv *= pDeltaHalf;
        qdivv *= qDeltaHalf;


        riemann_input.L.rho = P.density();
        riemann_input.R.rho = Q.density();
        riemann_input.L.v[0] = pvFrame[0];
        riemann_input.R.v[0] = qvFrame[0];
        riemann_input.L.v[1] = pvFrame[1];
        riemann_input.R.v[1] = qvFrame[1];
        riemann_input.L.v[2] = pvFrame[2];
        riemann_input.R.v[2] = qvFrame[2];
        riemann_input.L.p = psph.P;
        riemann_input.R.p = qsph.P;

//      printf("1) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("1) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

        // We add the gradients terms (from extrapolation and forward prediction)
        for (auto j = 0; j < 3; ++j) {
            riemann_input.L.rho += ( dr[j] - pDeltaHalf*pvFrame[j])*psph.gradRho[j];
            riemann_input.R.rho += (-dr[j] - qDeltaHalf*qvFrame[j])*qsph.gradRho[j];

            riemann_input.L.v[0] += ( dr[j]*psph.gradVx[j]);
            riemann_input.R.v[0] += (-dr[j]*qsph.gradVx[j]);

            riemann_input.L.v[1] += ( dr[j]*psph.gradVy[j]);
            riemann_input.R.v[1] += (-dr[j]*qsph.gradVy[j]);

            riemann_input.L.v[2] += ( dr[j]*psph.gradVz[j]);
            riemann_input.R.v[2] += (-dr[j]*qsph.gradVz[j]);

            riemann_input.L.p += ( dr[j] - pDeltaHalf*pvFrame[j])*psph.gradP[j];
            riemann_input.R.p += (-dr[j] - qDeltaHalf*qvFrame[j])*qsph.gradP[j];
        }
//      printf("2) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("2) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

        // Placing this here solved the convergence problem for the comoving soundwaves.
        //   This problem may be caused because we do not use the time extrapolated cell-centered states in
        //   this limiter
        /*
        genericPairwiseLimiter(p.density(), q.density(), &riemann_input.L.rho, &riemann_input.R.rho);
        genericPairwiseLimiter(psph.P, qsph.P, &riemann_input.L.p, &riemann_input.R.p);
        genericPairwiseLimiter(pvFrame[0], qvFrame[0], &riemann_input.L.v[0], &riemann_input.R.v[0]);
        genericPairwiseLimiter(pvFrame[1], qvFrame[1], &riemann_input.L.v[1], &riemann_input.R.v[1]);
        genericPairwiseLimiter(pvFrame[2], qvFrame[2], &riemann_input.L.v[2], &riemann_input.R.v[2]);
        */


        double temp;


        // Forward extrapolation of velocity
        for (auto j = 0; j < 3; ++j) {
            temp = pvFrame[j]*pdivv + psph.gradP[j]/pDensity*pDeltaHalf;
            riemann_input.L.v[j] -= temp;
            vFrame[j] -= 0.5*temp;

            temp = qvFrame[j]*qdivv + qsph.gradP[j]/Q.density()*qDeltaHalf;
            riemann_input.R.v[j] -= temp;
            vFrame[j] -= 0.5*temp;
        }

        for (auto j = 0; j < 3; ++j) {
            temp = psph.lastAcc[j]*pDeltaHalf*smf->a;
            riemann_input.L.v[j] += temp;
            vFrame[j] += 0.5*temp;

            temp = qsph.lastAcc[j]*qDeltaHalf*smf->a;
            riemann_input.R.v[j] += temp;
            vFrame[j] += 0.5*temp;
        }

        riemann_input.L.rho -= pDensity*pdivv;
        riemann_input.R.rho -= Q.density()*qdivv;
        riemann_input.L.p -= smf->dConstGamma*psph.P*pdivv;
        riemann_input.R.p -= smf->dConstGamma*qsph.P*qdivv;

        genericPairwiseLimiter(P.density(), Q.density(), &riemann_input.L.rho, &riemann_input.R.rho);
        genericPairwiseLimiter(psph.P, qsph.P, &riemann_input.L.p, &riemann_input.R.p);
        genericPairwiseLimiter(pvFrame[0], qvFrame[0], &riemann_input.L.v[0], &riemann_input.R.v[0]);
        genericPairwiseLimiter(pvFrame[1], qvFrame[1], &riemann_input.L.v[1], &riemann_input.R.v[1]);
        genericPairwiseLimiter(pvFrame[2], qvFrame[2], &riemann_input.L.v[2], &riemann_input.R.v[2]);

        if (pkd->csm->val.bComove) {

            for (auto j = 0; j < 3; ++j) {
                temp = smf->H * pDeltaHalf * smf->a * pvFrame[j];
                riemann_input.L.v[j] -= temp;
                vFrame[j] -= 0.5*temp;

                temp = smf->H * qDeltaHalf * smf->a * qvFrame[j];
                riemann_input.R.v[j] -= temp;
                vFrame[j] -= 0.5*temp;
            }

            riemann_input.L.p -= 3. * smf->H * pDeltaHalf * smf->a *
                                 (smf->dConstGamma - 1.) * psph.P;

            riemann_input.R.p -= 3. * smf->H * qDeltaHalf * smf->a *
                                 (smf->dConstGamma - 1.) * qsph.P;

        }

        // DEBUG: Tests for the riemann solver extracted from Toro (10.1007/b79761)
        // Test 1
//       riemann_input.L.rho = 1.0; riemann_input.L.p = 1.0; riemann_input.L.v[0] = 0.0;
//       riemann_input.L.rho = 0.125; riemann_input.L.p = 0.1; riemann_input.L.v[0] = 0.0;

        if (riemann_input.L.rho < 0) {
            riemann_input.L.rho = P.density();
            /* printf("WARNING, L.rho < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.R.rho < 0) {
            riemann_input.R.rho = Q.density();
            /* printf("WARNING, R.rho < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.L.p < 0) {
            riemann_input.L.p = psph.P;
            /* printf("WARNING, L.p < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.R.p < 0) {
            riemann_input.R.p = qsph.P;
            /* printf("WARNING, R.p < 0 : using first-order scheme \n");*/
        }

#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
        const double a_inv3 = 1./(smf->a * smf->a * smf->a);
        const double pLeEOS = eEOSPressureFloor(a_inv3, riemann_input.L.rho, ph,
                                                smf->dConstGamma, smf->eEOS);
        if (pLeEOS != NOT_IN_EEOS)
            riemann_input.L.p = std::max(riemann_input.L.p, pLeEOS);

        const double pReEOS = eEOSPressureFloor(a_inv3, riemann_input.R.rho, qh,
                                                smf->dConstGamma, smf->eEOS);
        if (pReEOS != NOT_IN_EEOS)
            riemann_input.R.p = std::max(riemann_input.R.p, pReEOS);
#endif

        //Riemann_solver(pkd, riemann_input, &riemann_output, face_unit, /*double press_tot_limiter TODO For now, just p>0: */ 0.0);
        double cs_L = sqrt(GAMMA * riemann_input.L.p / riemann_input.L.rho);
        double cs_R = sqrt(GAMMA * riemann_input.R.p / riemann_input.R.rho);
        riemann_input.L.u  = riemann_input.L.p / (GAMMA_MINUS1 * riemann_input.L.rho);
        riemann_input.R.u  = riemann_input.R.p / (GAMMA_MINUS1 * riemann_input.R.rho);
        double h_L = riemann_input.L.p/riemann_input.L.rho +
                     riemann_input.L.u +
                     0.5*(riemann_input.L.v[0]*riemann_input.L.v[0]+
                          riemann_input.L.v[1]*riemann_input.L.v[1]+
                          riemann_input.L.v[2]*riemann_input.L.v[2]);

        double h_R = riemann_input.R.p/riemann_input.R.rho +
                     riemann_input.R.u +
                     0.5*(riemann_input.R.v[0]*riemann_input.R.v[0]+
                          riemann_input.R.v[1]*riemann_input.R.v[1]+
                          riemann_input.R.v[2]*riemann_input.R.v[2]);

        double v_line_L = riemann_input.L.v[0]*face_unit[0] +
                          riemann_input.L.v[1]*face_unit[1] +
                          riemann_input.L.v[2]*face_unit[2];
        double v_line_R = riemann_input.R.v[0]*face_unit[0] +
                          riemann_input.R.v[1]*face_unit[1] +
                          riemann_input.R.v[2]*face_unit[2];

        // We just remove this to avoid the compiler from screaming.
        // This function is never called in this case.
#ifndef OPTIM_FLUX_VEC
        Riemann_solver_exact(smf, riemann_input, &riemann_output, face_unit, v_line_L, v_line_R, cs_L, cs_R, h_L, h_R);
#endif

#ifdef USE_MFM
        /*
        if (riemann_output.Fluxes.rho != 0){
           printf("Frho %e \n",riemann_output.Fluxes.rho);
           abort();
        }
        */
        riemann_output.Fluxes.rho = 0.;
        riemann_output.Fluxes.p = riemann_output.P_M * riemann_output.S_M;
        for (auto j = 0; j < 3; ++j)
            riemann_output.Fluxes.v[j] = riemann_output.P_M * face_unit[j];
#endif

        // Force 2D
#ifdef FORCE_1D
        riemann_output.Fluxes.v[2] = 0.;
        riemann_output.Fluxes.v[1] = 0.;
#endif
#ifdef FORCE_2D
        riemann_output.Fluxes.v[2] = 0.;
#endif


        // Check for NAN fluxes
        if (riemann_output.Fluxes.rho!=riemann_output.Fluxes.rho)
            riemann_output.Fluxes.rho = 0.;//abort();
        if (riemann_output.Fluxes.p!=riemann_output.Fluxes.p)
            riemann_output.Fluxes.p = 0.;//abort();


        if (pkd->csm->val.bComove)
            minDt /= smf->a; // 1/a term before \nabla



        // Now we de-boost the fluxes following Eq. A8 Hopkins 2015
        for (auto j = 0; j < 3; ++j) {
            riemann_output.Fluxes.p += vFrame[j] * riemann_output.Fluxes.v[j];
            riemann_output.Fluxes.p += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
        }

        // Now we just multiply by the face area
        riemann_output.Fluxes.p *= modApq;
        riemann_output.Fluxes.rho *= modApq;
        for (auto j = 0; j < 3; ++j) {
            riemann_output.Fluxes.v[j] *= modApq;
            riemann_output.Fluxes.v[j] += vFrame[j]*riemann_output.Fluxes.rho;
        }


        if (smf->dDelta > 0) {

#ifndef OPTIM_NO_REDUNDANT_FLUXES
            {
#else
            if ( (2.*qh < rpq) | !Q.is_active()) {
#endif


                Q.set_mass(Q.mass() + minDt * riemann_output.Fluxes.rho);

                qsph.mom[0] += minDt * riemann_output.Fluxes.v[0] ;
                qsph.mom[1] += minDt * riemann_output.Fluxes.v[1] ;
                qsph.mom[2] += minDt * riemann_output.Fluxes.v[2] ;

                qsph.E += minDt * riemann_output.Fluxes.p;

                qsph.Uint += minDt * ( riemann_output.Fluxes.p -
                                       riemann_output.Fluxes.v[0]*qv[0] -
                                       riemann_output.Fluxes.v[1]*qv[1] -
                                       riemann_output.Fluxes.v[2]*qv[2] +
                                       0.5*dot(qv,qv)*riemann_output.Fluxes.rho );
#ifndef USE_MFM
                qsph.drDotFrho[0] += minDt * riemann_output.Fluxes.rho * dx * smf->a;
                qsph.drDotFrho[1] += minDt * riemann_output.Fluxes.rho * dy * smf->a;
                qsph.drDotFrho[2] += minDt * riemann_output.Fluxes.rho * dz * smf->a;
#endif

                qsph.Frho -= riemann_output.Fluxes.rho;
                qsph.Fene -= riemann_output.Fluxes.p;
                for (auto j = 0; j < 3; ++j) {
                    qsph.Fmom[j] -= riemann_output.Fluxes.v[j];
                }
            }

            P.set_mass(P.mass() - minDt * riemann_output.Fluxes.rho);

            psph.mom[0] -= minDt * riemann_output.Fluxes.v[0];
            psph.mom[1] -= minDt * riemann_output.Fluxes.v[1];
            psph.mom[2] -= minDt * riemann_output.Fluxes.v[2];

            psph.E -= minDt * riemann_output.Fluxes.p;

            psph.Uint -= minDt * ( riemann_output.Fluxes.p -
                                   riemann_output.Fluxes.v[0]*pv[0] -
                                   riemann_output.Fluxes.v[1]*pv[1] -
                                   riemann_output.Fluxes.v[2]*pv[2] +
                                   0.5*dot(pv,pv)*riemann_output.Fluxes.rho );
#ifndef USE_MFM
            psph.drDotFrho[0] += minDt * riemann_output.Fluxes.rho * dx * smf->a;
            psph.drDotFrho[1] += minDt * riemann_output.Fluxes.rho * dy * smf->a;
            psph.drDotFrho[2] += minDt * riemann_output.Fluxes.rho * dz * smf->a;
#endif
        }
        /* Old fluxes update (see 15/04/19 )
         * TODO: This is not needed for the update of the conserved variables. Instead,
         * it is now only used for the acceleleration criteria. Memory-wise, this can be
         * substantially improved
         */
        // Own contribution is always added
        psph.Frho += riemann_output.Fluxes.rho;
        psph.Fene += riemann_output.Fluxes.p;
        for (auto j = 0; j < 3; ++j) {
            psph.Fmom[j] += riemann_output.Fluxes.v[j];
        }

    } // End of loop over neighbors

}



#ifdef OPTIM_FLUX_VEC
/* Vectorizable version of the riemann solver.
 *
 * There are a few differences with respect the previous one, the most importants:
 *   a) we use the input buffer directly,
 *      rather than accessing the particle data directly
 *   b) we omit all clauses that could terminate the loop (continue, abort, etc).
 *      If, for example, FORCE_2D is used, the loop may not be vectorized
 *
 * Now we have two hydroRiemann routines, which means that there is A LOT of
 * code duplication. This can cause bugs and deteriorate readability.
 *
 * At some point, one of them must be discontinued TODO
 */

/* When doing SIMD, a structure of arrays is used,
 * rather than an array of structures.
 *
 * To simplify the indexing of the elements, these enum should be
 * always used
 */
enum q_data {
    q_mass,
    q_ball,
    q_dx, q_dy, q_dz, q_dr,
    q_rung,
    q_rho,
    q_P,
#ifdef ENTROPY_SWITCH
    q_S,
#endif
    q_vx, q_vy, q_vz,
    q_gradRhoX, q_gradRhoY, q_gradRhoZ,
    q_gradPX, q_gradPY, q_gradPZ,
    q_gradVxX, q_gradVxY, q_gradVxZ,
    q_gradVyX, q_gradVyY, q_gradVyZ,
    q_gradVzX, q_gradVzY, q_gradVzZ,
    q_lastAccX, q_lastAccY, q_lastAccZ,
    q_lastUpdateTime,
    q_B_XX, q_B_YY, q_B_ZZ, q_B_XY, q_B_XZ, q_B_YZ,
    q_omega,
    q_last
};

enum FLUX_OUT {
    out_minDt,
    out_Frho,
    out_FmomX,out_FmomY,out_FmomZ,
    out_Fene,
#ifdef ENTROPY_SWITCH
    out_FS,
#endif
    out_last
};

// Simple macro to improve readability
#define q(X) input_buffer[q_##X][i]
void hydroRiemann_vec(PARTICLE *pIn,float fBall,int nSmooth,
                      my_real **restrict input_buffer,
                      my_real **restrict output_buffer, SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[pIn];

    const auto &pv = P.velocity();
    auto &psph = P.sph();
    const my_real ph = 0.5*P.ball();

    const my_real pDensity = P.density();
    const my_real p_omega = psph.omega;


#ifdef __INTEL_COMPILER
    __assume_aligned(input_buffer, 64);
    __assume_aligned(input_buffer[0], 64);
#pragma simd
#pragma vector aligned
#endif
#ifdef __GNUC__
//TODO Trick GCC into autovectorizing this!!
#endif
    for (auto i = 0; i < nSmooth; ++i) {

        const my_real qh = q(ball);

        const my_real dx = q(dx);
        const my_real dy = q(dy);
        const my_real dz = q(dz);

#ifdef FORCE_1D
        if (dz!=0) continue;
        if (dy!=0) continue;
#endif
#ifdef FORCE_2D
        if (dz!=0) continue;
#endif



        // Face where the riemann problem will be solved
        const my_real rpq = q(dr);


        /* We update the conservatives variables taking the minimum timestep
         * between the particles, as in AREPO
         */
        const my_real p_dt = smf->dDelta/(1<<P.rung());
        const my_real q_dt = q(rung);
        const my_real minDt = (p_dt > q_dt ? q_dt : p_dt) / smf->a;


        my_real qDeltaHalf=0.0, pDeltaHalf=0.0;
        if (smf->dDelta > 0) {
            pDeltaHalf = (smf->dTime - psph.lastUpdateTime + 0.5*p_dt)/smf->a;
            qDeltaHalf = (smf->dTime - q(lastUpdateTime) + 0.5*q_dt)/smf->a;
        }

        // DEBUG: Avoid temporal extrapolation
        //pDeltaHalf = 0.;
        //qDeltaHalf = 0.;

        const my_real omega_q = q(omega);

        // \tilde{\psi}_j (x_i)
        my_real psi = -cubicSplineKernel(rpq, ph)/p_omega;
        TinyVector<my_real,3> psiTilde_p, psiTilde_q;
        psiTilde_p[0] = (psph.B[XX]*dx + psph.B[XY]*dy + psph.B[XZ]*dz)*psi;
        psiTilde_p[1] = (psph.B[XY]*dx + psph.B[YY]*dy + psph.B[YZ]*dz)*psi;
        psiTilde_p[2] = (psph.B[XZ]*dx + psph.B[YZ]*dy + psph.B[ZZ]*dz)*psi;

        // \tilde{\psi}_i (x_j)
        psi = cubicSplineKernel(rpq, qh)/omega_q;
        psiTilde_q[0] = (q(B_XX)*dx + q(B_XY)*dy + q(B_XZ)*dz)*psi;
        psiTilde_q[1] = (q(B_XY)*dx + q(B_YY)*dy + q(B_YZ)*dz)*psi;
        psiTilde_q[2] = (q(B_XZ)*dx + q(B_YZ)*dy + q(B_ZZ)*dz)*psi;

        const TinyVector<my_real,3> Apq{psiTilde_p/p_omega - psiTilde_q/omega_q};
        const my_real modApq = sqrt(dot(Apq,Apq));

        /* DEBUG
        if (modApq<=0.0) {
           printf("dx %e \t dy %e \t dz %e \n", dx, dy, dz);
           printf("rpq %e hpq %e ratio %e Wpq %e \n", rpq, hpq, rpq/hpq, Wpq);
        }
        assert(modApq>0.0); // Area should be positive!
        */


        TinyVector<my_real,3> face_unit{0.0};
        if (modApq > 0.) {
            face_unit = Apq / modApq;
        }


        // Velocity of the quadrature mid-point
        my_real vFrame[3];
        vFrame[0] = 0.5 * (pv[0] + q(vx));
        vFrame[1] = 0.5 * (pv[1] + q(vy));
        vFrame[2] = 0.5 * (pv[2] + q(vz));

        my_real pvFrame[3], qvFrame[3];
        for (auto j = 0; j < 3; ++j) {
            // We boost to the reference of the p-q 'face'
            pvFrame[j] = pv[j] - vFrame[j];
        }

        qvFrame[0] = q(vx) - vFrame[0];
        qvFrame[1] = q(vy) - vFrame[1];
        qvFrame[2] = q(vz) - vFrame[2];

        // Mid-point rule
        const my_real dr[3] = {-0.5*dx, -0.5*dy, -0.5*dz};

        // DEBUG: Avoid spatial extrapolation
        //dr[0] = 0.0;
        //dr[1] = 0.0;
        //dr[2] = 0.0;

        // Divergence of the velocity field for the forward in time prediction
        my_real pdivv = (psph.gradVx[0] + psph.gradVy[1] + psph.gradVz[2])*pDeltaHalf;
        my_real qdivv = (q(gradVxX) + q(gradVyY) + q(gradVzZ))*qDeltaHalf;

        // At some point we should erase the need for this structs... FIXME
        struct Input_vec_Riemann riemann_input;
        struct Riemann_outputs riemann_output;

        riemann_input.L.rho = pDensity;
        riemann_input.R.rho = q(rho);
        riemann_input.L.v[0] = pvFrame[0];
        riemann_input.R.v[0] = qvFrame[0];
        riemann_input.L.v[1] = pvFrame[1];
        riemann_input.R.v[1] = qvFrame[1];
        riemann_input.L.v[2] = pvFrame[2];
        riemann_input.R.v[2] = qvFrame[2];
        riemann_input.L.p = psph.P;
        riemann_input.R.p = q(P);

//      printf("1) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("1) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

        // We add the gradients terms (from extrapolation and forward prediction)
        for (auto j = 0; j < 3; ++j) {
            riemann_input.L.rho += ( dr[j] - pDeltaHalf*pvFrame[j])*psph.gradRho[j];

            riemann_input.L.v[0] += ( dr[j]*psph.gradVx[j]);
            riemann_input.L.v[1] += ( dr[j]*psph.gradVy[j]);
            riemann_input.L.v[2] += ( dr[j]*psph.gradVz[j]);

            riemann_input.L.p += ( dr[j] - pDeltaHalf*pvFrame[j])*psph.gradP[j];
        }



        riemann_input.R.rho += (-dr[0] - qDeltaHalf*qvFrame[0])*q(gradRhoX);
        riemann_input.R.rho += (-dr[1] - qDeltaHalf*qvFrame[1])*q(gradRhoY);
        riemann_input.R.rho += (-dr[2] - qDeltaHalf*qvFrame[2])*q(gradRhoZ);

        riemann_input.R.v[0] += (-dr[0]*q(gradVxX));
        riemann_input.R.v[0] += (-dr[1]*q(gradVxY));
        riemann_input.R.v[0] += (-dr[2]*q(gradVxZ));

        riemann_input.R.v[1] += (-dr[0]*q(gradVyX));
        riemann_input.R.v[1] += (-dr[1]*q(gradVyY));
        riemann_input.R.v[1] += (-dr[2]*q(gradVyZ));

        riemann_input.R.v[2] += (-dr[0]*q(gradVzX));
        riemann_input.R.v[2] += (-dr[1]*q(gradVzY));
        riemann_input.R.v[2] += (-dr[2]*q(gradVzZ));

        riemann_input.R.p += (-dr[0] - qDeltaHalf*qvFrame[0])*q(gradPX);
        riemann_input.R.p += (-dr[1] - qDeltaHalf*qvFrame[1])*q(gradPY);
        riemann_input.R.p += (-dr[2] - qDeltaHalf*qvFrame[2])*q(gradPZ);



//      printf("2) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("2) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);



        my_real temp;



        for (auto j = 0; j < 3; ++j) { // Forward extrapolation of velocity
            temp = pvFrame[j]*pdivv + psph.gradP[j]/pDensity*pDeltaHalf;
            riemann_input.L.v[j] -= temp;
            vFrame[j] -= 0.5*temp;
        }

        temp = qvFrame[0]*qdivv + q(gradPX)/q(rho)*qDeltaHalf;
        riemann_input.R.v[0] -= temp;
        vFrame[0] -= 0.5*temp;

        temp = qvFrame[1]*qdivv + q(gradPY)/q(rho)*qDeltaHalf;
        riemann_input.R.v[1] -= temp;
        vFrame[1] -= 0.5*temp;

        temp = qvFrame[2]*qdivv + q(gradPZ)/q(rho)*qDeltaHalf;
        riemann_input.R.v[2] -= temp;
        vFrame[2] -= 0.5*temp;


        for (auto j = 0; j < 3; ++j) {
            temp = psph.lastAcc[j]*pDeltaHalf*smf->a;
            riemann_input.L.v[j] += temp;
            vFrame[j] += 0.5*temp;

        }
        temp = q(lastAccX)*qDeltaHalf*smf->a;
        riemann_input.R.v[0] += temp;
        vFrame[0] += 0.5*temp;

        temp = q(lastAccY)*qDeltaHalf*smf->a;
        riemann_input.R.v[1] += temp;
        vFrame[1] += 0.5*temp;

        temp = q(lastAccZ)*qDeltaHalf*smf->a;
        riemann_input.R.v[2] += temp;
        vFrame[2] += 0.5*temp;

        riemann_input.L.rho -= pDensity*pdivv;
        riemann_input.R.rho -= q(rho)*qdivv;
        riemann_input.L.p -= smf->dConstGamma*psph.P*pdivv;
        riemann_input.R.p -= smf->dConstGamma*q(P)*qdivv;

        genericPairwiseLimiter(pDensity, q(rho), &riemann_input.L.rho, &riemann_input.R.rho);
        genericPairwiseLimiter(psph.P, q(P), &riemann_input.L.p, &riemann_input.R.p);
        for (auto j = 0; j < 3; ++j) {
            genericPairwiseLimiter(pvFrame[j], qvFrame[j], &riemann_input.L.v[j], &riemann_input.R.v[j]);
        }

        if (pkd->csm->val.bComove) {

            for (auto j = 0; j < 3; ++j) {
                temp = smf->H * pDeltaHalf * smf->a * pvFrame[j];
                riemann_input.L.v[j] -= temp;
                vFrame[j] -= 0.5*temp;

                temp = smf->H * qDeltaHalf * smf->a * qvFrame[j];
                riemann_input.R.v[j] -= temp;
                vFrame[j] -= 0.5*temp;
            }

            riemann_input.L.p -= 3. * smf->H * pDeltaHalf * smf->a * (smf->dConstGamma - 1.) * psph.P;
            riemann_input.R.p -= 3. * smf->H * qDeltaHalf * smf->a * (smf->dConstGamma - 1.) * q(P);

        }

        // DEBUG: Tests for the riemann solver extracted from Toro (10.1007/b79761)
        // Test 1
//       riemann_input.L.rho = 1.0; riemann_input.L.p = 1.0; riemann_input.L.v[0] = 0.0;
//       riemann_input.L.rho = 0.125; riemann_input.L.p = 0.1; riemann_input.L.v[0] = 0.0;

        if (riemann_input.L.rho < 0) {
            riemann_input.L.rho = pDensity;
            /* printf("WARNING, L.rho < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.R.rho < 0) {
            riemann_input.R.rho = q(rho);
            /* printf("WARNING, R.rho < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.L.p < 0) {
            riemann_input.L.p = psph.P;
            /* printf("WARNING, L.p < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.R.p < 0) {
            riemann_input.R.p = q(P);
            /* printf("WARNING, R.p < 0 : using first-order scheme \n");*/
        }

#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
        const double a_inv3 = 1./(smf->a * smf->a * smf->a);
        const double pLeEOS = eEOSPressureFloor(a_inv3, riemann_input.L.rho, ph, //probably ball, not ph
                                                smf->dConstGamma, smf->eEOS);
        if (pLeEOS != NOT_IN_EEOS)
            riemann_input.L.p = std::max(riemann_input.L.p, pLeEOS);

        const double pReEOS = eEOSPressureFloor(a_inv3, riemann_input.R.rho, qh,
                                                smf->dConstGamma, smf->eEOS);
        if (pReEOS != NOT_IN_EEOS)
            riemann_input.R.p = std::max(riemann_input.R.p, pReEOS);
#endif

        double cs_L = sqrt(GAMMA * riemann_input.L.p / riemann_input.L.rho);
        double cs_R = sqrt(GAMMA * riemann_input.R.p / riemann_input.R.rho);
        riemann_input.L.u = riemann_input.L.p / (GAMMA_MINUS1 * riemann_input.L.rho);
        riemann_input.R.u = riemann_input.R.p / (GAMMA_MINUS1 * riemann_input.R.rho);
        double h_L = riemann_input.L.p/riemann_input.L.rho +
                     riemann_input.L.u +
                     0.5*(riemann_input.L.v[0]*riemann_input.L.v[0] +
                          riemann_input.L.v[1]*riemann_input.L.v[1] +
                          riemann_input.L.v[2]*riemann_input.L.v[2]);
        double h_R = riemann_input.R.p/riemann_input.R.rho +
                     riemann_input.R.u +
                     0.5*(riemann_input.R.v[0]*riemann_input.R.v[0] +
                          riemann_input.R.v[1]*riemann_input.R.v[1] +
                          riemann_input.R.v[2]*riemann_input.R.v[2]);

        double v_line_L = riemann_input.L.v[0]*face_unit[0] +
                          riemann_input.L.v[1]*face_unit[1] +
                          riemann_input.L.v[2]*face_unit[2];

        double v_line_R = riemann_input.R.v[0]*face_unit[0] +
                          riemann_input.R.v[1]*face_unit[1] +
                          riemann_input.R.v[2]*face_unit[2];

        int niter = Riemann_solver_exact(smf,
                                         riemann_input.R.rho, riemann_input.R.p, riemann_input.R.v,
                                         riemann_input.L.rho, riemann_input.L.p, riemann_input.L.v,
                                         &riemann_output.P_M, &riemann_output.S_M,
                                         &riemann_output.Fluxes.rho, &riemann_output.Fluxes.p, &riemann_output.Fluxes.v[0],
                                         face_unit.data(), v_line_L, v_line_R, cs_L, cs_R, h_L, h_R);





#ifdef ENTROPY_SWITCH

#ifdef USE_MFM
        // As we are in a truly lagrangian configuration,
        // there is no advection of entropy among particles.
        double fluxes_S = 0.;
#else
        // riemann_output.Fluxes contains now the face state given by the riemann solver.
        // We only need that for computing the entropy flux, and then can be overwritten
        double fluxes_S = 0.;
        for (auto j = 0; j < 3; ++j) fluxes_S += riemann_output.Fluxes.v[j]*face_unit[j];
        if (fluxes_S > 0) {
            // Maybe this values should be properly extrapolated to the faces..
            // but this is expensive!
            fluxes_S *= psph.S*pDensity/P.mass();
        }
        else {
            fluxes_S *= q(S)*q(rho)/q(mass);
        }
        fluxes_S *= riemann_output.Fluxes.rho*modApq;
        fluxes_S = 0.;

#endif //USE_MFM
#endif //ENTROPY_SWITCH

#ifdef USE_MFM
        riemann_output.Fluxes.rho = 0.;
        riemann_output.Fluxes.p = riemann_output.P_M * riemann_output.S_M;
        for (auto j = 0; j < 3; ++j)
            riemann_output.Fluxes.v[j] = riemann_output.P_M * face_unit[j];
#endif
        // End MFM

        // Force 2D
#ifdef FORCE_1D
        riemann_output.Fluxes.v[2] = 0.;
        riemann_output.Fluxes.v[1] = 0.;
#endif
#ifdef FORCE_2D
        riemann_output.Fluxes.v[2] = 0.;
#endif



        // DEBUG
//       abort();

        // Check for NAN fluxes
        if (riemann_output.Fluxes.rho!=riemann_output.Fluxes.rho)
            riemann_output.Fluxes.rho = 0.;//abort();
        if (riemann_output.Fluxes.p!=riemann_output.Fluxes.p)
            riemann_output.Fluxes.p = 0.;//abort();




        // Now we de-boost the fluxes following Eq. A8 Hopkins 2015
        for (auto j = 0; j < 3; ++j) {
            riemann_output.Fluxes.p += vFrame[j] * riemann_output.Fluxes.v[j];
            riemann_output.Fluxes.p += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
        }

        // Now we just multiply by the face area
        riemann_output.Fluxes.p *= modApq;
        riemann_output.Fluxes.rho *= modApq;
        for (auto j = 0; j < 3; ++j) {
            riemann_output.Fluxes.v[j] *= modApq;
            riemann_output.Fluxes.v[j] += vFrame[j]*riemann_output.Fluxes.rho;
        }

        // We fill the output buffer with the fluxes, which then
        // will be added to the corresponding particles
        output_buffer[out_Frho][i] = riemann_output.Fluxes.rho;
        output_buffer[out_Fene][i] = riemann_output.Fluxes.p;
        output_buffer[out_FmomX][i] = riemann_output.Fluxes.v[0];
        output_buffer[out_FmomY][i] = riemann_output.Fluxes.v[1];
        output_buffer[out_FmomZ][i] = riemann_output.Fluxes.v[2];
#ifdef ENTROPY_SWITCH
        output_buffer[out_FS][i] = fluxes_S;
#endif
        output_buffer[out_minDt][i] = minDt;

    } // End of loop over neighbors
}


void hydroFluxFillBuffer(my_real **buffer, PARTICLE *qIn, int i, double dr2,
                         TinyVector<double,3> dr, SMF *smf) {
    PKD pkd = smf->pkd;
    auto Q = pkd->particles[qIn];
    double dDelta = smf->dDelta;
    float qh = 0.5*Q.ball();
    auto &qsph = Q.sph();
    buffer[q_mass][i] = Q.mass();
    buffer[q_ball][i] = qh;
    buffer[q_dx][i] = dr[0];
    buffer[q_dy][i] = dr[1];
    buffer[q_dz][i] = dr[2];
    buffer[q_dr][i] = sqrt(dr2);
    buffer[q_rung][i] = dDelta/(1<<Q.rung());
    buffer[q_rho][i] = Q.density();
    buffer[q_P][i] = qsph.P;
#ifdef ENTROPY_SWITCH
    buffer[q_S][i] = qsph.S;
#endif
    const auto &qv = Q.velocity();
    buffer[q_vx][i] = qv[0];
    buffer[q_vy][i] = qv[1];
    buffer[q_vz][i] = qv[2];

    buffer[q_gradRhoX][i] = qsph.gradRho[0];
    buffer[q_gradRhoY][i] = qsph.gradRho[1];
    buffer[q_gradRhoZ][i] = qsph.gradRho[2];

    buffer[q_gradPX][i] = qsph.gradP[0];
    buffer[q_gradPY][i] = qsph.gradP[1];
    buffer[q_gradPZ][i] = qsph.gradP[2];

    buffer[q_gradVxX][i] = qsph.gradVx[0];
    buffer[q_gradVxY][i] = qsph.gradVx[1];
    buffer[q_gradVxZ][i] = qsph.gradVx[2];

    buffer[q_gradVyX][i] = qsph.gradVy[0];
    buffer[q_gradVyY][i] = qsph.gradVy[1];
    buffer[q_gradVyZ][i] = qsph.gradVy[2];

    buffer[q_gradVzX][i] = qsph.gradVz[0];
    buffer[q_gradVzY][i] = qsph.gradVz[1];
    buffer[q_gradVzZ][i] = qsph.gradVz[2];

    buffer[q_lastUpdateTime][i] = qsph.lastUpdateTime;
    buffer[q_lastAccX][i] = qsph.lastAcc[0];
    buffer[q_lastAccY][i] = qsph.lastAcc[1];
    buffer[q_lastAccZ][i] = qsph.lastAcc[2];
    buffer[q_B_XX][i] = qsph.B[XX];
    buffer[q_B_YY][i] = qsph.B[YY];
    buffer[q_B_ZZ][i] = qsph.B[ZZ];
    buffer[q_B_XY][i] = qsph.B[XY];
    buffer[q_B_XZ][i] = qsph.B[XZ];
    buffer[q_B_YZ][i] = qsph.B[YZ];
    buffer[q_omega][i] = qsph.omega;
}


void hydroFluxUpdateFromBuffer(my_real **out_buffer, my_real **in_buffer,
                               PARTICLE *pIn, PARTICLE *qIn, int i, SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[pIn];
    auto Q = pkd->particles[qIn];
    const auto &pv = P.velocity();
    const auto &qv = Q.velocity();
    auto &psph = P.sph();
    auto &qsph = Q.sph();
    const auto &dDelta = smf->dDelta;
    const auto &aFac = smf->a;
    if (dDelta>0) {
        P.set_mass(P.mass() - out_buffer[out_minDt][i] * out_buffer[out_Frho][i]);

        psph.mom[0] -= out_buffer[out_minDt][i] * out_buffer[out_FmomX][i];
        psph.mom[1] -= out_buffer[out_minDt][i] * out_buffer[out_FmomY][i];
        psph.mom[2] -= out_buffer[out_minDt][i] * out_buffer[out_FmomZ][i];

        psph.E -= out_buffer[out_minDt][i] * out_buffer[out_Fene][i];

        psph.Uint -= out_buffer[out_minDt][i] * ( out_buffer[out_Fene][i]
                     - out_buffer[out_FmomX][i]*pv[0]
                     - out_buffer[out_FmomY][i]*pv[1]
                     - out_buffer[out_FmomZ][i]*pv[2]
                     + 0.5*dot(pv,pv)*out_buffer[out_Frho][i] );

#ifdef ENTROPY_SWITCH
        psph.S -= out_buffer[out_minDt][i] * out_buffer[out_FS][i];
#endif

#ifndef USE_MFM
        psph.drDotFrho[0] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dx][i] * aFac;
        psph.drDotFrho[1] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dy][i] * aFac;
        psph.drDotFrho[2] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dz][i] * aFac;
#endif
        psph.Frho +=    out_buffer[out_Frho][i];
        psph.Fene +=    out_buffer[out_Fene][i];
        psph.Fmom[0] += out_buffer[out_FmomX][i];
        psph.Fmom[1] += out_buffer[out_FmomY][i];
        psph.Fmom[2] += out_buffer[out_FmomZ][i];
    }
    else {
        psph.Frho +=    out_buffer[out_Frho][i];
        psph.Fene +=    out_buffer[out_Fene][i];
        psph.Fmom[0] += out_buffer[out_FmomX][i];
        psph.Fmom[1] += out_buffer[out_FmomY][i];
        psph.Fmom[2] += out_buffer[out_FmomZ][i];
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
        assert( qsph.P == in_buffer[q_P][i] );
        if (dDelta>0) {
            Q.set_mass(Q.mass() + out_buffer[out_minDt][i] * out_buffer[out_Frho][i]);

            qsph.mom[0] += out_buffer[out_minDt][i] * out_buffer[out_FmomX][i];
            qsph.mom[1] += out_buffer[out_minDt][i] * out_buffer[out_FmomY][i];
            qsph.mom[2] += out_buffer[out_minDt][i] * out_buffer[out_FmomZ][i];

            qsph.E += out_buffer[out_minDt][i] * out_buffer[out_Fene][i];

            qsph.Uint += out_buffer[out_minDt][i] * ( out_buffer[out_Fene][i]
                         - out_buffer[out_FmomX][i]*qv[0]
                         - out_buffer[out_FmomY][i]*qv[1]
                         - out_buffer[out_FmomZ][i]*qv[2]
                         + 0.5*dot(qv,qv)*out_buffer[out_Frho][i] );
#ifdef ENTROPY_SWITCH
            qsph.S += out_buffer[out_minDt][i] * out_buffer[out_FS][i];
#endif

#ifndef USE_MFM
            qsph.drDotFrho[0] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dx][i] * aFac;
            qsph.drDotFrho[1] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dy][i] * aFac;
            qsph.drDotFrho[2] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dz][i] * aFac;
#endif
            qsph.Frho -=    out_buffer[out_Frho][i];
            qsph.Fene -=    out_buffer[out_Fene][i];
            qsph.Fmom[0] -= out_buffer[out_FmomX][i];
            qsph.Fmom[1] -= out_buffer[out_FmomY][i];
            qsph.Fmom[2] -= out_buffer[out_FmomZ][i];
        }
        else {
            qsph.Frho -=    out_buffer[out_Frho][i];
            qsph.Fene -=    out_buffer[out_Fene][i];
            qsph.Fmom[0] -= out_buffer[out_FmomX][i];
            qsph.Fmom[1] -= out_buffer[out_FmomY][i];
            qsph.Fmom[2] -= out_buffer[out_FmomZ][i];
        }

    } // q marked/active

}

void hydroFluxGetNvars(int *in, int *out) {
    *in = q_last;
    *out = out_last;
}

#endif // OPTIM_FLUX_VEC

