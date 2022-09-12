#include "hydro/hydro.h"
#include "master.h"
#include "eEOS/eEOS.h"
#include "hydro/limiters.h"
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
void hydroRiemann_old(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    //TODO Clean unused variables!
    PKD pkd = smf->pkd;
    auto P = pkd->particles[pIn];
    double minDt;
    double pvFrame[3], qvFrame[3], vFrame[3];
    double modApq, pDeltaHalf, qDeltaHalf;
    double pdivv, qdivv, psi;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3], dr[3];
    struct Input_vec_Riemann<double> riemann_input;
    struct Riemann_outputs<double> riemann_output;

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

        /*
        genericPairwiseLimiter(P.density(), Q.density(), &riemann_input.L.rho, &riemann_input.R.rho);
        genericPairwiseLimiter(psph.P, qsph.P, &riemann_input.L.p, &riemann_input.R.p);
        genericPairwiseLimiter(pvFrame[0], qvFrame[0], &riemann_input.L.v[0], &riemann_input.R.v[0]);
        genericPairwiseLimiter(pvFrame[1], qvFrame[1], &riemann_input.L.v[1], &riemann_input.R.v[1]);
        genericPairwiseLimiter(pvFrame[2], qvFrame[2], &riemann_input.L.v[2], &riemann_input.R.v[2]);
        */

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
        /*
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
        */

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


template <typename ftype=double>
static inline void extrapolateDensityInTime(ftype &rho, ftype rho0, ftype vx, ftype vy, ftype vz,
        ftype dt, ftype gradRhoX, ftype gradRhoY, ftype gradRhoZ, ftype divv) {
    rho -= dt*(vx*gradRhoX + vy*gradRhoY + vz*gradRhoZ + rho0*divv);
}

template <typename ftype=double>
static inline void extrapolateVelocityInTime(ftype &v, ftype &vFrame, ftype divv, ftype v0, ftype dt, ftype acc,
        ftype gradP, ftype rho0,
        ftype a) {
    ftype temp;
    temp = -(v0*divv + gradP/rho0)*dt;
    temp += acc*dt*a;
    v += temp;
    vFrame += 0.5*temp;
}

template <typename ftype=double>
static inline void extrapolatePressureInTime(ftype &p, ftype dt, ftype gradPx, ftype gradPy, ftype gradPz,
        ftype p0, ftype vx0, ftype vy0, ftype vz0, ftype divv, ftype dConstGamma) {
    p -= dt*( vx0*gradPx + vy0*gradPy + vz0*gradPz  + dConstGamma*p0*divv);
}

template <typename ftype=double>
static inline void extrapolateStateInTime(
    ftype &rho, ftype &vx, ftype &vy, ftype &vz, ftype &p,
    ftype &vFx, ftype &vFy, ftype &vFz,
    ftype rho0, ftype vx0, ftype vy0, ftype vz0, ftype p0,
    ftype dt,
    ftype gradRhoX, ftype gradRhoY, ftype gradRhoZ,
    ftype gradPX, ftype gradPY, ftype gradPZ,
    ftype accx, ftype accy, ftype accz, ftype divv,
    ftype dConstGamma, ftype a) {

    extrapolateDensityInTime( rho, rho0, vx, vy, vz, dt,
                              gradRhoX, gradRhoY, gradRhoZ, divv);

    extrapolateVelocityInTime( vx, vFx, divv, vx0, dt, accx, gradPX, rho0, a);
    extrapolateVelocityInTime( vy, vFy, divv, vy0, dt, accy, gradPY, rho0, a);
    extrapolateVelocityInTime( vz, vFz, divv, vz0, dt, accz, gradPZ, rho0, a);

    extrapolatePressureInTime( p, dt, gradPX, gradPY, gradPZ,
                               p0, vx0, vy0, vz0, divv, dConstGamma);
}

template <typename ftype=double>
static inline void extrapolateVariableInSpace(ftype &var, ftype dx, ftype dy, ftype dz,
        ftype gradx, ftype grady, ftype gradz) {
    var += dx*gradx + dy*grady + dz*gradz;
}


template <typename ftype=double>
static inline void extrapolateVelocityCosmology(ftype &v, ftype &vFrame, ftype v0, ftype dt, ftype H, ftype a) {
    ftype temp = H * dt * a * v;
    v -= temp;
    vFrame -= 0.5*temp;
}


template <typename ftype=double>
static inline void extrapolateCosmology(ftype &vx, ftype &vy, ftype &vz,
                                        ftype &vFramex, ftype &vFramey, ftype &vFramez, ftype &p,
                                        ftype dt, ftype vx0, ftype vy0, ftype vz0, ftype p0,
                                        ftype dConstGamma, ftype H, ftype a) {
    extrapolateVelocityCosmology( vx,  vFramex,  vx0,  dt,  H,  a);
    extrapolateVelocityCosmology( vy,  vFramey,  vy0,  dt,  H,  a);
    extrapolateVelocityCosmology( vz,  vFramez,  vz0,  dt,  H,  a);

    p -= 3. * H * dt * a * (dConstGamma - 1.) * p0;
}

template <typename ftype=double>
static inline void computeFace(ftype &modApq, std::array<ftype,3> &unit,
                               ftype rpq,  ftype dx, ftype dy, ftype dz,
                               ftype ph, ftype qh, ftype p_omega, ftype q_omega,
                               ftype pBxx, ftype pBxy, ftype pBxz,
                               ftype pByy, ftype pByz, ftype pBzz,
                               ftype qBxx, ftype qBxy, ftype qBxz,
                               ftype qByy, ftype qByz, ftype qBzz) {
    ftype psi;
    ftype psiTilde[3];
    ftype Apq[3] = {0.,0.,0.};

    // \tilde{\psi}_j (x_i)
    psi  = -cubicSplineKernel(rpq, ph)/p_omega;
    psiTilde[0] = (pBxx*dx + pBxy*dy + pBxz*dz)*psi;
    psiTilde[1] = (pBxy*dx + pByy*dy + pByz*dz)*psi;
    psiTilde[2] = (pBxz*dx + pByz*dy + pBzz*dz)*psi;
    for (auto j=0; j<3; j++) {
        Apq[j] += psiTilde[j]/p_omega;
    }

    // \tilde{\psi}_i (x_j)
    psi = cubicSplineKernel(rpq, qh)/q_omega;
    psiTilde[0] = (qBxx*dx + qBxy*dy + qBxz*dz)*psi;
    psiTilde[1] = (qBxy*dx + qByy*dy + qByz*dz)*psi;
    psiTilde[2] = (qBxz*dx + qByz*dy + qBzz*dz)*psi;
    for (auto j=0; j<3; j++) {
        Apq[j] -= psiTilde[j]/q_omega;
    }

    modApq = 0.0;
    for (auto j=0; j<3; j++) {
        modApq += Apq[j]*Apq[j];
    }
    modApq = sqrt(modApq);
    for (auto j=0; j<3; j++) {
        unit[j] = Apq[j]/modApq;
    }
}


template <typename ftype=double>
static inline void doSinglePPFlux(ftype &F_rho, std::array<ftype,3> &F_v, ftype &F_p, ftype &F_S, ftype &minDt,
                                  bool bComove, ftype dTime, ftype dDelta, ftype a, ftype H, ftype dConstGamma,
                                  ftype rpq, ftype dx, ftype dy, ftype dz,
                                  ftype pBall, ftype pLastUpdateTime, ftype pDt,
                                  ftype pOmega,
                                  ftype pBxx, ftype pBxy, ftype pBxz,
                                  ftype pByy, ftype pByz, ftype pBzz,
                                  ftype pDensity, ftype pVpredx, ftype pVpredy, ftype pVpredz, ftype pP, ftype pS,
                                  ftype pGradRhoX, ftype pGradRhoY, ftype pGradRhoZ,
                                  ftype pGradPX, ftype pGradPY, ftype pGradPZ,
                                  ftype pGradVxX, ftype pGradVxY, ftype pGradVxZ,
                                  ftype pGradVyX, ftype pGradVyY, ftype pGradVyZ,
                                  ftype pGradVzX, ftype pGradVzY, ftype pGradVzZ,
                                  ftype pLastAccX, ftype pLastAccY, ftype pLastAccZ,
                                  ftype qBall, ftype qLastUpdateTime, ftype qDt,
                                  ftype qOmega,
                                  ftype qBxx, ftype qBxy, ftype qBxz,
                                  ftype qByy, ftype qByz, ftype qBzz,
                                  ftype qDensity, ftype qVpredx, ftype qVpredy, ftype qVpredz, ftype qP, ftype qS,
                                  ftype qGradRhoX, ftype qGradRhoY, ftype qGradRhoZ,
                                  ftype qGradPX, ftype qGradPY, ftype qGradPZ,
                                  ftype qGradVxX, ftype qGradVxY, ftype qGradVxZ,
                                  ftype qGradVyX, ftype qGradVyY, ftype qGradVyZ,
                                  ftype qGradVzX, ftype qGradVzY, ftype qGradVzZ,
                                  ftype qLastAccX, ftype qLastAccY, ftype qLastAccZ) {


#ifdef FORCE_1D
    if (dz!=0) continue;
    if (dy!=0) continue;
#endif
#ifdef FORCE_2D
    if (dz!=0) continue;
#endif



    /* We update the conservatives variables taking the minimum timestep
     * between the particles, as in AREPO
     */
    minDt = min(pDt, qDt);
    minDt /=  a;


    ftype qDeltaHalf=0.0, pDeltaHalf=0.0;
    pDeltaHalf = (dTime - pLastUpdateTime + 0.5*pDt)/a;
    qDeltaHalf = (dTime - qLastUpdateTime + 0.5*qDt)/a;

    // DEBUG: Avoid temporal extrapolation
    //pDeltaHalf = 0.;
    //qDeltaHalf = 0.;

    ftype modApq;
    std::array<ftype, 3> face_unit;
    computeFace(modApq, face_unit,
                rpq,   dx,  dy,  dz,
                pBall,  qBall,  pOmega,  qOmega,
                pBxx,  pBxy,  pBxz,
                pByy,  pByz,  pBzz,
                qBxx, qBxy, qBxz,
                qByy, qByz, qBzz);


    // Velocity of the quadrature mid-point
    ftype vFrame[3];
    vFrame[0] = 0.5*(pVpredx+qVpredx);
    vFrame[1] = 0.5*(pVpredy+qVpredy);
    vFrame[2] = 0.5*(pVpredz+qVpredz);

    ftype pv[3], qv[3];
    // We boost to the reference of the p-q 'face'
    pv[0] = pVpredx - vFrame[0];
    qv[0] = qVpredx - vFrame[0];
    pv[1] = pVpredy - vFrame[1];
    qv[1] = qVpredy - vFrame[1];
    pv[2] = pVpredz - vFrame[2];
    qv[2] = qVpredz - vFrame[2];

    // Mid-point rule
    dx  = -0.5*dx;
    dy  = -0.5*dy;
    dz  = -0.5*dz;

    // DEBUG: Avoid spatial extrapolation
    //dx = 0.0;
    //dy = 0.0;
    //dz = 0.0;


    // Divergence of the velocity field for the forward in time prediction
    ftype pdivv = (pGradVxX + pGradVyY + pGradVzZ);
    ftype qdivv = (qGradVxX + qGradVyY + qGradVzZ);

    ftype L_v[3], R_v[3];
    ftype L_rho = pDensity;
    ftype R_rho = qDensity;
    ftype L_p = pP;
    ftype R_p = qP;
    L_v[0] = pv[0];
    R_v[0] = qv[0];
    L_v[1] = pv[1];
    R_v[1] = qv[1];
    L_v[2] = pv[2];
    R_v[2] = qv[2];

    extrapolateStateInTime(
        L_rho,
        L_v[0], L_v[1], L_v[2],
        L_p,
        vFrame[0],  vFrame[1], vFrame[2],
        pDensity, pv[0], pv[1], pv[2], pP,
        pDeltaHalf,
        pGradRhoX, pGradRhoY, pGradRhoZ,
        pGradPX, pGradPY, pGradPZ,
        pLastAccX, pLastAccY, pLastAccZ, pdivv,
        dConstGamma, a);

    extrapolateVariableInSpace( L_rho, dx,dy, dz,
                                pGradRhoX, pGradRhoY, pGradRhoZ);
    extrapolateVariableInSpace( L_p, dx,dy, dz,
                                pGradPX, pGradPY, pGradPZ);
    extrapolateVariableInSpace( L_v[0], dx,dy, dz,
                                pGradVxX, pGradVxY, pGradVxZ);
    extrapolateVariableInSpace( L_v[1], dx,dy, dz,
                                pGradVyX, pGradVyY, pGradVyZ);
    extrapolateVariableInSpace( L_v[2], dx,dy, dz,
                                pGradVzX, pGradVzY, pGradVzZ);



    dx = -dx;
    dy = -dy;
    dz = -dz;

    extrapolateStateInTime(
        R_rho,
        R_v[0], R_v[1], R_v[2],
        R_p,
        vFrame[0],  vFrame[1], vFrame[2],
        qDensity, qv[0], qv[1], qv[2], qP,
        qDeltaHalf,
        qGradRhoX, qGradRhoY, qGradRhoZ,
        qGradPX, qGradPY, qGradPZ,
        qLastAccX, qLastAccY, qLastAccZ, qdivv,
        dConstGamma, a);

    extrapolateVariableInSpace( R_rho, dx,dy, dz,
                                qGradRhoX, qGradRhoY, qGradRhoZ);
    extrapolateVariableInSpace( R_p, dx,dy, dz,
                                qGradPX, qGradPY, qGradPZ);
    extrapolateVariableInSpace( R_v[0], dx,dy, dz,
                                qGradVxX, qGradVxY, qGradVxZ);
    extrapolateVariableInSpace( R_v[1], dx,dy, dz,
                                qGradVyX, qGradVyY, qGradVyZ);
    extrapolateVariableInSpace( R_v[2], dx,dy, dz,
                                qGradVzX, qGradVzY, qGradVzZ);

    genericPairwiseLimiter<ftype>(pDensity, qDensity, &L_rho, &R_rho);
    genericPairwiseLimiter<ftype>(pP, qP, &L_p, &R_p);
    genericPairwiseLimiter<ftype>(pv[0], qv[0], &L_v[0], &R_v[0]);
    genericPairwiseLimiter<ftype>(pv[1], qv[1], &L_v[1], &R_v[1]);
    genericPairwiseLimiter<ftype>(pv[2], qv[2], &L_v[2], &R_v[2]);

    if (bComove) {
        extrapolateCosmology(
            R_v[0], R_v[1], R_v[2],
            vFrame[0],  vFrame[1], vFrame[2],
            R_p,
            qDeltaHalf,
            qv[0], qv[1], qv[2], qP,
            dConstGamma, H, a);

        extrapolateCosmology(
            L_v[0], L_v[1], L_v[2],
            vFrame[0],  vFrame[1], vFrame[2],
            L_p,
            pDeltaHalf,
            pv[0], pv[1], pv[2], pP,
            dConstGamma, H, a);
    }

    /*
    if (L_rho < 0) {
        L_rho = pDensity;
    }
    if (R_rho < 0) {
        R_rho = qDensity;
    }
    if (L_p < 0) {
        L_p = pP;
    }
    if (R_p < 0) {
        R_p = qP;
    }
    */

#ifdef EEOS_POLYTROPE
    const double pLpoly =
        polytropicPressureFloor(a_inv3, L_rho, smf->dConstGamma,
                                smf->dEOSPolyFloorIndex, smf->dEOSPolyFloorDen, smf->dEOSPolyFlooru);
    const double pRpoly =
        polytropicPressureFloor(a_inv3, R_rho, smf->dConstGamma,
                                smf->dEOSPolyFloorIndex, smf->dEOSPolyFloorDen, smf->dEOSPolyFlooru);
    L_p = MAX(L_p, pLpoly);
    R_p = MAX(R_p, pRpoly);
#endif
#ifdef EEOS_JEANS
    const double pLjeans =
        jeansPressureFloor(L_rho, ph, smf->dConstGamma, smf->dEOSNJeans);
    const double pRjeans =
        jeansPressureFloor(R_rho, q(ball), smf->dConstGamma, smf->dEOSNJeans);
    L_p = MAX(L_p, pLjeans);
    R_p = MAX(R_p, pRjeans);
#endif

    ftype P_M, S_M;
    Riemann_solver_exact(dConstGamma,
                         R_rho, R_p, R_v,
                         L_rho, L_p, L_v,
                         P_M, S_M,
                         &F_rho, &F_p, F_v.data(),
                         face_unit.data());



#ifdef ENTROPY_SWITCH

#ifdef USE_MFM
    // As we are in a truly lagrangian configuration,
    // there is no advection of entropy among particles.
#else
    F_S=0;
    for (auto j=0; j<3; j++) F_S += F_v[j]*face_unit[j];
    if (F_S > 0) {
        // Maybe this values should be properly extrapolated to the faces..
        // but this is expensive!
        F_S *= psph->S*pDensity/pkdMass(pkd,p);
    }
    else {
        F_S *= q(S)*q(rho)/q(mass);
    }
    F_S *= F_rho*modApq;
    F_S = 0.;

#endif //USE_MFM
#endif //ENTROPY_SWITCH

#ifdef USE_MFM
    F_rho = 0.;
    F_p = P_M * S_M;
    for (auto j=0; j<3; j++)
        F_v[j] = P_M * face_unit[j];
#endif
    // End MFM

    // Force 2D
#ifdef FORCE_1D
    F_v[2] = 0.;
    F_v[1] = 0.;
#endif
#ifdef FORCE_2D
    F_v[2] = 0.;
#endif


    /*
    // Check for NAN fluxes
    if (F_rho!=F_rho)
        F_rho = 0.;//abort();
    if (F_p!=F_p)
        F_p = 0.;//abort();
    */




    // Now we de-boost the fluxes following Eq. A8 Hopkins 2015
    for (auto j=0; j<3; j++) {
        F_p += vFrame[j] * F_v[j];
        F_p += (0.5*vFrame[j]*vFrame[j])*F_rho;
    }

    // Now we just multiply by the face area
    F_p *= modApq;
    F_rho *= modApq;
    for (auto j=0; j<3; j++) {
        F_v[j] *= modApq;
        F_v[j] += vFrame[j]*F_rho;
    }


    // Phew! Done ;)
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
#define q(X)    input_buffer[q_##X * nBuff + i]
#define qout(X) output_buffer[out_##X * nBuff + i]
void hydroRiemann_simd(PARTICLE *pIn,float fBall,int nSmooth, int nBuff,
                       my_real *restrict input_buffer,
                       my_real *restrict output_buffer, SMF *smf) {
    PKD pkd = smf->pkd;
    auto P = pkd->particles[pIn];

    const auto &pv = P.velocity();
    auto &psph = P.sph();

    const my_real pDensity = P.density();
    const my_real p_omega = psph.omega;
    double a_inv3 = 1./(smf->a * smf->a * smf->a);



    dvec dTime, dDelta,  a,  H,  dConstGamma;
    bool bComove = pkd->csm->val.bComove;
    dTime = smf->dTime;
    dDelta = smf->dDelta;
    a = smf->a;
    H = smf->H;
    dConstGamma = smf->dConstGamma;

    dvec pomega = psph.omega;
    dvec ph     = 0.5*P.ball();
    dvec plast  = psph.lastUpdateTime;
    dvec pDt    = smf->dDelta/(1<<P.rung());
    dvec pBXX   = psph.B[XX];
    dvec pBXY   = psph.B[XY];
    dvec pBXZ   = psph.B[XZ];
    dvec pBYY   = psph.B[YY];
    dvec pBYZ   = psph.B[YZ];
    dvec pBZZ   = psph.B[ZZ];
    dvec pDens  = pDensity;
    dvec pVpredx= P.velocity()[0];
    dvec pVpredy= P.velocity()[1];
    dvec pVpredz= P.velocity()[2];
    dvec pPres  = psph.P;
    dvec pS;
#ifdef ENTROPY_SWITCH
    pS    = psph.S;
#endif
    dvec pgradRhox = psph.gradRho[0];
    dvec pgradRhoy = psph.gradRho[1];
    dvec pgradRhoz = psph.gradRho[2];
    dvec pgradPx = psph.gradP[0];
    dvec pgradPy = psph.gradP[1];
    dvec pgradPz = psph.gradP[2];
    dvec pgradVxx = psph.gradVx[0];
    dvec pgradVxy = psph.gradVx[1];
    dvec pgradVxz = psph.gradVx[2];
    dvec pgradVyx = psph.gradVy[0];
    dvec pgradVyy = psph.gradVy[1];
    dvec pgradVyz = psph.gradVy[2];
    dvec pgradVzx = psph.gradVz[0];
    dvec pgradVzy = psph.gradVz[1];
    dvec pgradVzz = psph.gradVz[2];
    dvec plastAccx = psph.lastAcc[0];
    dvec plastAccy = psph.lastAcc[1];
    dvec plastAccz = psph.lastAcc[2];

    dvec qdr;
    dvec qdx;
    dvec qdy;
    dvec qdz;
    dvec qomega;
    dvec qh;
    dvec qlast;
    dvec qDt;
    dvec qBXX;
    dvec qBXY;
    dvec qBXZ;
    dvec qBYY;
    dvec qBYZ;
    dvec qBZZ;
    dvec qDens;
    dvec qvx;
    dvec qvy;
    dvec qvz;
    dvec qP;
    dvec qS;
    dvec qgradRhox;
    dvec qgradRhoy;
    dvec qgradRhoz;
    dvec qgradPx;
    dvec qgradPy;
    dvec qgradPz;
    dvec qgradVxx;
    dvec qgradVxy;
    dvec qgradVxz;
    dvec qgradVyx;
    dvec qgradVyy;
    dvec qgradVyz;
    dvec qgradVzx;
    dvec qgradVzy;
    dvec qgradVzz;
    dvec qlastAccx;
    dvec qlastAccy;
    dvec qlastAccz;

#ifdef __INTEL_COMPILER
    __assume_aligned(input_buffer, 64);
#pragma simd
#pragma vector aligned
#endif
#ifdef __GNUC__
//TODO Trick GCC into autovectorizing this!!
#endif
#pragma forceinline
//#pragma clang loop vectorize(assume_safety)
//#pragma clang loop vectorize(enable)
//TODO make sure about nSmooth limit and out-of-bounds
// This assert will be activated if nSmooth is close to nBuff, a permanent
// solution to this is needed
    assert(nBuff>(nSmooth+dvec::width()-1));
    for (auto i=0; i<nSmooth; i+=dvec::width()) {
        dvec F_rho;
        std::array<dvec,3> F_v;
        dvec F_P;
        dvec F_S;
        dvec minDt;

        qdr.load(       &q(dr));
        qdx.load(       &q(dx));
        qdy.load(       &q(dy));
        qdz.load(       &q(dz));
        qomega.load(    &q(omega));
        qh.load(        &q(ball));
        qlast.load(     &q(lastUpdateTime));
        qDt.load(       &q(rung));
        qBXX.load(      &q(B_XX));
        qBXY.load(      &q(B_XY));
        qBXZ.load(      &q(B_XZ));
        qBYY.load(      &q(B_YY));
        qBYZ.load(      &q(B_YZ));
        qBZZ.load(      &q(B_ZZ));
        qDens.load(     &q(rho));
        qvx.load(       &q(vx));
        qvy.load(       &q(vy));
        qvz.load(       &q(vz));
        qP.load(        &q(P));
#ifdef ENTROPY_SWITCH
        qS.load(        &q(S));
#endif
        qgradRhox.load( &q(gradRhoX));
        qgradRhoy.load( &q(gradRhoY));
        qgradRhoz.load( &q(gradRhoZ));
        qgradPx.load(   &q(gradPX));
        qgradPy.load(   &q(gradPY));
        qgradPz.load(   &q(gradPZ));
        qgradVxx.load(  &q(gradVxX));
        qgradVxy.load(  &q(gradVxY));
        qgradVxz.load(  &q(gradVxZ));
        qgradVyx.load(  &q(gradVyX));
        qgradVyy.load(  &q(gradVyY));
        qgradVyz.load(  &q(gradVyZ));
        qgradVzx.load(  &q(gradVzX));
        qgradVzy.load(  &q(gradVzY));
        qgradVzz.load(  &q(gradVzZ));
        qlastAccx.load( &q(lastAccX));
        qlastAccy.load( &q(lastAccY));
        qlastAccz.load( &q(lastAccZ));

        if (smf->dDelta <= 0.0)
            pDt = qDt = plast = qlast = dTime = 0;


        doSinglePPFlux<dvec>( F_rho, F_v, F_P, F_S, minDt,
                              bComove, dTime, dDelta,  a,  H,  dConstGamma,
                              qdr,  qdx,  qdy,  qdz,
                              ph,  plast,  pDt,
                              pomega,
                              pBXX,  pBXY,  pBXZ,
                              pBYY,  pBYZ,  pBZZ,
                              pDens,  pVpredx, pVpredy, pVpredz,  pPres,  pS,
                              pgradRhox, pgradRhoy, pgradRhoz,
                              pgradPx,   pgradPy,   pgradPz,
                              pgradVxx,  pgradVxy,  pgradVxz,
                              pgradVyx,  pgradVyy,  pgradVyz,
                              pgradVzx,  pgradVzy,  pgradVzz,
                              plastAccx, plastAccy, plastAccz,
                              qh,  qlast,  qDt,
                              qomega,
                              qBXX,  qBXY,  qBXZ,
                              qBYY,  qBYZ,  qBZZ,
                              qDens,  qvx, qvy, qvz,  qP,  qS,
                              qgradRhox,  qgradRhoy,  qgradRhoz,
                              qgradPx,    qgradPy,    qgradPz,
                              qgradVxx,   qgradVxy,   qgradVxz,
                              qgradVyx,   qgradVyy,   qgradVyz,
                              qgradVzx,   qgradVzy,   qgradVzz,
                              qlastAccx,  qlastAccy,  qlastAccz) ;

        // We fill the output buffer with the fluxes, which then
        // will be added to the corresponding particles
        F_rho.store(&output_buffer[out_Frho * nBuff + i]);
        F_P.store(&output_buffer[out_Fene * nBuff + i]);
        F_v[0].store(&output_buffer[out_FmomX * nBuff + i]);
        F_v[1].store(&output_buffer[out_FmomY * nBuff + i]);
        F_v[2].store(&output_buffer[out_FmomZ * nBuff + i]);
        minDt.store(&output_buffer[out_minDt * nBuff + i]);

        /* Repeat without using simd instruction for comparison
        double F_rho0;
        std::array<double,3> F_v0;
        double F_P0;
        double F_S0;
        double minDt0;

        double pS0, qS0;
        double ph0 = pkdBall(pkd,p);
        double pDt0 = smf->dDelta/(1<<p->uRung);

        for (auto k=0; k<SIMD_DWIDTH; k++){
           double qDt0 = q(rung);
        doSinglePPFlux( F_rho0, F_v0, F_P0, F_S0, minDt0,
                        bComove, smf->dTime, smf->dDelta,  smf->a,  smf->H,  smf->dConstGamma,
                        q(dr),  q(dx),  q(dy),  q(dz),
                        ph0,  psph->lastUpdateTime,  pDt0,
                        psph->omega,
                        psph->B[XX],  psph->B[XY],  psph->B[XZ],
                        psph->B[YY],  psph->B[YZ],  psph->B[ZZ],
                        pDensity,  psph->vPred[0], psph->vPred[1], psph->vPred[2], psph->P,  pS0,
                        psph->gradRho[0], psph->gradRho[1], psph->gradRho[2],
                        psph->gradP[0],   psph->gradP[1],   psph->gradP[2],
                        psph->gradVx[0],  psph->gradVx[1],  psph->gradVx[2],
                        psph->gradVy[0],  psph->gradVy[1],  psph->gradVy[2],
                        psph->gradVz[0],  psph->gradVz[1],  psph->gradVz[2],
                        psph->lastAcc[0], psph->lastAcc[1], psph->lastAcc[2],
                        q(ball),  q(lastUpdateTime),  qDt0,
                        q(omega),
                        q(B_XX),  q(B_XY),  q(B_XZ),
                        q(B_YY),  q(B_YZ),  q(B_ZZ),
                        q(rho),  q(vx), q(vy), q(vz),  q(P),  qS0,
                        q(gradRhoX),  q(gradRhoY),  q(gradRhoZ),
                        q(gradPX),    q(gradPY),    q(gradPZ),
                        q(gradVxX),   q(gradVxY),   q(gradVxZ),
                        q(gradVyX),   q(gradVyY),   q(gradVyZ),
                        q(gradVzX),   q(gradVzY),   q(gradVzZ),
                        q(lastAccX),  q(lastAccY),  q(lastAccZ) );
           i++;
        }
        i -= SIMD_DWIDTH;
        */

    } // End of loop over neighbors
}

void hydroRiemann(PARTICLE *pIn,float fBall,int nSmooth, int nBuff,
                  my_real *restrict input_buffer,
                  my_real *restrict output_buffer, SMF *smf) {
    PKD pkd = smf->pkd;

    auto P = pkd->particles[pIn];
    auto &psph = P.sph();
    const auto &pv = P.velocity();

    my_real pDensity = P.density();
    double a_inv3 = 1./(smf->a * smf->a * smf->a);



    bool bComove = pkd->csm->val.bComove;


#ifdef __INTEL_COMPILER
    __assume_aligned(input_buffer, 64);
#pragma simd
#pragma vector aligned
#endif
#ifdef __GNUC__
//TODO Trick GCC into autovectorizing this!!
#endif
#pragma forceinline
//#pragma clang loop vectorize(assume_safety)
//#pragma clang loop vectorize(enable)
    for (auto i=0; i<nSmooth; i++) {
        double F_rho;
        std::array<double,3> F_v;
        double F_P;
        double F_S;
        double minDt;

        double pS, qS;
        double ph = P.ball();
        double qDt = q(rung);
        double pDt = smf->dDelta/(1<<P.rung());

        doSinglePPFlux<double>( F_rho, F_v, F_P, F_S, minDt,
                        bComove, smf->dTime, smf->dDelta,  smf->a,  smf->H,  smf->dConstGamma,
                        q(dr),  q(dx),  q(dy),  q(dz),
                        ph,  psph.lastUpdateTime,  pDt,
                        psph.omega,
                        psph.B[XX],  psph.B[XY],  psph.B[XZ],
                        psph.B[YY],  psph.B[YZ],  psph.B[ZZ],
                        pDensity,  pv[0], pv[1], pv[2], psph.P,  pS,
                        psph.gradRho[0], psph.gradRho[1], psph.gradRho[2],
                        psph.gradP[0],   psph.gradP[1],   psph.gradP[2],
                        psph.gradVx[0],  psph.gradVx[1],  psph.gradVx[2],
                        psph.gradVy[0],  psph.gradVy[1],  psph.gradVy[2],
                        psph.gradVz[0],  psph.gradVz[1],  psph.gradVz[2],
                        psph.lastAcc[0], psph.lastAcc[1], psph.lastAcc[2],
                        q(ball),  q(lastUpdateTime),  qDt,
                        q(omega),
                        q(B_XX),  q(B_XY),  q(B_XZ),
                        q(B_YY),  q(B_YZ),  q(B_ZZ),
                        q(rho),  q(vx), q(vy), q(vz),  q(P),  qS,
                        q(gradRhoX),  q(gradRhoY),  q(gradRhoZ),
                        q(gradPX),    q(gradPY),    q(gradPZ),
                        q(gradVxX),   q(gradVxY),   q(gradVxZ),
                        q(gradVyX),   q(gradVyY),   q(gradVyZ),
                        q(gradVzX),   q(gradVzY),   q(gradVzZ),
                        q(lastAccX),  q(lastAccY),  q(lastAccZ) );
        //printf("%e %e\n", F_P, qout(Fene));
        qout(Frho) = F_rho;
        qout(Fene) = F_P;
        qout(FmomX) = F_v[0];
        qout(FmomY) = F_v[1];
        qout(FmomZ) = F_v[2];
#ifdef ENTROPY_SWITCH
        qout(FS) = F_S;
#endif
        qout(minDt) = minDt;

    } // End of loop over neighbors
}


void hydroFluxFillBuffer(my_real *input_buffer, PARTICLE *qIn, int i, int nBuff,
                         double dr2, TinyVector<double,3> dr  , SMF *smf) {
    PKD pkd = smf->pkd;
    auto Q = pkd->particles[qIn];
    double dDelta = smf->dDelta;
    float qh = 0.5*Q.ball();
    auto &qsph = Q.sph();
    q(mass) = Q.mass();
    q(ball) = qh;
    q(dx) = dr[0];
    q(dy) = dr[1];
    q(dz) = dr[2];
    q(dr) = sqrt(dr2);
    q(rung) = dDelta/(1<<Q.rung());
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
}


void hydroFluxUpdateFromBuffer(my_real *output_buffer, my_real *input_buffer,
                               PARTICLE *pIn, PARTICLE *qIn, int i, int nBuff, SMF *smf) {
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
        P.set_mass(P.mass() - qout(minDt) * qout(Frho));

        psph.mom[0] -= qout(minDt) * qout(FmomX);
        psph.mom[1] -= qout(minDt) * qout(FmomY);
        psph.mom[2] -= qout(minDt) * qout(FmomZ);

        psph.E -= qout(minDt) * qout(Fene);

        psph.Uint -= qout(minDt) * ( qout(Fene)
                     - qout(FmomX)*pv[0]
                     - qout(FmomY)*pv[1]
                     - qout(FmomZ)*pv[2]
                     + 0.5*dot(pv,pv)*qout(Frho) );

#ifdef ENTROPY_SWITCH
        psph.S -= qout(minDt)* qout(FS);
#endif

#ifndef USE_MFM
        psph.drDotFrho[0] += qout(minDt) * qout(Frho) * q(dx) * aFac;
        psph.drDotFrho[1] += qout(minDt) * qout(Frho) * q(dy) * aFac;
        psph.drDotFrho[2] += qout(minDt) * qout(Frho) * q(dz) * aFac;
#endif
        psph.Frho +=    qout(Frho);
        psph.Fene +=    qout(Fene);
        psph.Fmom[0] += qout(FmomX);
        psph.Fmom[1] += qout(FmomY);
        psph.Fmom[2] += qout(FmomZ);
    }
    else {
        psph.Frho +=    qout(Frho);
        psph.Fene +=    qout(Fene);
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
        assert( qsph.P == q(P) );
        if (dDelta>0) {
            Q.set_mass(Q.mass() + qout(minDt) * qout(Frho));

            qsph.mom[0] += qout(minDt) * qout(FmomX);
            qsph.mom[1] += qout(minDt) * qout(FmomY);
            qsph.mom[2] += qout(minDt) * qout(FmomZ);

            qsph.E += qout(minDt) * qout(Fene);

            qsph.Uint += qout(minDt) * ( qout(Fene)
                         - qout(FmomX)*qv[0]
                         - qout(FmomY)*qv[1]
                         - qout(FmomZ)*qv[2]
                         + 0.5*dot(qv,qv)*qout(Frho) );
#ifdef ENTROPY_SWITCH
            qsph.S += qout(minDt) * qout(FS);
#endif

#ifndef USE_MFM
            qsph.drDotFrho[0] += qout(minDt) * qout(Frho) * q(dx) * aFac;
            qsph.drDotFrho[1] += qout(minDt) * qout(Frho) * q(dy) * aFac;
            qsph.drDotFrho[2] += qout(minDt) * qout(Frho) * q(dz) * aFac;
#endif
            qsph.Frho -=    qout(Frho);
            qsph.Fene -=    qout(Fene);
            qsph.Fmom[0] -= qout(FmomX);
            qsph.Fmom[1] -= qout(FmomY);
            qsph.Fmom[2] -= qout(FmomZ);
        }
        else {
            qsph.Frho -=    qout(Frho);
            qsph.Fene -=    qout(Fene);
            qsph.Fmom[0] -= qout(FmomX);
            qsph.Fmom[1] -= qout(FmomY);
            qsph.Fmom[2] -= qout(FmomZ);
        }

    } // q marked/active

}

void hydroFluxGetNvars(int *in, int *out) {
    *in = q_last;
    *out = out_last;
}

#endif // OPTIM_FLUX_VEC

