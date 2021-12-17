#include "hydro/hydro.h"
#include "master.h"
#include "eEOS/eEOS.h"
#ifdef OPTIM_FLUX_VEC
#include "riemann_own.h"
#else
#include "riemann.h"
#endif


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

extern "C"
void pkdResetFluxes(PKD pkd, double dTime,double dDelta,double dDeltaVPred,double dDeltaTime) {
    PARTICLE *p;
    SPHFIELDS *psph;
    float* pmass;
    int i;
    int pLower, pUpper;

   pLower = 0;
   pUpper = pkdLocal(pkd);

    assert(pkd->oFieldOffset[oVelocity]);
    //assert(pkd->oMass);

    /*
    ** Add the computed flux to the conserved variables for each gas particle
    */
    assert(pkd->oFieldOffset[oSph]);
    for (i=pLower;i<pUpper;++i) { 
    p = pkdParticle(pkd,i);
       if (pkdIsGas(pkd,p)  && pkdIsActive(pkd,p)   ) {
          psph = pkdSph(pkd, p);
          psph->Frho = 0.0;
          psph->Fene = 0.0;
          psph->Fmom[0] = 0.0;
          psph->Fmom[1] = 0.0;
          psph->Fmom[2] = 0.0;
       }
    }

    }

void MSR::MeshlessFluxes(double dTime,double dDelta)
{
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

void initHydroFluxes(void *vpkd, void *vp)
{
}

/* Zero all the conserved quantities, which will be updated
 * during the hydro loop.
 *
 * Then those will be merged with the actual particle information inside
 * combThirdHydroLoop
 */
void initHydroFluxesCached(void *vpkd, void *vp)
{
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = (PARTICLE *) vp;
    assert(!pkd->bNoParticleOrder);
    // For the init*Cached and comb functions we still have to explicitly
    // check if we are handling gas particles even ifdef OPTIM_REORDER_IN_NODES
    // because these operations are done in a cache-line basis, and a line may
    // contain other particles that are not of interest!
    if (pkdIsGas(pkd,p)) {
        SPHFIELDS *psph = pkdSph(pkd,p);
        int i;

        float *pmass = (float *) pkdField(p,pkd->oFieldOffset[oMass]);
        *pmass = 0.0;
        psph->mom[0] = 0.;
        psph->mom[1] = 0.;
        psph->mom[2] = 0.;
        psph->E = 0.;
        psph->Uint = 0.;

        if (pkdIsActive(pkd,p)) {
            psph->Frho = 0.0;
            psph->Fene = 0.0;
            for (i=0; i<3; i++) {
                psph->Fmom[i] = 0.0;
            }
        }
    }
}


/* This version is deprecated. It may be removed in future versions of the code
 * without notice.
 *
 * The maintained Riemann solver is the vectorized version of this function,
 * which is activated wit the OPTIM_SMOOTH_NODE and OPTIM_FLUX_VEC flags
 */
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf)
{
    //TODO Clean unused variables!
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    float *pmass, *qmass;
    double minDt;
    double pv[3], qv[3], vFrame[3];
    double ph,qh, hpq, modApq, rpq, dx,dy,dz,Wpq,psi_p,psi_q,pDensity,pDeltaHalf, qDeltaHalf;
    double pdivv, qdivv, psi;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3], dr[3];
    struct Input_vec_Riemann riemann_input;
    struct Riemann_outputs riemann_output;
    int i,j;

    psph = pkdSph(pkd, p);
    ph = pkdBall(pkd, p);

    pDensity = pkdDensity(pkd,p);

    double a_inv3 = 1./(smf->a * smf->a * smf->a);

    for (i=0; i<nSmooth; ++i) {

        q = nnList[i].pPart;
        qsph = pkdSph(pkd, q);
        qh = pkdBall(pkd,q);

        dx = nnList[i].dx;
        dy = nnList[i].dy;
        dz = nnList[i].dz;
#ifdef FORCE_1D
        if (dz!=0) continue;
        if (dy!=0) continue;
#endif
#ifdef FORCE_2D
        if (dz!=0) continue;
#endif

        /* In the nnList there is a 'copy' of the own particle,
         * which we can omit as there are no fluxes to be computed here
         */
        if (dx==0 && dy==0 && dz==0) continue;

        hpq = ph;
        rpq = sqrt(nnList[i].fDist2);
        // We only compute the fluxes if both particles are within the kernel
        // of each other
        if (2.*qh < rpq) continue;

        Wpq = cubicSplineKernel(rpq, hpq);
        if (Wpq==0.0) {
            continue;
        }

        /* We update the conservatives variables taking the minimum timestep
         * between the particles, as in AREPO */
        if (!pkdIsActive(pkd,q)) {
            // If q is not active we now that p has the smallest dt
            minDt = smf->dDelta/(1<<p->uRung) ;
        } else {
            // Otherwise we need to explicitly check
            if (p->uRung > q->uRung) {
                minDt = smf->dDelta/(1<<p->uRung) ;
            } else {
                minDt = smf->dDelta/(1<<q->uRung) ;
            }
        }

        if (smf->dDelta > 0) {
            pDeltaHalf = smf->dTime - psph->lastUpdateTime + 0.5*smf->dDelta/(1<<p->uRung);
            qDeltaHalf = smf->dTime - qsph->lastUpdateTime + 0.5*smf->dDelta/(1<<q->uRung);
        } else {
            /* For the initialization step we do not extrapolate because we
             * dont have a reliable dDelta
             */
            qDeltaHalf = 0.0;
            pDeltaHalf = 0.0;
        }
        if(pkd->csm->val.bComove) {
            qDeltaHalf /= smf->a;
            pDeltaHalf /= smf->a;
        }

        // DEBUG: Avoid temporal extrapolation
        //pDeltaHalf = 0.;
        //qDeltaHalf = 0.;



        // \tilde{\psi}_j (x_i)
        psi = -cubicSplineKernel(rpq, ph)/psph->omega;
        psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
        psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
        psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

        // \tilde{\psi}_i (x_j)
        psi = cubicSplineKernel(rpq, qh)/qsph->omega;
        psiTilde_q[0] = (qsph->B[XX]*dx + qsph->B[XY]*dy + qsph->B[XZ]*dz)*psi;
        psiTilde_q[1] = (qsph->B[XY]*dx + qsph->B[YY]*dy + qsph->B[YZ]*dz)*psi;
        psiTilde_q[2] = (qsph->B[XZ]*dx + qsph->B[YZ]*dy + qsph->B[ZZ]*dz)*psi;

        modApq = 0.0;
        for (j=0; j<3; j++) {
            Apq[j] = psiTilde_p[j]/psph->omega - psiTilde_q[j]/qsph->omega;
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


        for (j=0; j<3; j++) {
            face_unit[j] = Apq[j]/modApq;
        }



        // Velocity of the quadrature mid-point
        for (j=0; j<3; j++) {
            vFrame[j] = 0.5*(psph->vPred[j]+qsph->vPred[j]);

            // We boost to the reference of the p-q 'face'
            pv[j] = psph->vPred[j] - vFrame[j];
            qv[j] = qsph->vPred[j] - vFrame[j];
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
        pdivv = psph->gradVx[0] + psph->gradVy[1] + psph->gradVz[2];
        qdivv = qsph->gradVx[0] + qsph->gradVy[1] + qsph->gradVz[2];

        pdivv *= pDeltaHalf;
        qdivv *= qDeltaHalf;


        riemann_input.L.rho = pkdDensity(pkd,p);
        riemann_input.R.rho = pkdDensity(pkd,q);
        riemann_input.L.v[0] = pv[0];
        riemann_input.R.v[0] = qv[0];
        riemann_input.L.v[1] = pv[1];
        riemann_input.R.v[1] = qv[1];
        riemann_input.L.v[2] = pv[2];
        riemann_input.R.v[2] = qv[2];
        riemann_input.L.p = psph->P;
        riemann_input.R.p = qsph->P;

//      printf("1) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("1) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

        // We add the gradients terms (from extrapolation and forward prediction)
        for (j=0; j<3; j++) {
            riemann_input.L.rho += ( dr[j] - pDeltaHalf*pv[j])*psph->gradRho[j];
            riemann_input.R.rho += (-dr[j] - qDeltaHalf*qv[j])*qsph->gradRho[j];

            riemann_input.L.v[0] += ( dr[j]*psph->gradVx[j]);
            riemann_input.R.v[0] += (-dr[j]*qsph->gradVx[j]);

            riemann_input.L.v[1] += ( dr[j]*psph->gradVy[j]);
            riemann_input.R.v[1] += (-dr[j]*qsph->gradVy[j]);

            riemann_input.L.v[2] += ( dr[j]*psph->gradVz[j]);
            riemann_input.R.v[2] += (-dr[j]*qsph->gradVz[j]);

            riemann_input.L.p += ( dr[j] - pDeltaHalf*pv[j])*psph->gradP[j];
            riemann_input.R.p += (-dr[j] - qDeltaHalf*qv[j])*qsph->gradP[j];
        }
//      printf("2) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("2) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

        // Placing this here solved the convergence problem for the comoving soundwaves.
        //   This problem may be caused because we do not use the time extrapolated cell-centered states in
        //   this limiter
        /*
        genericPairwiseLimiter(pkdDensity(pkd,p), pkdDensity(pkd,q), &riemann_input.L.rho, &riemann_input.R.rho);
        genericPairwiseLimiter(psph->P, qsph->P, &riemann_input.L.p, &riemann_input.R.p);
        genericPairwiseLimiter(pv[0], qv[0], &riemann_input.L.v[0], &riemann_input.R.v[0]);
        genericPairwiseLimiter(pv[1], qv[1], &riemann_input.L.v[1], &riemann_input.R.v[1]);
        genericPairwiseLimiter(pv[2], qv[2], &riemann_input.L.v[2], &riemann_input.R.v[2]);
        */


        double temp;


        // Forward extrapolation of velocity
        for (j=0; j<3; j++) {
            temp = pv[j]*pdivv + psph->gradP[j]/pDensity*pDeltaHalf;
            riemann_input.L.v[j] -= temp;
            vFrame[j] -= 0.5*temp;

            temp = qv[j]*qdivv + qsph->gradP[j]/pkdDensity(pkd,q)*qDeltaHalf;
            riemann_input.R.v[j] -= temp;
            vFrame[j] -= 0.5*temp;
        }

        for (j=0; j<3; j++) {
            temp = psph->lastAcc[j]*pDeltaHalf*smf->a;
            riemann_input.L.v[j] += temp;
            vFrame[j] += 0.5*temp;

            temp = qsph->lastAcc[j]*qDeltaHalf*smf->a;
            riemann_input.R.v[j] += temp;
            vFrame[j] += 0.5*temp;
        }

        riemann_input.L.rho -= pDensity*pdivv;
        riemann_input.R.rho -= pkdDensity(pkd,q)*qdivv;
        riemann_input.L.p -= smf->dConstGamma*psph->P*pdivv;
        riemann_input.R.p -= smf->dConstGamma*qsph->P*qdivv;

        genericPairwiseLimiter(pkdDensity(pkd,p), pkdDensity(pkd,q), &riemann_input.L.rho, &riemann_input.R.rho);
        genericPairwiseLimiter(psph->P, qsph->P, &riemann_input.L.p, &riemann_input.R.p);
        genericPairwiseLimiter(pv[0], qv[0], &riemann_input.L.v[0], &riemann_input.R.v[0]);
        genericPairwiseLimiter(pv[1], qv[1], &riemann_input.L.v[1], &riemann_input.R.v[1]);
        genericPairwiseLimiter(pv[2], qv[2], &riemann_input.L.v[2], &riemann_input.R.v[2]);

        if(pkd->csm->val.bComove) {

            for (j=0; j<3; j++) {
                temp = smf->H * pDeltaHalf * smf->a * pv[j];
                riemann_input.L.v[j] -= temp;
                vFrame[j] -= 0.5*temp;

                temp = smf->H * qDeltaHalf * smf->a * qv[j];
                riemann_input.R.v[j] -= temp;
                vFrame[j] -= 0.5*temp;
            }

            riemann_input.L.p -= 3. * smf->H * pDeltaHalf * smf->a *
                                 (smf->dConstGamma - 1.) * psph->P;

            riemann_input.R.p -= 3. * smf->H * qDeltaHalf * smf->a *
                                 (smf->dConstGamma - 1.) * qsph->P;

        }

        // DEBUG: Tests for the riemann solver extracted from Toro (10.1007/b79761)
        // Test 1
//       riemann_input.L.rho = 1.0; riemann_input.L.p = 1.0; riemann_input.L.v[0] = 0.0;
//       riemann_input.L.rho = 0.125; riemann_input.L.p = 0.1; riemann_input.L.v[0] = 0.0;

        if (riemann_input.L.rho < 0) {
            riemann_input.L.rho = pkdDensity(pkd,p);
            /* printf("WARNING, L.rho < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.R.rho < 0) {
            riemann_input.R.rho = pkdDensity(pkd,q);
            /* printf("WARNING, R.rho < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.L.p < 0) {
            riemann_input.L.p = psph->P;
            /* printf("WARNING, L.p < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.R.p < 0) {
            riemann_input.R.p = qsph->P;
            /* printf("WARNING, R.p < 0 : using first-order scheme \n");*/
        }

#ifdef EEOS_POLYTROPE
        const double pLpoly =
            polytropicPressureFloor(a_inv3, riemann_input.L.rho, smf->dConstGamma, 
                  smf->dEOSPolyFloorIndex, smf->dEOSPolyFloorDen, smf->dEOSPolyFlooru);
        const double pRpoly =
            polytropicPressureFloor(a_inv3, riemann_input.R.rho, smf->dConstGamma,
                  smf->dEOSPolyFloorIndex, smf->dEOSPolyFloorDen, smf->dEOSPolyFlooru);
        riemann_input.L.p = MAX(riemann_input.L.p, pLpoly);
        riemann_input.R.p = MAX(riemann_input.R.p, pRpoly);
#endif
#ifdef EEOS_JEANS
        const double pLjeans =
            jeansPressureFloor(riemann_input.L.rho, ph, smf->dConstGamma, smf->dEOSNJeans);
        const double pRjeans =
            jeansPressureFloor(riemann_input.R.rho, qh, smf->dConstGamma, smf->dEOSNJeans);
        riemann_input.L.p = MAX(riemann_input.L.p, pLjeans);
        riemann_input.R.p = MAX(riemann_input.R.p, pRjeans);
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
        for(j=0; j<3; j++)
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


        if(pkd->csm->val.bComove)
            minDt /= smf->a; // 1/a term before \nabla



        // Now we de-boost the fluxes following Eq. A8 Hopkins 2015
        for(j=0; j<3; j++) {
            riemann_output.Fluxes.p += vFrame[j] * riemann_output.Fluxes.v[j];
            riemann_output.Fluxes.p += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
        }

        // Now we just multiply by the face area
        riemann_output.Fluxes.p *= modApq;
        riemann_output.Fluxes.rho *= modApq;
        for (j=0; j<3; j++) {
            riemann_output.Fluxes.v[j] *= modApq;
            riemann_output.Fluxes.v[j] += vFrame[j]*riemann_output.Fluxes.rho;
        }


        if (smf->dDelta > 0) {
            pmass = (float *) pkdField(p,pkd->oFieldOffset[oMass]);
            qmass = (float *) pkdField(q,pkd->oFieldOffset[oMass]);



#ifndef OPTIM_NO_REDUNDANT_FLUXES
            {
#else
            if ( (2.*qh < rpq) | !pkdIsActive(pkd,q)) {
#endif


                *qmass += minDt * riemann_output.Fluxes.rho ;

                qsph->mom[0] += minDt * riemann_output.Fluxes.v[0] ;
                qsph->mom[1] += minDt * riemann_output.Fluxes.v[1] ;
                qsph->mom[2] += minDt * riemann_output.Fluxes.v[2] ;

                qsph->E += minDt * riemann_output.Fluxes.p;

                qsph->Uint += minDt * ( riemann_output.Fluxes.p -
                                        riemann_output.Fluxes.v[0]*qsph->vPred[0] -
                                        riemann_output.Fluxes.v[1]*qsph->vPred[1] -
                                        riemann_output.Fluxes.v[2]*qsph->vPred[2] +
                                        0.5*(qsph->vPred[0]*qsph->vPred[0] +
                                             qsph->vPred[1]*qsph->vPred[1] +
                                             qsph->vPred[2]*qsph->vPred[2]) *
                                        riemann_output.Fluxes.rho );
#ifndef USE_MFM
                qsph->drDotFrho[0] += minDt * riemann_output.Fluxes.rho * dx;
                qsph->drDotFrho[1] += minDt * riemann_output.Fluxes.rho * dy;
                qsph->drDotFrho[2] += minDt * riemann_output.Fluxes.rho * dz;
#endif

                qsph->Frho -= riemann_output.Fluxes.rho;
                qsph->Fene -= riemann_output.Fluxes.p;
                for(j=0; j<3; j++) {
                    qsph->Fmom[j] -= riemann_output.Fluxes.v[j];
                }
            }

            *pmass -= minDt * riemann_output.Fluxes.rho ;

            psph->mom[0] -= minDt * riemann_output.Fluxes.v[0] ;
            psph->mom[1] -= minDt * riemann_output.Fluxes.v[1] ;
            psph->mom[2] -= minDt * riemann_output.Fluxes.v[2] ;

            psph->E -= minDt * riemann_output.Fluxes.p;

            psph->Uint -= minDt * ( riemann_output.Fluxes.p -
                                    riemann_output.Fluxes.v[0]*psph->vPred[0] -
                                    riemann_output.Fluxes.v[1]*psph->vPred[1] -
                                    riemann_output.Fluxes.v[2]*psph->vPred[2] +
                                    0.5*(psph->vPred[0]*psph->vPred[0] +
                                         psph->vPred[1]*psph->vPred[1] +
                                         psph->vPred[2]*psph->vPred[2]) *
                                    riemann_output.Fluxes.rho );
#ifndef USE_MFM
            psph->drDotFrho[0] += minDt * riemann_output.Fluxes.rho * dx;
            psph->drDotFrho[1] += minDt * riemann_output.Fluxes.rho * dy;
            psph->drDotFrho[2] += minDt * riemann_output.Fluxes.rho * dz;
#endif
        }
        /* Old fluxes update (see 15/04/19 )
         * TODO: This is not needed for the update of the conserved variables. Instead,
         * it is now only used for the acceleleration criteria. Memory-wise, this can be
         * substantially improved
         */
        // Own contribution is always added
        psph->Frho += riemann_output.Fluxes.rho;
        psph->Fene += riemann_output.Fluxes.p;
        for(j=0; j<3; j++) {
            psph->Fmom[j] += riemann_output.Fluxes.v[j];
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
void hydroRiemann_vec(PARTICLE *p,float fBall,int nSmooth,
                      my_real** restrict input_buffer,
                      my_real** restrict output_buffer, SMF *smf)
{
    PKD pkd = smf->pkd;
    int i,j;

    SPHFIELDS* psph = pkdSph(pkd, p);
    my_real ph = pkdBall(pkd, p);

    my_real pDensity = pkdDensity(pkd,p);
    my_real p_omega = psph->omega;
    double a_inv3 = 1./(smf->a * smf->a * smf->a);


#ifdef __INTEL_COMPILER
    __assume_aligned(input_buffer, 64);
    __assume_aligned(input_buffer[0], 64);
#pragma simd
#pragma vector aligned
#endif
#ifdef __GNUC__
//TODO Trick GCC into autovectorizing this!!
#endif
    for (i=0; i<nSmooth; ++i) {

        my_real qh = q(ball);

        my_real dx = q(dx);
        my_real dy = q(dy);
        my_real dz = q(dz);

#ifdef FORCE_1D
        if (dz!=0) continue;
        if (dy!=0) continue;
#endif
#ifdef FORCE_2D
        if (dz!=0) continue;
#endif



        // Face where the riemann problem will be solved
        my_real rpq = q(dr);


        /* We update the conservatives variables taking the minimum timestep
         * between the particles, as in AREPO
         */
        my_real p_dt = smf->dDelta/(1<<p->uRung);
        my_real q_dt = q(rung);
        my_real minDt =  p_dt > q_dt ? q_dt : p_dt;
        minDt /=  smf->a;


        my_real qDeltaHalf=0.0, pDeltaHalf=0.0;
        if (smf->dDelta > 0) {
            pDeltaHalf = (smf->dTime - psph->lastUpdateTime + 0.5*p_dt)/smf->a;
            qDeltaHalf = (smf->dTime - q(lastUpdateTime) + 0.5*q_dt)/smf->a;
        }

        // DEBUG: Avoid temporal extrapolation
        //pDeltaHalf = 0.;
        //qDeltaHalf = 0.;

        my_real omega_q = q(omega);

        // \tilde{\psi}_j (x_i)
        my_real psi = -cubicSplineKernel(rpq, ph)/p_omega;
        my_real psiTilde_p[3], psiTilde_q[3];
        psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
        psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
        psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

        // \tilde{\psi}_i (x_j)
        psi = cubicSplineKernel(rpq, qh)/omega_q;
        psiTilde_q[0] = (q(B_XX)*dx + q(B_XY)*dy + q(B_XZ)*dz)*psi;
        psiTilde_q[1] = (q(B_XY)*dx + q(B_YY)*dy + q(B_YZ)*dz)*psi;
        psiTilde_q[2] = (q(B_XZ)*dx + q(B_YZ)*dy + q(B_ZZ)*dz)*psi;

        my_real modApq = 0.0;
        my_real Apq[3];
        for (j=0; j<3; j++) {
            Apq[j] = psiTilde_p[j]/p_omega - psiTilde_q[j]/omega_q;
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


        my_real face_unit[3];
        for (j=0; j<3; j++) {
            face_unit[j] = Apq[j]/modApq;
        }


        // Velocity of the quadrature mid-point
        my_real vFrame[3];
        vFrame[0] = 0.5*(psph->vPred[0]+q(vx));
        vFrame[1] = 0.5*(psph->vPred[1]+q(vy));
        vFrame[2] = 0.5*(psph->vPred[2]+q(vz));

        my_real pv[3], qv[3];
        for (j=0; j<3; j++) {
            // We boost to the reference of the p-q 'face'
            pv[j] = psph->vPred[j] - vFrame[j];
        }

        qv[0] = q(vx) - vFrame[0];
        qv[1] = q(vy) - vFrame[1];
        qv[2] = q(vz) - vFrame[2];

        // Mid-point rule
        my_real dr[3];
        dr[0] = -0.5*dx;
        dr[1] = -0.5*dy;
        dr[2] = -0.5*dz;

        // DEBUG: Avoid spatial extrapolation
        //dr[0] = 0.0;
        //dr[1] = 0.0;
        //dr[2] = 0.0;

        // Divergence of the velocity field for the forward in time prediction
        my_real pdivv = (psph->gradVx[0] + psph->gradVy[1] + psph->gradVz[2])*pDeltaHalf;
        my_real qdivv = (q(gradVxX) + q(gradVyY) + q(gradVzZ))*qDeltaHalf;

        // At some point we should erase the need for this structs... FIXME
        struct Input_vec_Riemann riemann_input;
        struct Riemann_outputs riemann_output;

        riemann_input.L.rho = pDensity;
        riemann_input.R.rho = q(rho);
        riemann_input.L.v[0] = pv[0];
        riemann_input.R.v[0] = qv[0];
        riemann_input.L.v[1] = pv[1];
        riemann_input.R.v[1] = qv[1];
        riemann_input.L.v[2] = pv[2];
        riemann_input.R.v[2] = qv[2];
        riemann_input.L.p = psph->P;
        riemann_input.R.p = q(P);

//      printf("1) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("1) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

        // We add the gradients terms (from extrapolation and forward prediction)
        for (j=0; j<3; j++) {
            riemann_input.L.rho += ( dr[j] - pDeltaHalf*pv[j])*psph->gradRho[j];

            riemann_input.L.v[0] += ( dr[j]*psph->gradVx[j]);
            riemann_input.L.v[1] += ( dr[j]*psph->gradVy[j]);
            riemann_input.L.v[2] += ( dr[j]*psph->gradVz[j]);

            riemann_input.L.p += ( dr[j] - pDeltaHalf*pv[j])*psph->gradP[j];
        }



        riemann_input.R.rho += (-dr[0] - qDeltaHalf*qv[0])*q(gradRhoX);
        riemann_input.R.rho += (-dr[1] - qDeltaHalf*qv[1])*q(gradRhoY);
        riemann_input.R.rho += (-dr[2] - qDeltaHalf*qv[2])*q(gradRhoZ);

        riemann_input.R.v[0] += (-dr[0]*q(gradVxX));
        riemann_input.R.v[0] += (-dr[1]*q(gradVxY));
        riemann_input.R.v[0] += (-dr[2]*q(gradVxZ));

        riemann_input.R.v[1] += (-dr[0]*q(gradVyX));
        riemann_input.R.v[1] += (-dr[1]*q(gradVyY));
        riemann_input.R.v[1] += (-dr[2]*q(gradVyZ));

        riemann_input.R.v[2] += (-dr[0]*q(gradVzX));
        riemann_input.R.v[2] += (-dr[1]*q(gradVzY));
        riemann_input.R.v[2] += (-dr[2]*q(gradVzZ));

        riemann_input.R.p += (-dr[0] - qDeltaHalf*qv[0])*q(gradPX);
        riemann_input.R.p += (-dr[1] - qDeltaHalf*qv[1])*q(gradPY);
        riemann_input.R.p += (-dr[2] - qDeltaHalf*qv[2])*q(gradPZ);



//      printf("2) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("2) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);



        my_real temp;



        for (j=0; j<3; j++) { // Forward extrapolation of velocity
            temp = pv[j]*pdivv + psph->gradP[j]/pDensity*pDeltaHalf;
            riemann_input.L.v[j] -= temp;
            vFrame[j] -= 0.5*temp;
        }

        temp = qv[0]*qdivv + q(gradPX)/q(rho)*qDeltaHalf;
        riemann_input.R.v[0] -= temp;
        vFrame[0] -= 0.5*temp;

        temp = qv[1]*qdivv + q(gradPY)/q(rho)*qDeltaHalf;
        riemann_input.R.v[1] -= temp;
        vFrame[1] -= 0.5*temp;

        temp = qv[2]*qdivv + q(gradPZ)/q(rho)*qDeltaHalf;
        riemann_input.R.v[2] -= temp;
        vFrame[2] -= 0.5*temp;


        for (j=0; j<3; j++) {
            temp = psph->lastAcc[j]*pDeltaHalf*smf->a;
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
        riemann_input.L.p -= smf->dConstGamma*psph->P*pdivv;
        riemann_input.R.p -= smf->dConstGamma*q(P)*qdivv;

        genericPairwiseLimiter(pDensity, q(rho), &riemann_input.L.rho, &riemann_input.R.rho);
        genericPairwiseLimiter(psph->P, q(P), &riemann_input.L.p, &riemann_input.R.p);
        for (j=0; j<3; j++) {
            genericPairwiseLimiter(pv[j], qv[j], &riemann_input.L.v[j], &riemann_input.R.v[j]);
        }

        if(pkd->csm->val.bComove) {

            for (j=0; j<3; j++) {
                temp = smf->H * pDeltaHalf * smf->a * pv[j];
                riemann_input.L.v[j] -= temp;
                vFrame[j] -= 0.5*temp;

                temp = smf->H * qDeltaHalf * smf->a * qv[j];
                riemann_input.R.v[j] -= temp;
                vFrame[j] -= 0.5*temp;
            }

            riemann_input.L.p -= 3. * smf->H * pDeltaHalf * smf->a * (smf->dConstGamma - 1.) * psph->P;
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
            riemann_input.L.p = psph->P;
            /* printf("WARNING, L.p < 0 : using first-order scheme \n");*/
        }
        if (riemann_input.R.p < 0) {
            riemann_input.R.p = q(P);
            /* printf("WARNING, R.p < 0 : using first-order scheme \n");*/
        }

#ifdef EEOS_POLYTROPE
        const double pLpoly =
            polytropicPressureFloor(a_inv3, riemann_input.L.rho, smf->dConstGamma, 
                  smf->dEOSPolyFloorIndex, smf->dEOSPolyFloorDen, smf->dEOSPolyFlooru);
        const double pRpoly =
            polytropicPressureFloor(a_inv3, riemann_input.R.rho, smf->dConstGamma,
                  smf->dEOSPolyFloorIndex, smf->dEOSPolyFloorDen, smf->dEOSPolyFlooru);
        riemann_input.L.p = MAX(riemann_input.L.p, pLpoly);
        riemann_input.R.p = MAX(riemann_input.R.p, pRpoly);
#endif
#ifdef EEOS_JEANS
        const double pLjeans =
            jeansPressureFloor(riemann_input.L.rho, ph, smf->dConstGamma, smf->dEOSNJeans);
        const double pRjeans =
            jeansPressureFloor(riemann_input.R.rho, q(ball), smf->dConstGamma, smf->dEOSNJeans);
        riemann_input.L.p = MAX(riemann_input.L.p, pLjeans);
        riemann_input.R.p = MAX(riemann_input.R.p, pRjeans);
#endif

        double cs_L = sqrt(GAMMA * riemann_input.L.p / riemann_input.L.rho);
        double cs_R = sqrt(GAMMA * riemann_input.R.p / riemann_input.R.rho);
        riemann_input.L.u  = riemann_input.L.p / (GAMMA_MINUS1 * riemann_input.L.rho);
        riemann_input.R.u  = riemann_input.R.p / (GAMMA_MINUS1 * riemann_input.R.rho);
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
                 face_unit, v_line_L, v_line_R, cs_L, cs_R, h_L, h_R);





#ifdef ENTROPY_SWITCH

#ifdef USE_MFM
        // As we are in a truly lagrangian configuration,
        // there is no advection of entropy among particles.
        double fluxes_S = 0.;
#else
        // riemann_output.Fluxes contains now the face state given by the riemann solver.
        // We only need that for computing the entropy flux, and then can be overwritten
        double fluxes_S = 0.;
        for (j=0; j<3; j++) fluxes_S += riemann_output.Fluxes.v[j]*face_unit[j];
        if (fluxes_S > 0) {
            // Maybe this values should be properly extrapolated to the faces..
            // but this is expensive!
            fluxes_S *= psph->S*pDensity/pkdMass(pkd,p);
        } else {
            fluxes_S *= q(S)*q(rho)/q(mass);
        }
        fluxes_S *= riemann_output.Fluxes.rho*modApq;
        fluxes_S = 0.;

#endif //USE_MFM
#endif //ENTROPY_SWITCH

#ifdef USE_MFM
        riemann_output.Fluxes.rho = 0.;
        riemann_output.Fluxes.p = riemann_output.P_M * riemann_output.S_M;
        for(j=0; j<3; j++)
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
        for(j=0; j<3; j++) {
            riemann_output.Fluxes.p += vFrame[j] * riemann_output.Fluxes.v[j];
            riemann_output.Fluxes.p += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
        }

        // Now we just multiply by the face area
        riemann_output.Fluxes.p *= modApq;
        riemann_output.Fluxes.rho *= modApq;
        for (j=0; j<3; j++) {
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


void hydroFluxFillBuffer(my_real **buffer, PARTICLE* q, int i, double dr2,
                                double dx, double dy, double dz, SMF *smf)
{
    PKD pkd = smf->pkd;
    double dDelta = smf->dDelta;
    float qh = pkdBall(pkd,q);
    SPHFIELDS* qsph = pkdSph(pkd,q);
    buffer[q_mass][i] = pkdMass(pkd,q);
    buffer[q_ball][i] = qh;
    buffer[q_dx][i] = dx;
    buffer[q_dy][i] = dy;
    buffer[q_dz][i] = dz;
    buffer[q_dr][i] = sqrt(dr2);
    buffer[q_rung][i] = dDelta/(1<<q->uRung);
    buffer[q_rho][i] = pkdDensity(pkd,q);
    buffer[q_P][i] = qsph->P;
#ifdef ENTROPY_SWITCH
    buffer[q_S][i] = qsph->S;
#endif
    buffer[q_vx][i] = qsph->vPred[0];
    buffer[q_vy][i] = qsph->vPred[1];
    buffer[q_vz][i] = qsph->vPred[2];

    buffer[q_gradRhoX][i] = qsph->gradRho[0];
    buffer[q_gradRhoY][i] = qsph->gradRho[1];
    buffer[q_gradRhoZ][i] = qsph->gradRho[2];

    buffer[q_gradPX][i] = qsph->gradP[0];
    buffer[q_gradPY][i] = qsph->gradP[1];
    buffer[q_gradPZ][i] = qsph->gradP[2];

    buffer[q_gradVxX][i] = qsph->gradVx[0];
    buffer[q_gradVxY][i] = qsph->gradVx[1];
    buffer[q_gradVxZ][i] = qsph->gradVx[2];

    buffer[q_gradVyX][i] = qsph->gradVy[0];
    buffer[q_gradVyY][i] = qsph->gradVy[1];
    buffer[q_gradVyZ][i] = qsph->gradVy[2];

    buffer[q_gradVzX][i] = qsph->gradVz[0];
    buffer[q_gradVzY][i] = qsph->gradVz[1];
    buffer[q_gradVzZ][i] = qsph->gradVz[2];

    buffer[q_lastUpdateTime][i] = qsph->lastUpdateTime;
    buffer[q_lastAccX][i] = qsph->lastAcc[0];
    buffer[q_lastAccY][i] = qsph->lastAcc[1];
    buffer[q_lastAccZ][i] = qsph->lastAcc[2];
    buffer[q_B_XX][i] = qsph->B[XX];
    buffer[q_B_YY][i] = qsph->B[YY];
    buffer[q_B_ZZ][i] = qsph->B[ZZ];
    buffer[q_B_XY][i] = qsph->B[XY];
    buffer[q_B_XZ][i] = qsph->B[XZ];
    buffer[q_B_YZ][i] = qsph->B[YZ];
    buffer[q_omega][i] = qsph->omega;
}


void hydroFluxUpdateFromBuffer(my_real **out_buffer, my_real **in_buffer,
                                      PARTICLE* p, PARTICLE* q, int i, SMF *smf)
{
    PKD pkd = smf->pkd;
    SPHFIELDS *psph = pkdSph(pkd,p);
    SPHFIELDS *qsph = pkdSph(pkd,q);
    double aFac = smf->a;
    double dDelta = smf->dDelta;
    float *qmass = (float*)pkdField(q,pkd->oFieldOffset[oMass]);
    float *pmass = (float*)pkdField(p,pkd->oFieldOffset[oMass]);
    if (dDelta>0) {
        *pmass -= out_buffer[out_minDt][i] * out_buffer[out_Frho][i] ;

        psph->mom[0] -= out_buffer[out_minDt][i] * out_buffer[out_FmomX][i] ;
        psph->mom[1] -= out_buffer[out_minDt][i] * out_buffer[out_FmomY][i] ;
        psph->mom[2] -= out_buffer[out_minDt][i] * out_buffer[out_FmomZ][i] ;

        psph->E -= out_buffer[out_minDt][i] * out_buffer[out_Fene][i];

        psph->Uint -= out_buffer[out_minDt][i] * ( out_buffer[out_Fene][i] - out_buffer[out_FmomX][i]*psph->vPred[0]
                      - out_buffer[out_FmomY][i]*psph->vPred[1]
                      - out_buffer[out_FmomZ][i]*psph->vPred[2]
                      + 0.5*(psph->vPred[0]*psph->vPred[0] + psph->vPred[1]*psph->vPred[1] + psph->vPred[2]*psph->vPred[2]) * out_buffer[out_Frho][i] );

#ifdef ENTROPY_SWITCH
        psph->S -= out_buffer[out_minDt][i] * out_buffer[out_FS][i];
#endif

#ifndef USE_MFM
        psph->drDotFrho[0] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dx][i] * aFac;
        psph->drDotFrho[1] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dy][i] * aFac;
        psph->drDotFrho[2] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dz][i] * aFac;
#endif
        psph->Frho +=      out_buffer[out_Frho][i] ;
        psph->Fene +=      out_buffer[out_Fene][i] ;
        psph->Fmom[0] +=   out_buffer[out_FmomX][i];
        psph->Fmom[1] +=   out_buffer[out_FmomY][i];
        psph->Fmom[2] +=   out_buffer[out_FmomZ][i];
    } else {
        psph->Frho +=      out_buffer[out_Frho][i] ;
        psph->Fene +=      out_buffer[out_Fene][i] ;
        psph->Fmom[0] +=   out_buffer[out_FmomX][i];
        psph->Fmom[1] +=   out_buffer[out_FmomY][i];
        psph->Fmom[2] +=   out_buffer[out_FmomZ][i];
    }

#ifndef OPTIM_NO_REDUNDANT_FLUXES
#ifdef OPTIM_AVOID_IS_ACTIVE
    if (!pkdIsActive(pkd,q))
#else
    if (!q->bMarked)
#endif
#endif
    {

        // If this is not the case, something VERY odd must have happened
        assert( qsph->P == in_buffer[q_P][i] );
        if (dDelta>0) {
            *qmass += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] ;

            qsph->mom[0] += out_buffer[out_minDt][i] * out_buffer[out_FmomX][i] ;
            qsph->mom[1] += out_buffer[out_minDt][i] * out_buffer[out_FmomY][i] ;
            qsph->mom[2] += out_buffer[out_minDt][i] * out_buffer[out_FmomZ][i] ;

            qsph->E += out_buffer[out_minDt][i] * out_buffer[out_Fene][i];

            qsph->Uint += out_buffer[out_minDt][i] * ( out_buffer[out_Fene][i] - out_buffer[out_FmomX][i]*qsph->vPred[0]
                          - out_buffer[out_FmomY][i]*qsph->vPred[1]
                          - out_buffer[out_FmomZ][i]*qsph->vPred[2]
                          + 0.5*(qsph->vPred[0]*qsph->vPred[0] + qsph->vPred[1]*qsph->vPred[1] + qsph->vPred[2]*qsph->vPred[2])*out_buffer[out_Frho][i]  );
#ifdef ENTROPY_SWITCH
            qsph->S += out_buffer[out_minDt][i] * out_buffer[out_FS][i];
#endif

#ifndef USE_MFM
            qsph->drDotFrho[0] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dx][i] ;
            qsph->drDotFrho[1] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dy][i] ;
            qsph->drDotFrho[2] += out_buffer[out_minDt][i] * out_buffer[out_Frho][i] * in_buffer[q_dz][i] ;
#endif
            qsph->Frho -=      out_buffer[out_Frho][i] ;
            qsph->Fene -=      out_buffer[out_Fene][i] ;
            qsph->Fmom[0] -=   out_buffer[out_FmomX][i];
            qsph->Fmom[1] -=   out_buffer[out_FmomY][i];
            qsph->Fmom[2] -=   out_buffer[out_FmomZ][i];
        } else {
            qsph->Frho -=       out_buffer[out_Frho][i];
            qsph->Fene -=       out_buffer[out_Fene][i];
            qsph->Fmom[0] -=    out_buffer[out_FmomX][i];
            qsph->Fmom[1] -=    out_buffer[out_FmomY][i];
            qsph->Fmom[2] -=    out_buffer[out_FmomZ][i];
        }

    } // q marked/active

}

void hydroFluxGetNvars(int *in, int *out)
{
    *in = q_last;
    *out = out_last;
}



#endif // OPTIM_FLUX_VEC


void combThirdHydroLoop(void *vpkd, void *v1,void *v2)
{
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = (PARTICLE *) v1;
    PARTICLE *p2 = (PARTICLE *) v2;

    assert(!pkd->bNoParticleOrder);
    if (pkdIsGas(pkd,p1) && pkdIsGas(pkd,p2)) {
        SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);
        int i;

        for (i=0; i<3; i++) {
            psph1->Fmom[i] += psph2->Fmom[i];
        }
        psph1->Frho += psph2->Frho;
        psph1->Fene += psph2->Fene;


        float *p1mass = (float *) pkdField(p1,pkd->oFieldOffset[oMass]);
        float *p2mass = (float *) pkdField(p2,pkd->oFieldOffset[oMass]);
        *p1mass += *p2mass;

        psph1->mom[0] += psph2->mom[0];
        psph1->mom[1] += psph2->mom[1];
        psph1->mom[2] += psph2->mom[2];
        psph1->E += psph2->E;
        psph1->Uint += psph2->Uint;
    }
}

