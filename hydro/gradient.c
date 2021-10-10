#include "hydro/hydro.h"


void msrMeshlessGradients(MSR msr,double dTime)
{
    double sec, dsec;
    printf("Computing gradients... ");

    msrTimerStart(msr, TIMER_GRADIENTS);
    if (msr->param.bConservativeReSmooth) {
#ifdef OPTIM_SMOOTH_NODE
#ifdef OPTIM_AVOID_IS_ACTIVE
        msrSelActive(msr);
#endif
        msrReSmoothNode(msr,dTime,SMX_HYDRO_GRADIENT,0,0);
#else
        msrReSmooth(msr,dTime,SMX_HYDRO_GRADIENT,0,0);
#endif
    } else {
        msrSmooth(msr,dTime,SMX_HYDRO_GRADIENT,0, msr->param.nSmooth);
    }

    msrTimerStop(msr, TIMER_GRADIENTS);
    dsec = msrTimerGet(msr, TIMER_GRADIENTS);
    printf("took %.5f seconds\n", dsec);
}



void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf)
{

    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double E[6];  //  We are assumming 3D here!
    double ph, rpq, hpq,Wpq, dx,dy,dz, diff;
    double rho_max, rho_min;
    double vx_max, vx_min, vy_max, vy_min, vz_max, vz_min;
    double p_max, p_min;
    double limRho, limVx, limVy, limVz, limP;
    double psi, psiTilde_p[3];
    float  *pv, *qv;
    int i, j;
    psph = pkdSph(pkd,p);
    ph = fBall;

    /* Compute the E matrix (Hopkins 2015, eq 14) */
    for (i=0; i<6; ++i) {
        E[i] = 0.0;
    }

    for (i=0; i<nSmooth; ++i) {
        q = nnList[i].pPart;

        dx = -nnList[i].dx;
        dy = -nnList[i].dy;
        dz = -nnList[i].dz;


        rpq = sqrt(nnList[i].fDist2);
        hpq = ph;

        Wpq = cubicSplineKernel(rpq, hpq);

        E[XX] += dx*dx*Wpq;
        E[YY] += dy*dy*Wpq;
        E[ZZ] += dz*dz*Wpq;

        E[XY] += dy*dx*Wpq;
        E[XZ] += dz*dx*Wpq;
        E[YZ] += dy*dz*Wpq;
    }

    /* Normalize the matrix */
    for (i=0; i<6; ++i) {
        E[i] /= psph->omega;
    }




    /* END of E matrix computation */
    //printf("nSmooth %d fBall %e E_q [XX] %e \t [XY] %e \t [XZ] %e \n
    //                       \t \t \t [YY] %e \t [YZ] %e \n
    //                \t\t\t \t \t \t [ZZ] %e \n",
    // nSmooth, fBall, E[XX], E[XY], E[XZ], E[YY], E[YZ], E[ZZ]);

    /* Now, we need to do the inverse */
    inverseMatrix(E, psph->B);
    psph->Ncond = conditionNumber(E, psph->B);

    // DEBUG: check if E^{-1} = B
    /*
    double *B = psph->B;
    double unityXX, unityYY, unityZZ, unityXY, unityXZ, unityYZ;
    unityXX = E[XX]*B[XX] + E[XY]*B[XY] + E[XZ]*B[XZ];
    unityYY = E[XY]*B[XY] + E[YY]*B[YY] + E[YZ]*B[YZ];
    unityZZ = E[XZ]*B[XZ] + E[YZ]*B[YZ] + E[ZZ]*B[ZZ];

    unityXY = E[XX]*B[XY] + E[XY]*B[YY] + E[XZ]*B[YZ];
    unityXZ = E[XX]*B[XZ] + E[XY]*B[YZ] + E[XZ]*B[ZZ];
    unityYZ = E[XY]*B[XZ] + E[YY]*B[YZ] + E[YZ]*B[ZZ];

    printf("XX %e \t YY %e \t ZZ %e \n", unityXX, unityYY, unityZZ);
    printf("XY %e \t XZ %e \t YZ %e \n", unityXY, unityXZ, unityYZ);
    */

    /* Now we can compute the gradients
     * This and the B matrix computation could be done in the first hydro loop
     * (where omega is computed) but at that time we do not know the densities
     * of *all* particles (because they depend on omega)
     */

    for (j=0; j<3; j++) {
        psph->gradRho[j] = 0.0;
        psph->gradVx[j] = 0.0;
        psph->gradVy[j] = 0.0;
        psph->gradVz[j] = 0.0;
        psph->gradP[j] = 0.0;
#if defined(MAKE_GLASS) || defined(REGULARIZE_MESH)
        psph->cellCM[j] = 0.0;
#endif
    }
    for (i=0; i<nSmooth; ++i) {

        q = nnList[i].pPart;
        qsph = pkdSph(pkd, q);

        dx = -nnList[i].dx;
        dy = -nnList[i].dy;
        dz = -nnList[i].dz;

        if (dx==0 && dy==0 && dz==0) continue;

        ph = fBall;
        rpq = sqrt(nnList[i].fDist2);
        hpq = ph;

        Wpq = cubicSplineKernel(rpq, hpq);
        psi = Wpq/psph->omega;

        psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
        psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
        psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

        diff = pkdDensity(pkd, q) - pkdDensity(pkd, p);
        for (j=0; j<3; j++) {
            psph->gradRho[j] += diff*psiTilde_p[j];
        }
        diff = qsph->vPred[0] - psph->vPred[0];
        for (j=0; j<3; j++) {
            psph->gradVx[j] += diff*psiTilde_p[j];
        }
        diff = qsph->vPred[1] - psph->vPred[1];
        for (j=0; j<3; j++) {
            psph->gradVy[j] += diff*psiTilde_p[j];
        }
        diff = qsph->vPred[2] - psph->vPred[2];
        for (j=0; j<3; j++) {
            psph->gradVz[j] += diff*psiTilde_p[j];
        }
        diff = qsph->P - psph->P;
        for (j=0; j<3; j++) {
            psph->gradP[j] += diff*psiTilde_p[j];
        }
#if defined(MAKE_GLASS) || defined(REGULARIZE_MESH)
        for (j=0; j<3; j++) {
            psph->cellCM[j] += psi*psiTilde_p[j];
        }
#endif
    }
#if defined(MAKE_GLASS) || defined(REGULARIZE_MESH)
    double CMfactor = - fBall*fBall*M_PI*fBall*fBall*fBall * psph->omega / 3.;
    for (j=0; j<3; j++) {
        psph->cellCM[j] *= CMfactor;
    }
#endif

    /* Now we can limit them */

    // First step, compute the maximum and minimum difference of each variable
    rho_min= vx_min= vy_min= vz_min= p_min =  HUGE_VAL;
    rho_max= vx_max= vy_max= vz_max= p_max = -HUGE_VAL;


    pv = pkdVel(pkd,p);
    for (i=0; i<nSmooth; ++i) {
//       if (nnList[i].dx==0 && nnList[i].dy==0 && nnList[i].dz==0) continue;
        q = nnList[i].pPart;
        qsph = pkdSph(pkd, q);
        qv = pkdVel(pkd,q);

        if (pkdDensity(pkd,q) < rho_min) rho_min = pkdDensity(pkd,q);
        if (pkdDensity(pkd,q) > rho_max) rho_max = pkdDensity(pkd,q);

        if (qv[0] < vx_min) vx_min = qv[0];
        if (qv[0] > vx_max) vx_max = qv[0];

        if (qv[1] < vy_min) vy_min = qv[1];
        if (qv[1] > vy_max) vy_max = qv[1];

        if (qv[2] < vz_min) vz_min = qv[2];
        if (qv[2] > vz_max) vz_max = qv[2];

        if (qsph->P < p_min) p_min = qsph->P;
        if (qsph->P > p_max) p_max = qsph->P;
    }


    limRho= limVx= limVy= limVz= limP = 1.;
#if defined(LIMITER_BARTH) || defined(LIMITER_CONDBARTH)
    for (i=0; i<nSmooth; ++i) {
        dx = -nnList[i].dx; //Vector from p to q
        dy = -nnList[i].dy;
        dz = -nnList[i].dz;
        if (dx==0 && dy==0 && dz==0) continue;
        q = nnList[i].pPart;
        qsph = pkdSph(pkd, q);

        // TODO: The differences could be computed outside of this loop
#ifdef LIMITER_BARTH
        BarthJespersenLimiter(&limRho, psph->gradRho, rho_max-pkdDensity(pkd,p), rho_min-pkdDensity(pkd,p), dx, dy, dz);
        BarthJespersenLimiter(&limVx, psph->gradVx, vx_max-pv[0], vx_min-pv[0], dx, dy, dz);
        BarthJespersenLimiter(&limVy, psph->gradVy, vy_max-pv[1], vy_min-pv[1], dx, dy, dz);
        BarthJespersenLimiter(&limVz, psph->gradVz, vz_max-pv[2], vz_min-pv[2], dx, dy, dz);
        BarthJespersenLimiter(&limP, psph->gradP, p_max-psph->P, p_min-psph->P, dx, dy, dz);
#endif

#ifdef LIMITER_CONDBARTH
        ConditionedBarthJespersenLimiter(&limRho, psph->gradRho, rho_max-pkdDensity(pkd,p), rho_min-pkdDensity(pkd,p), dx, dy, dz, 10., psph->Ncond);
        ConditionedBarthJespersenLimiter(&limVx, psph->gradVx, vx_max-pv[0], vx_min-pv[0], dx, dy, dz, 10., psph->Ncond);
        ConditionedBarthJespersenLimiter(&limVy, psph->gradVy, vy_max-pv[1], vy_min-pv[1], dx, dy, dz, 10., psph->Ncond);
        ConditionedBarthJespersenLimiter(&limVz, psph->gradVz, vz_max-pv[2], vz_min-pv[2], dx, dy, dz, 10., psph->Ncond);
        ConditionedBarthJespersenLimiter(&limP, psph->gradP, p_max-psph->P, p_min-psph->P, dx, dy, dz, 10., psph->Ncond);
#endif
    }
#endif

    for (j=0; j<3; j++) {
        psph->gradRho[j] *= limRho;
        psph->gradVx[j] *= limVx;
        psph->gradVy[j] *= limVy;
        psph->gradVz[j] *= limVz;
        psph->gradP[j] *= limP;
    }
    /* END OF LIMITER */
}



