#include <algorithm>
#include "hydro/hydro.h"
#include "master.h"
using blitz::TinyVector;
using blitz::dot;

void MSR::MeshlessGradients(double dTime, double dDelta) {
    double dsec;
    printf("Computing gradients... ");

    TimerStart(TIMER_GRADIENTS);
#ifdef OPTIM_SMOOTH_NODE
#ifdef OPTIM_AVOID_IS_ACTIVE
    SelActives();
#endif
    ReSmoothNode(dTime,dDelta,SMX_HYDRO_GRADIENT,0);
#else
    ReSmooth(dTime,dDelta,SMX_HYDRO_GRADIENT,0);
#endif

    TimerStop(TIMER_GRADIENTS);
    dsec = TimerGet(TIMER_GRADIENTS);
    printf("took %.5f seconds\n", dsec);
}



void hydroGradients(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {

    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    double rho_max, rho_min;
    double vx_max, vx_min, vy_max, vy_min, vz_max, vz_min;
    double p_max, p_min;
    double limRho, limVx, limVy, limVz, limP;

    auto &psph = p.sph();
    const double ph = 0.5*fBall;

#ifndef OPTIM_SMOOTH_NODE
    /* Compute the E matrix (Hopkins 2015, eq 14) */
    TinyVector<double,6> E{0.0};  //  We are assumming 3D here!

    for (auto i = 0; i < nSmooth; ++i) {

        const double rpq = sqrt(nnList[i].fDist2);
        const auto &hpq = ph;

        const double Wpq = cubicSplineKernel(rpq, hpq);
        const auto &dr = nnList[i].dr;

        E[XX] += dr[0]*dr[0]*Wpq;
        E[YY] += dr[1]*dr[1]*Wpq;
        E[ZZ] += dr[2]*dr[2]*Wpq;

        E[XY] += dr[1]*dr[0]*Wpq;
        E[XZ] += dr[2]*dr[0]*Wpq;
        E[YZ] += dr[1]*dr[2]*Wpq;
    }

    /* Normalize the matrix */
    E /= psph.omega;




    /* END of E matrix computation */
    //printf("nSmooth %d fBall %e E_q [XX] %e \t [XY] %e \t [XZ] %e \n
    //                       \t \t \t [YY] %e \t [YZ] %e \n
    //                \t\t\t \t \t \t [ZZ] %e \n",
    // nSmooth, fBall, E[XX], E[XY], E[XZ], E[YY], E[YZ], E[ZZ]);

    /* Now, we need to do the inverse */
    inverseMatrix(E.data(), psph.B.data());
    psph.Ncond = conditionNumber(E.data(), psph.B.data());

    // DEBUG: check if E^{-1} = B
    /*
    double *B = psph.B.data();
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
#endif

    /* Now we can compute the gradients
     * This and the B matrix computation could be done in the first hydro loop
     * (where omega is computed) but at that time we do not know the densities
     * of *all* particles (because they depend on omega)
     */

    psph.gradRho = 0.0;
    psph.gradVx = 0.0;
    psph.gradVy = 0.0;
    psph.gradVz = 0.0;
    psph.gradP = 0.0;
    for (auto i = 0; i < nSmooth; ++i) {
        if (pIn == nnList[i].pPart) continue;
        auto q = pkd->particles[nnList[i].pPart];
        auto &qsph = q.sph();

        const TinyVector<double,3> dr{-nnList[i].dr};

        const double rpq = sqrt(nnList[i].fDist2);
        const auto &hpq = ph;
        const double Wpq = cubicSplineKernel(rpq, hpq);
        const double psi = Wpq/psph.omega;

        TinyVector<double,3> psiTilde_p;
        psiTilde_p[0] = dot(dr, TinyVector<double,3> {psph.B[XX],psph.B[XY],psph.B[XZ]}) * psi;
        psiTilde_p[1] = dot(dr, TinyVector<double,3> {psph.B[XY],psph.B[YY],psph.B[YZ]}) * psi;
        psiTilde_p[2] = dot(dr, TinyVector<double,3> {psph.B[XZ],psph.B[YZ],psph.B[ZZ]}) * psi;

        psph.gradRho += psiTilde_p * (q.density() - p.density());
        psph.gradVx  += psiTilde_p * (qsph.vPred[0] - psph.vPred[0]);
        psph.gradVy  += psiTilde_p * (qsph.vPred[1] - psph.vPred[1]);
        psph.gradVz  += psiTilde_p * (qsph.vPred[2] - psph.vPred[2]);
        psph.gradP   += psiTilde_p * (qsph.P - psph.P);
    }

    /* Now we can limit them */

    // First step, compute the maximum and minimum difference of each variable
    rho_min= vx_min= vy_min= vz_min= p_min =  HUGE_VAL;
    rho_max= vx_max= vy_max= vz_max= p_max = -HUGE_VAL;


    const auto &pv = p.velocity();
    for (auto i = 0; i < nSmooth; ++i) {
//        if (pIn == nnList[i].pPart) continue;
        auto q = pkd->particles[nnList[i].pPart];
        auto &qsph = q.sph();
        const auto &qv = q.velocity();

        rho_min = std::min(rho_min, static_cast<double>(q.density()));
        rho_max = std::max(rho_max, static_cast<double>(q.density()));

        if (qv[0] < vx_min) vx_min = qv[0];
        if (qv[0] > vx_max) vx_max = qv[0];

        if (qv[1] < vy_min) vy_min = qv[1];
        if (qv[1] > vy_max) vy_max = qv[1];

        if (qv[2] < vz_min) vz_min = qv[2];
        if (qv[2] > vz_max) vz_max = qv[2];

        p_min = std::min(p_min, qsph.P);
        p_max = std::max(p_max, qsph.P);
    }


    limRho= limVx= limVy= limVz= limP = 1.;
#if defined(LIMITER_BARTH) || defined(LIMITER_CONDBARTH)
    for (auto i = 0; i < nSmooth; ++i) {
        if (pIn == nnList[i].pPart) continue;
        const TinyVector<double,3> dr{-nnList[i].dr}; //Vector from p to q

        // TODO: The differences could be computed outside of this loop
#ifdef LIMITER_BARTH
        BarthJespersenLimiter(&limRho, psph.gradRho, rho_max-p.density(), rho_min-p.density(), dr);
        BarthJespersenLimiter(&limVx, psph.gradVx, vx_max-pv[0], vx_min-pv[0], dr);
        BarthJespersenLimiter(&limVy, psph.gradVy, vy_max-pv[1], vy_min-pv[1], dr);
        BarthJespersenLimiter(&limVz, psph.gradVz, vz_max-pv[2], vz_min-pv[2], dr);
        BarthJespersenLimiter(&limP, psph.gradP, p_max-psph.P, p_min-psph.P, dr);
#endif

#ifdef LIMITER_CONDBARTH
        ConditionedBarthJespersenLimiter(&limRho, psph.gradRho, rho_max-p.density(), rho_min-p.density(), dr, 10., psph.Ncond);
        ConditionedBarthJespersenLimiter(&limVx, psph.gradVx, vx_max-pv[0], vx_min-pv[0], dr, 10., psph.Ncond);
        ConditionedBarthJespersenLimiter(&limVy, psph.gradVy, vy_max-pv[1], vy_min-pv[1], dr, 10., psph.Ncond);
        ConditionedBarthJespersenLimiter(&limVz, psph.gradVz, vz_max-pv[2], vz_min-pv[2], dr, 10., psph.Ncond);
        ConditionedBarthJespersenLimiter(&limP, psph.gradP, p_max-psph.P, p_min-psph.P, dr, 10., psph.Ncond);
#endif
    }
#endif

    psph.gradRho *= limRho;
    psph.gradVx *= limVx;
    psph.gradVy *= limVy;
    psph.gradVz *= limVz;
    psph.gradP *= limP;
    /* END OF LIMITER */
}

