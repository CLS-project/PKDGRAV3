#include <algorithm>
#include "hydro/hydro.h"
#include "hydro/limiters.h"
using blitz::dot;
using blitz::TinyVector;

auto ServiceHydroGradients::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->pkd();
    auto helper = std::make_shared<SmoothHelperHydroGradients>(pkd, in->bConditionedBarthJespersen, in->bStricterSlopeLimiter);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
#ifdef OPTIM_SMOOTH_NODE
    *out = smx.ReSmoothNode<smooth::policy::OnlyGasOptimised, smooth::policy::OnlyGas>();
#else
    *out = smx.ReSmooth<smooth::policy::OnlyGasOptimised>();
#endif
    return sizeof(output);
}

void SmoothHelperHydroGradients::pack(void* dst, const void* src) {
    auto p1 = static_cast<pack_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->velocity = p2.velocity();
        p1->P = p2.sph().P;
        p1->fBall = p2.ball();
        p1->fDensity = p2.density();
        p1->bMarked = p2.marked();
    }
}

void SmoothHelperHydroGradients::unpack(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const pack_type*>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_adj_decorator();

        p1.set_position(p2->position);
        p1.velocity() = p2->velocity;
        p1.sph().P = p2->P;
        p1.set_ball(p2->fBall);
        p1.set_density(p2->fDensity);
        p1.set_marked(p2->bMarked);
    }
}

void SmoothHelperHydroGradients::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    double rho_max, rho_min;
    double vx_max, vx_min, vy_max, vy_min, vz_max, vz_min;
    double p_max, p_min;
    double max_face_dist;

    const auto& pv = p.velocity();
    auto& psph = p.sph();
    const double pH = p.ball();

#ifndef OPTIM_SMOOTH_NODE
    /* Compute the E matrix (Hopkins 2015, eq 14) */
    TinyVector<double, 6> E{0.0}; //  We are assumming 3D here!

    for (auto& nn : nnList) {

        const double rpq = sqrt(nn.fDist2);
        const auto& Hpq = pH;

        const double Wpq = cubicSplineKernel(rpq, Hpq);
        const auto& dr = nn.dr;
        E[+MOI::XX] += dr[0] * dr[0] * Wpq;
        E[+MOI::YY] += dr[1] * dr[1] * Wpq;
        E[+MOI::ZZ] += dr[2] * dr[2] * Wpq;

        E[+MOI::XY] += dr[1] * dr[0] * Wpq;
        E[+MOI::XZ] += dr[2] * dr[0] * Wpq;
        E[+MOI::YZ] += dr[1] * dr[2] * Wpq;
    }

    /* Normalize the matrix */
    E /= psph.omega;

    /* END of E matrix computation */

    /* Now, we need to do the inverse */
    inverseMatrix(E.data(), psph.B.data());
    psph.Ncond = conditionNumber(E.data(), psph.B.data());

    // DEBUG: check if E^{-1} = B
    /*
    double *B = psph.B.data();
    double unityXX, unityYY, unityZZ, unityXY, unityXZ, unityYZ;
    unityXX = E[+MOI::XX]*B[+MOI::XX] + E[+MOI::XY]*B[+MOI::XY] + E[+MOI::XZ]*B[+MOI::XZ];
    unityYY = E[+MOI::XY]*B[+MOI::XY] + E[+MOI::YY]*B[+MOI::YY] + E[+MOI::YZ]*B[+MOI::YZ];
    unityZZ = E[+MOI::XZ]*B[+MOI::XZ] + E[+MOI::YZ]*B[+MOI::YZ] + E[+MOI::ZZ]*B[+MOI::ZZ];
    unityXY = E[+MOI::XX]*B[+MOI::XY] + E[+MOI::XY]*B[+MOI::YY] + E[+MOI::XZ]*B[+MOI::YZ];
    unityXZ = E[+MOI::XX]*B[+MOI::XZ] + E[+MOI::XY]*B[+MOI::YZ] + E[+MOI::XZ]*B[+MOI::ZZ];
    unityYZ = E[+MOI::XY]*B[+MOI::XZ] + E[+MOI::YY]*B[+MOI::YZ] + E[+MOI::YZ]*B[+MOI::ZZ];

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

    /* We also use this loop to get the primitive variables' local extrema */
    rho_min = rho_max = p.density();
    vx_min = vx_max = pv[0];
    vy_min = vy_max = pv[1];
    vz_min = vz_max = pv[2];
    p_min = p_max = psph.P;
    max_face_dist = 0.;

    for (auto& nn : nnList) {
        if (&p == nn.pPart) continue;
        auto q = pkd->particles[nn.pPart];
        const auto& qv = q.velocity();
        auto& qsph = q.sph();

        /* Vector from p to q; i.e., rq - rp */
        const TinyVector<double, 3> dr{-nn.dr};

        const double rpq = sqrt(nn.fDist2);
        const auto& Hpq = pH;
        const double Wpq = cubicSplineKernel(rpq, Hpq);
        const double psi = Wpq / psph.omega;

        TinyVector<double, 3> psiTilde_p;
        psiTilde_p[0] = dot(dr, TinyVector<double, 3>{psph.B[+MOI::XX], psph.B[+MOI::XY], psph.B[+MOI::XZ]}) * psi;
        psiTilde_p[1] = dot(dr, TinyVector<double, 3>{psph.B[+MOI::XY], psph.B[+MOI::YY], psph.B[+MOI::YZ]}) * psi;
        psiTilde_p[2] = dot(dr, TinyVector<double, 3>{psph.B[+MOI::XZ], psph.B[+MOI::YZ], psph.B[+MOI::ZZ]}) * psi;

        psph.gradRho += psiTilde_p * (q.density() - p.density());
        psph.gradVx += psiTilde_p * (qv[0] - pv[0]);
        psph.gradVy += psiTilde_p * (qv[1] - pv[1]);
        psph.gradVz += psiTilde_p * (qv[2] - pv[2]);
        psph.gradP += psiTilde_p * (qsph.P - psph.P);

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

        max_face_dist = std::max(max_face_dist, rpq);
    }

    /* Now we can limit the gradients */
    double beta = 1.;
    if (bConditionedBarthJespersen) {
        beta = std::max(1., 2. * std::min(1., 10. / psph.Ncond));
    }

    if (bStricterSlopeLimiter) {
        max_face_dist *= 0.5;
        double var_max_extrap, c_fac;

        var_max_extrap = sqrt(dot(psph.gradRho, psph.gradRho)) * max_face_dist;
        psph.gradRho *= BarthJespersenLimiter(rho_max - p.density(), p.density() - rho_min, var_max_extrap, beta);

        var_max_extrap = sqrt(dot(psph.gradVx, psph.gradVx)) * max_face_dist;
        psph.gradVx *= BarthJespersenLimiter(vx_max - pv[0], pv[0] - vx_min, var_max_extrap, beta);
        c_fac = psph.c / var_max_extrap;
        if (c_fac < 1.) psph.gradVx *= c_fac;

        var_max_extrap = sqrt(dot(psph.gradVy, psph.gradVy)) * max_face_dist;
        psph.gradVy *= BarthJespersenLimiter(vy_max - pv[1], pv[1] - vy_min, var_max_extrap, beta);
        c_fac = psph.c / var_max_extrap;
        if (c_fac < 1.) psph.gradVy *= c_fac;

        var_max_extrap = sqrt(dot(psph.gradVz, psph.gradVz)) * max_face_dist;
        psph.gradVz *= BarthJespersenLimiter(vz_max - pv[2], pv[2] - vz_min, var_max_extrap, beta);
        c_fac = psph.c / var_max_extrap;
        if (c_fac < 1.) psph.gradVz *= c_fac;

        var_max_extrap = sqrt(dot(psph.gradP, psph.gradP)) * max_face_dist;
        psph.gradP *= BarthJespersenLimiter(p_max - psph.P, psph.P - p_min, var_max_extrap, beta);
    }
    else {
        /* We need an extra loop to compute the maximum and minimum extrapolation
         * to this cell's interfaces */
        double rho_max_extrap, rho_min_extrap, p_max_extrap, p_min_extrap;
        double vx_max_extrap, vx_min_extrap, vy_max_extrap, vy_min_extrap, vz_max_extrap, vz_min_extrap;
        rho_max_extrap = vx_max_extrap = vy_max_extrap = vz_max_extrap = p_max_extrap = -FLT_MAX;
        rho_min_extrap = vx_min_extrap = vy_min_extrap = vz_min_extrap = p_min_extrap = FLT_MAX;

        for (auto& nn : nnList) {
            /* Vector from p to the pq face */
            const TinyVector<double, 3> dr{-0.5 * nn.dr};

            const auto rho_q_extrap = dot(dr, psph.gradRho);
            rho_min_extrap = std::min(rho_min_extrap, rho_q_extrap);
            rho_max_extrap = std::max(rho_max_extrap, rho_q_extrap);

            const auto vx_q_extrap = dot(dr, psph.gradVx);
            vx_min_extrap = std::min(vx_min_extrap, vx_q_extrap);
            vx_max_extrap = std::max(vx_max_extrap, vx_q_extrap);

            const auto vy_q_extrap = dot(dr, psph.gradVy);
            vy_min_extrap = std::min(vy_min_extrap, vy_q_extrap);
            vy_max_extrap = std::max(vy_max_extrap, vy_q_extrap);

            const auto vz_q_extrap = dot(dr, psph.gradVz);
            vz_min_extrap = std::min(vz_min_extrap, vz_q_extrap);
            vz_max_extrap = std::max(vz_max_extrap, vz_q_extrap);

            const auto p_q_extrap = dot(dr, psph.gradP);
            p_min_extrap = std::min(p_min_extrap, p_q_extrap);
            p_max_extrap = std::max(p_max_extrap, p_q_extrap);
        }

        psph.gradRho *= BarthJespersenLimiter(rho_max - p.density(), p.density() - rho_min, rho_max_extrap, -rho_min_extrap, beta);
        psph.gradVx *= BarthJespersenLimiter(vx_max - pv[0], pv[0] - vx_min, vx_max_extrap, -vx_min_extrap, beta);
        psph.gradVy *= BarthJespersenLimiter(vy_max - pv[1], pv[1] - vy_min, vy_max_extrap, -vy_min_extrap, beta);
        psph.gradVz *= BarthJespersenLimiter(vz_max - pv[2], pv[2] - vz_min, vz_max_extrap, -vz_min_extrap, beta);
        psph.gradP *= BarthJespersenLimiter(p_max - psph.P, psph.P - p_min, p_max_extrap, -p_min_extrap, beta);

        double c_fac;
        c_fac = psph.c / std::max(std::abs(vx_min_extrap), std::abs(vx_max_extrap));
        if (c_fac < 1.) psph.gradVx *= c_fac;

        c_fac = psph.c / std::max(std::abs(vy_min_extrap), std::abs(vy_max_extrap));
        if (c_fac < 1.) psph.gradVy *= c_fac;

        c_fac = psph.c / std::max(std::abs(vz_min_extrap), std::abs(vz_max_extrap));
        if (c_fac < 1.) psph.gradVz *= c_fac;
    }
    /* END OF LIMITER */
}
