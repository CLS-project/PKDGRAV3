#include "hydro/hydro.h"
#include "master.h"
#include "eEOS/eEOS.h"
#include "hydro/limiters.h"
#include "riemann.h"
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

template <typename dtype=dvec, typename mtype=dmask>
class MeshlessHydroSolver {

private:

    inline void extrapolateDensityInTime(dtype &rho, dtype rho0, dtype vx, dtype vy, dtype vz,
                                         dtype dt, dtype gradRhoX, dtype gradRhoY, dtype gradRhoZ, dtype divv) {
        rho -= dt*(vx*gradRhoX + vy*gradRhoY + vz*gradRhoZ + rho0*divv);
    }

    inline void extrapolateVelocityInTime(dtype &v, dtype &vFrame, dtype divv, dtype v0, dtype dt, dtype acc,
                                          dtype gradP, dtype rho0,
                                          dtype a) {
        dtype temp;
        temp = -(v0*divv + gradP/rho0)*dt;
        temp += acc*dt*a;
        v += temp;
        vFrame += 0.5*temp;
    }

    inline void extrapolatePressureInTime(dtype &p, dtype dt, dtype gradPx, dtype gradPy, dtype gradPz,
                                          dtype p0, dtype vx0, dtype vy0, dtype vz0, dtype divv, dtype dConstGamma) {
        p -= dt*( vx0*gradPx + vy0*gradPy + vz0*gradPz  + dConstGamma*p0*divv);
    }

    inline void extrapolateStateInTime(
        dtype &rho, dtype &vx, dtype &vy, dtype &vz, dtype &p,
        dtype &vFx, dtype &vFy, dtype &vFz,
        dtype rho0, dtype vx0, dtype vy0, dtype vz0, dtype p0,
        dtype dt,
        dtype gradRhoX, dtype gradRhoY, dtype gradRhoZ,
        dtype gradPX, dtype gradPY, dtype gradPZ,
        dtype accx, dtype accy, dtype accz, dtype divv,
        dtype dConstGamma, dtype a) {

        extrapolateDensityInTime( rho, rho0, vx, vy, vz, dt,
                                  gradRhoX, gradRhoY, gradRhoZ, divv);

        extrapolateVelocityInTime( vx, vFx, divv, vx0, dt, accx, gradPX, rho0, a);
        extrapolateVelocityInTime( vy, vFy, divv, vy0, dt, accy, gradPY, rho0, a);
        extrapolateVelocityInTime( vz, vFz, divv, vz0, dt, accz, gradPZ, rho0, a);

        extrapolatePressureInTime( p, dt, gradPX, gradPY, gradPZ,
                                   p0, vx0, vy0, vz0, divv, dConstGamma);
    }

    inline void extrapolateVariableInSpace(dtype &var, dtype dx,dtype dy, dtype dz,
                                           dtype gradx, dtype grady, dtype gradz) {
        var += dx*gradx + dy*grady + dz*gradz;
    }

    inline void extrapolateVelocityCosmology(dtype &v, dtype &vFrame, dtype v0, dtype dt, dtype H, dtype a) {
        dtype temp = H * dt * a * v;
        v -= temp;
        vFrame -= 0.5*temp;
    }

    inline void extrapolateCosmology(dtype &vx, dtype &vy, dtype &vz,
                                     dtype &vFramex, dtype &vFramey, dtype &vFramez, dtype &p,
                                     dtype dt, dtype vx0, dtype vy0, dtype vz0, dtype p0,
                                     dtype dConstGamma, dtype H, dtype a) {
        extrapolateVelocityCosmology( vx,  vFramex,  vx0,  dt,  H,  a);
        extrapolateVelocityCosmology( vy,  vFramey,  vy0,  dt,  H,  a);
        extrapolateVelocityCosmology( vz,  vFramez,  vz0,  dt,  H,  a);

        p -= 3. * H * dt * a * (dConstGamma - 1.) * p0;
    }

    inline void computeFace(dtype &modApq, std::array<dtype,3> &unit,
                            dtype rpq,  dtype dx, dtype dy, dtype dz,
                            dtype ph, dtype qh, dtype p_omega, dtype q_omega,
                            dtype pBxx, dtype pBxy, dtype pBxz,
                            dtype pByy, dtype pByz, dtype pBzz,
                            dtype qBxx, dtype qBxy, dtype qBxz,
                            dtype qByy, dtype qByz, dtype qBzz) {
        dtype psi;
        dtype psiTilde[3];
        dtype Apq[3] = {0.,0.,0.};

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
        dtype zero = 0.;
        mtype cond = modApq>zero;
        for (auto j=0; j<3; j++) {
            unit[j] = Apq[j]/modApq;
        }
        if (testz(cond)) { // Some of the elements are zero
            for (auto j=0; j<3; j++)
                unit[j] = 0.0;
        }
    }

    inline void low_limit(dtype &var, dtype def, dtype min) {
        mtype cond = var<min;
        var = mask_mov(var, cond, def);
    }

    inline void doSinglePPFlux(dtype &F_rho, std::array<dtype,3> &F_v, dtype &F_p, dtype &F_S, dtype &minDt,
                               bool bComove, dtype dTime, dtype dDelta, dtype a, dtype H, dtype dConstGamma,
                               dtype rpq, dtype dx, dtype dy, dtype dz,
                               dtype pBall, dtype pLastUpdateTime, dtype pDt,
                               dtype pOmega,
                               dtype pBxx, dtype pBxy, dtype pBxz,
                               dtype pByy, dtype pByz, dtype pBzz,
                               dtype pDensity, dtype pVpredx, dtype pVpredy, dtype pVpredz, dtype pP, dtype pS,
                               dtype pGradRhoX, dtype pGradRhoY, dtype pGradRhoZ,
                               dtype pGradPX, dtype pGradPY, dtype pGradPZ,
                               dtype pGradVxX, dtype pGradVxY, dtype pGradVxZ,
                               dtype pGradVyX, dtype pGradVyY, dtype pGradVyZ,
                               dtype pGradVzX, dtype pGradVzY, dtype pGradVzZ,
                               dtype pLastAccX, dtype pLastAccY, dtype pLastAccZ,
                               dtype qBall, dtype qLastUpdateTime, dtype qDt,
                               dtype qOmega,
                               dtype qBxx, dtype qBxy, dtype qBxz,
                               dtype qByy, dtype qByz, dtype qBzz,
                               dtype qDensity, dtype qVpredx, dtype qVpredy, dtype qVpredz, dtype qP, dtype qS,
                               dtype qGradRhoX, dtype qGradRhoY, dtype qGradRhoZ,
                               dtype qGradPX, dtype qGradPY, dtype qGradPZ,
                               dtype qGradVxX, dtype qGradVxY, dtype qGradVxZ,
                               dtype qGradVyX, dtype qGradVyY, dtype qGradVyZ,
                               dtype qGradVzX, dtype qGradVzY, dtype qGradVzZ,
                               dtype qLastAccX, dtype qLastAccY, dtype qLastAccZ,
                               struct eEOSparam &eEOS) {

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

        dtype qDeltaHalf=0.0, pDeltaHalf=0.0;
        pDeltaHalf = (dTime - pLastUpdateTime + 0.5*pDt)/a;
        qDeltaHalf = (dTime - qLastUpdateTime + 0.5*qDt)/a;

        // DEBUG: Avoid temporal extrapolation
        //pDeltaHalf = 0.;
        //qDeltaHalf = 0.;

        dtype modApq;
        std::array<dtype, 3> face_unit;
        computeFace(modApq, face_unit,
                    rpq,   dx,  dy,  dz,
                    pBall,  qBall,  pOmega,  qOmega,
                    pBxx,  pBxy,  pBxz,
                    pByy,  pByz,  pBzz,
                    qBxx, qBxy, qBxz,
                    qByy, qByz, qBzz);

        // Velocity of the quadrature mid-point
        dtype vFrame[3];
        vFrame[0] = 0.5*(pVpredx+qVpredx);
        vFrame[1] = 0.5*(pVpredy+qVpredy);
        vFrame[2] = 0.5*(pVpredz+qVpredz);

        dtype pv[3], qv[3];
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
        dtype pdivv = (pGradVxX + pGradVyY + pGradVzZ);
        dtype qdivv = (qGradVxX + qGradVyY + qGradVzZ);

        dtype L_v[3], R_v[3];
        dtype L_rho = pDensity;
        dtype R_rho = qDensity;
        dtype L_p = pP;
        dtype R_p = qP;
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

        genericPairwiseLimiter(pDensity, qDensity, L_rho, R_rho);
        genericPairwiseLimiter(pP, qP, L_p, R_p);
        genericPairwiseLimiter(pv[0], qv[0], L_v[0], R_v[0]);
        genericPairwiseLimiter(pv[1], qv[1], L_v[1], R_v[1]);
        genericPairwiseLimiter(pv[2], qv[2], L_v[2], R_v[2]);

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

        dtype zero = 0.0;
        low_limit(L_rho, pDensity, zero);
        low_limit(R_rho, qDensity, zero);
        low_limit(L_p, pP, zero);
        low_limit(R_p, qP, zero);

#ifdef EEOS_POLYTROPE
        dtype a_inv3 = 1.0/(a*a*a);
        const dtype pLpoly =
            polytropicPressureFloor<dtype,mtype>(a_inv3, L_rho, dConstGamma,eEOS);
        const dtype pRpoly =
            polytropicPressureFloor<dtype,mtype>(a_inv3, R_rho, dConstGamma,eEOS);
        L_p = max(L_p, pLpoly);
        R_p = max(R_p, pRpoly);
#endif
        /*
        #ifdef EEOS_JEANS
        const double pLjeans =
            jeansPressureFloor(L_rho, ph, smf->dConstGamma, smf->dEOSNJeans);
        const double pRjeans =
            jeansPressureFloor(R_rho, q(ball), smf->dConstGamma, smf->dEOSNJeans);
        L_p = max(L_p, pLjeans);
        R_p = max(R_p, pRjeans);
        #endif
        */

        dtype P_M, S_M;

        RiemannSolverExact<dtype,mtype> riemann(dConstGamma);
        int niter = riemann.solve(
                        R_rho, R_p, R_v,
                        L_rho, L_p, L_v,
                        P_M, S_M,
                        &F_rho, &F_p, F_v.data(),
                        face_unit.data());

        /*
        int nan;
        // Only works if compiling with -fno-finite-math-only !!
        nan = nan_guard(S_M, zero);
        nan = nan_guard(P_M, zero);
        if (nan){
           printf("-----\n");
           dump(S_M);
           dump(P_M);
           dump(R_rho);
           dump(L_rho);
           dump(R_p);
           dump(L_p);
           dump(R_v[0]);
           dump(L_v[0]);
           dump(R_v[1]);
           dump(L_v[1]);
           dump(R_v[2]);
           dump(L_v[2]);
           dump(face_unit[0]);
           dump(face_unit[1]);
           dump(face_unit[2]);
        }
        */

#ifdef ENTROPY_SWITCH

#ifdef USE_MFM
        // As we are in a truly lagrangian configuration,
        // there is no adtypetion of entropy among particles.
        F_S=0.;
#else
        /*
        for (auto j=0; j<3; j++) F_S += F_v[j]*face_unit[j];
        if (F_S > 0) {
            // Maybe this values should be properly extrapolated to the faces..
            // but this is expensive!
            F_S *= pS*pOmega;
        }
        else {
            F_S *= qS*qOmega;
        }
        F_S *= F_rho*modApq;
        */
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

        /*
        printf("----\n");
        printf("%d\n", niter);
        dump(modApq);
        dump(F_p);
        dump(F_rho);
        dump(F_v[0]);
        dump(F_v[1]);
        dump(F_v[2]);
        */

        /*
        assert(!nan_guard(minDt, zero));
        assert(!nan_guard(F_rho, zero));
        assert(!nan_guard(F_p, zero));
        assert(!nan_guard(F_v[0], zero));
        assert(!nan_guard(F_v[1], zero));
        assert(!nan_guard(F_v[2], zero));
        */

        // Phew! Done ;)
    }

#ifdef OPTIM_FLUX_VEC

public:

// Simple macro to improve readability
#define q(X)    input_buffer[q_##X * nBuff + i]
#define qout(X) output_buffer[out_##X * nBuff + i]
    void hydroRiemann(PARTICLE *pIn,float fBall,int nSmooth, int nBuff,
                      my_real *restrict input_buffer,
                      my_real *restrict output_buffer, SMF *smf) {
        PKD pkd = smf->pkd;
        auto P = pkd->particles[pIn];

        const auto &pv = P.velocity();
        auto &psph = P.sph();

        const my_real pDensity = P.density();
        const my_real p_omega = psph.omega;

        bool bComove = pkd->csm->val.bComove;
        dtype dTime = smf->dTime;
        dtype dDelta = smf->dDelta;
        dtype a = smf->a;
        dtype H = smf->H;
        dtype dConstGamma = smf->dConstGamma;

        dtype pomega = psph.omega;
        dtype ph     = 0.5*P.ball();
        dtype plast  = psph.lastUpdateTime;
        dtype pDt    = smf->dDelta/(1<<P.rung());
        dtype pBXX   = psph.B[XX];
        dtype pBXY   = psph.B[XY];
        dtype pBXZ   = psph.B[XZ];
        dtype pBYY   = psph.B[YY];
        dtype pBYZ   = psph.B[YZ];
        dtype pBZZ   = psph.B[ZZ];
        dtype pDens  = pDensity;
        dtype pVpredx= P.velocity()[0];
        dtype pVpredy= P.velocity()[1];
        dtype pVpredz= P.velocity()[2];
        dtype pPres  = psph.P;
        dtype pS;
#ifdef ENTROPY_SWITCH
        pS    = psph.S;
#endif
        dtype pgradRhox = psph.gradRho[0];
        dtype pgradRhoy = psph.gradRho[1];
        dtype pgradRhoz = psph.gradRho[2];
        dtype pgradPx = psph.gradP[0];
        dtype pgradPy = psph.gradP[1];
        dtype pgradPz = psph.gradP[2];
        dtype pgradVxx = psph.gradVx[0];
        dtype pgradVxy = psph.gradVx[1];
        dtype pgradVxz = psph.gradVx[2];
        dtype pgradVyx = psph.gradVy[0];
        dtype pgradVyy = psph.gradVy[1];
        dtype pgradVyz = psph.gradVy[2];
        dtype pgradVzx = psph.gradVz[0];
        dtype pgradVzy = psph.gradVz[1];
        dtype pgradVzz = psph.gradVz[2];
        dtype plastAccx = psph.lastAcc[0];
        dtype plastAccy = psph.lastAcc[1];
        dtype plastAccz = psph.lastAcc[2];

        dtype qdr;
        dtype qdx;
        dtype qdy;
        dtype qdz;
        dtype qomega;
        dtype qh;
        dtype qlast;
        dtype qDt;
        dtype qBXX;
        dtype qBXY;
        dtype qBXZ;
        dtype qBYY;
        dtype qBYZ;
        dtype qBZZ;
        dtype qDens;
        dtype qvx;
        dtype qvy;
        dtype qvz;
        dtype qP;
        dtype qS;
        dtype qgradRhox;
        dtype qgradRhoy;
        dtype qgradRhoz;
        dtype qgradPx;
        dtype qgradPy;
        dtype qgradPz;
        dtype qgradVxx;
        dtype qgradVxy;
        dtype qgradVxz;
        dtype qgradVyx;
        dtype qgradVyy;
        dtype qgradVyz;
        dtype qgradVzx;
        dtype qgradVzy;
        dtype qgradVzz;
        dtype qlastAccx;
        dtype qlastAccy;
        dtype qlastAccz;

#ifdef __INTEL_COMPILER
//     __assume_aligned(input_buffer, 64);
// #pragma simd
// #pragma vector aligned
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
        assert(nBuff>(nSmooth+dtype::width()-1));
        for (auto i=0; i<nSmooth; i+=dtype::width()) {
            dtype F_rho;
            std::array<dtype,3> F_v;
            dtype F_P;
            dtype F_S;
            dtype minDt;

            qdr.load(       &q(dr));
            qdx.load(       &q(dx));
            qdy.load(       &q(dy));
            qdz.load(       &q(dz));
            qomega.load(    &q(omega));
            // Change name to h, not ball
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

            if (smf->dDelta <= 0.0) {
                pDt = 0.;
                qDt = 0.;
                plast = 0.;
                qlast = 0.;
                dTime = 0.;
            }

            doSinglePPFlux( F_rho, F_v, F_P, F_S, minDt,
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
                            qlastAccx,  qlastAccy,  qlastAccz,
                            smf->eEOS );

            // We fill the output buffer with the fluxes, which then
            // will be added to the corresponding particles
            F_rho.store(&output_buffer[out_Frho * nBuff + i]);
            F_P.store(&output_buffer[out_Fene * nBuff + i]);
            F_v[0].store(&output_buffer[out_FmomX * nBuff + i]);
            F_v[1].store(&output_buffer[out_FmomY * nBuff + i]);
            F_v[2].store(&output_buffer[out_FmomZ * nBuff + i]);
#ifdef ENTROPY_SWITCH
            F_S.store(&output_buffer[out_FS * nBuff + i]);
#endif
            //dump(F_v[0]);
            minDt.store(&output_buffer[out_minDt * nBuff + i]);

        }
    }
};

void hydroFluxFillBuffer(my_real *input_buffer, PARTICLE *qIn, int i, int nBuff,
                         double dr2, blitz::TinyVector<double,3> dr, SMF *smf) {
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

void hydroFluxGetBufferInfo(int *in, int *out) {
    *in = q_last;
    *out = out_last;
}

#endif // OPTIM_FLUX_VEC

void hydroRiemann_wrapper(PARTICLE *p,float fBall,int nSmooth, int nBuff,
                          my_real *restrict input_buffer,
                          my_real *restrict output_buffer, SMF *smf) {

#if defined(USE_SIMD_FLUX) && ( defined(HAVE_MM_POW) || defined(HAVE_MM256_POW) || defined(HAVE_MM512_POW) )
    MeshlessHydroSolver<dvec,dmask> solver;
#else
    MeshlessHydroSolver<vec<double,double>,mmask<bool>> solver;
#endif
    solver.hydroRiemann(p,fBall,nSmooth, nBuff,
                        input_buffer,
                        output_buffer, smf) ;
}
