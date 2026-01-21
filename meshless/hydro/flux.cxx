#include <new>
#include "blitz/tinyvec2.h"
#include "hydro/hydro.h"
#include "eEOS/eEOS.h"
#include "hydro/limiters.h"
#include "riemann.h"
#include "smooth/smooth.h"
#include "smooth/smoothfcn.h"
using blitz::dot;
using blitz::TinyVector;

using namespace meshless::FLUX; // IFLD and OFLD

void SmoothHelperHydroFluxes::pack(void* dst, const void* src) {
    auto p1 = static_cast<pack_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        const auto& sph = p2.sph();

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

void SmoothHelperHydroFluxes::unpack(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const pack_type*>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_adj_decorator();

        auto& sph = p1.sph();

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
void SmoothHelperHydroFluxes::init(void* dst) {
    auto p = pkd->particles[static_cast<PARTICLE*>(dst)];
    // For the init*Cached and comb functions we still have to explicitly
    // check if we are handling gas particles even ifdef OPTIM_REORDER_IN_NODES
    // because these operations are done in a cache-line basis, and a line may
    // contain other particles that are not of interest!
    if (p.is_gas()) {
        // assert((int)((char*)p.decorator_field()-(char*)dst) == pkd->particles.ParticleSize());
        auto& sph = p.sph();

        sph.Frho = 0.0;
        sph.Fmom = 0.0;
        sph.Fene = 0.0;
        sph.mom = 0.0;
        sph.E = 0.0;
        sph.Uint = 0.0;
        if (p.have_mfv()) {
            auto& mfv = p.mfv();
            mfv.drDotFrho = 0.0;
            p.set_mass(0.0);
        }
    }
}

void SmoothHelperHydroFluxes::flush(void* dst, const void* src) {
    auto p1 = static_cast<flush_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    if (p2.is_gas()) {
        const auto& sph = p2.sph();
        // assert((int)((char*)p2.decorator_field()-(char*)src) == pkd->particles.ParticleSize());

        p1->Frho = sph.Frho;
        p1->Fmom = sph.Fmom;
        p1->Fene = sph.Fene;
        p1->mom = sph.mom;
        p1->E = sph.E;
        p1->Uint = sph.Uint;
        if (p2.have_mfv()) {
            const auto& mfv = p2.mfv();
            p1->drDotFrho = mfv.drDotFrho;
            p1->fMass = p2.mass();
        }
    }
}

void SmoothHelperHydroFluxes::combine(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const flush_type*>(src);

    if (p1.is_gas()) {
        auto& sph = p1.sph();

        sph.Frho += p2->Frho;
        sph.Fmom += p2->Fmom;
        sph.Fene += p2->Fene;
        sph.mom += p2->mom;
        sph.E += p2->E;
        sph.Uint += p2->Uint;
        if (p1.have_mfv()) {
            auto& mfv = p1.mfv();
            mfv.drDotFrho += p2->drDotFrho;
            p1.set_mass(p1.mass() + p2->fMass);
        }
    }
}

template <typename dtype = dvec, typename mtype = dmask>
class MeshlessHydroSolver {
    bool bMFV;

  public:
    explicit MeshlessHydroSolver(bool bMFV) : bMFV(bMFV) {}

  private:
    inline void extrapolateDensityInTime(dtype& rho, dtype rho0, dtype vx, dtype vy, dtype vz, dtype dt, dtype gradRhoX, dtype gradRhoY, dtype gradRhoZ,
                                         dtype divv) {
        rho -= dt * (vx * gradRhoX + vy * gradRhoY + vz * gradRhoZ + rho0 * divv);
    }

    inline void extrapolateVelocityInTime(dtype& v, dtype& vFrame, dtype divv, dtype v0, dtype dt, dtype acc, dtype gradP, dtype rho0, dtype a) {
        dtype temp;
        temp = -(v0 * divv + gradP / rho0) * dt;
        temp += acc * dt * a;
        v += temp;
        vFrame += 0.5 * temp;
    }

    inline void extrapolatePressureInTime(dtype& p, dtype dt, dtype gradPx, dtype gradPy, dtype gradPz, dtype p0, dtype vx0, dtype vy0, dtype vz0, dtype divv,
                                          dtype dConstGamma) {
        p -= dt * (vx0 * gradPx + vy0 * gradPy + vz0 * gradPz + dConstGamma * p0 * divv);
    }

    inline void extrapolateStateInTime(dtype& rho, dtype& vx, dtype& vy, dtype& vz, dtype& p, dtype& vFx, dtype& vFy, dtype& vFz, dtype rho0, dtype vx0,
                                       dtype vy0, dtype vz0, dtype p0, dtype dt, dtype gradRhoX, dtype gradRhoY, dtype gradRhoZ, dtype gradPX, dtype gradPY,
                                       dtype gradPZ, dtype accx, dtype accy, dtype accz, dtype divv, dtype dConstGamma, dtype a) {

        extrapolateDensityInTime(rho, rho0, vx, vy, vz, dt, gradRhoX, gradRhoY, gradRhoZ, divv);

        extrapolateVelocityInTime(vx, vFx, divv, vx0, dt, accx, gradPX, rho0, a);
        extrapolateVelocityInTime(vy, vFy, divv, vy0, dt, accy, gradPY, rho0, a);
        extrapolateVelocityInTime(vz, vFz, divv, vz0, dt, accz, gradPZ, rho0, a);

        extrapolatePressureInTime(p, dt, gradPX, gradPY, gradPZ, p0, vx0, vy0, vz0, divv, dConstGamma);
    }

    inline void extrapolateVariableInSpace(dtype& var, dtype dx, dtype dy, dtype dz, dtype gradx, dtype grady, dtype gradz) {
        var += dx * gradx + dy * grady + dz * gradz;
    }

    inline void extrapolateVelocityCosmology(dtype& v, dtype& vFrame, dtype v0, dtype dt, dtype H, dtype a) {
        dtype temp = H * dt * a * v;
        v -= temp;
        vFrame -= 0.5 * temp;
    }

    inline void extrapolateCosmology(dtype& vx, dtype& vy, dtype& vz, dtype& vFramex, dtype& vFramey, dtype& vFramez, dtype& p, dtype dt, dtype vx0, dtype vy0,
                                     dtype vz0, dtype p0, dtype dConstGamma, dtype H, dtype a) {
        extrapolateVelocityCosmology(vx, vFramex, vx0, dt, H, a);
        extrapolateVelocityCosmology(vy, vFramey, vy0, dt, H, a);
        extrapolateVelocityCosmology(vz, vFramez, vz0, dt, H, a);

        p -= 3. * H * dt * a * (dConstGamma - 1.) * p0;
    }

    inline void computeFace(dtype& modApq, std::array<dtype, 3>& unit, dtype rpq, dtype dx, dtype dy, dtype dz, dtype pH, dtype qH, dtype p_omega,
                            dtype q_omega, dtype pBxx, dtype pBxy, dtype pBxz, dtype pByy, dtype pByz, dtype pBzz, dtype qBxx, dtype qBxy, dtype qBxz,
                            dtype qByy, dtype qByz, dtype qBzz) {
        dtype psi;
        dtype psiTilde[3];
        dtype Apq[3] = {0., 0., 0.};

        // \tilde{\psi}_j (x_i)
        psi = -cubicSplineKernel(rpq, pH) / p_omega;
        psiTilde[0] = (pBxx * dx + pBxy * dy + pBxz * dz) * psi;
        psiTilde[1] = (pBxy * dx + pByy * dy + pByz * dz) * psi;
        psiTilde[2] = (pBxz * dx + pByz * dy + pBzz * dz) * psi;
        for (auto j = 0; j < 3; j++) {
            Apq[j] += psiTilde[j] / p_omega;
        }

        // \tilde{\psi}_i (x_j)
        psi = cubicSplineKernel(rpq, qH) / q_omega;
        psiTilde[0] = (qBxx * dx + qBxy * dy + qBxz * dz) * psi;
        psiTilde[1] = (qBxy * dx + qByy * dy + qByz * dz) * psi;
        psiTilde[2] = (qBxz * dx + qByz * dy + qBzz * dz) * psi;
        for (auto j = 0; j < 3; j++) {
            Apq[j] -= psiTilde[j] / q_omega;
        }

        modApq = 0.0;
        for (auto j = 0; j < 3; j++) {
            modApq += Apq[j] * Apq[j];
        }
        modApq = sqrt(modApq);
        for (auto j = 0; j < 3; j++) {
            unit[j] = Apq[j] / modApq;
        }
        // Avoid NaNs when the face has zero area
        dtype zero = 0.;
        mtype cond = modApq == zero;
        if (movemask(cond)) {
            for (auto j = 0; j < 3; j++) unit[j] = mask_mov(unit[j], cond, zero);
        }
    }

    inline void low_limit(dtype& var, dtype def, dtype min) {
        mtype cond = var < min;
        var = mask_mov(var, cond, def);
    }

    inline void doSinglePPFlux(const mtype mask, dtype& F_rho, std::array<dtype, 3>& F_v, dtype& F_p, dtype& F_S, dtype& minDt, bool bComove, dtype dTime,
                               dtype dDelta, dtype a, dtype H, dtype dConstGamma, dtype rpq, dtype dx, dtype dy, dtype dz, dtype pBall, dtype pLastUpdateTime,
                               dtype pDt, dtype pOmega, dtype pSPHEoM, dtype pBxx, dtype pBxy, dtype pBxz, dtype pByy, dtype pByz, dtype pBzz, dtype pDensity,
                               dtype pVpredx, dtype pVpredy, dtype pVpredz, dtype pP, dtype pS, dtype pGradRhoX, dtype pGradRhoY, dtype pGradRhoZ,
                               dtype pGradPX, dtype pGradPY, dtype pGradPZ, dtype pGradVxX, dtype pGradVxY, dtype pGradVxZ, dtype pGradVyX, dtype pGradVyY,
                               dtype pGradVyZ, dtype pGradVzX, dtype pGradVzY, dtype pGradVzZ, dtype pLastAccX, dtype pLastAccY, dtype pLastAccZ, dtype qBall,
                               dtype qLastUpdateTime, dtype qDt, dtype qOmega, dtype qSPHEoM, dtype qBxx, dtype qBxy, dtype qBxz, dtype qByy, dtype qByz,
                               dtype qBzz, dtype qDensity, dtype qVpredx, dtype qVpredy, dtype qVpredz, dtype qP, dtype qS, dtype qGradRhoX, dtype qGradRhoY,
                               dtype qGradRhoZ, dtype qGradPX, dtype qGradPY, dtype qGradPZ, dtype qGradVxX, dtype qGradVxY, dtype qGradVxZ, dtype qGradVyX,
                               dtype qGradVyY, dtype qGradVyZ, dtype qGradVzX, dtype qGradVzY, dtype qGradVzZ, dtype qLastAccX, dtype qLastAccY,
                               dtype qLastAccZ, const eEOSparam& eEOS) {

#ifdef FORCE_1D
        if (dz != 0) continue;
        if (dy != 0) continue;
#endif
#ifdef FORCE_2D
        if (dz != 0) continue;
#endif

        /* We update the conservatives variables taking the minimum timestep
         * between the particles, as in AREPO
         */
        minDt = min(pDt, qDt);
        minDt /= a;

        dtype qDeltaHalf = 0.0, pDeltaHalf = 0.0;
        pDeltaHalf = (dTime - pLastUpdateTime + 0.5 * pDt) / a;
        qDeltaHalf = (dTime - qLastUpdateTime + 0.5 * qDt) / a;

        // DEBUG: Avoid temporal extrapolation
        // pDeltaHalf = 0.;
        // qDeltaHalf = 0.;

        dtype modApq;
        std::array<dtype, 3> face_unit;
        computeFace(modApq, face_unit, rpq, dx, dy, dz, pBall, qBall, pOmega, qOmega, pBxx, pBxy, pBxz, pByy, pByz, pBzz, qBxx, qBxy, qBxz, qByy, qByz, qBzz);

        // Velocity of the quadrature mid-point
        dtype vFrame[3];
        vFrame[0] = 0.5 * (pVpredx + qVpredx);
        vFrame[1] = 0.5 * (pVpredy + qVpredy);
        vFrame[2] = 0.5 * (pVpredz + qVpredz);

        dtype pv[3], qv[3];
        // We boost to the reference of the p-q 'face'
        pv[0] = pVpredx - vFrame[0];
        qv[0] = qVpredx - vFrame[0];
        pv[1] = pVpredy - vFrame[1];
        qv[1] = qVpredy - vFrame[1];
        pv[2] = pVpredz - vFrame[2];
        qv[2] = qVpredz - vFrame[2];

        // Mid-point rule
        dx = -0.5 * dx;
        dy = -0.5 * dy;
        dz = -0.5 * dz;

        // DEBUG: Avoid spatial extrapolation
        // dx = 0.0;
        // dy = 0.0;
        // dz = 0.0;

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

        extrapolateStateInTime(L_rho, L_v[0], L_v[1], L_v[2], L_p, vFrame[0], vFrame[1], vFrame[2], pDensity, pv[0], pv[1], pv[2], pP, pDeltaHalf, pGradRhoX,
                               pGradRhoY, pGradRhoZ, pGradPX, pGradPY, pGradPZ, pLastAccX, pLastAccY, pLastAccZ, pdivv, dConstGamma, a);

        extrapolateVariableInSpace(L_rho, dx, dy, dz, pGradRhoX, pGradRhoY, pGradRhoZ);
        extrapolateVariableInSpace(L_p, dx, dy, dz, pGradPX, pGradPY, pGradPZ);
        extrapolateVariableInSpace(L_v[0], dx, dy, dz, pGradVxX, pGradVxY, pGradVxZ);
        extrapolateVariableInSpace(L_v[1], dx, dy, dz, pGradVyX, pGradVyY, pGradVyZ);
        extrapolateVariableInSpace(L_v[2], dx, dy, dz, pGradVzX, pGradVzY, pGradVzZ);

        dx = -dx;
        dy = -dy;
        dz = -dz;

        extrapolateStateInTime(R_rho, R_v[0], R_v[1], R_v[2], R_p, vFrame[0], vFrame[1], vFrame[2], qDensity, qv[0], qv[1], qv[2], qP, qDeltaHalf, qGradRhoX,
                               qGradRhoY, qGradRhoZ, qGradPX, qGradPY, qGradPZ, qLastAccX, qLastAccY, qLastAccZ, qdivv, dConstGamma, a);

        extrapolateVariableInSpace(R_rho, dx, dy, dz, qGradRhoX, qGradRhoY, qGradRhoZ);
        extrapolateVariableInSpace(R_p, dx, dy, dz, qGradPX, qGradPY, qGradPZ);
        extrapolateVariableInSpace(R_v[0], dx, dy, dz, qGradVxX, qGradVxY, qGradVxZ);
        extrapolateVariableInSpace(R_v[1], dx, dy, dz, qGradVyX, qGradVyY, qGradVyZ);
        extrapolateVariableInSpace(R_v[2], dx, dy, dz, qGradVzX, qGradVzY, qGradVzZ);

        genericPairwiseLimiter(pDensity, qDensity, L_rho, R_rho);
        genericPairwiseLimiter(pP, qP, L_p, R_p);
        genericPairwiseLimiter(pv[0], qv[0], L_v[0], R_v[0]);
        genericPairwiseLimiter(pv[1], qv[1], L_v[1], R_v[1]);
        genericPairwiseLimiter(pv[2], qv[2], L_v[2], R_v[2]);

        if (bComove) {
            extrapolateCosmology(R_v[0], R_v[1], R_v[2], vFrame[0], vFrame[1], vFrame[2], R_p, qDeltaHalf, qv[0], qv[1], qv[2], qP, dConstGamma, H, a);

            extrapolateCosmology(L_v[0], L_v[1], L_v[2], vFrame[0], vFrame[1], vFrame[2], L_p, pDeltaHalf, pv[0], pv[1], pv[2], pP, dConstGamma, H, a);
        }

        dtype zero = 0.0;
        low_limit(L_rho, pDensity, zero);
        low_limit(R_rho, qDensity, zero);
        low_limit(L_p, pP, zero);
        low_limit(R_p, qP, zero);

        if (eEOS.bDoEOS) {
            dtype a_inv3 = 1.0 / (a * a * a);
            const dtype pLpoly = polytropicPressureFloor<dtype, mtype>(a_inv3, L_rho, dConstGamma, eEOS);
            const dtype pRpoly = polytropicPressureFloor<dtype, mtype>(a_inv3, R_rho, dConstGamma, eEOS);
            L_p = max(L_p, pLpoly);
            R_p = max(R_p, pRpoly);
        }
        /*
        #ifdef EEOS_JEANS
        const double pLjeans =
            jeansPressureFloor(L_rho, pH, dConstGamma, dEOSNJeans);
        const double pRjeans =
            jeansPressureFloor(R_rho, bi(IFLD::ball,i), dConstGamma, dEOSNJeans);
        L_p = max(L_p, pLjeans);
        R_p = max(R_p, pRjeans);
        #endif
        */

        dtype P_M, S_M;

        RiemannSolverExact<dtype, mtype> riemann(dConstGamma, mask, bMFV);
        int niter = riemann.solve(R_rho, R_p, R_v, L_rho, L_p, L_v, P_M, S_M, &F_rho, &F_p, F_v.data(), face_unit.data());

        /*
        int nan;
        // Only works if compiling with -fno-finite-math-only !!
        auto aa = maskz_mov(mask, P_M);
        auto bb = maskz_mov(mask, S_M);
        auto cc = maskz_mov(mask, modApq);
        nan = nan_guard(aa, zero) | nan_guard(bb, zero) | nan_guard(cc, zero);
        if (nan) {
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
            dump(modApq);
        }
        */

#ifdef ENTROPY_SWITCH

        if (!bMFV) {
            // As we are in a truly lagrangian configuration,
            // there is no adtypetion of entropy among particles.
            F_S = 0.;
        }
        else {
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
        }
#endif // ENTROPY_SWITCH

        if (!bMFV) { // MFM
            F_rho = 0.;
            F_p = P_M * S_M;
            for (auto j = 0; j < 3; j++) F_v[j] = P_M * face_unit[j];
        }

        // Force 2D
#ifdef FORCE_1D
        F_v[2] = 0.;
        F_v[1] = 0.;
#endif
#ifdef FORCE_2D
        F_v[2] = 0.;
#endif

        // Now we de-boost the fluxes following Eq. A8 Hopkins 2015
        for (auto j = 0; j < 3; j++) {
            F_p += vFrame[j] * F_v[j];
            F_p += (0.5 * vFrame[j] * vFrame[j]) * F_rho;
        }

        // Now we just multiply by the face area
        F_p *= modApq;
        F_rho *= modApq;
        for (auto j = 0; j < 3; j++) {
            F_v[j] *= modApq;
            F_v[j] += vFrame[j] * F_rho;
        }

#ifdef MFM_FLUX_CORRECTION
        if (!bMFV) { // MFM
            const dtype cs_p{sqrt(dConstGamma * pP / pDensity)};
            const dtype cs_q{sqrt(dConstGamma * qP / qDensity)};
            mtype cond{abs(S_M / min(cs_p, cs_q)) < 0.6};
            if (movemask(cond)) {
                const dtype dU_mfm_p = pVpredx * F_v[0] + pVpredy * F_v[1] + pVpredz * F_v[2] - F_p;
                const dtype dU_mfm_q = F_p - qVpredx * F_v[0] - qVpredy * F_v[1] - qVpredz * F_v[2];

                const dtype dU_sph_fac = 2. * P_M * ((qVpredx - pVpredx) * dx + (qVpredy - pVpredy) * dy + (qVpredz - pVpredz) * dz) / rpq;
                const dtype dU_sph_p = dU_sph_fac * pSPHEoM * cubicSplineKernelDeriv(rpq, pBall) / pOmega;
                const dtype dU_sph_q = dU_sph_fac * qSPHEoM * cubicSplineKernelDeriv(rpq, qBall) / qOmega;

                mtype pIsCooler{cs_p < cs_q};
                mtype pDoCorr{cond & pIsCooler & (dU_mfm_p < 0.0) & (dU_sph_p > dU_mfm_p)};
                mtype qDoCorr{cond & ~pIsCooler & (dU_mfm_q < 0.0) & (dU_sph_q > dU_mfm_q)};

                F_p = mask_mov(F_p, pDoCorr, F_p + dU_mfm_p - dU_sph_p);
                F_p = mask_mov(F_p, qDoCorr, F_p - dU_mfm_q + dU_sph_q);
            }
        }
#endif

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
    void hydroRiemann(PKD pkd, PARTICLE* pIn, float fBall, int nSmooth, EnumArray<meshless::myreal, IFLD>& bi, EnumArray<meshless::myreal, OFLD>& bo,
                      double dTime, double dDelta, double a, double H, double dConstGamma, const eEOSparam& eEOS) {
        auto P = pkd->particles[pIn];

        const auto& pv = P.velocity();
        auto& psph = P.sph();

        bool bComove = pkd->csm->val.bComove;

        dtype pomega = psph.omega;
        dtype pspheom = psph.spheom;
        dtype pH = P.ball();
        dtype plast = psph.lastUpdateTime;
        dtype pDt = dDelta / (1 << P.rung());
        dtype pBXX = psph.B[+MOI::XX];
        dtype pBXY = psph.B[+MOI::XY];
        dtype pBXZ = psph.B[+MOI::XZ];
        dtype pBYY = psph.B[+MOI::YY];
        dtype pBYZ = psph.B[+MOI::YZ];
        dtype pBZZ = psph.B[+MOI::ZZ];
        dtype pDens = P.density();
        dtype pVpredx = P.velocity()[0];
        dtype pVpredy = P.velocity()[1];
        dtype pVpredz = P.velocity()[2];
        dtype pPres = psph.P;
        dtype pS;
    #ifdef ENTROPY_SWITCH
        pS = psph.S;
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
        dtype qspheom;
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

        double index[dtype::width()];
        for (auto i = 0; i < dtype::width(); i++) index[i] = i;
        dtype index_v;
        index_v.load(index);

        assert(bi.extent(1) >= ((nSmooth + dtype::mask()) & ~dtype::mask()));
        for (auto i = 0; i < nSmooth; i += dtype::width()) {
            dtype F_rho;
            std::array<dtype, 3> F_v;
            dtype F_P;
            dtype F_S;
            dtype minDt;

            const mtype mask = index_v < (double)nSmooth;
            index_v += (double)dtype::width();

            qdr.load(&bi(IFLD::dr, i));
            qdx.load(&bi(IFLD::dx, i));
            qdy.load(&bi(IFLD::dy, i));
            qdz.load(&bi(IFLD::dz, i));
            qomega.load(&bi(IFLD::omega, i));
            qspheom.load(&bi(IFLD::spheom, i));
            // Change name to h, not ball
            qh.load(&bi(IFLD::ball, i));
            qlast.load(&bi(IFLD::lastUpdateTime, i));
            qDt.load(&bi(IFLD::rung, i));
            qBXX.load(&bi(IFLD::B_XX, i));
            qBXY.load(&bi(IFLD::B_XY, i));
            qBXZ.load(&bi(IFLD::B_XZ, i));
            qBYY.load(&bi(IFLD::B_YY, i));
            qBYZ.load(&bi(IFLD::B_YZ, i));
            qBZZ.load(&bi(IFLD::B_ZZ, i));
            qDens.load(&bi(IFLD::rho, i));
            qvx.load(&bi(IFLD::vx, i));
            qvy.load(&bi(IFLD::vy, i));
            qvz.load(&bi(IFLD::vz, i));
            qP.load(&bi(IFLD::P, i));
    #ifdef ENTROPY_SWITCH
            qS.load(&bi(IFLD::S, i));
    #endif
            qgradRhox.load(&bi(IFLD::gradRhoX, i));
            qgradRhoy.load(&bi(IFLD::gradRhoY, i));
            qgradRhoz.load(&bi(IFLD::gradRhoZ, i));
            qgradPx.load(&bi(IFLD::gradPX, i));
            qgradPy.load(&bi(IFLD::gradPY, i));
            qgradPz.load(&bi(IFLD::gradPZ, i));
            qgradVxx.load(&bi(IFLD::gradVxX, i));
            qgradVxy.load(&bi(IFLD::gradVxY, i));
            qgradVxz.load(&bi(IFLD::gradVxZ, i));
            qgradVyx.load(&bi(IFLD::gradVyX, i));
            qgradVyy.load(&bi(IFLD::gradVyY, i));
            qgradVyz.load(&bi(IFLD::gradVyZ, i));
            qgradVzx.load(&bi(IFLD::gradVzX, i));
            qgradVzy.load(&bi(IFLD::gradVzY, i));
            qgradVzz.load(&bi(IFLD::gradVzZ, i));
            qlastAccx.load(&bi(IFLD::lastAccX, i));
            qlastAccy.load(&bi(IFLD::lastAccY, i));
            qlastAccz.load(&bi(IFLD::lastAccZ, i));

            if (dDelta <= 0.0) {
                pDt = 0.;
                qDt = 0.;
                plast = 0.;
                qlast = 0.;
                dTime = 0.;
            }

            doSinglePPFlux(mask, F_rho, F_v, F_P, F_S, minDt, bComove, dTime, dDelta, a, H, dConstGamma, qdr, qdx, qdy, qdz, pH, plast, pDt, pomega, pspheom,
                           pBXX, pBXY, pBXZ, pBYY, pBYZ, pBZZ, pDens, pVpredx, pVpredy, pVpredz, pPres, pS, pgradRhox, pgradRhoy, pgradRhoz, pgradPx, pgradPy,
                           pgradPz, pgradVxx, pgradVxy, pgradVxz, pgradVyx, pgradVyy, pgradVyz, pgradVzx, pgradVzy, pgradVzz, plastAccx, plastAccy, plastAccz,
                           qh, qlast, qDt, qomega, qspheom, qBXX, qBXY, qBXZ, qBYY, qBYZ, qBZZ, qDens, qvx, qvy, qvz, qP, qS, qgradRhox, qgradRhoy, qgradRhoz,
                           qgradPx, qgradPy, qgradPz, qgradVxx, qgradVxy, qgradVxz, qgradVyx, qgradVyy, qgradVyz, qgradVzx, qgradVzy, qgradVzz, qlastAccx,
                           qlastAccy, qlastAccz, eEOS);

            // We fill the output buffer with the fluxes, which then
            // will be added to the corresponding particles
            F_rho.store(&bo(OFLD::Frho, i));
            F_P.store(&bo(OFLD::Fene, i));
            F_v[0].store(&bo(OFLD::FmomX, i));
            F_v[1].store(&bo(OFLD::FmomY, i));
            F_v[2].store(&bo(OFLD::FmomZ, i));
    #ifdef ENTROPY_SWITCH
            F_S.store(&bo(OFLD::FS, i));
    #endif
            // dump(F_v[0]);
            minDt.store(&bo(OFLD::minDt, i));
        }
    }
};

template <typename T, typename E>
void hydroFluxFillBuffer(PKD pkd, EnumArray<T, E>& bi, const NN& nn, int i, double dDelta) {
    auto qIn = nn.pPart;
    auto dr2 = nn.fDist2;
    auto dr = nn.dr;
    auto Q = pkd->particles[qIn];
    auto& qsph = Q.sph();
    bi(IFLD::ball, i) = Q.ball();
    bi(IFLD::dx, i) = dr[0];
    bi(IFLD::dy, i) = dr[1];
    bi(IFLD::dz, i) = dr[2];
    bi(IFLD::dr, i) = sqrt(dr2);
    bi(IFLD::rung, i) = dDelta / (1 << Q.rung());
    bi(IFLD::rho, i) = Q.density();
    bi(IFLD::P, i) = qsph.P;
    #ifdef ENTROPY_SWITCH
    bi(IFLD::S, i) = qsph.S;
    #endif
    const auto& qv = Q.velocity();
    bi(IFLD::vx, i) = qv[0];
    bi(IFLD::vy, i) = qv[1];
    bi(IFLD::vz, i) = qv[2];

    bi(IFLD::gradRhoX, i) = qsph.gradRho[0];
    bi(IFLD::gradRhoY, i) = qsph.gradRho[1];
    bi(IFLD::gradRhoZ, i) = qsph.gradRho[2];

    bi(IFLD::gradPX, i) = qsph.gradP[0];
    bi(IFLD::gradPY, i) = qsph.gradP[1];
    bi(IFLD::gradPZ, i) = qsph.gradP[2];

    bi(IFLD::gradVxX, i) = qsph.gradVx[0];
    bi(IFLD::gradVxY, i) = qsph.gradVx[1];
    bi(IFLD::gradVxZ, i) = qsph.gradVx[2];

    bi(IFLD::gradVyX, i) = qsph.gradVy[0];
    bi(IFLD::gradVyY, i) = qsph.gradVy[1];
    bi(IFLD::gradVyZ, i) = qsph.gradVy[2];

    bi(IFLD::gradVzX, i) = qsph.gradVz[0];
    bi(IFLD::gradVzY, i) = qsph.gradVz[1];
    bi(IFLD::gradVzZ, i) = qsph.gradVz[2];

    bi(IFLD::lastUpdateTime, i) = qsph.lastUpdateTime;
    bi(IFLD::lastAccX, i) = qsph.lastAcc[0];
    bi(IFLD::lastAccY, i) = qsph.lastAcc[1];
    bi(IFLD::lastAccZ, i) = qsph.lastAcc[2];
    bi(IFLD::B_XX, i) = qsph.B[+MOI::XX];
    bi(IFLD::B_YY, i) = qsph.B[+MOI::YY];
    bi(IFLD::B_ZZ, i) = qsph.B[+MOI::ZZ];
    bi(IFLD::B_XY, i) = qsph.B[+MOI::XY];
    bi(IFLD::B_XZ, i) = qsph.B[+MOI::XZ];
    bi(IFLD::B_YZ, i) = qsph.B[+MOI::YZ];
    bi(IFLD::omega, i) = qsph.omega;
    bi(IFLD::spheom, i) = qsph.spheom;
}

template <typename T, typename E>
void hydroFluxUpdateFromBuffer(PKD pkd, PARTICLE* pIn, EnumArray<T, E>& bo, const NN& nn, int i, double dDelta, double aFac) {
    auto qIn = nn.pPart;
    auto dr = nn.dr;
    auto P = pkd->particles[pIn];
    auto Q = pkd->particles[qIn];
    const auto& pv = P.velocity();
    const auto& qv = Q.velocity();
    auto& psph = P.sph();
    auto& qsph = Q.sph();
    if (dDelta > 0) {
        if (P.have_mfv()) {
            P.set_mass(P.mass() - bo(OFLD::minDt, i) * bo(OFLD::Frho, i));
        }

        psph.mom[0] -= bo(OFLD::minDt, i) * bo(OFLD::FmomX, i);
        psph.mom[1] -= bo(OFLD::minDt, i) * bo(OFLD::FmomY, i);
        psph.mom[2] -= bo(OFLD::minDt, i) * bo(OFLD::FmomZ, i);

        psph.E -= bo(OFLD::minDt, i) * bo(OFLD::Fene, i);

        psph.Uint -= bo(OFLD::minDt, i) * (bo(OFLD::Fene, i) - bo(OFLD::FmomX, i) * pv[0] - bo(OFLD::FmomY, i) * pv[1] - bo(OFLD::FmomZ, i) * pv[2] +
                                           0.5 * dot(pv, pv) * bo(OFLD::Frho, i));
    #ifdef ENTROPY_SWITCH
        psph.S -= bo(OFLD::minDt, i) * bo(OFLD::FS, i);
    #endif

        if (P.have_mfv()) {
            auto& mfv = P.mfv();
            mfv.drDotFrho[0] += bo(OFLD::minDt, i) * bo(OFLD::Frho, i) * dr[0] * aFac;
            mfv.drDotFrho[1] += bo(OFLD::minDt, i) * bo(OFLD::Frho, i) * dr[1] * aFac;
            mfv.drDotFrho[2] += bo(OFLD::minDt, i) * bo(OFLD::Frho, i) * dr[2] * aFac;
        }
        psph.Frho += bo(OFLD::Frho, i);
        psph.Fene += bo(OFLD::Fene, i);
        psph.Fmom[0] += bo(OFLD::FmomX, i);
        psph.Fmom[1] += bo(OFLD::FmomY, i);
        psph.Fmom[2] += bo(OFLD::FmomZ, i);
    }
    else {
        psph.Frho += bo(OFLD::Frho, i);
        psph.Fene += bo(OFLD::Fene, i);
        psph.Fmom[0] += bo(OFLD::FmomX, i);
        psph.Fmom[1] += bo(OFLD::FmomY, i);
        psph.Fmom[2] += bo(OFLD::FmomZ, i);
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
        // assert(qsph.P == bi(IFLD::P, i));
        if (dDelta > 0) {
            if (Q.have_mfv()) {
                Q.set_mass(Q.mass() + bo(OFLD::minDt, i) * bo(OFLD::Frho, i));
            }

            qsph.mom[0] += bo(OFLD::minDt, i) * bo(OFLD::FmomX, i);
            qsph.mom[1] += bo(OFLD::minDt, i) * bo(OFLD::FmomY, i);
            qsph.mom[2] += bo(OFLD::minDt, i) * bo(OFLD::FmomZ, i);

            qsph.E += bo(OFLD::minDt, i) * bo(OFLD::Fene, i);

            qsph.Uint += bo(OFLD::minDt, i) * (bo(OFLD::Fene, i) - bo(OFLD::FmomX, i) * qv[0] - bo(OFLD::FmomY, i) * qv[1] - bo(OFLD::FmomZ, i) * qv[2] +
                                               0.5 * dot(qv, qv) * bo(OFLD::Frho, i));
    #ifdef ENTROPY_SWITCH
            qsph.S += bo(OFLD::minDt, i) * bo(OFLD::FS, i);
    #endif

            if (Q.have_mfv()) {
                auto& mfv = P.mfv();
                mfv.drDotFrho[0] += bo(OFLD::minDt, i) * bo(OFLD::Frho, i) * dr[0] * aFac;
                mfv.drDotFrho[1] += bo(OFLD::minDt, i) * bo(OFLD::Frho, i) * dr[1] * aFac;
                mfv.drDotFrho[2] += bo(OFLD::minDt, i) * bo(OFLD::Frho, i) * dr[2] * aFac;
            }
            qsph.Frho -= bo(OFLD::Frho, i);
            qsph.Fene -= bo(OFLD::Fene, i);
            qsph.Fmom[0] -= bo(OFLD::FmomX, i);
            qsph.Fmom[1] -= bo(OFLD::FmomY, i);
            qsph.Fmom[2] -= bo(OFLD::FmomZ, i);
        }
        else {
            qsph.Frho -= bo(OFLD::Frho, i);
            qsph.Fene -= bo(OFLD::Fene, i);
            qsph.Fmom[0] -= bo(OFLD::FmomX, i);
            qsph.Fmom[1] -= bo(OFLD::FmomY, i);
            qsph.Fmom[2] -= bo(OFLD::FmomZ, i);
        }

    } // q marked/active
}

#endif // OPTIM_FLUX_VEC

struct SmoothBallOptimisedFluxes : smooth::policy::ball {
    static auto inline accept(PKD pkd, PARTICLE* P, double fDist2, const TinyVector<double, 3>& dr_node, const NN& nn) -> bool {
#ifndef OPTIM_AVOID_IS_ACTIVE
        abort(); // not supported
#endif
        const auto& qBall = nn.fBall;
        if (qBall * qBall < fDist2) return false;
#ifdef OPTIM_NO_REDUNDANT_FLUXES
        if (nn.bMarked) { // OPTIM_AVOID_IS_ACTIVE
            if (pkd->Self() == nn.iPid) {
                if (P < nn.pPart) return false;
            }
            else if (pkd->Self() < nn.iPid) return false;
        }
#endif
        return true;
    }
};

auto ServiceHydroFluxes::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->lcl()->pkd;
    assert(pkd->decorators.present(DECORATOR_FIELD::oMFM)); /* Validate memory model */
    auto helper = std::make_shared<SmoothHelperHydroFluxes>(pkd, in->dTime, in->dDelta, in->a, in->H, in->dConstGamma, in->eEOS);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
#ifdef OPTIM_SMOOTH_NODE
    *out = smx.ReSmoothNode<smooth::policy::OnlyGasOptimised, smooth::policy::OnlyGas, SmoothBallOptimisedFluxes>();
#else
        *out = smx.ReSmooth<smooth::policy::OnlyGasOptimised>();
#endif
    return sizeof(output);
}

void SmoothHelperHydroFluxes::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
#if defined(USE_SIMD_FLUX)
    using DVEC = dvec;
    using MASK = dmask;
#else
        using DVEC = vec<double, double>;
        using MASK = mmask<bool>;
#endif

    // Allocate or extend the buffers if required
    auto nBuff = (nnList.size() + DVEC::mask()) & ~DVEC::mask();
    if (ia.extent(1) < nBuff) {
        auto nByte = std::max(static_cast<size_t>(256), 2 * nBuff) * static_cast<size_t>(IFLD::last);
        EnumArray na(new (std::align_val_t(sizeof(DVEC))) meshless::myreal[nByte], IFLD::last, nBuff);
        ia.reference(na);
    }
    if (oa.extent(1) < nBuff) {
        auto nByte = std::max(static_cast<size_t>(256), 2 * nBuff) * static_cast<size_t>(OFLD::last);
        EnumArray na(new (std::align_val_t(sizeof(DVEC))) meshless::myreal[nByte], OFLD::last, nBuff);
        oa.reference(na);
    }

    // Vectorize the quantities from the neighbors into the input buffer
    for (auto i = 0; i < nnList.size(); i++) {
        hydroFluxFillBuffer(pkd, ia, nnList[i], i, dDelta);
    }

    MeshlessHydroSolver<DVEC, MASK> solver(p.have_mfv());
    solver.hydroRiemann(pkd, &p, fBall, nnList.size(), ia, oa, dTime, dDelta, a_, H_, dConstGamma, eEOS);

    // Gather the results from the Riemann solver
    for (auto i = 0; i < nnList.size(); i++) {
        hydroFluxUpdateFromBuffer(pkd, &p, oa, nnList[i], i, dDelta, a_);
    }
}
