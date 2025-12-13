#ifndef SOLVER_H
#define SOLVER_H
#include "core/element.h"
#include "core/vmath.h"
#include "riemann.h"
#include "limiters.h"
#include "kernels.h"

#ifdef __CUDACC__
    #define PP_CUDA_BOTH __host__ __device__
#else
    #define PP_CUDA_BOTH
#endif

enum q_data {
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
    q_spheom,
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

#define XX 0
#define YY 3
#define ZZ 5
#define XY 1
#define XZ 2
#define YZ 4

template<class F>
struct ResultMFM {
    F minDt, Frho, Fmom[3], Fene, FS;
    PP_CUDA_BOTH void zero() {
        Frho=Fene=FS=0.0f;
        Fmom[0] = 0.0f;
        Fmom[1] = 0.0f;
        Fmom[2] = 0.0f;
        minDt=0.0f;
    }
    PP_CUDA_BOTH ResultMFM<F> operator+=(const ResultMFM<F> rhs) {
        Frho += rhs.Frho;
        Fmom[0] += rhs.Fmom[0];
        Fmom[1] += rhs.Fmom[1];
        Fmom[2] += rhs.Fmom[2];
        Fene += rhs.Fene;
        FS += rhs.FS;
        return *this;
    }
};

template <typename dtype=dvec, typename mtype=dmask>
class MeshlessHydroSolver {
    bool bMFV;
public:
    PP_CUDA_BOTH explicit MeshlessHydroSolver(bool bMFV) : bMFV(bMFV) {}
private:


    PP_CUDA_BOTH inline void extrapolateDensityInTime(dtype &rho, dtype rho0, dtype vx, dtype vy, dtype vz,
            dtype dt, dtype gradRhoX, dtype gradRhoY, dtype gradRhoZ, dtype divv) {
        rho -= dt*(vx*gradRhoX + vy*gradRhoY + vz*gradRhoZ + rho0*divv);
    }

    PP_CUDA_BOTH inline void extrapolateVelocityInTime(dtype &v, dtype &vFrame, dtype divv, dtype v0, dtype dt, dtype acc,
            dtype gradP, dtype rho0,
            dtype a) {
        dtype temp;
        temp = -(v0*divv + gradP/rho0)*dt;
        temp += acc*dt*a;
        v += temp;
        vFrame += 0.5*temp;
    }

    PP_CUDA_BOTH inline void extrapolatePressureInTime(dtype &p, dtype dt, dtype gradPx, dtype gradPy, dtype gradPz,
            dtype p0, dtype vx0, dtype vy0, dtype vz0, dtype divv, dtype dConstGamma) {
        p -= dt*( vx0*gradPx + vy0*gradPy + vz0*gradPz  + dConstGamma*p0*divv);
    }

    PP_CUDA_BOTH inline void extrapolateStateInTime(
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

    PP_CUDA_BOTH inline void extrapolateVariableInSpace(dtype &var, dtype dx,dtype dy, dtype dz,
            dtype gradx, dtype grady, dtype gradz) {
        var += dx*gradx + dy*grady + dz*gradz;
    }

    PP_CUDA_BOTH inline void extrapolateVelocityCosmology(dtype &v, dtype &vFrame, dtype v0, dtype dt, dtype H, dtype a) {
        dtype temp = H * dt * a * v;
        v -= temp;
        vFrame -= 0.5*temp;
    }

    PP_CUDA_BOTH inline void extrapolateCosmology(dtype &vx, dtype &vy, dtype &vz,
            dtype &vFramex, dtype &vFramey, dtype &vFramez, dtype &p,
            dtype dt, dtype vx0, dtype vy0, dtype vz0, dtype p0,
            dtype dConstGamma, dtype H, dtype a) {
        extrapolateVelocityCosmology( vx,  vFramex,  vx0,  dt,  H,  a);
        extrapolateVelocityCosmology( vy,  vFramey,  vy0,  dt,  H,  a);
        extrapolateVelocityCosmology( vz,  vFramez,  vz0,  dt,  H,  a);

        p -= 3. * H * dt * a * (dConstGamma - 1.) * p0;
    }

    PP_CUDA_BOTH inline void computeFace(dtype &modApq, std::array<dtype,3> &unit,
                                         dtype rpq,  dtype dx, dtype dy, dtype dz,
                                         dtype pH, dtype qH, dtype p_omega, dtype q_omega,
                                         dtype pBxx, dtype pBxy, dtype pBxz,
                                         dtype pByy, dtype pByz, dtype pBzz,
                                         dtype qBxx, dtype qBxy, dtype qBxz,
                                         dtype qByy, dtype qByz, dtype qBzz) {
        dtype psi;
        dtype psiTilde[3];
        dtype Apq[3] = {0.,0.,0.};

        // \tilde{\psi}_j (x_i)
        psi  = -cubicSplineKernel<dtype,mtype>(rpq, pH)/p_omega;
        psiTilde[0] = (pBxx*dx + pBxy*dy + pBxz*dz)*psi;
        psiTilde[1] = (pBxy*dx + pByy*dy + pByz*dz)*psi;
        psiTilde[2] = (pBxz*dx + pByz*dy + pBzz*dz)*psi;
        for (auto j=0; j<3; j++) {
            Apq[j] += psiTilde[j]/p_omega;
        }

        // \tilde{\psi}_i (x_j)
        psi = cubicSplineKernel<dtype,mtype>(rpq, qH)/q_omega;
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
        // Avoid NaNs when the face has zero area
        dtype zero = 0.;
        mtype cond = modApq==zero;
        if (movemask(cond)) {
            for (auto j=0; j<3; j++)
                unit[j] = mask_mov(unit[j], cond, zero);
        }
    }

    PP_CUDA_BOTH inline void low_limit(dtype &var, dtype def, dtype min) {
        mtype cond = var<min;
        var = mask_mov(var, cond, def);
    }

public:

    PP_CUDA_BOTH inline ResultMFM<dtype> doSinglePPFlux(const mtype mask,
            bool bComove, dtype dTime, dtype dDelta, dtype a, dtype H, dtype dConstGamma,
            dtype rpq, dtype dx, dtype dy, dtype dz,
            dtype pBall, dtype pLastUpdateTime, dtype pDt,
            dtype pOmega, dtype pSPHEoM,
            dtype pBxx, dtype pBxy, dtype pBxz,
            dtype pByy, dtype pByz, dtype pBzz,
            dtype pDensity, dtype pVpredx, dtype pVpredy, dtype pVpredz, dtype pP, dtype pS,
            // dtype pGradRhoX, dtype pGradRhoY, dtype pGradRhoZ,
            // dtype pGradPX, dtype pGradPY, dtype pGradPZ,
            // dtype pGradVxX, dtype pGradVxY, dtype pGradVxZ,
            // dtype pGradVyX, dtype pGradVyY, dtype pGradVyZ,
            // dtype pGradVzX, dtype pGradVzY, dtype pGradVzZ,
            // dtype pLastAccX, dtype pLastAccY, dtype pLastAccZ,
            dtype qBall, dtype qLastUpdateTime, dtype qDt,
            dtype qOmega, dtype qSPHEoM,
            dtype qBxx, dtype qBxy, dtype qBxz,
            dtype qByy, dtype qByz, dtype qBzz,
            dtype qDensity, dtype qVpredx, dtype qVpredy, dtype qVpredz, dtype qP, dtype qS,
            // dtype qGradRhoX, dtype qGradRhoY, dtype qGradRhoZ,
            // dtype qGradPX, dtype qGradPY, dtype qGradPZ,
            // dtype qGradVxX, dtype qGradVxY, dtype qGradVxZ,
            // dtype qGradVyX, dtype qGradVyY, dtype qGradVyZ,
            // dtype qGradVzX, dtype qGradVzY, dtype qGradVzZ,
            // dtype qLastAccX, dtype qLastAccY, dtype qLastAccZ,
            struct eEOSparam &eEOS) {

        ResultMFM<dtype> out;
        dtype z{0.0};
        out.minDt = z;
        out.Frho = z;
        out.Fmom[0] = z;
        out.Fmom[1] = z;
        out.Fmom[2] = z;
        out.Fene = z;
        out.FS = z;

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
        out.minDt = min(pDt, qDt);
        out.minDt /=  a;
        // printf("CC \t");
        // dump(out.minDt);
        // dump(pDt);
        // dump(qDt);
        // printf("---\n");

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


        /*
        // Divergence of the velocity field for the forward in time prediction
        dtype pdivv = (pGradVxX + pGradVyY + pGradVzZ);
        dtype qdivv = (qGradVxX + qGradVyY + qGradVzZ);

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
        */

        dx = -dx;
        dy = -dy;
        dz = -dz;

        /*
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
        */

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
            jeansPressureFloor(L_rho, pH, smf->dConstGamma, smf->dEOSNJeans);
        const double pRjeans =
            jeansPressureFloor(R_rho, q(ball), smf->dConstGamma, smf->dEOSNJeans);
        L_p = max(L_p, pLjeans);
        R_p = max(R_p, pRjeans);
        #endif
        */

        dtype P_M, S_M;

        RiemannSolverExact<dtype,mtype> riemann(dConstGamma, mask, bMFV);
        // int niter = riemann.hllc(
        //                 R_rho, R_p, R_v,
        //                 L_rho, L_p, L_v,
        //                 P_M, S_M,
        //                 &out.Frho, &out.Fene, out.Fmom,
        //                 face_unit.data());
        int niter = riemann.solve(
                        R_rho, R_p, R_v,
                        L_rho, L_p, L_v,
                        P_M, S_M,
                        &out.Frho, &out.Fene, out.Fmom,
                        face_unit.data());

        /*
        int nan;
        // Only works if compiling with -fno-finite-math-only !!
        auto aa = maskz_mov(mask, P_M);
        auto bb = maskz_mov(mask, S_M);
        auto cc = maskz_mov(mask, modApq);
        nan = nan_guard(aa, zero) | nan_guard(bb, zero) | nan_guard(cc, zero);
        if (nan) {
        */
        // printf("-----\n");
        // dump(S_M); // nan
        // dump(P_M);
        // dump(R_rho);
        // dump(L_rho);
        // dump(R_p);
        // dump(L_p);
        // dump(R_v[0]); // nan
        // dump(L_v[0]);
        // dump(R_v[1]); // nan
        // dump(L_v[1]);
        // dump(R_v[2]); // nan
        // dump(L_v[2]);
        // dump(face_unit[0]);
        // dump(face_unit[1]);
        // dump(face_unit[2]);
        // dump(modApq);
        //}

#ifdef ENTROPY_SWITCH

        if (!bMFV) {
            // As we are in a truly lagrangian configuration,
            // there is no adtypetion of entropy among particles.
            F_S=0.;
        }
        else {
            /*
            for (auto j=0; j<3; j++) F_S += Fmom[j]*face_unit[j];
            if (F_S > 0) {
                // Maybe this values should be properly extrapolated to the faces..
                // but this is expensive!
                F_S *= pS*pOmega;
            }
            else {
                F_S *= qS*qOmega;
            }
            F_S *= Frho*modApq;
            */
            F_S = 0.;
        }
#endif //ENTROPY_SWITCH

        if (!bMFV)  { // MFM
            out.Frho = 0.;
            out.Fene = P_M * S_M;
            for (auto j=0; j<3; j++)
                out.Fmom[j] = P_M * face_unit[j];
        }

        // Force 2D
#ifdef FORCE_1D
        out.Fmom[2] = 0.;
        out.Fmom[1] = 0.;
#endif
#ifdef FORCE_2D
        out.Fmom[2] = 0.;
#endif

        // Now we de-boost the fluxes following Eq. A8 Hopkins 2015
        for (auto j=0; j<3; j++) {
            out.Fene += vFrame[j] * out.Fmom[j];
            out.Fene += (0.5*vFrame[j]*vFrame[j])*out.Frho;
        }

        // Now we just multiply by the face area
        out.Fene *= modApq;
        out.Frho *= modApq;
        for (auto j=0; j<3; j++) {
            out.Fmom[j] *= modApq;
            out.Fmom[j] += vFrame[j]*out.Frho;
        }

#ifdef MFM_FLUX_CORRECTION
        if (!bMFV)  { // MFM
            const dtype cs_p{sqrt(dConstGamma *pP/pDensity)};
            const dtype cs_q{sqrt(dConstGamma *qP/qDensity)};
            mtype cond{abs(S_M/min(cs_p,cs_q)) < 0.6};
            if (movemask(cond)) {
                const dtype dU_mfm_p = pVpredx*Fmom[0] + pVpredy*Fmom[1] + pVpredz*Fmom[2] - out.Fene;
                const dtype dU_mfm_q = out.Fene - qVpredx*Fmom[0] - qVpredy*Fmom[1] - qVpredz*Fmom[2];

                const dtype dU_sph_fac = 2. * P_M * ((qVpredx-pVpredx)*dx + (qVpredy-pVpredy)*dy + (qVpredz-pVpredz)*dz) / rpq;
                const dtype dU_sph_p = dU_sph_fac * pSPHEoM * cubicSplineKernelDeriv(rpq, pBall) / pOmega;
                const dtype dU_sph_q = dU_sph_fac * qSPHEoM * cubicSplineKernelDeriv(rpq, qBall) / qOmega;

                mtype pIsCooler{cs_p < cs_q};
                mtype pDoCorr{cond &pIsCooler &(dU_mfm_p < 0.0) &(dU_sph_p > dU_mfm_p)};
                mtype qDoCorr{cond &~pIsCooler &(dU_mfm_q < 0.0) &(dU_sph_q > dU_mfm_q)};

                out.Fene = mask_mov(out.Fene, pDoCorr, out.Fene + dU_mfm_p - dU_sph_p);
                out.Fene = mask_mov(out.Fene, qDoCorr, out.Fene - dU_mfm_q + dU_sph_q);
            }
        }
#endif


        // If I uncomment this there are no nan..
        //printf("----\n");
        //printf("%d\n", niter);
        // dump(modApq);
        // dump(rpq);
        // dump(out.Fene);
        // dump(out.Frho);
        // dump(out.Fmom[0]);
        // dump(out.Fmom[1]);
        // dump(out.Fmom[2]);

        /*
        assert(!nan_guard(minDt, zero));
        assert(!nan_guard(Frho, zero));
        assert(!nan_guard(Fene, zero));
        assert(!nan_guard(Fmom[0], zero));
        assert(!nan_guard(Fmom[1], zero));
        assert(!nan_guard(Fmom[2], zero));
        */

        // Phew! Done ;)
        return out;
    }


// Simple macro to improve readability
#define q(X)    input_buffer[q_##X * nBuff + i]
#define qout(X) output_buffer[out_##X * nBuff + i]
    void hydroRiemann(particleStore::ParticleReference &P,float fBall,int nSmooth, int nBuff,
                      decorators::myreal *restrict input_buffer,
                      decorators::myreal *restrict output_buffer,
                      bool bComove, double dTime, double dDelta, double a, double H,
                      double dConstGamma, struct eEOSparam &eEOS
                     );
};

#endif
