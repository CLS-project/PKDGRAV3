#include <stdio.h>
#include "core/particle.h"
#include "solver.h"


#ifdef OPTIM_FLUX_VEC


// Simple macro to improve readability
#define q(X)    input_buffer[q_##X * nBuff + i]
#define qout(X) output_buffer[out_##X * nBuff + i]
template <typename dtype, typename mtype>
void MeshlessHydroSolver<dtype, mtype>::hydroRiemann(particleStore::ParticleReference &P,float fBall,int nSmooth, int nBuff,
        decorators::myreal *restrict input_buffer,
        decorators::myreal *restrict output_buffer,
        bool bComove, double dTime, double dDelta, double a, double H,
        double dConstGamma, struct eEOSparam &eEOS) {
// TODO: this only compiles correctly when using dvec, because it is used in the
// "old" implementation. In the new use we use fvec and thus we can conditionally
// omit this function.
    if constexpr (std::is_same_v<dtype, dvec>) {

        const auto &pv = P.velocity();
        auto &psph = P.sph();

        dtype pomega = psph.omega;
        dtype pspheom= psph.spheom;
        dtype pH     = P.ball();
        dtype plast  = psph.lastUpdateTime;
        dtype pDt    = dDelta/(1<<P.rung());
        dtype pBXX   = psph.B[XX];
        dtype pBXY   = psph.B[XY];
        dtype pBXZ   = psph.B[XZ];
        dtype pBYY   = psph.B[YY];
        dtype pBYZ   = psph.B[YZ];
        dtype pBZZ   = psph.B[ZZ];
        dtype pDens  = P.density();
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
        for (auto i=0; i<dtype::width(); i++)
            index[i] = i;
        dtype index_v;
        index_v.load(index);

#pragma forceinline
        assert(nBuff>(nSmooth+dtype::width()-1));
        for (auto i=0; i<nSmooth; i+=dtype::width()) {
            dtype F_rho;
            dtype F_v[3];
            dtype F_P;
            dtype F_S;
            dtype minDt;

            const mtype mask = index_v<(double)nSmooth;
            index_v += (double)dtype::width();

            qdr.load(       &q(dr));
            qdx.load(       &q(dx));
            qdy.load(       &q(dy));
            qdz.load(       &q(dz));
            qomega.load(    &q(omega));
            qspheom.load(   &q(spheom));
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

            if (dDelta <= 0.0) {
                pDt = 0.;
                qDt = 0.;
                plast = 0.;
                qlast = 0.;
                dTime = 0.;
            }

            auto out = doSinglePPFlux( mask,
                                       bComove, dTime, dDelta,  a,  H,  dConstGamma,
                                       qdr,  qdx,  qdy,  qdz,
                                       pH, plast, pDt,
                                       pomega, pspheom,
                                       pBXX,  pBXY,  pBXZ,
                                       pBYY,  pBYZ,  pBZZ,
                                       pDens,  pVpredx, pVpredy, pVpredz,  pPres,  pS,
                                       // pgradRhox, pgradRhoy, pgradRhoz,
                                       // pgradPx,   pgradPy,   pgradPz,
                                       // pgradVxx,  pgradVxy,  pgradVxz,
                                       // pgradVyx,  pgradVyy,  pgradVyz,
                                       // pgradVzx,  pgradVzy,  pgradVzz,
                                       // plastAccx, plastAccy, plastAccz,
                                       qh, qlast, qDt,
                                       qomega, qspheom,
                                       qBXX,  qBXY,  qBXZ,
                                       qBYY,  qBYZ,  qBZZ,
                                       qDens,  qvx, qvy, qvz,  qP,  qS,
                                       // qgradRhox,  qgradRhoy,  qgradRhoz,
                                       // qgradPx,    qgradPy,    qgradPz,
                                       // qgradVxx,   qgradVxy,   qgradVxz,
                                       // qgradVyx,   qgradVyy,   qgradVyz,
                                       // qgradVzx,   qgradVzy,   qgradVzz,
                                       // qlastAccx,  qlastAccy,  qlastAccz,
                                       eEOS );

            // We fill the output buffer with the fluxes, which then
            // will be added to the corresponding particles
            out.Frho.store(&output_buffer[out_Frho * nBuff + i]);
            out.Fene.store(&output_buffer[out_Fene * nBuff + i]);
            out.Fmom[0].store(&output_buffer[out_FmomX * nBuff + i]);
            out.Fmom[1].store(&output_buffer[out_FmomY * nBuff + i]);
            out.Fmom[2].store(&output_buffer[out_FmomZ * nBuff + i]);
#ifdef ENTROPY_SWITCH
            out.F_S.store(&output_buffer[out_FS * nBuff + i]);
#endif
            out.minDt.store(&output_buffer[out_minDt * nBuff + i]);

        }
    }
}
#endif // OPTIM_FLUX_VEC

// Explicit instantiation of the classes for the different vector types
template class MeshlessHydroSolver<dvec, dmask>;
template class MeshlessHydroSolver<fvec, fmask>;
