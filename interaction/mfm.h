#ifndef INTERACTION_HYDRO_H
#define INTERACTION_HYDRO_H
#include <stack>
#include "baseinteraction.h"
#include "basetype.h"
#include "core/simd.h"
#include "core/vmath.h"
#include "interaction/ilhydro.h"
#include "meshless/hydro/solver.h"
#include "meshless/hydro/hydro.h"
#include "pkd.h"
#include "pkd_config.h"
#include "walk/stats.h"

#define XX 0
#define YY 3
#define ZZ 5
#define XY 1
#define XZ 2
#define YZ 4

/*
 * INTERACTION TRAITS
 */
class HydroInteractionTraits {
public:
    struct name {
        static constexpr std::string_view value = "HydroInteraction";
    };

    struct Options {
        double a;
        double H;
        double gamma;
        double dTime;
        double dDelta;

        int bPeriodic;
        int bGPU;

        bool doExtensiveTest;
        bool bComove;

        uint8_t uRungLo;
        uint8_t uRungHi;
    };

    typedef struct PINFOIN {
        blitz::TinyVector<float, 3> r;
        blitz::TinyVector<float, 3> a;
        blitz::TinyVector<float, 3> v;
        float Omega;
        float fBall;
        float isGas;
        float rho;
        float P;
        float fSmooth2;
        float fDensity;
        float spheom;
        float H;
        float last;
        float Dt;
        float BXX;
        float BXY;
        float BXZ;
        float BYY;
        float BYZ;
        float BZZ;
        float S;
        float gradRhox;
        float gradRhoy;
        float gradRhoz;
        float gradPx;
        float gradPy;
        float gradPz;
        float gradVxx;
        float gradVxy;
        float gradVxz;
        float gradVyx;
        float gradVyy;
        float gradVyz;
        float gradVzx;
        float gradVzy;
        float gradVzz;
        float lastAccx;
        float lastAccy;
        float lastAccz;
        float isActive;
        //float isMarked; // This is only for testing, can be removed
        void init(particleStore::ParticlePointer p) {
            isGas = p->is_gas();
            if (!isGas) return;
            const auto &psph = p->sph();
            a = p->acceleration();
            fSmooth2 = fBall * fBall;
            fDensity = p->density();
            fBall = p->ball();
            isActive = p->is_active() ? 1.0f : 0.0f;
            //isMarked = p->is_marked() ? 1.0f : 0.0f;
            Omega = psph.omega;
            rho = p->density();
            spheom = psph.spheom;
            last = psph.lastUpdateTime;
            BXX = psph.B[XX];
            BXY = psph.B[XY];
            BXZ = psph.B[XZ];
            BYY = psph.B[YY];
            BYZ = psph.B[YZ];
            BZZ = psph.B[ZZ];
            v = p->velocity();
            P = psph.P;
            S = 0.0f; // psph.S;
            gradRhox = psph.gradRho[0];
            gradRhoy = psph.gradRho[1];
            gradRhoz = psph.gradRho[2];
            gradPx = psph.gradP[0];
            gradPy = psph.gradP[1];
            gradPz = psph.gradP[2];
            gradVxx = psph.gradVx[0];
            gradVxy = psph.gradVx[1];
            gradVxz = psph.gradVx[2];
            gradVyx = psph.gradVy[0];
            gradVyy = psph.gradVy[1];
            gradVyz = psph.gradVy[2];
            gradVzx = psph.gradVz[0];
            gradVzy = psph.gradVz[1];
            gradVzz = psph.gradVz[2];
            lastAccx = psph.lastAcc[0];
            lastAccy = psph.lastAcc[1];
            lastAccz = psph.lastAcc[2];
        };
    } PINFOIN;

    typedef struct PINFOOUT {
        float Frho, Fmom[3], Fene, FS;
        float minDt;

        void init(particleStore::ParticlePointer p) {
            Frho = 0.0f;
            Fmom[0] = 0.0f;
            Fmom[1] = 0.0;
            Fmom[2] = 0.0;
            Fene = 0.0f;
            FS = 0.0;
            minDt = 0.0f;
        }
    } PINFOOUT;

    using IL = ilhydroList;

    // In this case, the result class is defined outside, in a header
    // that defines the operations of this kernel, but we maintain
    // the typename for consistency
    template <class F>
    using ResultInteraction = struct ResultMFM<F>;
};

/*
 * INTERACTION EVALUATION
 */

template <class F, class M>
/*PP_CUDA_BOTH*/ HydroInteractionTraits::ResultInteraction<F>
EvalMFM(F fx, F fy, F fz, F Ph, // Particle
        F Plast, F PDt, F Pomega, F Pspheom, F PBXX, F PBXY, F PBXZ, F PBYY,
        F PBYZ, F PBZZ, F PDens, F Pvx, F Pvy, F Pvz, F PP, F PS,
        // F PgradRhox, F PgradRhoy, F PgradRhoz, F PgradPx,
        // F PgradPy, F PgradPz, F PgradVxx, F PgradVxy,
        // F PgradVxz, F PgradVyx, F PgradVyy, F PgradVyz,
        // F PgradVzx, F PgradVzy, F PgradVzz, F PlastAccx,
        // F PlastAccy, F PlastAccz,
        F PisGas, F PisActive, F Idx, F Idy, F Idz,
        F Ih, // Neighbour
        F Ilast, F IDt, F Iomega, F Ispheom, F IBXX, F IBXY, F IBXZ, F IBYY,
        F IBYZ, F IBZZ, F IDens, F Ivx, F Ivy, F Ivz, F IP, F IS,
        // F IgradRhox, F IgradRhoy, F IgradRhoz, F IgradPx,
        // F IgradPy, F IgradPz, F IgradVxx, F IgradVxy,
        // F IgradVxz, F IgradVyx, F IgradVyy, F IgradVyz,
        // F IgradVzx, F IgradVzy, F IgradVzz, F IlastAccx,
        // F IlastAccy, F IlastAccz,
        F IisGas, F IisActive, bool bComove, float dTime,
        float dDelta, // Options
        float a, float H, float gamma) {
    MeshlessHydroSolver<F, M> solver(false);
    M isGasInteraction = PisGas > 0 & IisGas > 0;

    F dx = fx + Idx;
    F dy = fy + Idy;
    F dz = fz + Idz;
    F dr2 = dx * dx + dy * dy + dz * dz;
    F q = dr2 / (Ph * Ph);
    M mask = isGasInteraction & q < 1 & dr2 / (Ih * Ih) < 1;
    // For non-active particles (that are marked) we only compute the
    // interaction with active particles
    mask &= (PisActive > 0 | IisActive > 0);
    mask &= dr2 > 0;
    //printf("mask %d \t %d %d \t %d %d %d\n", movemask(mask), movemask(PisGas>0),  movemask(IisGas>0), movemask(q < 1), movemask(dr2/(Ih*Ih) <1), movemask(dr2>0));
    if (movemask(mask) == 0) {
        HydroInteractionTraits::ResultInteraction<F> z;
        z.zero();
        return z;
    }
    // Make sure that we would eventually compute the other side of this
    // interaction (isMarked is not really needed anywhere else)
    // M marked = IisMarked > 0;
    // if (movemask(mask & ~marked)>0) {
    //     printf("--");
    //     dump(marked);
    //     dump(mask);
    //     assert(0);
    // }
    F dr = sqrt(dr2);
    struct eEOSparam placeholder;
    auto iresult = solver.doSinglePPFlux(
                       mask, bComove, dTime, dDelta, a, H, gamma, dr, dx, dy, dz, Ph, Plast,
                       PDt, Pomega, Pspheom, PBXX, PBXY, PBXZ, PBYY, PBYZ, PBZZ, PDens, Pvx,
                       Pvy, Pvz, PP, PS,
                       // PgradRhox, PgradRhoy,
                       // PgradRhoz, PgradPx, PgradPy, PgradPz, PgradVxx, PgradVxy,
                       // PgradVxz, PgradVyx, PgradVyy, PgradVyz, PgradVzx, PgradVzy,
                       // PgradVzz, PlastAccx, PlastAccy, PlastAccz,
                       Ih, Ilast, IDt, Iomega, Ispheom, IBXX, IBXY, IBXZ, IBYY, IBYZ, IBZZ,
                       IDens, Ivx, Ivy, Ivz, IP, IS,
                       // IgradRhox, IgradRhoy, IgradRhoz,
                       // IgradPx, IgradPy, IgradPz, IgradVxx, IgradVxy, IgradVxz,
                       // IgradVyx, IgradVyy, IgradVyz, IgradVzx, IgradVzy, IgradVzz,
                       // IlastAccx, IlastAccy, IlastAccz,
                       placeholder);

    iresult.Frho = maskz_mov(mask, iresult.Frho);
    iresult.Fmom[0] = maskz_mov(mask, iresult.Fmom[0]);
    iresult.Fmom[1] = maskz_mov(mask, iresult.Fmom[1]);
    iresult.Fmom[2] = maskz_mov(mask, iresult.Fmom[2]);
    iresult.Fene = maskz_mov(mask, iresult.Fene);

    iresult.Frho *= iresult.minDt;
    iresult.Fmom[0] *= iresult.minDt;
    iresult.Fmom[1] *= iresult.minDt;
    iresult.Fmom[2] *= iresult.minDt;
    iresult.Fene *= iresult.minDt;
    iresult.FS *= iresult.minDt;

    // printf("----\n");
    // dump(iresult.Fene);
    // dump(iresult.Frho);
    // dump(iresult.Fmom[0]);
    // dump(iresult.Fmom[1]);
    // dump(iresult.Fmom[2]);
    // dump(iresult.minDt);

    // for (int index = 0; index < iresult.minDt.width(); index++) {
    //     if (!isfinite(iresult.Frho[index]) ||
    //     !isfinite(iresult.Fmom[0][index]) ||
    //     !isfinite(iresult.Fmom[1][index]) ||
    //     !isfinite(iresult.Fmom[2][index]) || !isfinite(iresult.Fene[index])
    //     || !isfinite(iresult.minDt[index]) ) {
    //         printf("An error occured in doSinglePPFlux:\n\
    //         Frho = %.5e, Fmom[0] = %.5e, Fmom[1] = %.5e, Fmom[2] = %.5e, Fene
    //         = %.5e\n\
    //         dr = %.5e, pDt = %.5e, qDt = %.5e\n\
    //         pP = %.5e, qP = %.5e, pDens = %.5e, qDens = %.5e\n"
    //         ,iresult.Frho[index],iresult.Fmom[0][index],iresult.Fmom[1][index],iresult.Fmom[2][index],iresult.Fene[index],
    //         dr[index], pDt[index], IDt.v[i][index], pP[index],
    //         IP.v[i][index], pDens[index], IDens.v[i][index]
    //         );
    //     }
    // }
    return iresult;
}

template <typename BLOCK>
struct ilist::EvalBlock<HydroInteractionTraits::ResultInteraction<fvec>,
           BLOCK> {
    typedef HydroInteractionTraits::ResultInteraction<fvec> result_type;
    const fvec fx, fy, fz, ph, plast, pDt, pomega, pspheom, pBXX, pBXY, pBXZ,
          pBYY, pBYZ, pBZZ, pDens, pvx, pvy, pvz, pP, pS, pgradRhox, pgradRhoy,
          pgradRhoz, pgradPx, pgradPy, pgradPz, pgradVxx, pgradVxy, pgradVxz,
          pgradVyx, pgradVyy, pgradVyz, pgradVzx, pgradVzy, pgradVzz, plastAccx,
          plastAccy, plastAccz, pisGas, pisActive;
    const HydroInteractionTraits::Options *const options;
    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy, fvec fz, fvec ph, fvec plast, fvec pDt,
              fvec pomega, fvec pspheom, fvec pBXX, fvec pBXY, fvec pBXZ,
              fvec pBYY, fvec pBYZ, fvec pBZZ, fvec pDens, fvec pvx, fvec pvy,
              fvec pvz, fvec pP, fvec pS, fvec pgradRhox, fvec pgradRhoy,
              fvec pgradRhoz, fvec pgradPx, fvec pgradPy, fvec pgradPz,
              fvec pgradVxx, fvec pgradVxy, fvec pgradVxz, fvec pgradVyx,
              fvec pgradVyy, fvec pgradVyz, fvec pgradVzx, fvec pgradVzy,
              fvec pgradVzz, fvec plastAccx, fvec plastAccy, fvec plastAccz,
              fvec pisGas, fvec pisActive,
              HydroInteractionTraits::Options *options)
        : fx(fx), fy(fy), fz(fz), ph(ph), plast(plast), pDt(pDt),
          pomega(pomega), pspheom(pspheom), pBXX(pBXX), pBXY(pBXY), pBXZ(pBXZ),
          pBYY(pBYY), pBYZ(pBYZ), pBZZ(pBZZ), pDens(pDens), pvx(pvx), pvy(pvy),
          pvz(pvz), pP(pP), pS(pS), pgradRhox(pgradRhox), pgradRhoy(pgradRhoy),
          pgradRhoz(pgradRhoz), pgradPx(pgradPx), pgradPy(pgradPy),
          pgradPz(pgradPz), pgradVxx(pgradVxx), pgradVxy(pgradVxy),
          pgradVxz(pgradVxz), pgradVyx(pgradVyx), pgradVyy(pgradVyy),
          pgradVyz(pgradVyz), pgradVzx(pgradVzx), pgradVzy(pgradVzy),
          pgradVzz(pgradVzz), plastAccx(plastAccx), plastAccy(plastAccy),
          plastAccz(plastAccz), pisGas(pisGas), pisActive(pisActive),
          options(options) {}

    result_type operator()(int n, BLOCK &blk) {
        // Sentinel values
        while (n & fvec::mask()) {
            // This is a little trick to speed up the calculation.
            // We set this values such that interactions at the end of simd
            // vector result in zero net forces
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e14f;
            blk.isGas.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i = 0; i < n; ++i) {

            result += EvalMFM<fvec, fmask>(
                          fx, fy, fz, ph, plast, pDt, pomega, pspheom, // Particle
                          pBXX, pBXY, pBXZ, pBYY, pBYZ, pBZZ, pDens, pvx, pvy, pvz, pP,
                          pS,
                          // pgradRhox, pgradRhoy, pgradRhoz, pgradPx, pgradPy, pgradPz,
                          // pgradVxx, pgradVxy, pgradVxz, pgradVyx, pgradVyy, pgradVyz,
                          // pgradVzx, pgradVzy, pgradVzz, plastAccx, plastAccy,
                          // plastAccz,
                          pisGas, pisActive, blk.dx.v[i], blk.dy.v[i],
                          blk.dz.v[i], blk.qh.v[i], blk.qlast.v[i], blk.qDt.v[i],
                          blk.qomega.v[i], blk.qspheom.v[i], // Particle
                          blk.qBXX.v[i], blk.qBXY.v[i], blk.qBXZ.v[i], blk.qBYY.v[i],
                          blk.qBYZ.v[i], blk.qBZZ.v[i], blk.qDens.v[i], blk.qvx.v[i],
                          blk.qvy.v[i], blk.qvz.v[i], blk.qP.v[i], blk.qS.v[i],
                          // blk.qgradRhox.v[i], blk.qgradRhoy.v[i], blk.qgradRhoz.v[i],
                          // blk.qgradPx.v[i], blk.qgradPy.v[i], blk.qgradPz.v[i],
                          // blk.qgradVxx.v[i], blk.qgradVxy.v[i], blk.qgradVxz.v[i],
                          // blk.qgradVyx.v[i], blk.qgradVyy.v[i], blk.qgradVyz.v[i],
                          // blk.qgradVzx.v[i], blk.qgradVzy.v[i], blk.qgradVzz.v[i],
                          // blk.qlastAccx.v[i], blk.qlastAccy.v[i], blk.qlastAccz.v[i],
                          blk.isGas.v[i], blk.isActive.v[i],
                          options->bComove, options->dTime, options->dDelta, options->a,
                          options->H, options->gamma);
        }
        return result;
    }
};

/*
 * INTERACTION CLASS
 */
template <typename Derived>
class HydroInteractionImpl
    : public BaseInteraction<HydroInteractionTraits, Derived> {
    using base = BaseInteraction<HydroInteractionTraits, Derived>;
    using impl = BaseInteractionInterface<HydroInteractionTraits, Derived>;
    friend impl;

public:
    using IL = typename HydroInteractionTraits::IL;
    using ILtile = typename IL::tile_type;
    using Options = typename HydroInteractionTraits::Options;
    using PINFOIN = typename HydroInteractionTraits::PINFOIN;
    using PINFOOUT = typename HydroInteractionTraits::PINFOOUT;
    using name = typename HydroInteractionTraits::name;
    using ResultInteraction =
        typename HydroInteractionTraits::template ResultInteraction<fvec>;
    using workParticle = typename impl::interWorkParticle;

    HydroInteractionImpl() : base() { this->il = new ilhydroList; };

    ~HydroInteractionImpl() { delete this->il; }

    void addParticleILImpl(double x, double y, double z, float fMass,
                           float fourh2, particleStore::ParticlePointer p) {
        if (p->is_gas()) {
            auto &psph = p->sph();

            float qomega = psph.omega;
            float qspheom = psph.spheom;
            float qH = p->ball();
            float qlast = psph.lastUpdateTime;
            float qDt = this->options->dDelta / (1 << p->rung());
            float qBXX = psph.B[XX];
            float qBXY = psph.B[XY];
            float qBXZ = psph.B[XZ];
            float qBYY = psph.B[YY];
            float qBYZ = psph.B[YZ];
            float qBZZ = psph.B[ZZ];
            float qDens = p->density();
            float qvx = p->velocity()[0];
            float qvy = p->velocity()[1];
            float qvz = p->velocity()[2];
            float qP = psph.P;
            float qS;
#ifdef ENTROPY_SWITCH
            qS = psph.S;
#endif
            float qgradRhox = psph.gradRho[0];
            float qgradRhoy = psph.gradRho[1];
            float qgradRhoz = psph.gradRho[2];
            float qgradPx = psph.gradP[0];
            float qgradPy = psph.gradP[1];
            float qgradPz = psph.gradP[2];
            float qgradVxx = psph.gradVx[0];
            float qgradVxy = psph.gradVx[1];
            float qgradVxz = psph.gradVx[2];
            float qgradVyx = psph.gradVy[0];
            float qgradVyy = psph.gradVy[1];
            float qgradVyz = psph.gradVy[2];
            float qgradVzx = psph.gradVz[0];
            float qgradVzy = psph.gradVz[1];
            float qgradVzz = psph.gradVz[2];
            float qlastAccx = psph.lastAcc[0];
            float qlastAccy = psph.lastAcc[1];
            float qlastAccz = psph.lastAcc[2];
            this->il->append(
                x, y, z, qH, qlast, qDt, qomega, qspheom, qBXX, qBXY, qBXZ,
                qBYY, qBYZ, qBZZ, qDens, qvx, qvy, qvz, qP, qS,
                // qgradRhox, qgradRhoy, qgradRhoz, qgradPx,
                // qgradPy, qgradPz, qgradVxx, qgradVxy, qgradVxz,
                // qgradVyx, qgradVyy, qgradVyz, qgradVzx, qgradVzy,
                // qgradVzz, qlastAccx, qlastAccy, qlastAccz,
                1.0f, p->is_active() ? 1.0f : 0.0f,
                p->is_marked() ? 1.0f : 0.0f);
        }
    }

private:
    void initializeInfoInImpl(workParticle *wp, int nP,
                              particleStore::ParticlePointer p) {
        // We overload this because the positions are now relative
        // to the IL center and floats, we could easily change
        // the constructor of PINFOIN to take the reference frame.
        base::initializeInfoInImpl(wp, nP, p);
        auto r = p->position();
        wp->pInfoIn[nP].r = this->il->get_dr(r);
        wp->pInfoIn[nP].Dt = this->options->dDelta / (1 << p->rung());
    }
    void setupCacheImpl() {
        this->pkd->mdl->CacheInitialize(CID_PARTICLE,std::make_shared<SmoothHelperHydroFluxes>(this->pkd, false));
    }

    bool particleSkipConditionImpl(particleStore::ParticlePointer p) {
        return !p->is_marked();
    }

    void particleWorkDoneImpl(workParticle *wp) {

        // First fill general statistics
        base::particleWorkDoneImpl(wp);

        // Now actually update the particles with the
        // information gathered during this interaction
        auto p = this->pkd->particles[this->pkd->Local()];
        for (auto i = 0; i < wp->nP; i++) {
            p = this->pkd->particles[wp->pPart[i]];
            auto &sph = p.sph();
            auto v = p.velocity();
            assert(isfinite(sph.mom[0]));
            assert(isfinite(wp->pInfoOut[i].Fmom[0]));
            sph.mom[0] -= wp->pInfoOut[i].Fmom[0];
            sph.mom[1] -= wp->pInfoOut[i].Fmom[1];
            sph.mom[2] -= wp->pInfoOut[i].Fmom[2];
            sph.E -= wp->pInfoOut[i].Fene;
            sph.Uint -= (wp->pInfoOut[i].Fene - v[0] * wp->pInfoOut[i].Fmom[0] -
                         v[1] * wp->pInfoOut[i].Fmom[1] -
                         v[2] * wp->pInfoOut[i].Fmom[2]);
            // IA: I do not understand this, but I will keep it
            auto q =
                this->pkd->particles[this->pkd->mdl->template AcquireWrite<PARTICLE>(CID_PARTICLE, wp->iPart[i])];
            q = p;
            this->pkd->mdl->ReleaseWrite(CID_PARTICLE, &q);
        }
    }


    void evalInteractionImpl(const PINFOIN &Part, ILtile &tile, PINFOOUT &Out) {
        ilist::EvalBlock<HydroInteractionTraits::ResultInteraction<fvec>,
              BlockHYDRO<ILHYDRO_PART_PER_BLK>>
              eval(Part.r[0], Part.r[1], Part.r[2], Part.fBall, Part.last,
                   Part.Dt, Part.Omega, Part.spheom, Part.BXX, Part.BXY, Part.BXZ,
                   Part.BYY, Part.BYZ, Part.BZZ, Part.rho, Part.v[0], Part.v[1],
                   Part.v[2], Part.P, Part.S, Part.gradRhox, Part.gradRhoy,
                   Part.gradRhoz, Part.gradPx, Part.gradPy, Part.gradPz,
                   Part.gradVxx, Part.gradVxy, Part.gradVxz, Part.gradVyx,
                   Part.gradVyy, Part.gradVyz, Part.gradVzx, Part.gradVzy,
                   Part.gradVzz, Part.lastAccx, Part.lastAccy, Part.lastAccz,
                   Part.isGas, Part.isActive, this->options);

        auto result = EvalTile(tile, eval);
        Out.Frho += hadd(result.Frho);
        Out.Fmom[0] += hadd(result.Fmom[0]);
        Out.Fmom[1] += hadd(result.Fmom[1]);
        Out.Fmom[2] += hadd(result.Fmom[2]);
        Out.Fene += hadd(result.Fene);
        Out.FS += hadd(result.FS);
    }

    void extensiveILTestImpl(workParticle *wp) {
        std::stack<std::pair<int,int>> cellStack;
        bool dofBallFactor = false;
        float fBallFactor = 1.0;
        float ballSizeLimit = INFINITY;

        // Add the toptree cell corresponding to the ROOT cell to the stack
        cellStack.push(std::make_pair(this->pkd->iTopTree[ROOT],this->pkd->Self()));
        int nParticles = 0;

        while (!cellStack.empty()) {
            auto [iCell, id] = cellStack.top();
            cellStack.pop();
            auto c = (id == this->pkd->Self()) ?
                     this->pkd->tree[iCell] :
                     this->pkd->tree[this->pkd->mdl->template Fetch<KDN>(CID_CELL,iCell,id)];

            // Walking the top tree
            if (c->is_top_tree()) {
                auto [iLower, idLower, iUpper, idUpper] = c->get_child_cells(id);
                cellStack.push(std::make_pair(iUpper, idUpper));
                cellStack.push(std::make_pair(iLower, idLower));
                continue;
            }

            /* We are either on
            ** the toptree cell that corresponds to my local cell
            ** or on one of the 2 children for cell corresponding to remote threads id != pkd->Self()
            ** so we can loop through the particles between c->lower() and c->upper()
            ** and mdlFetch them from that thread
            */
            for (auto pj=c->lower(); pj<=c->upper(); ++pj) {
                auto p = (id == this->pkd->Self()) ?
                         this->pkd->particles[pj] :
                         this->pkd->particles[this->pkd->mdl->template Fetch<PARTICLE>(CID_PARTICLE,pj,id)];

                ++nParticles;

                // Get position and ball size of particle p
                auto pr = p.position();
                float fBallFactor = (dofBallFactor) ? fBallFactor : 1.0f;
                float pBall2 = std::min(ballSizeLimit, p.ball() * fBallFactor);
                pBall2 *= pBall2;

                int needsCheck = 0;
                // Loop over all particles in the wp
                for (auto i=0; i<wp->nP; ++i) {
                    auto q = this->pkd->particles[wp->iPart[i]];

                    // Get position and ball size of particle q
                    auto qr = q.position();
                    float qBall2 = std::min(ballSizeLimit, q.ball() * fBallFactor);
                    qBall2 *= qBall2;

                    // Calculate distance squared between particle p and q
                    auto dist = pr - qr;
                    float dist2 = blitz::dot(dist,dist);

                    // Do check for gather
                    if (dist2 < qBall2) {
                        needsCheck = 1;
                        // printf("%ld: %e %e %e \t %ld: %e %e %e\n", p.ParticleID(), pr[0], pr[1], pr[2], q.ParticleID(), qr[0], qr[1], qr[2]);
                        break;
                    }
                    // Do check for scatter
                    if (dist2 < pBall2) {
                        needsCheck = 1;
                        break;
                    }
                }
                // Now we need to check if particle p is on the ILP list
                if (needsCheck) {

                    // Get distance relative to the reference of the il
                    fvec dx = (float)(this->il->getReference(0) - pr[0]);
                    fvec dy = (float)(this->il->getReference(1) - pr[1]);
                    fvec dz = (float)(this->il->getReference(2) - pr[2]);
                    // if (p.ParticleID() == 128) {
                    //     printf("Looking for\n");
                    //     dump(dx);
                    //     dump(dy);
                    //     dump(dz);
                    //     printf("~~~~");
                    // }

                    fvec occurrence = 0.0f;
                    fmask mask;

                    // Compare every entry of the il to the positions and count the matches
                    for ( auto &tile : *this->il ) {
                        auto nBlocks = tile.count() / tile.width;
                        for (auto iBlock=0; iBlock<=nBlocks; ++iBlock) {
                            int n = (iBlock == nBlocks ?  tile.count() - nBlocks *tile.width : tile.width);
                            auto &blk = tile[iBlock];
                            while (n & fvec::mask()) {
                                blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e14f;
                                ++n;
                            }
                            n /= fvec::width();
                            for (auto i=0; i<n; ++i) {
                                mask = ((blk.dx.v[i] == dx) & (blk.dy.v[i] == dy) & (blk.dz.v[i] == dz));
                                // if (p.ParticleID() == 128) {
                                //     printf("---\n");
                                //     dump(blk.dx.v[i]);
                                //     dump(blk.dy.v[i]);
                                //     dump(blk.dz.v[i]);
                                // }
                                occurrence += maskz_mov(mask,fvec(1.0f));
                            }
                        }
                    }
                    // Sum up the matches and compare to 1
                    float occurrences = hadd(occurrence);
                    // printf("%ld: %e %d\n", p.ParticleID(), occurrences, p.is_marked());
                    assert(occurrences == 1.0f);
                }
            }
        }
        assert(nParticles == this->pkd->nGas);
    }

};
class HydroInteraction : public HydroInteractionImpl<HydroInteraction> {
public:
    explicit HydroInteraction() : HydroInteractionImpl<HydroInteraction>() {}
};
#endif
