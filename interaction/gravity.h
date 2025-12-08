#ifndef INTERACTION_GRAVITY_H
#define INTERACTION_GRAVITY_H
#include "baseinteraction.h"
#include "basetype.h"
#include "core/simd.h"
#include "core/vmath.h"
#include "gravity/qeval.h"
#include "pkd.h"
#include "pkd_config.h"
#include "smooth/smooth.h"
#include "walk/stats.h"

/*
 * INTERACTION TRAITS
 */
class GravityInteractionTraits {
public:
    struct name {
        static constexpr std::string_view value = "GravityInteraction";
    };

    struct Options {
        struct pkdKickParameters kick;
        struct pkdLightconeParameters lc;
        struct pkdTimestepParameters ts;
        double dTime;
        double dEwCut;
        double dEwhCut;
        double dTheta;
        int nReps;
        int bPeriodic;
        int bEwald;
        int bGPU;
        int iRoot1;
        int iRoot2;
        // SPHOptions SPHoptions;
        bool doExtensiveTest;
        uint8_t uRungLo;
        uint8_t uRungHi;
    };

    typedef struct PINFOIN {
        blitz::TinyVector<float, 3> r;
        blitz::TinyVector<float, 3> a;
        float fSmooth2;
        float fDensity;
        void init(particleStore::ParticlePointer p) {};
    } PINFOIN;

    typedef struct PINFOOUT {
        blitz::TinyVector<float, 3> a;
        float fPot;
        float dirsum, normsum;
        float rhopmax;

        void init(particleStore::ParticlePointer p, float dirLsum,
                  float normLsum) {
            a = 0.0f;
            fPot = 0.0f;
            dirsum = dirLsum;
            normsum = normLsum;
            rhopmax = 0.0f;
        }
    } PINFOOUT;

    using IL = ilpList;
    using ILC = ilcList;

    template <class F> struct ResultInteraction {};

    template <class F = float> struct ResultPC {
        F ax, ay, az, pot;
        F ir, norm;
        /* PC_CUDA_BOTH */ void zero() { ax = ay = az = pot = ir = norm = 0; }
        /* PC_CUDA_BOTH */ ResultPC<F> operator+=(const ResultPC<F> rhs) {
            ax += rhs.ax;
            ay += rhs.ay;
            az += rhs.az;
            pot += rhs.pot;
            ir += rhs.ir;
            norm += rhs.norm;
            return *this;
        }
    };

    template <class F = float> struct ResultPP {
        typedef F type;
        F ax, ay, az, pot;
        F ir, norm;
        /*PP_CUDA_BOTH*/ void zero() { ax = ay = az = pot = ir = norm = 0; }
        /*PP_CUDA_BOTH*/ ResultPP<F> operator+=(const ResultPP<F> rhs) {
            ax += rhs.ax;
            ay += rhs.ay;
            az += rhs.az;
            pot += rhs.pot;
            ir += rhs.ir;
            norm += rhs.norm;
            return *this;
        }
    };
};

/*
 * INTERACTION EVALUATION
 */
template <class F, class E, class M>
static int evalEwald(const E &ew, const M &mom, F &ax, F &ay, F &az, F &fPot,
                     const F &x, const F &y, const F &z, const F &g0,
                     const F &g1, const F &g2, const F &g3, const F &g4,
                     const F &g5) {
    F onethird = 1.0 / 3.0;

    F xx = 0.5 * x * x;
    F xxx = onethird * xx * x;
    F xxy = xx * y;
    F xxz = xx * z;
    F yy = 0.5 * y * y;
    F yyy = onethird * yy * y;
    F xyy = yy * x;
    F yyz = yy * z;
    F zz = 0.5 * z * z;
    F zzz = onethird * zz * z;
    F xzz = zz * x;
    F yzz = zz * y;
    F xy = x * y;
    F xyz = xy * z;
    F xz = x * z;
    F yz = y * z;
    F Q2mirx = F(mom.xx) * x + F(mom.xy) * y + F(mom.xz) * z;
    F Q2miry = F(mom.xy) * x + F(mom.yy) * y + F(mom.yz) * z;
    F Q2mirz = F(mom.xz) * x + F(mom.yz) * y + F(mom.zz) * z;
    F Q3mirx = F(mom.xxx) * xx + F(mom.xxy) * xy + F(mom.xxz) * xz +
               F(mom.xyy) * yy + F(mom.xyz) * yz + F(mom.xzz) * zz;
    F Q3miry = F(mom.xxy) * xx + F(mom.xyy) * xy + F(mom.xyz) * xz +
               F(mom.yyy) * yy + F(mom.yyz) * yz + F(mom.yzz) * zz;
    F Q3mirz = F(mom.xxz) * xx + F(mom.xyz) * xy + F(mom.xzz) * xz +
               F(mom.yyz) * yy + F(mom.yzz) * yz + F(mom.zzz) * zz;
    F Q4mirx = F(mom.xxxx) * xxx + F(mom.xxxy) * xxy + F(mom.xxxz) * xxz +
               F(mom.xxyy) * xyy + F(mom.xxyz) * xyz + F(mom.xxzz) * xzz +
               F(mom.xyyy) * yyy + F(mom.xyyz) * yyz + F(mom.xyzz) * yzz +
               F(mom.xzzz) * zzz;
    F Q4miry = F(mom.xxxy) * xxx + F(mom.xxyy) * xxy + F(mom.xxyz) * xxz +
               F(mom.xyyy) * xyy + F(mom.xyyz) * xyz + F(mom.xyzz) * xzz +
               F(mom.yyyy) * yyy + F(mom.yyyz) * yyz + F(mom.yyzz) * yzz +
               F(mom.yzzz) * zzz;
    F Q4mirz = F(mom.xxxz) * xxx + F(mom.xxyz) * xxy + F(mom.xxzz) * xxz +
               F(mom.xyyz) * xyy + F(mom.xyzz) * xyz + F(mom.xzzz) * xzz +
               F(mom.yyyz) * yyy + F(mom.yyzz) * yyz + F(mom.yzzz) * yzz +
               F(mom.zzzz) * zzz;
    F Q4x = F(ew.Q4xx) * x + F(ew.Q4xy) * y + F(ew.Q4xz) * z;
    F Q4y = F(ew.Q4xy) * x + F(ew.Q4yy) * y + F(ew.Q4yz) * z;
    F Q4z = F(ew.Q4xz) * x + F(ew.Q4yz) * y + F(ew.Q4zz) * z;
    F Q2mir = 0.5 * (Q2mirx * x + Q2miry * y + Q2mirz * z) -
              (F(ew.Q3x) * x + F(ew.Q3y) * y + F(ew.Q3z) * z) + F(ew.Q4);
    F Q3mir = onethird * (Q3mirx * x + Q3miry * y + Q3mirz * z) -
              0.5 * (Q4x * x + Q4y * y + Q4z * z);
    F Q4mir = 0.25 * (Q4mirx * x + Q4miry * y + Q4mirz * z);
    F Qta =
        g1 * F(mom.m) - g2 * F(ew.Q2) + g3 * Q2mir + g4 * Q3mir + g5 * Q4mir;
    fPot -=
        g0 * F(mom.m) - g1 * F(ew.Q2) + g2 * Q2mir + g3 * Q3mir + g4 * Q4mir;
    ax +=
        g2 * (Q2mirx - F(ew.Q3x)) + g3 * (Q3mirx - Q4x) + g4 * Q4mirx - x * Qta;
    ay +=
        g2 * (Q2miry - F(ew.Q3y)) + g3 * (Q3miry - Q4y) + g4 * Q4miry - y * Qta;
    az +=
        g2 * (Q2mirz - F(ew.Q3z)) + g3 * (Q3mirz - Q4z) + g4 * Q4mirz - z * Qta;

    return COST_FLOP_EWALD;
}

#if defined(USE_SIMD_EWALD) && (defined(__SSE2__) || defined(__ARM_NEON))
// We changed the name to avoid collision with the previous definiton (although
// they are the same)
inline double evalEwaldSIMD2(PKD pkd, ewaldSIMD *ews, dvec &ax, dvec &ay,
                             dvec &az, dvec &dPot, dvec x, dvec y, dvec z,
                             dvec r2, const dmask &doerfc) {
    dvec dir, dir2, a, g0, g1, g2, g3, g4, g5, alphan;
    dvec rerf, rerfc, ex2;
    dvec alpha2x2 = 2.0 * dvec(ews->ewp.alpha2);

    dir = rsqrt(r2);
    dir2 = dir * dir;
    ex2 = exp(-r2 * dvec(ews->ewp.alpha2));
    a = ex2 * dvec(ews->ewp.ka) * dir2;

    verf(dvec(ews->ewp.alpha) * r2 * dir, dvec(ews->ewp.ialpha) * dir, ex2,
         rerf, rerfc);

    g0 = dir * mask_mov(-rerf, doerfc, rerfc);
    g1 = g0 * dir2 + a;
    alphan = alpha2x2;
    g2 = 3.0 * g1 * dir2 + alphan * a;
    alphan *= alpha2x2;
    g3 = 5.0 * g2 * dir2 + alphan * a;
    alphan *= alpha2x2;
    g4 = 7.0 * g3 * dir2 + alphan * a;
    alphan *= alpha2x2;
    g5 = 9.0 * g4 * dir2 + alphan * a;

    evalEwald<dvec, struct ewaldSIMD::PEWALDVARS, struct ewaldSIMD::PMOMC>(
            ews->ewp, ews->ewm, ax, ay, az, dPot, x, y, z, g0, g1, g2, g3, g4, g5);

    return COST_FLOP_EWALD * dvec::width();
}
#endif

template <class F, class M>
/*PP_CUDA_BOTH*/ typename GravityInteractionTraits::template ResultPP<F>
EvalPP(F Pdx, F Pdy, F Pdz, F Psmooth2,       // Particle
       F Idx, F Idy, F Idz, F fourh2, F Im) { // Interaction(s)
    GravityInteractionTraits::ResultPP<F> result;
    constexpr float minSoftening = 1e-18f;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx * dx + dy * dy + dz * dz;
    M vcmp = d2 < fourh2;
    F td2 = max(minSoftening, max(d2, fourh2));
    F dir = rsqrt(td2);
    F dir2 = dir * dir;

    td2 = dir2 * d2; /* for SOFTENED */
    dir2 = dir2 * dir;

    /* dir and dir2 are valid now for both softened and unsoftened particles */
    /* Now we apply the fix to softened particles only */
    if (!testz(vcmp)) {
        td2 = maskz_mov(vcmp, F(1.0f - td2));
        dir *=
            1.0f + td2 * (0.5f + td2 * (3.0f / 8.0f + td2 * (45.0f / 32.0f)));
        dir2 *= 1.0f + td2 * (1.5f + td2 * (135.0f / 16.0f));
    }
    dir2 *= -Im;
    result.ax = dx * dir2;
    result.ay = dy * dir2;
    result.az = dz * dir2;
    result.pot = -Im * dir;
    result.ir = dir;
    result.norm = d2;
    return result;
}
// Calculate additional terms for GravStep
template <class F, class M>
/*PP_CUDA_BOTH*/ GravityInteractionTraits::ResultPP<F>
EvalPP(F Pdx, F Pdy, F Pdz, F Psmooth2,     // Particle
       F Idx, F Idy, F Idz, F fourh2, F Im, // Interaction(s)
       F Pax, F Pay, F Paz, F imaga) {
    GravityInteractionTraits::ResultPP<F> result =
        EvalPP<F, M>(Pdx, Pdy, Pdz, Psmooth2, Idx, Idy, Idz, fourh2, Im);
    F adotai = Pax * result.ax + Pay * result.ay + Paz * result.az;
    adotai =
        maskz_mov((adotai > 0.0f) & (result.norm >= Psmooth2), adotai) * imaga;
    result.norm = adotai * adotai;
    result.ir *= result.norm;
    return result;
}

template <typename BLOCK>
struct ilist::EvalBlock<
    typename GravityInteractionTraits::template ResultPP<fvec>, BLOCK> {
    typedef GravityInteractionTraits::ResultPP<fvec> result_type;
    const fvec fx, fy, fz, pSmooth2, Pax, Pay, Paz, imaga;

    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy, fvec fz, fvec pSmooth2, fvec Pax, fvec Pay,
              fvec Paz, fvec imaga)
        : fx(fx), fy(fy), fz(fz), pSmooth2(pSmooth2), Pax(Pax), Pay(Pay),
          Paz(Paz), imaga(imaga) {}

    result_type operator()(int n, BLOCK &blk) {
        // Sentinal values
        while (n & fvec::mask()) {
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e18f;
            blk.m.s[n] = 0.0f;
            blk.fourh2.s[n] = 1e-18f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i = 0; i < n; ++i) {
            result += EvalPP<fvec, fmask>(
                          fx, fy, fz, pSmooth2, blk.dx.v[i], blk.dy.v[i], blk.dz.v[i],
                          blk.fourh2.v[i], blk.m.v[i], Pax, Pay, Paz, imaga);
        }
        return result;
    }
};

template <class F, class M, bool bGravStep>
/*PC_CUDA_BOTH*/ GravityInteractionTraits::ResultPC<F>
EvalPC(F Pdx, F Pdy, F Pdz, F Psmooth2, // Particle
       F Idx, F Idy, F Idz, F Im, F Iu, // Interaction(s)
       F Ixxxx, F Ixxxy, F Ixxxz, F Ixxyz, F Ixxyy, F Iyyyz, F Ixyyz, F Ixyyy,
       F Iyyyy, F Ixxx, F Ixyy, F Ixxy, F Iyyy, F Ixxz, F Iyyz, F Ixyz, F Ixx,
       F Ixy, F Ixz, F Iyy, F Iyz,
#ifdef USE_DIAPOLE
    F Ix, F Iy, F Iz,
#endif
       F Pax, F Pay, F Paz, F imaga) {
    GravityInteractionTraits::ResultPC<F> result;
    const F onethird = 1.0f / 3.0f;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx * dx + dy * dy + dz * dz;
    F dir = rsqrt(d2);
    F x = dx * dir;
    F y = dy * dir;
    F z = dz * dir;
    {
        // Order 0
        F g = -dir * Im;
        result.pot = g;
        result.ax = x * g;
        result.ay = y * g;
        result.az = z * g;
    }
    F u = Iu * dir;
    F g1 = dir * u;
    F g2 = 3.0f * g1 * u;
    F g3 = 5.0f * g2 * u;
    F g4 = 7.0f * g3 * u;

    F xx = 0.5f * x * x;
    F yy = 0.5f * y * y;
    F zz = 0.5f * z * z;

    // Now calculate the interactions up to Hexadecapole order.
    {
        // Order 4
        F xxx = x * (onethird * xx - zz);
        F tx = Ixxxx * xxx;
        F ty = Ixxxy * xxx;
        F tz = Ixxxz * xxx;
        F xxz = z * (xx - onethird * zz);
        tx += Ixxxz * xxz;
        ty += Ixxyz * xxz;
        tz -= Ixxxx * xxz + Ixxyy * xxz;
        F yyz = z * (yy - onethird * zz);
        tx += Ixyyz * yyz;
        ty += Iyyyz * yyz;
        tz -= Ixxyy * yyz + Iyyyy * yyz;
        F yyy = y * (onethird * yy - zz);
        tx += Ixyyy * yyy;
        ty += Iyyyy * yyy;
        tz += Iyyyz * yyy;

        xx -= zz;
        yy -= zz;

        F xyy = x * yy;
        tx += Ixxyy * xyy;
        ty += Ixyyy * xyy;
        tz += Ixyyz * xyy;
        F xxy = y * xx;
        tx += Ixxxy * xxy;
        ty += Ixxyy * xxy;
        tz += Ixxyz * xxy;
        F xyz = x * y * z;
        tx += Ixxyz * xyz;
        ty += Ixyyz * xyz;
        tz -= (Ixyyy + Ixxxy) * xyz;

        result.ax += (tx *= g4);
        result.ay += (ty *= g4);
        result.az += (tz *= g4);
        F g = 0.25f * (tx * x + ty * y + tz * z);
        result.pot -= g;
        result.ax -= x * 9.0f * g;
        result.ay -= y * 9.0f * g;
        result.az -= z * 9.0f * g;
    }
    {
        // Order 3
        F tx = Ixxx * xx + Ixyy * yy;
        F ty = Ixxy * xx + Iyyy * yy;
        F tz = Ixxz * xx + Iyyz * yy;
        F xy = x * y;
        tx += Ixxy * xy;
        ty += Ixyy * xy;
        tz += Ixyz * xy;
        F xz = x * z;
        tx += Ixxz * xz;
        ty += Ixyz * xz;
        tz -= (Ixxx + Ixyy) * xz;
        F yz = y * z;
        tx += Ixyz * yz;
        ty += Iyyz * yz;
        tz -= (Ixxy + Iyyy) * yz;
        result.ax += (tx *= g3); // xxx
        result.ay += (ty *= g3); // xxy
        result.az += (tz *= g3); // xxz

        F g = onethird * (tx * x + ty * y + tz * z);
        result.pot -= g;
        result.ax -= x * 7.0f * g;
        result.ay -= y * 7.0f * g;
        result.az -= z * 7.0f * g;
    }
    {
        // Order 2
        F t, g;
        result.ax += (t = g2 * (Ixx * x + Ixy * y + Ixz * z)); // xx
        g = t * x;
        result.ay += (t = g2 * (Iyy * y + Ixy * x + Iyz * z)); // xy
        g += t * y;
        result.az += (t = g2 * (-(Ixx + Iyy) * z + Ixz * x + Iyz * y)); // xz
        g += t * z;
        g *= 0.5f;
        result.pot -= g;
        result.ax -= x * 5.0f * g;
        result.ay -= y * 5.0f * g;
        result.az -= z * 5.0f * g;
    }

#ifdef USE_DIAPOLE
    {
        F yy = g1 * Ix;
        result.ax += yy;
        F yz = g1 * Iy;
        result.ay += yz;
        F zz = g1 * Iz;
        result.az += zz;
        F g = yy * x + yz * y + zz * z;
        result.pot -= g;
        result.ax -= x * 3.0f * g;
        result.ay -= y * 3.0f * g;
        result.az -= z * 3.0f * g;
    }
#endif
    result.ax *= dir;
    result.ay *= dir;
    result.az *= dir;

    /* Calculations for determining the timestep. */
    if (bGravStep) {
        F adotai = Pax * result.ax + Pay * result.ay + Paz * result.az;
        adotai = maskz_mov((adotai > 0.0f) & (d2 > Psmooth2), adotai) * imaga;
        result.norm = adotai * adotai;
        result.ir = dir * result.norm;
    }
    return result;
}

template <typename BLOCK>
struct ilist::EvalBlock<GravityInteractionTraits::ResultPC<fvec>, BLOCK> {
    typedef GravityInteractionTraits::ResultPC<fvec> result_type;
    const fvec fx, fy, fz, pSmooth2, Pax, Pay, Paz, imaga;

    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy, fvec fz, fvec pSmooth2, fvec Pax, fvec Pay,
              fvec Paz, fvec imaga)
        : fx(fx), fy(fy), fz(fz), pSmooth2(pSmooth2), Pax(Pax), Pay(Pay),
          Paz(Paz), imaga(imaga) {}

    result_type operator()(int n, BLOCK &blk) {
        // Sentinal values
        while (n & fvec::mask()) {
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e18f;
            blk.m.s[n] = 0.0f;
            blk.u.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        GravityInteractionTraits::ResultPC<fvec> result;
        result.zero();
        for (auto i = 0; i < n; ++i) {
            result += EvalPC<fvec, fmask, true>(
                          fx, fy, fz, pSmooth2, blk.dx.v[i], blk.dy.v[i], blk.dz.v[i],
                          blk.m.v[i], blk.u.v[i], blk.xxxx.v[i], blk.xxxy.v[i],
                          blk.xxxz.v[i], blk.xxyz.v[i], blk.xxyy.v[i], blk.yyyz.v[i],
                          blk.xyyz.v[i], blk.xyyy.v[i], blk.yyyy.v[i], blk.xxx.v[i],
                          blk.xyy.v[i], blk.xxy.v[i], blk.yyy.v[i], blk.xxz.v[i],
                          blk.yyz.v[i], blk.xyz.v[i], blk.xx.v[i], blk.xy.v[i],
                          blk.xz.v[i], blk.yy.v[i], blk.yz.v[i],
#ifdef USE_DIAPOLE
                          blk.x.v[i], blk.y.v[i], blk.z.v[i],
#endif
                          Pax, Pay, Paz, imaga);
        }
        return result;
    }
};

/*
 * INTERACTION CLASS
 */
template <typename Derived>
class GravityInteractionImpl
    : public BaseInteraction<GravityInteractionTraits, Derived> {
    using base = BaseInteraction<GravityInteractionTraits, Derived>;
    using impl = BaseInteractionInterface<GravityInteractionTraits, Derived>;
    friend impl;

public:
    using IL = typename GravityInteractionTraits::IL;
    using Options = typename GravityInteractionTraits::Options;
    using PINFOIN = typename GravityInteractionTraits::PINFOIN;
    using PINFOOUT = typename GravityInteractionTraits::PINFOOUT;
    using name = typename GravityInteractionTraits::name;
    using ResultInteraction =
        typename GravityInteractionTraits::template ResultInteraction<fvec>;
    using workParticle = typename impl::workParticle;
    // Internal state for the interaction, that should be
    // setup during the the walk and may vary between calls
    // to interact
    float dirLsum;
    float normLsum;
    LOCR *pLoc;
    SMX smx;
    SMF smf;

    GravityInteractionImpl() : base() {};

    void setupImpl(PKD pkd, Options *options, WalkStats *walkStats) {
        base::setupImpl(pkd, options, walkStats);
        // We just use the PKD ones, but could be
        // separated in the future
        this->il = &this->pkd->ilp;
        this->ilc = &this->pkd->ilc;
        this->ill = &this->pkd->ill;

        if (this->options->bEwald)
            pkdEwaldInit(this->pkd, this->options->nReps, this->options->dEwCut,
                         this->options->dEwhCut,
                         this->options->bGPU); /* ignored in Flop count! */
    }

    void clearILImpl() {
        this->il->clear();
        ilc->clear();
        // ill->clear();
    };
    void clearILL() { ill->clear(); }
    void setReferenceILImpl(blitz::TinyVector<double, 3> ref) {
        this->il->setReference(ref);
        ilc->setReference(ref);
    };

    // These two functions are required in the special case of the gravity
    // because there are some IL that are filled during the opening criteria
    // of the tree walk
    void addLocalIL(blitz::TinyVector<double, 3> dx, FMOMR *mom, double bMax) {
        ill->append((float)dx[0], (float)dx[1], (float)dx[2], mom, bMax);
    };

    void addCellIL(blitz::TinyVector<double, 3> r, FMOMR *mom, double bMax) {
        ilc->append(r[0], r[1], r[2], mom, bMax);
    }

private:
    ilcList *ilc;
    ilcList *ill;

    void initializeInfoInImpl(workParticle *wp, int nP,
                              particleStore::ParticlePointer p) {
        auto r = p->position();
        wp->pInfoIn[nP].r = this->pkd->ilc.get_dr(r);
        if (this->options->ts.bGravStep) {
            /*
            ** Calculate local density using smooth; this is fast because the
            *particles are
            ** likely to be cached already because they will be on the P-P list.
            */
            smf.pfDensity = &wp->pInfoIn[nP].fDensity;
            auto fBall = smSmoothSingle(smx, &smf, *p, this->options->iRoot1,
                                        this->options->iRoot2);
            wp->pInfoIn[nP].fSmooth2 = fBall * fBall;
        }
        else {
            /*
            ** We are not using GravStep!
            */
            wp->pInfoIn[nP].fSmooth2 = 0.0;
        }
    };

    void initializeInfoOutImpl(workParticle *wp, int nP,
                               particleStore::ParticlePointer p) {
        wp->pInfoOut[nP].init(p, dirLsum, normLsum);
    };

    void queuePC(PKD pkd, workParticle *wp, ilcList &ilc, int bGravStep,
                 bool bGPU = true) {
        for (auto &tile : ilc) {
            ++wp->nTilesTotal;
            if (bGPU) {
                // #ifdef USE_CUDA
                //                 if (pkd->cudaClient->queuePC(wp, tile,
                //                 bGravStep))
                //                     continue;
                // #endif
                // #ifdef USE_METAL
                //                 if (pkd->metalClient->queuePC(wp, tile,
                //                 bGravStep))
                //                     continue;
                // #endif
            }
            ++wp->nTilesCPU;
            for (auto i = 0; i < wp->nP; ++i) {
                pkdGravEvalPC(wp->pInfoIn[i], tile, wp->pInfoOut[i]);
                wp->dFlopSingleCPU += COST_FLOP_PC * tile.size();
            }
        }
    }

    void pkdGravEvalPC(const PINFOIN &Part, ilcTile &tile, PINFOOUT &Out) {
        float a2 = blitz::dot(Part.a, Part.a);
        fvec imaga = a2 > 0.0f ? 1.0f / sqrtf(a2) : 0.0f;
        ilist::EvalBlock<GravityInteractionTraits::ResultPC<fvec>, ilcBlock>
        eval(Part.r[0], Part.r[1], Part.r[2], Part.fSmooth2, Part.a[0],
             Part.a[1], Part.a[2], imaga);
        ;
        auto result = EvalTile(tile, eval);
        Out.a[0] += hadd(result.ax);
        Out.a[1] += hadd(result.ay);
        Out.a[2] += hadd(result.az);
        Out.fPot += hadd(result.pot);
        Out.dirsum += hadd(result.ir);
        Out.normsum += hadd(result.norm);
    }

    void pkdGravEvalPP(const PINFOIN &Part, ilpTile &tile, PINFOOUT &Out) {
        float a2 = blitz::dot(Part.a, Part.a);
        fvec imaga = a2 > 0.0f ? 1.0f / sqrtf(a2) : 0.0f;

        ilist::EvalBlock<GravityInteractionTraits::ResultPP<fvec>, ilpBlock>
        eval(Part.r[0], Part.r[1], Part.r[2], Part.fSmooth2, Part.a[0],
             Part.a[1], Part.a[2], imaga);
        ;
        auto result = EvalTile(tile, eval);
        Out.a[0] += hadd(result.ax);
        Out.a[1] += hadd(result.ay);
        Out.a[2] += hadd(result.az);
        Out.fPot += hadd(result.pot);
        Out.dirsum += hadd(result.ir);
        Out.normsum += hadd(result.norm);
    }

    void queuePP(PKD pkd, workParticle *wp, ilpList &ilp, int bGravStep,
                 bool bGPU = true) {
        for (auto &tile : ilp) {
            ++wp->nTilesTotal;
            // if (bGPU) {
            // #ifdef USE_CUDA
            //     if (pkd->cudaClient->queuePP(wp, tile, bGravStep))
            //         continue;
            // #endif
            // #ifdef USE_METAL
            //     if (pkd->metalClient->queuePP(wp, tile, bGravStep))
            //         continue;
            //#endif
            //}
            ++wp->nTilesCPU;
            for (auto i = 0; i < wp->nP; ++i) {
                pkdGravEvalPP(wp->pInfoIn[i], tile, wp->pInfoOut[i]);
                wp->dFlopSingleCPU += COST_FLOP_PP * tile.size();
            }
        }
    }

    void queueEwald(PKD pkd, workParticle *wp, bool bGPU = true) {
        int i;
        int nQueued = 0;
        if (bGPU) {
            // #ifdef USE_CUDA
            //         nQueued = pkd->cudaClient->queueEwald(wp);
            // #endif
        }
        ++wp->nRefs;
        for (i = nQueued; i < wp->nP; ++i) {
            PINFOIN *in = &wp->pInfoIn[i];
            PINFOOUT *out = &wp->pInfoOut[i];
            double r[3];
            r[0] = wp->c[0] + in->r[0];
            r[1] = wp->c[1] + in->r[1];
            r[2] = wp->c[2] + in->r[2];
            wp->dFlop +=
                pkdParticleEwald(pkd, r, out->a.data(), &out->fPot,
                                 &wp->dFlopSingleCPU, &wp->dFlopDoubleCPU);
        }
        this->particleWorkDone(wp);
    }

    double pkdParticleEwald(PKD pkd, double *r, float *pa, float *pPot,
                            double *pdFlopSingle, double *pdFlopDouble) {
        struct EwaldVariables &ew = pkd->ew;
        EwaldTable *ewt = &pkd->ewt;
        const MOMC &restrict mom = ew.mom;
        double L, Pot, ax, ay, az, dx, dy, dz, x, y, z, r2;
#ifdef USE_SIMD_EWALD
        dvec dPot, dax, day, daz;
        fvec fPot, fax, fay, faz, fx, fy, fz;
        alignas(dvec) dvec::array_t px, py, pz, pr2, pInHole;
        int nSIMD = 0;
#endif
        int ix, iy, iz;
        int bInHole, bInHolex, bInHolexy;
        double dFlopSingle = 0;
        double dFlopDouble = 0;
        int nLoop = 0;

        L = ew.Lbox;
        dx = r[0] - ew.r[0];
        dy = r[1] - ew.r[1];
        dz = r[2] - ew.r[2];

        ax = ay = az = 0.0;
        Pot = mom.m * ew.k1;
#ifdef USE_SIMD_EWALD
        dPot.zero();
        dax.zero();
        day.zero();
        daz.zero();
#endif
        for (ix = -ew.nEwReps; ix <= ew.nEwReps; ++ix) {
            bInHolex = (abs(ix) <= ew.nReps);
            x = dx + ix * L;
            for (iy = -ew.nEwReps; iy <= ew.nEwReps; ++iy) {
                bInHolexy = (bInHolex && abs(iy) <= ew.nReps);
                y = dy + iy * L;
                for (iz = -ew.nEwReps; iz <= ew.nEwReps; ++iz) {
                    bInHole = (bInHolexy && abs(iz) <= ew.nReps);
                    z = dz + iz * L;
                    r2 = x * x + y * y + z * z;
                    if (r2 > ew.fEwCut2 && !bInHole)
                        continue;
                    if (r2 < ew.fInner2) {
                        double g0, g1, g2, g3, g4, g5, alphan;
                        /*
                         * For small r, series expand about
                         * the origin to avoid errors caused
                         * by cancellation of large terms.
                         */
                        alphan = ew.ka;
                        r2 *= ew.alpha2;
                        g0 = alphan * ((1.0 / 3.0) * r2 - 1.0);
                        alphan *= 2 * ew.alpha2;
                        g1 = alphan * ((1.0 / 5.0) * r2 - (1.0 / 3.0));
                        alphan *= 2 * ew.alpha2;
                        g2 = alphan * ((1.0 / 7.0) * r2 - (1.0 / 5.0));
                        alphan *= 2 * ew.alpha2;
                        g3 = alphan * ((1.0 / 9.0) * r2 - (1.0 / 7.0));
                        alphan *= 2 * ew.alpha2;
                        g4 = alphan * ((1.0 / 11.0) * r2 - (1.0 / 9.0));
                        alphan *= 2 * ew.alpha2;
                        g5 = alphan * ((1.0 / 13.0) * r2 - (1.0 / 11.0));

                        dFlopDouble +=
                            evalEwald<double, struct EwaldVariables, MOMC>(
                                ew, ew.mom, ax, ay, az, Pot, x, y, z, g0, g1,
                                g2, g3, g4, g5);
                    }
                    else {
#if defined(USE_SIMD_EWALD)
                        px[nSIMD] = x;
                        py[nSIMD] = y;
                        pz[nSIMD] = z;
                        pr2[nSIMD] = r2;
                        pInHole[nSIMD] = bInHole;
                        //          doerfc.i[nSIMD] = bInHole ? 0 : UINT64_MAX;
                        if (++nSIMD == dvec::width()) {
                            dFlopDouble += evalEwaldSIMD2(
                                               pkd, &pkd->es, dax, day, daz, dPot, dvec(px),
                                               dvec(py), dvec(pz), dvec(pr2),
                                               dvec(pInHole) == 0.0);
                            nSIMD = 0;
                        }
#else
                        double dir, dir2, a;
                        double g0, g1, g2, g3, g4, g5, alphan;
                        dir = 1 / sqrt(r2);
                        dir2 = dir * dir;
                        a = exp(-r2 * ew.alpha2);
                        a *= ew.ka * dir2;

                        if (bInHole) {
                            g0 = -erf(ew.alpha * r2 * dir);
                        }
                        else {
                            g0 = erfc(ew.alpha * r2 * dir);
                        }
                        g0 *= dir;
                        g1 = g0 * dir2 + a;
                        alphan = 2 * ew.alpha2;
                        g2 = 3 * g1 * dir2 + alphan * a;
                        alphan *= 2 * ew.alpha2;
                        g3 = 5 * g2 * dir2 + alphan * a;
                        alphan *= 2 * ew.alpha2;
                        g4 = 7 * g3 * dir2 + alphan * a;
                        alphan *= 2 * ew.alpha2;
                        g5 = 9 * g4 * dir2 + alphan * a;
                        dFlopDouble +=
                            evalEwald<double, struct EwaldVariables, MOMC>(
                                ew, ew.mom, ax, ay, az, Pot, x, y, z, g0, g1,
                                g2, g3, g4, g5);
                        // dFlopDouble +=
                        // evalEwald(ew,&ax,&ay,&az,&Pot,x,y,z,g0,g1,g2,g3,g4,g5);
#endif
                    }
                    ++nLoop;
                }
            }
        }
#if defined(USE_SIMD_EWALD)
        /* Finish remaining SIMD operations if necessary */
        if (nSIMD) { /* nSIMD can be 0 through 7 */
#define M (~0LL)
#if defined(__AVX512F__)
            static const i64v::array_t keepmask[] = {
                {0, 0, 0, 0, 0, 0, 0, 0}, {M, 0, 0, 0, 0, 0, 0, 0},
                {M, M, 0, 0, 0, 0, 0, 0}, {M, M, M, 0, 0, 0, 0, 0},
                {M, M, M, M, 0, 0, 0, 0}, {M, M, M, M, M, 0, 0, 0},
                {M, M, M, M, M, M, 0, 0}, { M, M, M, M, M, M, M, 0 }
            };
#elif defined(__AVX__)
            static const i64v::array_t keepmask[] = {{0, 0, 0, 0},
                {M, 0, 0, 0},
                {M, M, 0, 0},
                {
                    M,
                    M,
                    M,
                    0
                }
            };
#else
            static const i64v::array_t keepmask[] = {{0, 0}, { M, 0 }};
#endif
#undef M
            dvec t, tax = 0, tay = 0, taz = 0, tpot = 0;
            evalEwaldSIMD2(pkd, &pkd->es, tax, tay, taz, tpot, dvec(px),
                           dvec(py), dvec(pz), dvec(pr2), dvec(pInHole) == 0.0);
            dFlopDouble += COST_FLOP_EWALD * nSIMD;
            t = cast_dvec(i64v(keepmask[nSIMD]));
            tax = tax & t;
            tay = tay & t;
            taz = taz & t;
            tpot = tpot & t;
            dax += tax;
            day += tay;
            daz += taz;
            dPot -= tpot;
            nSIMD = 0;
        }
#endif

#ifdef USE_SIMD_EWALD
        /* h-loop is done in float precision */
        fax = cvt_fvec(dax);
        fay = cvt_fvec(day);
        faz = cvt_fvec(daz);
        fPot = cvt_fvec(dPot);

        fx = dx;
        fy = dy;
        fz = dz;

        nLoop = (ew.nEwhLoop + fvec::mask()) / fvec::width();
        int i = 0;
        do {
            fvec hdotx, s, c, t;
            hdotx = fvec(ewt->hx.p[i]) * fx + fvec(ewt->hy.p[i]) * fy +
                    fvec(ewt->hz.p[i]) * fz;
            fvec svec, cvec;

            sincosf(fvec(hdotx), svec, cvec);
            s = svec;
            c = cvec;

            fPot += fvec(ewt->hSfac.p[i]) * s + fvec(ewt->hCfac.p[i]) * c;
            s *= fvec(ewt->hCfac.p[i]);
            c *= fvec(ewt->hSfac.p[i]);
            t = s - c;
            fax += fvec(ewt->hx.p[i]) * t;
            fay += fvec(ewt->hy.p[i]) * t;
            faz += fvec(ewt->hz.p[i]) * t;
        } while (++i < nLoop);
        dFlopSingle += ew.nEwhLoop * COST_FLOP_HLOOP;

        ax += hadd(fax);
        ay += hadd(fay);
        az += hadd(faz);
        Pot += hadd(fPot);
#else
        /*
        ** Scoring for the h-loop (+,*)
        **  Without trig = (10,14)
        **      Trig est.    = 2*(6,11)  same as 1/sqrt scoring.
        **      Total        = (22,36)
        **                   = 58
        */
        for (auto i = 0; i < ew.nEwhLoop; ++i) {
            double hdotx, s, c, t;
            hdotx = ewt->hx.f[i] * dx + ewt->hy.f[i] * dy + ewt->hz.f[i] * dz;
            c = cos(hdotx);
            s = sin(hdotx);
            Pot += ewt->hCfac.f[i] * c + ewt->hSfac.f[i] * s;
            t = ewt->hCfac.f[i] * s - ewt->hSfac.f[i] * c;
            ax += ewt->hx.f[i] * t;
            ay += ewt->hy.f[i] * t;
            az += ewt->hz.f[i] * t;
        }
        dFlopDouble += ew.nEwhLoop * COST_FLOP_HLOOP;
#endif
        pa[0] += ax;
        pa[1] += ay;
        pa[2] += az;
        *pPot += Pot;

        *pdFlopSingle += dFlopSingle;
        *pdFlopDouble += dFlopDouble;
        return dFlopDouble + dFlopSingle;
    }

    void pkdEwaldInit(PKD pkd, int nReps, double fEwCut, double fhCut,
                      bool bGPU) {
        struct EwaldVariables *const ew = &pkd->ew;
        EwaldTable *const ewt = &pkd->ewt;
        const MOMC *restrict mom = &ew->mom;
        int i, hReps, hx, hy, hz, h2;
        double k4, L;
        double gam[6], mfacc, mfacs;
        double ax, ay, az;
        const int iOrder = 4;

        L = pkd->fPeriod[0];
        ew->Lbox = L;
        /*
        ** Create SIMD versions of the moments.
        */
#if defined(USE_SIMD_EWALD) && (defined(__SSE2__) || defined(__ARM_NEON))
        pkd->es.ewm.m = dvec(mom->m);
        pkd->es.ewm.xx = dvec(mom->xx);
        pkd->es.ewm.yy = dvec(mom->yy);
        pkd->es.ewm.xy = dvec(mom->xy);
        pkd->es.ewm.xz = dvec(mom->xz);
        pkd->es.ewm.yz = dvec(mom->yz);
        pkd->es.ewm.xxx = dvec(mom->xxx);
        pkd->es.ewm.xyy = dvec(mom->xyy);
        pkd->es.ewm.xxy = dvec(mom->xxy);
        pkd->es.ewm.yyy = dvec(mom->yyy);
        pkd->es.ewm.xxz = dvec(mom->xxz);
        pkd->es.ewm.yyz = dvec(mom->yyz);
        pkd->es.ewm.xyz = dvec(mom->xyz);
        pkd->es.ewm.xxxx = dvec(mom->xxxx);
        pkd->es.ewm.xyyy = dvec(mom->xyyy);
        pkd->es.ewm.xxxy = dvec(mom->xxxy);
        pkd->es.ewm.yyyy = dvec(mom->yyyy);
        pkd->es.ewm.xxxz = dvec(mom->xxxz);
        pkd->es.ewm.yyyz = dvec(mom->yyyz);
        pkd->es.ewm.xxyy = dvec(mom->xxyy);
        pkd->es.ewm.xxyz = dvec(mom->xxyz);
        pkd->es.ewm.xyyz = dvec(mom->xyyz);
        pkd->es.ewm.zz = dvec(mom->zz);
        pkd->es.ewm.xzz = dvec(mom->xzz);
        pkd->es.ewm.yzz = dvec(mom->yzz);
        pkd->es.ewm.zzz = dvec(mom->zzz);
        pkd->es.ewm.xxzz = dvec(mom->xxzz);
        pkd->es.ewm.xyzz = dvec(mom->xyzz);
        pkd->es.ewm.xzzz = dvec(mom->xzzz);
        pkd->es.ewm.yyzz = dvec(mom->yyzz);
        pkd->es.ewm.yzzz = dvec(mom->yzzz);
        pkd->es.ewm.zzzz = dvec(mom->zzzz);
#endif

        /*
        ** Set up traces of the complete multipole moments.
        */
        ew->Q4xx = 0.5 * (mom->xxxx + mom->xxyy + mom->xxzz);
        ew->Q4xy = 0.5 * (mom->xxxy + mom->xyyy + mom->xyzz);
        ew->Q4xz = 0.5 * (mom->xxxz + mom->xyyz + mom->xzzz);
        ew->Q4yy = 0.5 * (mom->xxyy + mom->yyyy + mom->yyzz);
        ew->Q4yz = 0.5 * (mom->xxyz + mom->yyyz + mom->yzzz);
        ew->Q4zz = 0.5 * (mom->xxzz + mom->yyzz + mom->zzzz);
        ew->Q4 = 0.25 * (ew->Q4xx + ew->Q4yy + ew->Q4zz);
        ew->Q3x = 0.5 * (mom->xxx + mom->xyy + mom->xzz);
        ew->Q3y = 0.5 * (mom->xxy + mom->yyy + mom->yzz);
        ew->Q3z = 0.5 * (mom->xxz + mom->yyz + mom->zzz);
        ew->Q2 = 0.5 * (mom->xx + mom->yy + mom->zz);
        ew->nReps = nReps;
        ew->nEwReps = ceil(fEwCut);
        ew->nEwReps = ew->nEwReps > nReps ? ew->nEwReps : nReps;
        ew->fEwCut2 = fEwCut * fEwCut * L * L;
        ew->fInner2 = 1.2e-3 * L * L;
        ew->alpha = 2.0 / L;
        ew->ialpha = 0.5 * L;
        ew->alpha2 = ew->alpha * ew->alpha;
        ew->k1 = M_PI / (ew->alpha2 * L * L * L);
        ew->ka = 2.0 * ew->alpha / sqrt(M_PI);
#if defined(USE_SIMD_EWALD) && (defined(__SSE2__) || defined(__ARM_NEON))
        pkd->es.ewp.Q4xx = dvec(ew->Q4xx);
        pkd->es.ewp.Q4xy = dvec(ew->Q4xy);
        pkd->es.ewp.Q4xz = dvec(ew->Q4xz);
        pkd->es.ewp.Q4yy = dvec(ew->Q4yy);
        pkd->es.ewp.Q4yz = dvec(ew->Q4yz);
        pkd->es.ewp.Q4zz = dvec(ew->Q4zz);
        pkd->es.ewp.Q4 = dvec(ew->Q4);
        pkd->es.ewp.Q3x = dvec(ew->Q3x);
        pkd->es.ewp.Q3y = dvec(ew->Q3y);
        pkd->es.ewp.Q3z = dvec(ew->Q3z);
        pkd->es.ewp.Q2 = dvec(ew->Q2);
        pkd->es.ewp.fEwCut2 = dvec(ew->fEwCut2);
        pkd->es.ewp.fInner2 = dvec(ew->fInner2);
        pkd->es.ewp.alpha = dvec(ew->alpha);
        pkd->es.ewp.ialpha = dvec(ew->ialpha);
        pkd->es.ewp.alpha2 = dvec(ew->alpha2);
        pkd->es.ewp.k1 = dvec(ew->k1);
        pkd->es.ewp.ka = dvec(ew->ka);
#endif

        /*
        ** Now setup stuff for the h-loop.
        */
        hReps = ceil(fhCut);
        k4 = M_PI * M_PI / (ew->alpha * ew->alpha * L * L);

        i = (int)pow(1 + 2 * hReps, 3);
#if defined(USE_SIMD_EWALD) && defined(__SSE__)
        i = (i + fvec::mask()) & ~fvec::mask();
#endif
        if (i > ew->nMaxEwhLoop) {
            ew->nMaxEwhLoop = i;
            ewt->hx.f = new (std::align_val_t(sizeof(fvec)))
            ewaldFloatType[ew->nMaxEwhLoop];
            assert(ewt->hx.f != NULL);
            ewt->hy.f = new (std::align_val_t(sizeof(fvec)))
            ewaldFloatType[ew->nMaxEwhLoop];
            assert(ewt->hy.f != NULL);
            ewt->hz.f = new (std::align_val_t(sizeof(fvec)))
            ewaldFloatType[ew->nMaxEwhLoop];
            assert(ewt->hz.f != NULL);
            ewt->hCfac.f = new (std::align_val_t(sizeof(fvec)))
            ewaldFloatType[ew->nMaxEwhLoop];
            assert(ewt->hCfac.f != NULL);
            ewt->hSfac.f = new (std::align_val_t(sizeof(fvec)))
            ewaldFloatType[ew->nMaxEwhLoop];
            assert(ewt->hSfac.f != NULL);
        }
        ew->nEwhLoop = i;
        i = (int)pow(1 + 2 * ew->nEwReps, 3);
#if defined(USE_SIMD_EWALD) && (defined(__SSE2__) || defined(__ARM_NEON))
        i = (i + fvec::mask()) & ~fvec::mask();
#endif
        i = 0;
        for (hx = -hReps; hx <= hReps; ++hx) {
            for (hy = -hReps; hy <= hReps; ++hy) {
                for (hz = -hReps; hz <= hReps; ++hz) {
                    h2 = hx * hx + hy * hy + hz * hz;
                    if (h2 == 0)
                        continue;
                    if (h2 > fhCut * fhCut)
                        continue;
                    assert(i < ew->nMaxEwhLoop);
                    gam[0] = exp(-k4 * h2) / (M_PI * h2 * L);
                    gam[1] = 2 * M_PI / L * gam[0];
                    gam[2] = -2 * M_PI / L * gam[1];
                    gam[3] = 2 * M_PI / L * gam[2];
                    gam[4] = -2 * M_PI / L * gam[3];
                    gam[5] = 2 * M_PI / L * gam[4];
                    gam[1] = 0.0;
                    gam[3] = 0.0;
                    gam[5] = 0.0;
                    ax = 0.0;
                    ay = 0.0;
                    az = 0.0;
                    mfacc = 0.0;
                    QEVAL(iOrder, ew->mom, gam, hx, hy, hz, ax, ay, az, mfacc);
                    gam[0] = exp(-k4 * h2) / (M_PI * h2 * L);
                    gam[1] = 2 * M_PI / L * gam[0];
                    gam[2] = -2 * M_PI / L * gam[1];
                    gam[3] = 2 * M_PI / L * gam[2];
                    gam[4] = -2 * M_PI / L * gam[3];
                    gam[5] = 2 * M_PI / L * gam[4];
                    gam[0] = 0.0;
                    gam[2] = 0.0;
                    gam[4] = 0.0;
                    ax = 0.0;
                    ay = 0.0;
                    az = 0.0;
                    mfacs = 0.0;
                    QEVAL(iOrder, ew->mom, gam, hx, hy, hz, ax, ay, az, mfacs);
                    ewt->hx.f[i] = 2 * M_PI / L * hx;
                    ewt->hy.f[i] = 2 * M_PI / L * hy;
                    ewt->hz.f[i] = 2 * M_PI / L * hz;
                    ewt->hCfac.f[i] = mfacc;
                    ewt->hSfac.f[i] = mfacs;
                    ++i;
                }
            }
        }
        ew->nEwhLoop = i;
        while (i < ew->nMaxEwhLoop) {
            ewt->hx.f[i] = 0;
            ewt->hy.f[i] = 0;
            ewt->hz.f[i] = 0;
            ewt->hCfac.f[i] = 0;
            ewt->hSfac.f[i] = 0;
            ++i;
        }
        // if (bGPU) {
        // #ifdef USE_CL
        //     clEwaldInit(pkd->mdl->clCtx, ew, ewt);
        //     mdlThreadBarrier(pkd->mdl);
        // #endif
        // #ifdef USE_CUDA
        //     // Only one thread needs to transfer the tables to the GPU
        //     if (pkd->mdl->Core() == 0)
        //         pkd->cudaClient->setupEwald(ew, ewt);
        //     pkd->mdl->ThreadBarrier();
        // #endif
        // }
    }
    void particleWorkDoneImpl(workParticle *wp) {
        base::particleWorkDoneImpl(wp);
        // We overwrite the dFlop, in the case of the gravity
        // wp->dFlop only accounts for the Ewald, the rest is
        // accounted for manually
        this->walkStats->dFlop =
            wp->dFlop + wp->nP * (this->pkd->ilp.count() * COST_FLOP_PP +
                                  this->pkd->ilc.count() * COST_FLOP_PC);
    }

    void processWorkParticleImpl(workParticle *wp) {
        /*
        ** Evaluate the local expansion.
        */
        if (pLoc) {
            for (auto i = 0; i < wp->nP; i++) {
                momFloat ax, ay, az, dPot;
                blitz::TinyVector<double, 3> r = wp->c + wp->pInfoIn[i].r;

                /*
                ** Evaluate local expansion.
                */
                blitz::TinyVector<double, 3> dr = r - this->kdn_r;
                dPot = 0;
                ax = 0;
                ay = 0;
                az = 0;

                momEvalLocr(pLoc, dr[0], dr[1], dr[2], &dPot, &ax, &ay, &az);

                wp->pInfoOut[i].fPot += dPot;
                wp->pInfoOut[i].a[0] += ax;
                wp->pInfoOut[i].a[1] += ay;
                wp->pInfoOut[i].a[2] += az;
            }
        }

        queuePC(this->pkd, wp, *this->ilc, this->options->ts.bGravStep,
                this->options->bGPU);

        queuePP(this->pkd, wp, *this->il, this->options->ts.bGravStep,
                this->options->bGPU);

        if (this->options->bEwald)
            queueEwald(this->pkd, wp, this->options->bGPU);
    }
};
class GravityInteraction : public GravityInteractionImpl<GravityInteraction> {
public:
    explicit GravityInteraction()
        : GravityInteractionImpl<GravityInteraction>() {}
};
#endif
