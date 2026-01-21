#ifndef HYDRO_H
#define HYDRO_H
#include "TraversePST.h"
#include "pkd.h"
#include "smooth/smooth.h"
#include <vector>
#include "blitz/array.h"
#include "blitz/memblock.h"

class ServiceHydroDensity : public TraverseCount<ServiceHydroDensity> {
  public:
    static constexpr std::string_view name = "HydroDensity";
    struct input {
        double dNeighborsStd;
        double dhMinOverSoft;
        int nSmooth;
        bool bPeriodic;
        bool bIterativeSmoothingLength;
    };
    using output = uint64_t;
    explicit ServiceHydroDensity(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceHydroDensityFinal : public TraverseCount<ServiceHydroDensityFinal> {
  public:
    static constexpr std::string_view name = "HydroDensityFinal";
    struct input {
        int nSmooth;
        bool bPeriodic;
    };
    using output = uint64_t;
    explicit ServiceHydroDensityFinal(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceBHMerger : public TraverseCount<ServiceBHMerger> {
  public:
    static constexpr std::string_view name = "BHMerger";
    struct input {
        int nSmooth;
        bool bPeriodic;
    };
    using output = uint64_t;
    explicit ServiceBHMerger(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceHydroFluxes : public TraverseCount<ServiceHydroFluxes> {
  public:
    static constexpr std::string_view name = "HydroFluxes";
    struct input {
        int nSmooth;
        bool bPeriodic;
        double dTime;
        double dDelta;
        double a;
        double H;
        double dConstGamma;
        struct eEOSparam eEOS;
    };
    using output = uint64_t;
    explicit ServiceHydroFluxes(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceBHGasPin : public TraversePST {
  public:
    static constexpr std::string_view name = "BHGasPin";
    struct input {
        int nSmooth;
        bool bPeriodic;
        double a;
    };
    using output = void;
    explicit ServiceBHGasPin(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceBHStep : public TraversePST {
  public:
    static constexpr std::string_view name = "BHStep";
    struct input {
        int nSmooth;
        bool bPeriodic;
    };
    using output = void;
    explicit ServiceBHStep(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceHydroGradients : public TraverseCount<ServiceHydroGradients> {
  public:
    static constexpr std::string_view name = "HydroGradients";
    struct input {
        int nSmooth;
        bool bPeriodic;
        bool bConditionedBarthJespersen;
        bool bStricterSlopeLimiter;
    };
    using output = uint64_t;
    explicit ServiceHydroGradients(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

class ServiceBHEvolve : public TraversePST {
  public:
    static constexpr std::string_view name = "BHEvolve";
    struct input {
        int nSmooth;
        bool bPeriodic;
        double a;
        double dBHAccretionAlpha;
        double dBHAccretionCvisc;
        double dBHAccretionEddFac;
        double dConstGamma;
        double dBHFBEff;
        double dBHFBEcrit;
        bool bBHAccretion;
        bool bBHFeedback;
    };
    using output = void;
    explicit ServiceBHEvolve(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

template <typename T, typename Enum>
class EnumArray {
    static_assert(std::is_enum_v<Enum>, "EnumArray requires an enum type");

  public:
    using value_type = T;

    // Raw-buffer constructor
    EnumArray(T* data, Enum last, int nBuff) : A_(data, blitz::shape(static_cast<int>(last), nBuff), blitz::deleteDataWhenDone) {}

    // Owning: allocate internal storage
    EnumArray(int nBuff) : A_(blitz::shape(static_cast<int>(Enum::last), nBuff)) {}

    // Construct from existing Blitz array
    explicit EnumArray(const blitz::Array<std::remove_const_t<T>, 2>& array) : A_(array) {}

    EnumArray() = default;

    // Rebind to another EnumArray (no copy)
    void reference(EnumArray const& other) { A_.reference(other.A_); }

    [[nodiscard]] auto extent(int i) const noexcept { return A_.extent(i); }

    // Resize the second dimension
    void resize(int nBuff) { A_.resize(blitz::shape(A_.extent(0), nBuff)); }

    // Indexing (constness propagates automatically)
    template <typename... I>
    auto operator()(Enum e, I... i) -> decltype(auto) {
        return A_(static_cast<int>(e), i...);
    }

    template <typename... I>
    auto operator()(Enum e, I... i) const -> decltype(auto) {
        return A_(static_cast<int>(e), i...);
    }

    [[nodiscard]] auto array() const noexcept -> const blitz::Array<std::remove_const_t<T>, 2>& { return A_; }

    // Only enable mutable access when T is non-const
    template <typename U = T>
    [[nodiscard]] auto array() noexcept -> blitz::Array<U, 2>&
        requires(!std::is_const_v<U>)
    {
        return A_;
    }

  private:
    blitz::Array<T, 2> A_;
};
template <typename T, typename Enum>
EnumArray(T const*, Enum, int) -> EnumArray<const T, Enum>;
template <typename T, typename Enum>
EnumArray(T*, Enum, int) -> EnumArray<T, Enum>;

// Field orders for moment of inertia. Operator+ is overloaded to convert enum to size_t
enum class MOI : int { XX = 0, XY = 1, XZ = 2, YY = 3, YZ = 4, ZZ = 5 };
constexpr auto operator+(MOI m) noexcept -> std::size_t {
    return static_cast<std::size_t>(m);
}
namespace meshless {
using myreal = double;
} // namespace meshless

namespace meshless::FLUX {
/* When doing SIMD, a structure of arrays is used,
 * rather than an array of structures.
 *
 * To simplify the indexing of the elements, these enum should be
 * always used
 */
enum class IFLD : int {
    ball,
    dx,
    dy,
    dz,
    dr,
    rung,
    rho,
    P,
#ifdef ENTROPY_SWITCH
    S,
#endif
    vx,
    vy,
    vz,
    gradRhoX,
    gradRhoY,
    gradRhoZ,
    gradPX,
    gradPY,
    gradPZ,
    gradVxX,
    gradVxY,
    gradVxZ,
    gradVyX,
    gradVyY,
    gradVyZ,
    gradVzX,
    gradVzY,
    gradVzZ,
    lastAccX,
    lastAccY,
    lastAccZ,
    lastUpdateTime,
    B_XX,
    B_YY,
    B_ZZ,
    B_XY,
    B_XZ,
    B_YZ,
    omega,
    spheom,
    last
};

enum class OFLD : int {
    minDt,
    Frho,
    FmomX,
    FmomY,
    FmomZ,
    Fene,
#ifdef ENTROPY_SWITCH
    FS,
#endif
    last
};
} // namespace meshless::FLUX

/* -----------------
 * MAIN FUNCTIONS AND CLASSES
 * -----------------
 */

// Cubic spline kernel as defined in Dehnen & Aly (2012)
//
// (1 − q)^3 − 4(1/2 − q)^3  for   0 < q <= 1/2
//                (1 − q)^3  for 1/2 < q <= 1
//                        0  for       q > 1
//
// where q = r/H.
template <typename dtype = double, typename mtype = bool>
inline dtype cubicSplineKernel(dtype r, dtype H) {
    const dtype dNormFac = 16 * M_1_PI;
    const dtype dHInv = 1.0 / H;
    const dtype q = r * dHInv;
    const dtype dNorm = dNormFac * dHInv * dHInv * dHInv;

    if (q < 0.5) {
        return dNorm * (0.5 - 3 * q * q * (1.0 - q));
    }
    else if (q < 1.0) {
        return dNorm * (1.0 - q) * (1.0 - q) * (1.0 - q);
    }
    else {
        return 0.0;
    }
}

#ifdef SIMD_H
template <>
inline dvec cubicSplineKernel(dvec r, dvec H) {
    const dvec dNormFac = 16 * M_1_PI;
    const dvec dHInv = 1.0 / H;
    const dvec q = r * dHInv;
    const dvec dNorm = dNormFac * dHInv * dHInv * dHInv;
    const dvec oneMq = 1.0 - q;
    dvec out = 0.0;
    dvec t;
    t = dNorm * oneMq * oneMq * oneMq;
    out = mask_mov(out, q < 1.0, t);
    t = dNorm * (0.5 - 3 * q * q * oneMq);
    out = mask_mov(out, q < 0.5, t);
    return out;
}
#endif

// First derivative (i.e. dW / dq) of the cubic spline kernel
template <typename dtype = double, typename mtype = bool>
inline dtype cubicSplineKernelDeriv(dtype r, dtype H) {
    const dtype dHInv = 1.0 / H;
    const dtype q = r * dHInv;
    const dtype dNorm = 48. * M_1_PI * dHInv * dHInv * dHInv;

    if (q < 0.5) {
        return dNorm * q * (3. * q - 2.);
    }
    else if (q < 1.0) {
        return -dNorm * (1. - q) * (1. - q);
    }
    else {
        return 0.0;
    }
}

#ifdef SIMD_H
template <>
inline dvec cubicSplineKernelDeriv(dvec r, dvec H) {
    const dvec dHInv = 1.0 / H;
    const dvec q = r * dHInv;
    const dvec dNorm = 48. * M_1_PI * dHInv * dHInv * dHInv;
    dvec out{0.}, t;

    t = -dNorm * (1. - q) * (1. - q);
    out = mask_mov(out, q < 1.0, t);
    t = dNorm * q * (3. * q - 2.);
    out = mask_mov(out, q < 0.5, t);
    return out;
}
#endif

/* Density loop */
struct hydroDensityPack {
    blitz::TinyVector<double, 3> position;
    float fBall;
    uint8_t iClass;
    bool bMarked;
};

class SmoothHelperHydroDensityBase : public SmoothHelperPacked<hydroDensityPack> {
    void pack(void* dst, const void* src) final;
    void unpack(void* dst, const void* src, const void* key) final;
    void init_shared(void* dst) final { init_shared_gas_decorator(pkd, dst); }

  protected:
    SmoothHelperHydroDensityBase(pkdContext* pkd) : SmoothHelperPacked<hydroDensityPack>(pkd, false) {}
};

class SmoothHelperHydroDensity : public SmoothHelperHydroDensityBase {
    double dNeighborsStd;
    bool bIterativeSmoothingLength;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;

  public:
    SmoothHelperHydroDensity(pkdContext* pkd, double dNeighborsStd, bool bIterativeSmoothingLength)
        : SmoothHelperHydroDensityBase(pkd), dNeighborsStd(dNeighborsStd), bIterativeSmoothingLength(bIterativeSmoothingLength) {}
};

class SmoothHelperHydroDensityFinal : public SmoothHelperHydroDensityBase {
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;

  public:
    SmoothHelperHydroDensityFinal(pkdContext* pkd) : SmoothHelperHydroDensityBase(pkd) {}
};

/* Gradient loop */
struct hydroGradientsPack {
    blitz::TinyVector<double, 3> position;
    blitz::TinyVector<double, 3> velocity;
    double P;
    float fBall;
    float fDensity;
    uint8_t iClass;
    bool bMarked;
};

class SmoothHelperHydroGradients : public SmoothHelperPacked<hydroGradientsPack> {
    bool bConditionedBarthJespersen;
    bool bStricterSlopeLimiter;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    void pack(void* dst, const void* src) final;
    void unpack(void* dst, const void* src, const void* key) final;
    void init_shared(void* dst) final { init_shared_gas_decorator(pkd, dst); }

  public:
    SmoothHelperHydroGradients(pkdContext* pkd, bool bConditionedBarthJespersen, bool bStricterSlopeLimiter)
        : SmoothHelperPacked<hydroGradientsPack>(pkd, false), bConditionedBarthJespersen(bConditionedBarthJespersen),
          bStricterSlopeLimiter(bStricterSlopeLimiter) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};

/* Flux loop */
struct hydroFluxesPack {
    blitz::TinyVector<double, 3> position;
    blitz::TinyVector<double, 3> velocity;
    blitz::TinyVector<double, 6> B;
    blitz::TinyVector<meshless::myreal, 3> gradRho, gradVx, gradVy, gradVz, gradP;
    meshless::myreal lastUpdateTime;
    blitz::TinyVector<meshless::myreal, 3> lastAcc;
    double omega;
    double spheom;
    double P;
    float fBall;
    float fDensity;
    uint8_t uRung;
    uint8_t iClass;
    bool bMarked;
};

struct hydroFluxesFlush {
    meshless::myreal Frho;
    blitz::TinyVector<meshless::myreal, 3> Fmom;
    meshless::myreal Fene;
    blitz::TinyVector<double, 3> mom;
    double E;
    double Uint;
    // Used for the MFV method, ignored otherwise
    blitz::TinyVector<double, 3> drDotFrho;
    float fMass;
};

class SmoothHelperHydroFluxes : public SmoothHelperPacked<hydroFluxesPack, hydroFluxesFlush> {
    EnumArray<meshless::myreal, meshless::FLUX::IFLD> ia;
    EnumArray<meshless::myreal, meshless::FLUX::OFLD> oa;
    const double dTime, dDelta, a_, H_, dConstGamma;
    eEOSparam eEOS;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;
    virtual void init(void* dst) override final;
    virtual void flush(void* dst, const void* src) override final;
    virtual void combine(void* dst, const void* src, const void* key) override final;
    virtual void init_shared(void* dst) override final { init_shared_gas_decorator(pkd, dst); }

  public:
    SmoothHelperHydroFluxes(pkdContext* pkd, double dTime, double dDelta, double a, double H, double dConstGamma, const eEOSparam& eEOS)
        : SmoothHelperPacked<hydroFluxesPack, hydroFluxesFlush>(pkd, true), dTime(dTime), dDelta(dDelta), a_(a), H_(H), dConstGamma(dConstGamma), eEOS(eEOS) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};

/* Time step loop */
struct hydroStepPack {
    blitz::TinyVector<double, 3> position;
    blitz::TinyVector<double, 3> velocity;
    float c;
    float fBall;
    uint8_t uRung;
    uint8_t uNewRung;
    uint8_t iClass;
    bool bMarked;
};

struct hydroStepFlush {
    uint8_t uNewRung;
    uint8_t uWake;
};

class SmoothHelperHydroStep : public SmoothHelperPacked<hydroStepPack, hydroStepFlush> {
    double dEtaCourant, a_, dCFLacc, dDelta;
    bool bDoGravity;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;
    virtual void init(void* dst) override final;
    virtual void flush(void* dst, const void* src) override final;
    virtual void combine(void* dst, const void* src, const void* key) override final;
    virtual void init_shared(void* dst) override final { init_shared_gas_decorator(pkd, dst); }

  public:
    SmoothHelperHydroStep(pkdContext* pkd, double dEtaCourant, double a, double dCFLacc, double dDelta, bool bDoGravity)
        : SmoothHelperPacked<hydroStepPack, hydroStepFlush>(pkd, true), dEtaCourant(dEtaCourant), a_(a), dCFLacc(dCFLacc), dDelta(dDelta),
          bDoGravity(bDoGravity) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
};

void pkdWakeParticles(PKD pkd, int iRoot, double dTime, double dDelta);

/* Source terms */
void hydroSourceGravity(PKD pkd, particleStore::ParticleReference& p, decorators::MFM* psph, double pDelta, blitz::TinyVector<double, 3>& pa,
                        double dScaleFactor, int bComove);
void hydroSourceExpansion(PKD pkd, particleStore::ParticleReference& p, decorators::MFM* psph, double pDelta, double dScaleFactor, double dHubble, int bComove,
                          double dConstGamma);
void hydroSyncEnergies(PKD pkd, particleStore::ParticleReference& p, decorators::MFM* psph, const blitz::TinyVector<double, 3>& pa, double dConstGamma);
void hydroSetPrimitives(PKD pkd, particleStore::ParticleReference& p, decorators::MFM* psph, double dTuFac, double dConstGamma);
void hydroSetLastVars(PKD pkd, particleStore::ParticleReference& p, decorators::MFM* psph, const blitz::TinyVector<double, 3>& pa, double dScaleFactor,
                      double dTime, double dDelta, double dConstGamma);
void hydroResetFluxes(decorators::MFM* psph);

/* -----------------
 * HELPERS
 * -----------------
 */
void inverseMatrix(double* E, double* B);
double conditionNumber(double* E, double* B);
void compute_Ustar(double rho_K, double S_K, double v_K, double p_K, double h_K, double S_s, double* rho_sK, double* rhov_sK, double* e_sK);

#endif
