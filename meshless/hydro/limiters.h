#ifndef LIMITER_H
#define LIMITER_H
namespace meshless {
template<typename T>
int sign(T v) {return (v>0) - (v<0);}
}

#ifdef __CUDACC__
    #define PP_CUDA_BOTH __host__ __device__
#else
    #define PP_CUDA_BOTH
#endif


#ifndef PASS
    #define PASS
    static int equal = 0;
    static int pass1 = 0;
    static int pass2 = 0;
    static int pass3 = 0;
    static int pass4 = 0;
    static int pass5 = 0;
    static int pass6 = 0;
#endif
#define psi1 0.5
#define psi2 0.25
template <typename dtype=dvec>
PP_CUDA_BOTH static inline void genericPairwiseLimiter(dtype Lstate, dtype Rstate,
        dtype &Lstate_face, dtype &Rstate_face) {
#ifdef DEBUG_FLUX_NOLIMITER
    return;
#endif
    dtype phi_max, phi_min, d1, d2, phi_mean, phi_p, phi_m;

    auto equal = Lstate == Rstate;

    d1 = psi1*abs(Lstate - Rstate);
    d2 = psi2*abs(Lstate - Rstate);

    phi_mean = 0.5*(Lstate+Rstate);

    phi_min = min(Lstate, Rstate);
    phi_max = max(Lstate, Rstate);

    {
        // Probably there is a more intelligent way of doing this using
        // dtype::sign_mask()
        auto ncond = ((phi_min-d1) > 0) ^ (phi_min > 0);
        auto cond = ~ncond;
        phi_m = mask_mov(phi_m, cond, phi_min-d1);
        phi_m = mask_mov(phi_m, ncond, phi_min/(1. + d1/abs(phi_min)));
    }
    /*
    if (meshless::sign(phi_min - d1) == meshless::sign(phi_min) ) {
        phi_m = phi_min - d1;
    }
    else {
        phi_m = phi_min/(1. + d1/abs(phi_min));
    }
    */

    {
        auto ncond = (phi_max+d1 > 0) ^ (phi_max > 0);
        auto cond = ~ncond;
        phi_p = mask_mov(phi_p, cond, phi_max+d1);
        phi_p = mask_mov(phi_p, ncond, phi_max/(1. + d1/abs(phi_max)));
    }
    /*
    if (meshless::sign(phi_max + d1) == meshless::sign(phi_max) ) {
        phi_p = phi_max + d1;
    }
    else {
        phi_p = phi_max/(1. + d1/abs(phi_max));
    }
    */

    {
        auto cond = Lstate < Rstate;
        auto ncond = ~cond;
        Lstate_face = mask_mov( Lstate_face, cond, max(phi_m, min(phi_mean+d2, Lstate_face)));
        Rstate_face = mask_mov( Rstate_face, cond, min(phi_p, max(phi_mean-d2, Rstate_face)));
        Rstate_face = mask_mov( Rstate_face,ncond, max(phi_m, min(phi_mean+d2, Rstate_face)));
        Lstate_face = mask_mov( Lstate_face,ncond, min(phi_p, max(phi_mean-d2, Lstate_face)));
    }
    /*
    if (Lstate < Rstate) {
        *Lstate_face = max(phi_m, min(phi_mean+d2, *Lstate_face));
        *Rstate_face = min(phi_p, max(phi_mean-d2, *Rstate_face));
    }
    else {
        *Rstate_face = max(phi_m, min(phi_mean+d2, *Rstate_face));
        *Lstate_face = min(phi_p, max(phi_mean-d2, *Lstate_face));
    }
    */

    Lstate_face = mask_mov(Lstate_face, equal, Lstate);
    Rstate_face = mask_mov(Rstate_face, equal, Rstate);
}

#if !defined(SIMD_H) || defined(__CUDACC__)
template <>
PP_CUDA_BOTH inline void genericPairwiseLimiter(double Lstate, double Rstate,
        double &Lstate_face, double &Rstate_face) {
#ifdef DEBUG_FLUX_NOLIMITER
    return;
#endif
    double phi_max, phi_min, d1, d2, phi_mean, phi_p, phi_m;

    if (Lstate == Rstate) {
        Lstate_face = Lstate;
        Rstate_face = Rstate;
        //equal++;
    }
    else {

        d1 = psi1*fabs(Lstate - Rstate);
        d2 = psi2*fabs(Lstate - Rstate);

        phi_mean = 0.5*(Lstate+Rstate);

        phi_min = min(Lstate, Rstate);
        phi_max = max(Lstate, Rstate);

        if (meshless::sign(phi_min - d1) == meshless::sign(phi_min) ) {
            phi_m = phi_min - d1;
            //pass1++;
        }
        else {
            phi_m = phi_min/(1. + d1/fabs(phi_min));
            //pass2++;
        }

        if (meshless::sign(phi_max + d1) == meshless::sign(phi_max) ) {
            phi_p = phi_max + d1;
            //pass3++;
        }
        else {
            phi_p = phi_max/(1. + d1/fabs(phi_max));
            //pass4++;
        }

        if (Lstate < Rstate) {
            Lstate_face = max(phi_m, min(phi_mean+d2, Lstate_face));
            Rstate_face = min(phi_p, max(phi_mean-d2, Rstate_face));
            //pass5++;
        }
        else {
            Rstate_face = max(phi_m, min(phi_mean+d2, Rstate_face));
            Lstate_face = min(phi_p, max(phi_mean-d2, Lstate_face));
            //pass6++;
        }

    }
}
#endif

static inline double BarthJespersenLimiter(double max_diff, double min_diff,
        double max_extrap, double min_extrap, double beta) {
#ifdef DEBUG_FLUX_NOLIMITER
    return 1.;
#endif
    min_diff = std::max(0., std::min(max_diff/max_extrap, min_diff/min_extrap));
    return std::min(1., beta * min_diff);
}

static inline double BarthJespersenLimiter(double max_diff, double min_diff, double max_extrap, double beta) {
#ifdef DEBUG_FLUX_NOLIMITER
    return 1.;
#endif
    min_diff = std::max(0., std::min(max_diff, min_diff));
    return std::min(1., beta * min_diff / max_extrap);
}
#endif

