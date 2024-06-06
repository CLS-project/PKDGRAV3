#ifndef LIMITER_H
#define LIMITER_H
namespace messhless {
template<typename T>
int sign(T v) {return (v>0) - (v<0);}
}

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
template <typename dtype=double>
static inline void genericPairwiseLimiter(dtype Lstate, dtype Rstate,
        dtype &Lstate_face, dtype &Rstate_face) {
#ifdef DEBUG_FLUX_NOLIMITER
    return;
#endif
    dtype phi_max, phi_min, d1, d2, phi_mean, phi_p, phi_m;

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

        if (messhless::sign(phi_min - d1) == messhless::sign(phi_min) ) {
            phi_m = phi_min - d1;
            //pass1++;
        }
        else {
            phi_m = phi_min/(1. + d1/fabs(phi_min));
            //pass2++;
        }

        if (messhless::sign(phi_max + d1) == messhless::sign(phi_max) ) {
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

#ifdef SIMD_H
template <>
inline void genericPairwiseLimiter(dvec Lstate, dvec Rstate,
                                   dvec &Lstate_face, dvec &Rstate_face) {
#ifdef DEBUG_FLUX_NOLIMITER
    return;
#endif
    dvec phi_max, phi_min, d1, d2, phi_mean, phi_p, phi_m;

    auto equal = Lstate == Rstate;

    d1 = psi1*abs(Lstate - Rstate);
    d2 = psi2*abs(Lstate - Rstate);

    phi_mean = 0.5*(Lstate+Rstate);

    phi_min = min(Lstate, Rstate);
    phi_max = max(Lstate, Rstate);

    {
        // Probably there is a more intelligent way of doing this using
        // dvec::sign_mask()
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
#endif

static inline void BarthJespersenLimiter(double *limVar, double *gradVar,
        double var_max, double var_min,
        double dx, double dy, double dz) {
#ifdef DEBUG_FLUX_NOLIMITER
    *limVar = 1;
    return;
#endif
    double diff, lim;

    diff = (gradVar[0]*dx + gradVar[1]*dy + gradVar[2]*dz);
    if (var_min > 0) { var_min=0; } //IA: Can happen due to machine precision
    if (var_max < 0) { var_max=0; } //IA: Can happen due to machine precision
    if (diff > 0.) {
        lim = var_max/diff;
    }
    else if (diff < 0.) {
        lim = var_min/diff;
    }
    else {
        lim = 1.;
    }
    if (lim > 1.) lim = 1.; // min(1,lim)
    if (lim < 0.) lim = 0.;
    if (lim < (*limVar)) { *limVar = lim;}
    // FIXME IA: Option to avoid extrapolation or limiter
//    *limVar = 1.0;
//    *limVar = 0.0;
}

/* IA: In this version we take into account the condition number,
 * which give us an idea about how 'well aligned' are the particles
 */
static inline void ConditionedBarthJespersenLimiter(double *limVar, double *gradVar,
        double var_max, double var_min,
        double dx, double dy, double dz,
        double Ncrit, double Ncond) {
#ifdef DEBUG_FLUX_NOLIMITER
    *limVar = 1;
    return;
#endif
    double diff, lim, beta;

    diff = Ncrit/Ncond;
    diff = diff < 1. ? diff : 1.;
    diff *= 2.;
    beta = (1. < diff) ? diff : 1.;

    diff = (gradVar[0]*dx + gradVar[1]*dy + gradVar[2]*dz);
    if (var_min > 0) { var_min=0; } //IA: Can happen due to machine precision
    if (var_max < 0) { var_max=0; } //IA: Can happen due to machine precision
    if (diff > 0.) {
        lim = var_max/diff;
    }
    else if (diff < 0.) {
        lim = var_min/diff;
    }
    else {
        lim = 1.;
    }
    lim *= beta;
    if (lim > 1.) lim = 1.; // min(1,lim)
    if (lim < 0.) lim = 0.;
    if (lim < (*limVar)) { *limVar = lim;}
    // FIXME IA: Option to avoid extrapolation or limiter
//    *limVar = 1.0;
//    *limVar = 0.0;
}
#endif
