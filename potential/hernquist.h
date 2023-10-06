#ifndef HERNQUIST_HINCLUDED
#define HERNQUIST_HINCLUDED

static inline std::tuple<blitz::TinyVector<double,3>, double>  hernquist(blitz::TinyVector<double,3> pos) {
    // Hard coded just for the isolated galaxy test
    const double const_reduced_hubble_cgs = 3.2407789e-18;
    const double dMsolUnit = 1e10;
    const double dKpcUnit = 1.0;
    const double dKmPerSecUnit = sqrt(GCGS*dMsolUnit*MSOLG
                                      /(dKpcUnit*KPCCM))/1e5;
    const double H0 = 70.4/ dKmPerSecUnit * ( dKpcUnit / 1e3);

    const double concentration = 9.0;
    const double M200 = 135.28423603962767;
    const double V200 = cbrt(10.*M200*H0);
    const double R200 = cbrt(M200/(100.*H0*H0));
    const double RS = R200 / concentration;

    const double al = RS * sqrt(2. * (log(1. + concentration) -
                                      concentration / (1. + concentration)));

    const double mass = M200*(1.-0.041);
    const double sqrtgm_inv = 1.f / sqrt(mass);
    const double epsilon =  0.2/dKpcUnit;
    const double epsilon2 = epsilon*epsilon;

    /* Calculate the acceleration (assume centred in [0,0,0])*/
    const double rr = sqrtf(blitz::dot(pos,pos) + epsilon2);
    const double r_plus_a_inv = 1.f / (rr + al);
    const double r_plus_a_inv2 = r_plus_a_inv * r_plus_a_inv;
    const double term = -mass * r_plus_a_inv2 / rr;

    /* Calculate the circular orbital period */
    const double period = 2.f * M_PI * sqrtf(rr) * al *
                          (1 + rr / al) * sqrtgm_inv;

    auto a = term * pos;
    const double time_step = 0.01 * period;

    return { a, time_step };
}

#endif
