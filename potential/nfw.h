#ifndef NFW_HINCLUDED
#define NFW_HINCLUDED

/* This implementation follows the notation and formulae of
 * Mo, van den Bosch and White, Galaxy Formation and Evolution, 2010
 * Section 7.5.1
 *
 * Remember, G=1
 */

static inline std::tuple<blitz::TinyVector<double,3>, double>  nfw(blitz::TinyVector<double,3> pos) {

    const double dMsolUnit = 1e10;
    const double dKpcUnit = 1.0;
    const double H0 = 70.4;
    const double Omega_m = 0.3;

    const double c = 5;
    const double Delta_h = 200;
    const double M_h = 1e12; // Msol
    const double eps = 5;  // kpc

    const double dKmPerSecUnit = sqrt(GCGS*dMsolUnit*MSOLG
                                      /(dKpcUnit*KPCCM))/1e5;
    const double H0code = H0/ dKmPerSecUnit * ( dKpcUnit / 1e3);
    const double rho_crit = 3. * H0code * H0code / ( 8. * M_PI);
    const double rho_mean = rho_crit * Omega_m;

    /* Add a softnening to avoid divergence at r->0 */
    const double epsilon =  eps/dKpcUnit;
    const double epsilon2 = epsilon*epsilon;

    /* Assume centred in [0,0,0] */
    const double rr = sqrtf(blitz::dot(pos,pos) + epsilon2);

    /* Eq. 7.136 */
    const double rho_h = Delta_h * rho_crit * Omega_m;
    const double r_h = cbrt( 3. * M_h/dMsolUnit / (4. * M_PI * rho_h ) );
    /* Eq. 7.140 */
    const double r_s = r_h / c;

    /* Eq. 7.141 */
    const double den = log(1.+c) - c/(1.+c);
    const double delta_char = Delta_h / 3. * ( c * c * c )/den;

    /* Eq. 7.139 */
    const double x = rr/r_h;
    const double cx = c * x;
    const double enclosed_mass = 4. * M_PI * rho_mean * delta_char *
                                 r_s * r_s * r_s * ( log(1 + cx) - cx/(1+cx) );

    const double term = - enclosed_mass  / (rr*rr*rr);

    /* Calculate the circular orbital period */
    /* Eq. 11.24 */
    const double V_c = sqrt( enclosed_mass / rr ) ;
    const double period = 2.f * M_PI * rr / V_c;

    auto a = term * pos;
    const double time_step = 0.01 * period;

    return { a, time_step };
}

#endif
