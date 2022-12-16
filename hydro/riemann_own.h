#include "core/vmath.h"

#define GAMMA (dConstGamma)
#define GAMMA_MINUS1 ((GAMMA-1.0))
/* END IA */
#define GAMMA_G1 ((GAMMA-1.0)/(2.0*GAMMA))
#define GAMMA_G2 ((GAMMA+1.0)/(2.0*GAMMA))
#define GAMMA_G3 ((2.0*GAMMA/(GAMMA-1.0)))
#define GAMMA_G4 (2.0/(GAMMA-1.0))
#define GAMMA_G5 (2.0/(GAMMA+1.0))
#define GAMMA_G6 ((GAMMA-1.0)/(GAMMA+1.0))
#define GAMMA_G7 (0.5*(GAMMA-1.0))
#define GAMMA_G8 (1.0/GAMMA)
#define GAMMA_G9 (GAMMA-1.0)

#define TOL_ITER 1.e-6
#define NMAX_ITER 100

template <typename ftype=double>
static inline void dump(ftype a) {
    printf("%e\n", a);
}
template <>
static inline void dump(dvec a) {
    for (auto i=0; i<dvec::width(); i++)
        printf(" %e ", a[i]);
    printf("\n");
}

template <typename ftype=double>
static inline int nan_guard(ftype &var, ftype def) {
    if (isnan(var)) {
        var=def;
        return 1;
    }
    else {
        return 0;
    }
}

template <>
inline int nan_guard(dvec &var, dvec def) {
    int nan = 0;
    var = mask_mov(var, var!=var, def);
    for (auto j=0; j<dvec::width(); j++) {
        if (isnan(var[j])) {
            nan = 1;
        }
    }
    return nan;
}

template <typename ftype=double>
struct Conserved_var_Riemann { //IA: TODO: change the name of this struct..!!
    ftype rho;
    ftype p;  // IA: Pressure is a primitive variable, not a conserved one...
    ftype v[3]; // IA: same...
    ftype u;
    ftype cs;
};

template <typename ftype=double>
struct Input_vec_Riemann {
    struct Conserved_var_Riemann<ftype> L;
    struct Conserved_var_Riemann<ftype> R;
};
template <typename ftype=double>
struct Riemann_outputs {
    ftype P_M;
    ftype S_M;
    struct Conserved_var_Riemann<ftype> Fluxes;
};

template <typename ftype=double>
static inline ftype DMAX(ftype a, ftype b) { return (a > b) ? a : b; }
template <typename ftype=double>
static inline ftype DMIN(ftype a, ftype b) { return (a < b) ? a : b; }

static inline double min(double a, double b) {
    return (a<b)? a : b;
}
static inline double max(double a, double b) {
    return (a>b)? a : b;
}

template <typename ftype=double>
static inline void limit_decrease(ftype &Pg, ftype Pg_prev) {
    if (Pg < 0.1 * Pg_prev)
        Pg = 0.1 * Pg_prev;
}
template <>
inline void limit_decrease(dvec &Pg, dvec Pg_prev) {
    dvec limit = 0.1*Pg_prev;
    Pg = mask_mov(Pg, Pg<limit, limit);
}


template <typename ftype=double>
static inline ftype guess_two_rarefaction(ftype dConstGamma,
        ftype R_rho,ftype R_p,ftype L_rho,ftype L_p,
        ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    ftype pnu = (cs_L+cs_R) - GAMMA_G7 * (v_line_R - v_line_L);
    ftype pde = cs_L * pow(L_p, -GAMMA_G1) + cs_R * pow(R_p, -GAMMA_G1);
    return pow(pnu / pde, GAMMA_G3);
}


template <typename ftype=double>
static inline ftype guess_two_shock(ftype dConstGamma, ftype pv,
                                    ftype R_rho,ftype R_p,ftype L_rho,ftype L_p,
                                    ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    // two-shock approximation
    ftype gel = sqrt((GAMMA_G5 / L_rho) / (GAMMA_G6 * L_p + pv));
    ftype ger = sqrt((GAMMA_G5 / R_rho) / (GAMMA_G6 * R_p + pv));
    ftype x = (gel * L_p + ger * R_p - (v_line_R - v_line_L)) / (gel + ger);
    return x;
}


template <typename ftype=double>
static inline ftype guess_for_pressure(ftype dConstGamma,
                                       ftype R_rho,ftype R_p,ftype L_rho,ftype L_p,
                                       ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    ftype pmin, pmax;
    // start with the usual lowest-order guess for the contact wave pressure
    ftype pv = 0.5*(L_p+R_p) - 0.125*(v_line_R-v_line_L)*(L_p+R_p)*(cs_L+cs_R);

    pmin = min(L_p,R_p);
    pmax = max(L_p,R_p);

    // if one side is vacuum, guess half the mean
    if (pmin<=0.0)
        return 0.5*(pmin+pmax);

    // if the two are sufficiently close, and pv is between both values, return it
    ftype qrat = pmax / pmin;
    if (qrat <= 2.0 && (pmin <= pv && pv <= pmax)) {
        //printf("close\n");
        return pv;
    }

    if (pv < pmin) {
        // use two-rarefaction solution
        pv = guess_two_rarefaction(dConstGamma, R_rho, R_p, L_rho, L_p,
                                   v_line_L, v_line_R, cs_L, cs_R);
        //printf("rarefaction\n");
        return pv;
    }
    else {
        // two-shock approximation
        pv = guess_two_shock(dConstGamma, pv, R_rho, R_p, L_rho, L_p,
                             v_line_L, v_line_R, cs_L, cs_R);
        //printf("twoshock\n");
        if (pv < pmin || pv > pmax) {
            pv = pmin;
            //printf("twoshock 2\n");
        }
        return pv;
    }
}

template <>
inline dvec guess_for_pressure(dvec dConstGamma,
                               dvec R_rho,dvec R_p,dvec L_rho,dvec L_p,
                               dvec v_line_L, dvec v_line_R, dvec cs_L, dvec cs_R) {
    dvec pmin, pmax;
    // start with the usual lowest-order guess for the contact wave pressure
    dvec pv = 0.5*(L_p+R_p) - 0.125*(v_line_R-v_line_L)*(L_p+R_p)*(cs_L+cs_R);

    pmin = min(L_p,R_p);
    pmax = max(L_p,R_p);

    // if one side is vacuum, guess half the mean
    dvec zero = 0.0;
    auto vac = pmin<=zero;
    auto nvac = ~vac;
    pv = mask_mov(pv, vac, 0.5*(pmin+pmax));
    //dump(pv);

    // if the two are sufficiently close, and pv is between both values, return it
    dvec qrat = pmax / pmin;
    dvec two  = 2.0;
    auto done = nvac & (qrat <= two & (pv >= pmin & pv <= pmax));
    auto ndone = ~done & nvac;
    if (movemask(ndone)) {

        auto cond = pv<pmin;
        auto ncond = ~cond;
        if (movemask(cond)) {
            //printf("Rarefaction\n");
            pv = mask_mov(pv, cond&ndone,
                          guess_two_rarefaction(dConstGamma, R_rho, R_p, L_rho, L_p,
                                                v_line_L, v_line_R, cs_L, cs_R) );
            //dump(pv);
        }

        if (movemask(ncond)) {
            //printf("Two shock\n");
            // two-shock approximation
            dvec pshock = guess_two_shock(dConstGamma, pv, R_rho, R_p, L_rho, L_p,
                                          v_line_L, v_line_R, cs_L, cs_R);
            auto cond2 = pshock<pmin | pshock>pmax;
            pshock = mask_mov( pshock, ncond&cond2&ndone, pmin);
            //pshock = mask_mov( pshock, cond2, pmin); // This solve the problem but I do not get why...
            //dump(pshock);
            pv = mask_mov( pv, ndone&ncond, pshock);
            //dump(pv);
        }
    }
    return pv;
}

/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/*  This is the "normal" Riemann fan, with no vacuum on L or R state! */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
#ifndef USE_MFM
template <typename ftype=double>
static inline void sample_reimann_standard(ftype dConstGamma, ftype S,
        ftype R_rho,ftype R_p, ftype R_v[3],ftype L_rho,ftype L_p, ftype L_v[3],
        ftype P_M, ftype S_M, ftype *rho_f, ftype *p_f, ftype *v_f,
        ftype n_unit[3], ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    ftype C_eff,S_eff;
    if (S <= S_M) { /* sample point is left of contact discontinuity */
        if (P_M <= L_p) { /* left fan (rarefaction) */
            ftype S_check_L = v_line_L - cs_L;
            if (S <= S_check_L) { /* left data state */
                *p_f = L_p;
                *rho_f = L_rho;
                for (auto k=0; k<3; k++)
                    v_f[k] = L_v[k];
                return;
            }
            else {
                ftype C_eff_L = cs_L * pow(P_M / L_p, GAMMA_G1);
                ftype S_tmp_L = S_M - C_eff_L;

                if (S > S_tmp_L) { /* middle left state */
                    *rho_f = L_rho * pow(P_M / L_p, GAMMA_G8);
                    *p_f = P_M;
                    for (auto k=0; k<3; k++)
                        v_f[k] = L_v[k] + (S_M-v_line_L)*n_unit[k];
                    return;
                }
                else {      /* left state inside fan */
                    S_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
                    C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
                    *rho_f = L_rho * pow(C_eff / cs_L, GAMMA_G4);
                    *p_f = L_p * pow(C_eff / cs_L, GAMMA_G3);
                    for (auto k=0; k<3; k++)
                        v_f[k] = L_v[k] + (S_eff-v_line_L)*n_unit[k];
                    return;
                }
            }
        }
        else {          /* left shock */
            if (L_p > 0) {
                ftype pml = P_M / L_p;
                ftype S_L = v_line_L - cs_L * sqrt(GAMMA_G2 * pml + GAMMA_G1);

                if (S <= S_L) { /* left data state */
                    *p_f = L_p;
                    *rho_f = L_rho;
                    for (auto k=0; k<3; k++)
                        v_f[k] = L_v[k];
                    return;
                }
                else {      /* middle left state behind shock */
                    *rho_f = L_rho * (pml + GAMMA_G6) / (pml * GAMMA_G6 + 1.0);
                    *p_f = P_M;
                    for (auto k=0; k<3; k++)
                        v_f[k] = L_v[k] + (S_M-v_line_L)*n_unit[k];
                    return;
                }
            }
            else {
                *rho_f = L_rho / GAMMA_G6;
                *p_f = P_M;
                for (auto k=0; k<3; k++)
                    v_f[k] = L_v[k] + (S_M-v_line_L)*n_unit[k];
                return;
            }
        }
    }
    else {  /* sample point is right of contact discontinuity */
        if (P_M > R_p) { /* right shock */
            if (R_p > 0) {
                ftype pmr = P_M / R_p;
                ftype S_R = v_line_R + cs_R * sqrt(GAMMA_G2 * pmr + GAMMA_G1);

                if (S >= S_R) { /* right data state */
                    *p_f = R_p;
                    *rho_f = R_rho;
                    for (auto k=0; k<3; k++)
                        v_f[k] = R_v[k];
                    return;
                }
                else {      /* middle right state behind shock */
                    *rho_f = R_rho * (pmr + GAMMA_G6) / (pmr * GAMMA_G6 + 1.0);
                    *p_f = P_M;
                    for (auto k=0; k<3; k++)
                        v_f[k] = R_v[k] + (S_M-v_line_R)*n_unit[k];
                    return;
                }
            }
            else {
                *rho_f = R_rho / GAMMA_G6;
                *p_f = P_M;
                for (auto k=0; k<3; k++)
                    v_f[k] = R_v[k] + (S_M-v_line_R)*n_unit[k];
                return;
            }
        }
        else {          /* right fan */
            ftype S_check_R = v_line_R + cs_R;
            if (S >= S_check_R) {   /* right data state */
                *p_f = R_p;
                *rho_f = R_rho;
                for (auto k=0; k<3; k++)
                    v_f[k] = R_v[k];
                return;
            }
            else {
                ftype C_eff_R = cs_R * pow(P_M / R_p, GAMMA_G1);
                ftype S_tmp_R = S_M + C_eff_R;

                if (S <= S_tmp_R) { /* middle right state */
                    *rho_f = R_rho * pow(P_M / R_p, GAMMA_G8);
                    *p_f = P_M;
                    for (auto k=0; k<3; k++)
                        v_f[k] = R_v[k] + (S_M-v_line_R)*n_unit[k];
                    return;
                }
                else {      /* fan right state */
                    S_eff = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
                    C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
                    *rho_f = R_rho * pow(C_eff / cs_R, GAMMA_G4);
                    *p_f = R_p * pow(C_eff / cs_R, GAMMA_G3);
                    for (auto k=0; k<3; k++)
                        v_f[k] = R_v[k] + (S_eff-v_line_R)*n_unit[k];
                    return;
                }
            }
        }
    }
}
// TODO: otherwise, MFV can not use simd
template <>
inline void sample_reimann_standard(dvec dConstGamma, dvec S,
                                    dvec R_rho,dvec R_p, dvec R_v[3],dvec L_rho,dvec L_p, dvec L_v[3],
                                    dvec P_M, dvec S_M, dvec *rho_f, dvec *p_f, dvec *v_f,
                                    dvec n_unit[3], dvec v_line_L, dvec v_line_R, dvec cs_L, dvec cs_R) {}
#endif //!USE_MFM


/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/* right state is a vacuum, but left state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
template <typename ftype=double>
static inline void sample_reimann_vaccum_right(ftype dConstGamma, ftype S,
        ftype R_rho,ftype R_p, ftype R_v[3],ftype L_rho,ftype L_p, ftype L_v[3],
        ftype &P_M, ftype &S_M, ftype *rho_f, ftype *p_f, ftype *v_f,
        ftype n_unit[3], ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    //ftype S_L = v_line_L - GAMMA_G4 * cs_L;
    ftype S_L = v_line_L + GAMMA_G4 * cs_L; // above line was a sign error, caught by Bert Vandenbroucke
    //printf("Vaccuum right\n");
#ifdef USE_MFM
    /* in this code mode, we are -always- moving with the contact discontinuity so density flux = 0;
     this constrains where we reside in the solution fan */
    P_M = 0;
    S_M = S_L;
    return;
#endif

    if (S_L < S) {
        /* vacuum */
        P_M = 0;
        S_M = S_L;
        *rho_f = 0;
    }
    else {
        /* left fan */
        ftype S_L_check = v_line_L - cs_L;
        if (S_L_check < S) {
            /* rarefaction fan left state */
            ftype C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
            P_M = L_p * pow(C_eff / cs_L, GAMMA_G3);
            S_M = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
            *rho_f = L_rho * pow(C_eff / cs_L, GAMMA_G4);
        }
        else {
            /* left data state */
            P_M = L_p;
            S_M = v_line_L;
            *rho_f = L_rho;
        }
    }
    *p_f = P_M;
    for (auto k=0; k<3; k++)
        v_f[k] = L_v[k] + (S_M - v_line_L) * n_unit[k];
    return;
}



/* --------------------------------------------------------------------------------- */
/* part of exact Riemann solver: */
/* left state is a vacuum, but right state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
template <typename ftype=double>
static inline void sample_reimann_vaccum_left(ftype dConstGamma, ftype S,
        ftype R_rho,ftype R_p, ftype R_v[3],ftype L_rho,ftype L_p, ftype L_v[3],
        ftype &P_M, ftype &S_M, ftype *rho_f, ftype *p_f, ftype *v_f,
        ftype n_unit[3], ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    //printf("Vaccuum left\n");
    ftype S_R = v_line_R - GAMMA_G4 * cs_R;
#ifdef USE_MFM
    /* in this code mode, we are -always- moving with the contact discontinuity so density flux = 0;
     this constrains where we reside in the solution fan */
    P_M = 0;
    S_M = S_R;
    return;
#endif

    if (S_R > S) {
        /* vacuum */
        P_M = 0;
        S_M = S_R;
        *rho_f = 0;
    }
    else {
        /* right fan */
        ftype S_R_check = v_line_R + cs_R;
        if (S_R_check > S) {
            /* rarefaction fan right state */
            ftype C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
            P_M = R_p * pow(C_eff / cs_R, GAMMA_G3);
            S_M = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
            *rho_f = R_rho * pow(C_eff / cs_R, GAMMA_G4);
        }
        else {
            /* right data state */
            P_M = R_p;
            S_M = v_line_R;
            *rho_f = R_rho;
        }
    }
    *p_f = P_M;
    for (auto k=0; k<3; k++)
        v_f[k] = R_v[k] + (S_M - v_line_R) * n_unit[k];
    return;
}


template <typename ftype=double>
static inline void get_shock_f_fp(ftype dConstGamma, ftype &f, ftype &fp,
                                  ftype pg, ftype p, ftype r, ftype cs) {
    ftype a0, a1, a2;
    a0 = GAMMA_G5 / r;
    a1 = GAMMA_G6 * p;
    a2 = sqrt(a0 / (pg+a1));

    f = (pg-p) * a2;
    fp = a2 * (1.0 - 0.5*(pg-p)/(a1+pg));
}


template <typename ftype=double>
static inline void get_rarefaction_f_fp(ftype dConstGamma, ftype &f, ftype &fp,
                                        ftype pg, ftype p, ftype r, ftype cs) {
    ftype pratio, a0;
    pratio = pg / p;
    a0 = pow(pratio, GAMMA_G1);

    f = GAMMA_G4 * cs * (a0-1.);
    fp = a0 / (r*cs*pratio);
}

template <typename ftype=double>
static inline void get_f_fp(ftype dConstGamma, ftype &f, ftype &fp,
                            ftype pg, ftype p, ftype r, ftype cs) {
    if (pg>p) {
        /* shock wave */
        get_shock_f_fp(dConstGamma, f, fp, pg, p, r, cs);
    }
    else {
        /* rarefaction wave */
        get_rarefaction_f_fp(dConstGamma, f, fp, pg, p, r, cs);
    }
}

template <>
static inline void get_f_fp(dvec dConstGamma, dvec &f, dvec &fp,
                            dvec pg, dvec p, dvec r, dvec cs) {
    auto cond = pg>p;
    auto ncond = pg<=p; //~cond;
    dvec f1, fp1;

    //if (movemask(cond)) {
    get_shock_f_fp(dConstGamma, f1, fp1, pg, p, r, cs);
    f = mask_mov(f, cond, f1);
    fp = mask_mov(fp, cond, fp1);
    //}

    //if (movemask(ncond)) {
    get_rarefaction_f_fp(dConstGamma, f1, fp1, pg, p, r, cs);
    f = mask_mov(f, ncond, f1);
    fp = mask_mov(fp, ncond, fp1);
    //}
}

template <typename ftype=double>
static inline bool isConverged(const ftype &tol, const int &itt) {
    if (itt<NMAX_ITER) {
        if (tol>TOL_ITER) {
            return false;
        }
        else {
            return true;
        }
    }
    else {
        return true;
    }
}

template <>
static inline bool isConverged(const dvec &tol, const int &itt) {
    if (itt<NMAX_ITER) {
        auto notconverged = tol > TOL_ITER;
        if (movemask(notconverged)) {
            return false;
        }
        else {
            return true;
        }
    }
    else {
        return true;
    }
}

template <typename ftype=double>
static inline int newrap_solver(ftype dConstGamma, ftype &Pg, ftype &W_L, ftype &W_R,
                                ftype dvel,
                                ftype L_p, ftype L_rho, ftype cs_L,
                                ftype R_p, ftype R_rho, ftype cs_R) {
    int itt = 0;
    ftype tol = 100.;
    ftype Pg_prev, Z_L, Z_R;
    //while ((tol>TOL_ITER)&&(itt<NMAX_ITER)) {
//TODO: For now just a fixed iteration count is used
    //for (auto k=0; k<4; k++) {
    while ( !isConverged(tol, itt) ) {
        Pg_prev=Pg;
        get_f_fp(dConstGamma, W_L, Z_L, Pg, L_p, L_rho, cs_L);
        get_f_fp(dConstGamma, W_R, Z_R, Pg, R_p, R_rho, cs_R);

        if (itt < 10)
            Pg -= (W_L + W_R + dvel) / (Z_L + Z_R);
        else
            Pg -= 0.5 * (W_L + W_R + dvel) / (Z_L + Z_R);

        limit_decrease(Pg, Pg_prev);

        tol = 2.0 * abs((Pg-Pg_prev)/(Pg+Pg_prev));
        itt += 1;
    }
    return itt;
}

template <typename ftype=double>
static inline int iterative_Riemann_solver(ftype dConstGamma,
        ftype R_rho,ftype R_p,ftype L_rho,ftype L_p,
        ftype &P_M, ftype &S_M,
        ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    /* before going on, let's compare this to an exact Riemann solution calculated iteratively */
    ftype Pg, W_L, W_R;
    int niter_Riemann;
    ftype dvel;
    dvel = v_line_R - v_line_L;

    Pg = guess_for_pressure(dConstGamma, R_rho, R_p, L_rho, L_p, v_line_L, v_line_R, cs_L, cs_R);
    //dump(Pg);

    niter_Riemann = newrap_solver(dConstGamma, Pg, W_L, W_R, dvel,
                                  L_p, L_rho, cs_L, R_p, R_rho, cs_R);
    //if (niter_Riemann<NMAX_ITER)
    {
        P_M = Pg;
        S_M = 0.5*(v_line_L+v_line_R) + 0.5*(W_R-W_L);

        return niter_Riemann;
    }
    //else {
    //    assert(0);
    //}
}




template <typename ftype=double>
static inline void sample_reimann_vaccum_internal(ftype dConstGamma, ftype S,
        ftype R_rho,ftype R_p, ftype R_v[3],ftype L_rho,ftype L_p, ftype L_v[3],
        ftype &P_M, ftype &S_M, ftype *rho_f, ftype *p_f, ftype *v_f,
        ftype n_unit[3], ftype v_line_L, ftype v_line_R, ftype cs_L, ftype cs_R) {
    ftype dvel = v_line_R - v_line_L;
    ftype check_vel = GAMMA_G4 * (cs_R + cs_L) - dvel;
    if (check_vel >= 0) return;
    //printf("VACUUM\n");

    ftype S_L = v_line_L + GAMMA_G4 * cs_L;
    ftype S_R = v_line_R - GAMMA_G4 * cs_R;
    if (S <= S_L) {
        /* left fan */
        sample_reimann_vaccum_right(dConstGamma,  S,
                                    R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                                    P_M,  S_M,  rho_f,  p_f,  v_f,
                                    n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
    }
    else if (S >= S_R) {
        /* right fan */
        sample_reimann_vaccum_left(dConstGamma,  S,
                                   R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                                   P_M,  S_M,  rho_f,  p_f,  v_f,
                                   n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
    }
    else {
        /* vacuum in between */
        P_M = 0;
        S_M = S;
#ifndef USE_MFM
        *rho_f = 0;
        *p_f = P_M;
        for (auto k=0; k<3; k++)
            v_f[k] = (L_v[k] + (R_v[k]-L_v[k]) * (S-S_L)/(S_R-S_L)) *
                     (1-n_unit[k]) + S * n_unit[k];
#endif
    }
    //printf("%e %e\n", P_M, S_M);
}

// Template specialization correct only if using MFM
template <>
inline void sample_reimann_vaccum_internal(dvec dConstGamma, dvec S,
        dvec R_rho,dvec R_p, dvec R_v[3],dvec L_rho,dvec L_p, dvec L_v[3],
        dvec &P_M, dvec &S_M, dvec *rho_f, dvec *p_f, dvec *v_f,
        dvec n_unit[3], dvec v_line_L, dvec v_line_R, dvec cs_L, dvec cs_R) {
#ifndef USE_MFM
    assert(0);
#endif
    dvec dvel = v_line_R - v_line_L;
    dvec check_vel = GAMMA_G4 * (cs_R + cs_L) - dvel;
    auto vac = check_vel < 0;
    if (movemask(vac)) {
        //printf("vaccum internal\n");
        dvec S_L = v_line_L + GAMMA_G4 * cs_L;
        dvec S_R = v_line_R - GAMMA_G4 * cs_R;

        auto lfan = vac&(S<=S_L);
        auto rfan = vac&(S>=S_R);
        dvec zeros = 0.;

        P_M = mask_mov(P_M, vac, zeros);
        S_M = mask_mov(S_M, vac, zeros);
        S_M = mask_mov(S_M, lfan, S_L);
        S_M = mask_mov(S_M, rfan, S_R);
    }
}



/* -------------------------------------------------------------------------------------------------------------- */
/*  Part of exact Riemann solver: */
/*    take the face state we have calculated from the exact Riemann solution and get the corresponding fluxes */
/*   (written by V. Springel for AREPO, with minor modifications) */
/* -------------------------------------------------------------------------------------------------------------- */
template <typename ftype=double>
static void convert_face_to_flux(ftype dConstGamma,
                                 ftype *rho_f, ftype *p_f, ftype *v_f,
                                 ftype n_unit[3]) {
    ftype rho, P, v[3], v_line=0, v_frame=0, h=0;
    rho = *rho_f;
    P = *p_f;
    for (auto k=0; k<3; k++) {
        v[k] = v_f[k];
        v_line += v[k] * n_unit[k];
        h += v[k] * v[k];
    }
    v_line -= v_frame;
    h *= 0.5 * rho; /* h is the kinetic energy density */
    h += (GAMMA/GAMMA_MINUS1) * P; /* now h is the enthalpy */
    /* now we just compute the standard fluxes for a given face state */
    *p_f = h * v_line;
    *rho_f = rho * v_line;
    for (auto k=0; k<3; k++)
        v_f[k] = (*rho_f) * v[k] + P * n_unit[k];
    return;
}

template<typename ftype=double>
static inline int handle_input_vacuum(ftype R_rho, ftype R_p,
                                      ftype L_rho, ftype L_p,
                                      ftype &P_M, ftype &S_M, ftype *rho_f,
                                      ftype *p_f, ftype *v_f) {
    if ((L_p == 0 && R_p == 0) || (L_rho==0 && R_rho==0)) {
        //printf("Vacuum!!\n");
        // we're in a Vaccuum!
        P_M = 0.;
        S_M = 0.;
#ifndef USE_MFM
        // *rho_f = *p_f = v_f[0] = v_f[1] = v_f[2] = 0;
        *rho_f = 0.;
        *p_f = 0.;
        v_f[0] = 0.;
        v_f[1] = 0.;
        v_f[2] = 0.;
#endif
        return 0;
    }
    else {
        return 1;
    }
}

template<>
inline int handle_input_vacuum(dvec R_rho, dvec R_p,
                               dvec L_rho, dvec L_p,
                               dvec &P_M, dvec &S_M, dvec *rho_f,
                               dvec *p_f, dvec *v_f) {
    dvec zeros = 0.0;
    auto vac = ((L_p <= zeros) & (R_p <= zeros)) | ((L_rho<=zeros) & (R_rho<=zeros));
    //dump(L_p);
    P_M = mask_mov(P_M, vac, zeros);
    S_M = mask_mov(S_M, vac, zeros);

    int notvac = movemask(~vac);
    if (movemask(vac)) {
        //printf("vaccum input\n");
        //dump(P_M);
        //dump(S_M);
    }

    if (notvac) {    // At least one non-vacuum
        return 1;
    }
    else {           // All vacuum, the iterative solver can be skipped
        return 0;
    }
}

template <typename ftype=double>
static inline int Riemann_solver_exact(ftype dConstGamma,
                                       ftype R_rho,ftype R_p, ftype R_v[3], ftype L_rho,ftype L_p, ftype L_v[3],
                                       ftype &P_M, ftype &S_M, ftype *rho_f, ftype *p_f, ftype *v_f,
                                       ftype n_unit[3]) {

    int niter;
    ftype cs_L = sqrt(GAMMA * L_p / L_rho);
    ftype cs_R = sqrt(GAMMA * R_p / R_rho);

    ftype v_line_L = L_v[0]*n_unit[0] +
                     L_v[1]*n_unit[1] +
                     L_v[2]*n_unit[2];

    ftype v_line_R = R_v[0]*n_unit[0] +
                     R_v[1]*n_unit[1] +
                     R_v[2]*n_unit[2];
    //printf("Going for the exact Riemann Solver \n");
    /* first, we need to check for all the special/exceptional cases that will cause things to go haywire */
    //ftype niter = 0;

    /* the usual situation is here:: */
    //if ((L_rho > 0) && (R_rho > 0))
    {
        niter=iterative_Riemann_solver(dConstGamma, R_rho, R_p, L_rho, L_p, P_M, S_M, v_line_L, v_line_R, cs_L, cs_R);
        //if (niter)
        {
//           printf("Normal Riemann solution %e \n", Riemann_out->Fluxes.rho);
            /* this is the 'normal' Reimann solution */

#ifndef USE_MFM
            ftype S = 0;
            sample_reimann_standard(dConstGamma, S,
                                    R_rho, R_p, R_v,
                                    L_rho, L_p, L_v,
                                    P_M, S_M,
                                    rho_f, p_f, v_f,
                                    n_unit,v_line_L,v_line_R,cs_L,cs_R);
#endif
            //
//           printf("Normal Riemann solution %e \n", Riemann_out->Fluxes.rho);
            //}
            //else {
            // ICs lead to vacuum, need to sample vacuum solution
            //sample_reimann_vaccum_internal(dConstGamma, 0.0,R_rho, R_p, R_v, L_rho, L_p, L_v, P_M, S_M, rho_f, p_f, v_f, n_unit,v_line_L,v_line_R,cs_L,cs_R);
        }
        ftype zero = 0.;
        sample_reimann_vaccum_internal(dConstGamma, zero, R_rho, R_p, R_v, L_rho, L_p, L_v, P_M, S_M, rho_f, p_f, v_f, n_unit,v_line_L,v_line_R,cs_L,cs_R);
        handle_input_vacuum(R_rho, R_p, L_rho, L_p, P_M, S_M, rho_f, p_f, v_f);
        //if (!proceed) return 0;
    }
    /*
    else {
        // one of the densities is zero or negative
        if ((L_rho<0)||(R_rho<0))
            //assert(0);
            if (L_rho>0) {
                sample_reimann_vaccum_right(dConstGamma,  0.0,
                                            R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                                            P_M,  S_M,  rho_f,  p_f,  v_f,
                                            n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
            }
        if (R_rho>0) {
            sample_reimann_vaccum_left(dConstGamma, 0.0,
                                       R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                                       P_M,  S_M,  rho_f,  p_f,  v_f,
                                       n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
        }
    }
    */
#ifndef USE_MFM
    /* if we got a valid solution, this solver returns face states: need to convert these to fluxes */
    convert_face_to_flux(dConstGamma, rho_f, p_f, v_f, n_unit);
#endif
    return niter;
}
