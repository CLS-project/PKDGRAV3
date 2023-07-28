#include "core/vmath.h"

#define TOL_ITER 1.e-6
#define NMAX_ITER 100

template <typename dtype=double>
static inline void dump(dtype a) {
    printf("%e\n", a);
}
template <>
static inline void dump(dvec a) {
    for (auto i=0; i<dvec::width(); i++)
        printf(" %e ", a[i]);
    printf("\n");
}

template <typename dtype>
static inline int nan_guard(dtype &var, dtype def) {
    int nan = 0;
    var = mask_mov(var, var!=var, def);
    for (auto j=0; j<dtype::width(); j++) {
        if (isnan(var[j])) {
            nan = 1;
        }
    }
    return nan;
}

template <>
inline int nan_guard(double &var, double def) {
    if (isnan(var)) {
        var=def;
        return 1;
    }
    else {
        return 0;
    }
}

template <typename dtype=dvec, typename mtype=dmask>
class RiemannSolverExact {
public:

    RiemannSolverExact(dtype gamma) : gamma(gamma) {
        G1 = (gamma-1.0)/(2.0*gamma);
        G2 = (gamma+1.0)/(2.0*gamma);
        G3 = (2.0*gamma/(gamma-1.0));
        G4 = 2.0/(gamma-1.0);
        G5 = 2.0/(gamma+1.0);
        G6 = (gamma-1.0)/(gamma+1.0);
        G7 = 0.5*(gamma-1.0);
        G8 = 1.0/gamma;
        G9 = gamma-1.0;
    };

private:
    dtype gamma;

    dtype G1;
    dtype G2;
    dtype G3;
    dtype G4;
    dtype G5;
    dtype G6;
    dtype G7;
    dtype G8;
    dtype G9;

    struct Conserved_var_Riemann { //IA: TODO: change the name of this struct..!!
        dtype rho;
        dtype p;  // IA: Pressure is a primitive variable, not a conserved one...
        dtype v[3]; // IA: same...
        dtype u;
        dtype cs;
    };

    struct Input_vec_Riemann {
        struct Conserved_var_Riemann L;
        struct Conserved_var_Riemann R;
    };

    struct Riemann_outputs {
        dtype P_M;
        dtype S_M;
        struct Conserved_var_Riemann Fluxes;
    };

    inline void limit_decrease(dtype &Pg, dtype Pg_prev) {
        dtype limit = 0.1*Pg_prev;
        Pg = mask_mov(Pg, Pg<limit, limit);
    }

    inline dtype guess_two_rarefaction(
        dtype R_rho,dtype R_p,dtype L_rho,dtype L_p,
        dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
        dtype pnu = (cs_L+cs_R) - G7 * (v_line_R - v_line_L);
        dtype pde = cs_L * pow(L_p, -G1) + cs_R * pow(R_p, -G1);
        return pow(pnu / pde, G3);
    }

    inline dtype guess_two_shock( dtype pv,
                                  dtype R_rho,dtype R_p,dtype L_rho,dtype L_p,
                                  dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
        // two-shock approximation
        dtype gel = sqrt((G5 / L_rho) / (G6 * L_p + pv));
        dtype ger = sqrt((G5 / R_rho) / (G6 * R_p + pv));
        dtype x = (gel * L_p + ger * R_p - (v_line_R - v_line_L)) / (gel + ger);
        return x;
    }

    inline dtype guess_for_pressure(
        dtype R_rho,dtype R_p,dtype L_rho,dtype L_p,
        dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
        dtype pmin, pmax;
        // start with the usual lowest-order guess for the contact wave pressure
        dtype pv = 0.5*(L_p+R_p) - 0.125*(v_line_R-v_line_L)*(L_p+R_p)*(cs_L+cs_R);

        pmin = min(L_p,R_p);
        pmax = max(L_p,R_p);

        // if one side is vacuum, guess half the mean
        dtype zero = 0.0;
        mtype vac = pmin<=zero;
        mtype nvac = ~vac;
        pv = mask_mov(pv, vac, 0.5*(pmin+pmax));
        //dump(pv);

        // if the two are sufficiently close, and pv is between both values, return it
        dtype qrat = pmax / pmin;
        dtype two  = 2.0;
        // All this casting can be made static or defining custom and/or/not functions for vector
        // Or overloading when one of the arguments is an int
        mtype done = nvac & (mtype)((mtype)(qrat <= two) & (mtype)((mtype)(pv >= pmin) & (mtype)(pv <= pmax)));
        mtype ndone = ~done & nvac;
        if (movemask(ndone)) {

            mtype cond = pv<pmin;
            mtype ncond = ~cond;
            if (movemask(cond)) {
                //printf("Rarefaction\n");
                pv = mask_mov(pv, cond&ndone,
                              guess_two_rarefaction( R_rho, R_p, L_rho, L_p,
                                                     v_line_L, v_line_R, cs_L, cs_R) );
                //dump(pv);
            }

            if (movemask(ncond)) {
                //printf("Two shock\n");
                // two-shock approximation
                dtype pshock = guess_two_shock( pv, R_rho, R_p, L_rho, L_p,
                                                v_line_L, v_line_R, cs_L, cs_R);
                mtype cond2 = (mtype)(pshock<pmin) | (mtype)(pshock>pmax);
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

    inline void sample_reimann_standard( dtype S,
                                         dtype R_rho,dtype R_p, dtype R_v[3],dtype L_rho,dtype L_p, dtype L_v[3],
                                         dtype P_M, dtype S_M, dtype *rho_f_out, dtype *p_f_out, dtype *v_f_out,
                                         dtype n_unit[3], dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
        // This is a very inefficient simd implementation until someone bothers to
        // write it properly (including testing for all the branches!)
        typename dtype::scalar_t rho_f[dtype::width()];
        typename dtype::scalar_t p_f[dtype::width()];
        typename dtype::scalar_t v_f[3][dtype::width()];
        for (auto v=0; v < dtype::width(); v++) {
            typename dtype::scalar_t C_eff,S_eff;
            if (S[v] <= S_M[v]) { /* sample point is left of contact discontinuity */
                if (P_M[v] <= L_p[v]) { /* left fan (rarefaction) */
                    typename dtype::scalar_t S_check_L = v_line_L[v] - cs_L[v];
                    if (S[v] <= S_check_L) { /* left data state */
                        p_f[v] = L_p[v];
                        rho_f[v] = L_rho[v];
                        for (auto k=0; k<3; k++)
                            v_f[k][v] = L_v[k][v];
                        continue;
                    }
                    else {
                        typename dtype::scalar_t C_eff_L = cs_L[v] * pow(P_M[v] / L_p[v], G1[v]);
                        typename dtype::scalar_t S_tmp_L = S_M[v] - C_eff_L;

                        if (S[v] > S_tmp_L) { /* middle left state */
                            rho_f[v] = L_rho[v] * pow(P_M[v] / L_p[v], G8[v]);
                            p_f[v] = P_M[v];
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = L_v[k][v] + (S_M[v]-v_line_L[v])*n_unit[k][v];
                            continue;
                        }
                        else {      /* left state inside fan */
                            S_eff = G5[v] * (cs_L[v] + G7[v] * v_line_L[v] + S[v]);
                            C_eff = G5[v] * (cs_L[v] + G7[v] * (v_line_L[v] - S[v]));
                            rho_f[v] = L_rho[v] * pow(C_eff / cs_L[v], G4[v]);
                            p_f[v] = L_p[v] * pow(C_eff / cs_L[v], G3[v]);
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = L_v[k][v] + (S_eff-v_line_L[v])*n_unit[k][v];
                            continue;
                        }
                    }
                }
                else {          /* left shock */
                    if (L_p[v] > 0) {
                        typename dtype::scalar_t pml = P_M[v] / L_p[v];
                        typename dtype::scalar_t S_L = v_line_L[v] - cs_L[v] * sqrt(G2[v] * pml + G1[v]);

                        if (S[v] <= S_L) { /* left data state */
                            p_f[v] = L_p[v];
                            rho_f[v] = L_rho[v];
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = L_v[k][v];
                            continue;
                        }
                        else {      /* middle left state behind shock */
                            rho_f[v] = L_rho[v] * (pml + G6[v]) / (pml * G6[v] + 1.0);
                            p_f[v] = P_M[v];
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = L_v[k][v] + (S_M[v]-v_line_L[v])*n_unit[k][v];
                            continue;
                        }
                    }
                    else {
                        rho_f[v] = L_rho[v] / G6[v];
                        p_f[v] = P_M[v];
                        for (auto k=0; k<3; k++)
                            v_f[k][v] = L_v[k][v] + (S_M[v]-v_line_L[v])*n_unit[k][v];
                        continue;
                    }
                }
            }
            else {  /* sample point is right of contact discontinuity */
                if (P_M[v] > R_p[v]) { /* right shock */
                    if (R_p[v] > 0) {
                        typename dtype::scalar_t pmr = P_M[v] / R_p[v];
                        typename dtype::scalar_t S_R = v_line_R[v] + cs_R[v] * sqrt(G2[v] * pmr + G1[v]);

                        if (S[v] >= S_R) { /* right data state */
                            p_f[v] = R_p[v];
                            rho_f[v] = R_rho[v];
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = R_v[k][v];
                            continue;
                        }
                        else {      /* middle right state behind shock */
                            rho_f[v] = R_rho[v] * (pmr + G6[v]) / (pmr * G6[v] + 1.0);
                            p_f[v] = P_M[v];
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = R_v[k][v] + (S_M[v]-v_line_R[v])*n_unit[k][v];
                            continue;
                        }
                    }
                    else {
                        rho_f[v] = R_rho[v] / G6[v];
                        p_f[v] = P_M[v];
                        for (auto k=0; k<3; k++)
                            v_f[k][v] = R_v[k][v] + (S_M[v]-v_line_R[v])*n_unit[k][v];
                        continue;
                    }
                }
                else {          /* right fan */
                    typename dtype::scalar_t S_check_R = v_line_R[v] + cs_R[v];
                    if (S[v] >= S_check_R) {   /* right data state */
                        p_f[v] = R_p[v];
                        rho_f[v] = R_rho[v];
                        for (auto k=0; k<3; k++)
                            v_f[k][v] = R_v[k][v];
                        continue;
                    }
                    else {
                        typename dtype::scalar_t C_eff_R = cs_R[v] * pow(P_M[v] / R_p[v], G1[v]);
                        typename dtype::scalar_t S_tmp_R = S_M[v] + C_eff_R;

                        if (S[v] <= S_tmp_R) { /* middle right state */
                            rho_f[v] = R_rho[v] * pow(P_M[v] / R_p[v], G8[v]);
                            p_f[v] = P_M[v];
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = R_v[k][v] + (S_M[v]-v_line_R[v])*n_unit[k][v];
                            continue;
                        }
                        else {      /* fan right state */
                            S_eff = G5[v] * (-cs_R[v] + G7[v] * v_line_R[v] + S[v]);
                            C_eff = G5[v] * (cs_R[v] - G7[v] * (v_line_R [v]- S[v]));
                            rho_f[v] = R_rho[v] * pow(C_eff / cs_R[v], G4[v]);
                            p_f[v] = R_p[v] * pow(C_eff / cs_R[v], G3[v]);
                            for (auto k=0; k<3; k++)
                                v_f[k][v] = R_v[k][v] + (S_eff-v_line_R[v])*n_unit[k][v];
                            continue;
                        }
                    }
                }
            }
        }
        rho_f_out->load(rho_f);
        p_f_out->load(p_f);
        for (auto k=0; k<3; k++)
            v_f_out[k].load(v_f[k]);

    }
// TODO: otherwise, MFV can not use simd

#else

    inline void sample_reimann_standard( dtype S,
                                         dtype R_rho,dtype R_p, dtype R_v[3],dtype L_rho,dtype L_p, dtype L_v[3],
                                         dtype P_M, dtype S_M, dtype *rho_f, dtype *p_f, dtype *v_f,
                                         dtype n_unit[3], dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {}
#endif //!USE_MFM

    /* --------------------------------------------------------------------------------- */
    /* Part of exact Riemann solver: */
    /* right state is a vacuum, but left state is not: sample the fan appropriately */
    /*  (written by V. Springel for AREPO) */
    /* --------------------------------------------------------------------------------- */

    inline void sample_reimann_vaccum_right( dtype S,
            dtype R_rho,dtype R_p, dtype R_v[3],dtype L_rho,dtype L_p, dtype L_v[3],
            dtype &P_M, dtype &S_M, dtype *rho_f, dtype *p_f, dtype *v_f,
            dtype n_unit[3], dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
        //dtype S_L = v_line_L - G4 * cs_L;
        dtype S_L = v_line_L + G4 * cs_L; // above line was a sign error, caught by Bert Vandenbroucke
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
            dtype S_L_check = v_line_L - cs_L;
            if (S_L_check < S) {
                /* rarefaction fan left state */
                dtype C_eff = G5 * (cs_L + G7 * (v_line_L - S));
                P_M = L_p * pow(C_eff / cs_L, G3);
                S_M = G5 * (cs_L + G7 * v_line_L + S);
                *rho_f = L_rho * pow(C_eff / cs_L, G4);
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
    inline void sample_reimann_vaccum_left( dtype S,
                                            dtype R_rho,dtype R_p, dtype R_v[3],dtype L_rho,dtype L_p, dtype L_v[3],
                                            dtype &P_M, dtype &S_M, dtype *rho_f, dtype *p_f, dtype *v_f,
                                            dtype n_unit[3], dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
        //printf("Vaccuum left\n");
        dtype S_R = v_line_R - G4 * cs_R;
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
            dtype S_R_check = v_line_R + cs_R;
            if (S_R_check > S) {
                /* rarefaction fan right state */
                dtype C_eff = G5 * (cs_R - G7 * (v_line_R - S));
                P_M = R_p * pow(C_eff / cs_R, G3);
                S_M = G5 * (-cs_R + G7 * v_line_R + S);
                *rho_f = R_rho * pow(C_eff / cs_R, G4);
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

    inline void get_shock_f_fp( dtype &f, dtype &fp,
                                dtype pg, dtype p, dtype r, dtype cs) {
        dtype a0, a1, a2;
        a0 = G5 / r;
        a1 = G6 * p;
        a2 = sqrt(a0 / (pg+a1));

        f = (pg-p) * a2;
        fp = a2 * (1.0 - 0.5*(pg-p)/(a1+pg));
    }

    inline void get_rarefaction_f_fp( dtype &f, dtype &fp,
                                      dtype pg, dtype p, dtype r, dtype cs) {
        dtype pratio, a0;
        pratio = pg / p;
        a0 = pow(pratio, G1);

        f = G4 * cs * (a0-1.);
        fp = a0 / (r*cs*pratio);
    }

    /*

     inline void get_f_fp( dtype &f, dtype &fp,
                                dtype pg, dtype p, dtype r, dtype cs) {
        if (pg>p) {
            // shock wave
            get_shock_f_fp( f, fp, pg, p, r, cs);
        }
        else {
            // rarefaction wave
            get_rarefaction_f_fp( f, fp, pg, p, r, cs);
        }
    }
    */

    inline void get_f_fp( dtype &f, dtype &fp,
                          dtype pg, dtype p, dtype r, dtype cs) {
        mtype cond = pg>p;
        mtype ncond = pg<=p; //~cond;
        dtype f1, fp1;

        //if (movemask(cond)) {
        get_shock_f_fp( f1, fp1, pg, p, r, cs);
        f = mask_mov(f, cond, f1);
        fp = mask_mov(fp, cond, fp1);
        //}

        //if (movemask(ncond)) {
        get_rarefaction_f_fp( f1, fp1, pg, p, r, cs);
        f = mask_mov(f, ncond, f1);
        fp = mask_mov(fp, ncond, fp1);
        //}
    }

    inline bool isConverged(const dtype &tol, const int &itt) {
        if (itt<NMAX_ITER) {
            mtype notconverged = tol > TOL_ITER;
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

    inline int newrap_solver( dtype &Pg, dtype &W_L, dtype &W_R,
                              dtype dvel,
                              dtype L_p, dtype L_rho, dtype cs_L,
                              dtype R_p, dtype R_rho, dtype cs_R) {
        int itt = 0;
        dtype tol = 100.;
        dtype Pg_prev, Z_L, Z_R;
        //while ((tol>TOL_ITER)&&(itt<NMAX_ITER)) {
//TODO: For now just a fixed iteration count is used
        //for (auto k=0; k<4; k++) {
        while ( !isConverged(tol, itt) ) {
            Pg_prev=Pg;
            get_f_fp( W_L, Z_L, Pg, L_p, L_rho, cs_L);
            get_f_fp( W_R, Z_R, Pg, R_p, R_rho, cs_R);

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

    inline int iterative_Riemann_solver(
        dtype R_rho,dtype R_p,dtype L_rho,dtype L_p,
        dtype &P_M, dtype &S_M,
        dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
        /* before going on, let's compare this to an exact Riemann solution calculated iteratively */
        dtype Pg, W_L, W_R;
        int niter_Riemann;
        dtype dvel;
        dvel = v_line_R - v_line_L;

        Pg = this->guess_for_pressure( R_rho, R_p, L_rho, L_p, v_line_L, v_line_R, cs_L, cs_R);
        //dump(Pg);

        niter_Riemann = newrap_solver( Pg, W_L, W_R, dvel,
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

    inline void sample_reimann_vaccum_internal( dtype S,
            dtype R_rho,dtype R_p, dtype R_v[3],dtype L_rho,dtype L_p, dtype L_v[3],
            dtype &P_M, dtype &S_M, dtype *rho_f, dtype *p_f, dtype *v_f,
            dtype n_unit[3], dtype v_line_L, dtype v_line_R, dtype cs_L, dtype cs_R) {
#ifndef USE_MFM
        //assert(0);
#endif
        dtype dvel = v_line_R - v_line_L;
        dtype check_vel = G4 * (cs_R + cs_L) - dvel;
        mtype vac = check_vel < 0;
        if (movemask(vac)) {
            //printf("vaccum internal\n");
            dtype S_L = v_line_L + G4 * cs_L;
            dtype S_R = v_line_R - G4 * cs_R;

            mtype lfan = vac&(S<=S_L);
            mtype rfan = vac&(S>=S_R);
            dtype zeros = 0.;

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

    void convert_face_to_flux(
        dtype *rho_f, dtype *p_f, dtype *v_f,
        dtype n_unit[3]) {
        dtype rho, P, v[3], v_line=0, v_frame=0, h=0;
        rho = *rho_f;
        P = *p_f;
        for (auto k=0; k<3; k++) {
            v[k] = v_f[k];
            v_line += v[k] * n_unit[k];
            h += v[k] * v[k];
        }
        v_line -= v_frame;
        h *= 0.5 * rho; /* h is the kinetic energy density */
        h += (gamma/G9) * P; /* now h is the enthalpy */
        /* now we just compute the standard fluxes for a given face state */
        *p_f = h * v_line;
        *rho_f = rho * v_line;
        for (auto k=0; k<3; k++)
            v_f[k] = (*rho_f) * v[k] + P * n_unit[k];
        return;
    }

    inline int handle_input_vacuum(dtype R_rho, dtype R_p,
                                   dtype L_rho, dtype L_p,
                                   dtype &P_M, dtype &S_M, dtype *rho_f,
                                   dtype *p_f, dtype *v_f) {
        dtype zeros = 0.0;
        mtype vac = (mtype)((mtype)(L_p <= zeros) & (mtype)(R_p <= zeros)) | (mtype)((mtype)(L_rho<=zeros) & (mtype)(R_rho<=zeros));
        //dump(L_p);
        P_M = mask_mov(P_M, vac, zeros);
        S_M = mask_mov(S_M, vac, zeros);

        if (movemask(vac)) {
            //printf("vaccum input\n");
            //dump(P_M);
            //dump(S_M);
        }

        mtype nvac = ~vac;
        if (movemask(nvac)) {    // At least one non-vacuum
            return 1;
        }
        else {           // All vacuum, the iterative solver can be skipped
            return 0;
        }
    }

public:

    int solve(
        dtype R_rho,dtype R_p, dtype R_v[3], dtype L_rho,dtype L_p, dtype L_v[3],
        dtype &P_M, dtype &S_M, dtype *rho_f, dtype *p_f, dtype *v_f,
        dtype n_unit[3]) {
        int niter;
        dtype cs_L = sqrt(gamma * L_p / L_rho);
        dtype cs_R = sqrt(gamma * R_p / R_rho);

        dtype v_line_L = L_v[0]*n_unit[0] +
                         L_v[1]*n_unit[1] +
                         L_v[2]*n_unit[2];

        dtype v_line_R = R_v[0]*n_unit[0] +
                         R_v[1]*n_unit[1] +
                         R_v[2]*n_unit[2];
        //printf("Going for the exact Riemann Solver \n");
        /* first, we need to check for all the special/exceptional cases that will cause things to go haywire */
        //dtype niter = 0;

        /* the usual situation is here:: */
        //if ((L_rho > 0) && (R_rho > 0))
        {
            niter=iterative_Riemann_solver( R_rho, R_p, L_rho, L_p, P_M, S_M, v_line_L, v_line_R, cs_L, cs_R);
            //if (niter)
            {
//           printf("Normal Riemann solution %e \n", Riemann_out->Fluxes.rho);
                /* this is the 'normal' Reimann solution */

#ifndef USE_MFM
                dtype S = 0;
                sample_reimann_standard( S,
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
                //sample_reimann_vaccum_internal( 0.0,R_rho, R_p, R_v, L_rho, L_p, L_v, P_M, S_M, rho_f, p_f, v_f, n_unit,v_line_L,v_line_R,cs_L,cs_R);
            }
            dtype zero = 0.;
            sample_reimann_vaccum_internal( zero, R_rho, R_p, R_v, L_rho, L_p, L_v, P_M, S_M, rho_f, p_f, v_f, n_unit,v_line_L,v_line_R,cs_L,cs_R);
            handle_input_vacuum(R_rho, R_p, L_rho, L_p, P_M, S_M, rho_f, p_f, v_f);
            //if (!proceed) return 0;
        }
        /*
        else {
            // one of the densities is zero or negative
            if ((L_rho<0)||(R_rho<0))
                //assert(0);
                if (L_rho>0) {
                    sample_reimann_vaccum_right(  0.0,
                                                R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                                                P_M,  S_M,  rho_f,  p_f,  v_f,
                                                n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
                }
            if (R_rho>0) {
                sample_reimann_vaccum_left( 0.0,
                                           R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                                           P_M,  S_M,  rho_f,  p_f,  v_f,
                                           n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
            }
        }
        */
#ifndef USE_MFM
        /* if we got a valid solution, this solver returns face states: need to convert these to fluxes */
        convert_face_to_flux( rho_f, p_f, v_f, n_unit);
#endif
        return niter;
    }

};
