#define GAMMA (pkd->param.dConstGamma)//1.6666666666666666666 // 1.4 //1.666666666666
#define GAMMA_MINUS1 (pkd->param.dConstGamma-1.0) //2./3. //0.4 //0.66666666
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
#define NMAX_ITER 1000

struct Conserved_var_Riemann //IA: TODO: change the name of this struct..!!
{
    double rho;
    double p;  // IA: Pressure is a primitive variable, not a conserved one...
    double v[3]; // IA: same...
    double u;
    double cs;
};

struct Input_vec_Riemann
{
    struct Conserved_var_Riemann L;
    struct Conserved_var_Riemann R;
};
struct Riemann_outputs
{
    double P_M;
    double S_M;
    struct Conserved_var_Riemann Fluxes;
};

INLINE double DMAX(double a, double b) { return (a > b) ? a : b; }
INLINE double DMIN(double a, double b) { return (a < b) ? a : b; }
INLINE int Riemann_solver_exact(PKD pkd, 
      double R_rho, double R_p, double R_v[3], double L_rho,double L_p, double L_v[3],
      double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f,
      double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R);

INLINE int iterative_Riemann_solver(PKD pkd, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double v_line_L, double v_line_R, double cs_L, double cs_R);

INLINE void convert_face_to_flux(PKD pkd, 
      double *rho_f, double *p_f, double *v_f, 
      double n_unit[3]);


INLINE double guess_for_pressure(PKD pkd, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double v_line_L, double v_line_R, double cs_L, double cs_R);

INLINE void sample_reimann_vaccum_internal(PKD pkd, double S, 
      double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
      double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f,
      double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);

INLINE void sample_reimann_standard(PKD pkd, double S, 
      double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
      double P_M, double S_M, double *rho_f, double *p_f, double *v_f, 
      double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);

INLINE void sample_reimann_vaccum_right(PKD pkd, double S, 
                                double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
                                double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f, 
                                double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);

INLINE void sample_reimann_vaccum_left(PKD pkd, double S, 
                                double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
                                double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f, 
                                double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);

inline int Riemann_solver_exact(PKD pkd, 
            double R_rho,double R_p, double R_v[3], double L_rho,double L_p, double L_v[3],
            double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f,
            double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R)
{
    //printf("Going for the exact Riemann Solver \n");
    /* first, we need to check for all the special/exceptional cases that will cause things to go haywire */
   int niter = 0;
    if((L_p == 0 && L_p == 0) || (L_rho==0 && R_rho==0))
    {
        /* we're in a Vaccuum! */
        *P_M = 0.;
        *S_M = 0.;
#ifndef USE_MFM
        //*rho_f = *p_f = v_f[0] = v_f[1] = v_f[2] = 0;
        *rho_f = 0.;
        *p_f = 0.;
        v_f[0] = 0.;
        v_f[1] = 0.;
        v_f[2] = 0.;
#endif
        //printf("VACCUUM!\n");
        return 0;
    }
    /* the usual situation is here:: */
    if((L_rho > 0) && (R_rho > 0))
    {
        niter=iterative_Riemann_solver(pkd, R_rho, R_p, L_rho, L_p, P_M, S_M, v_line_L, v_line_R, cs_L, cs_R);
        if(niter)
        {
//           printf("Normal Riemann solution %e \n", Riemann_out->Fluxes.rho);
            /* this is the 'normal' Reimann solution */

#ifndef USE_MFM
            sample_reimann_standard(pkd, 0.0,
                  R_rho, R_p, R_v,
                  L_rho, L_p, L_v,
                  *P_M, *S_M,
                  rho_f, p_f, v_f,
                  n_unit,v_line_L,v_line_R,cs_L,cs_R);
#endif
            //
//           printf("Normal Riemann solution %e \n", Riemann_out->Fluxes.rho);
        }
        else
        {
            /* ICs lead to vacuum, need to sample vacuum solution */
            sample_reimann_vaccum_internal(pkd, 0.0,R_rho, R_p, R_v, L_rho, L_p, L_v, P_M, S_M, rho_f, p_f, v_f, n_unit,v_line_L,v_line_R,cs_L,cs_R);
        }
    } else {
        /* one of the densities is zero or negative */
        if((L_rho<0)||(R_rho<0))
            assert(0);
        if(L_rho>0){
             sample_reimann_vaccum_right(pkd,  0.0,
                      R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                      P_M,  S_M,  rho_f,  p_f,  v_f,
                      n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
        }
        if(R_rho>0){
             sample_reimann_vaccum_left(pkd, 0.0,
                      R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                      P_M,  S_M,  rho_f,  p_f,  v_f,
                      n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
        }
    }
#ifndef USE_MFM
    /* if we got a valid solution, this solver returns face states: need to convert these to fluxes */
    convert_face_to_flux(pkd, rho_f, p_f, v_f, n_unit);
#endif
    return niter;
}


inline int iterative_Riemann_solver(PKD pkd, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    /* before going on, let's compare this to an exact Riemann solution calculated iteratively */
    double Pg,Pg_prev,W_L,W_R,Z_L,Z_R,tol,pratio; int niter_Riemann=0;
    double a0,a1,a2,dvel,check_vel;
    dvel = v_line_R - v_line_L;
    check_vel = GAMMA_G4 * (cs_R + cs_L) - dvel;
    /* if check_vel<0, this will produce a vacuum: need to use vacuum-specific subroutine */
    if(check_vel < 0) return 0;
    
    tol=100.0;
    Pg = guess_for_pressure(pkd, R_rho, R_p, L_rho, L_p, P_M, S_M, v_line_L, v_line_R, cs_L, cs_R);
    while((tol>TOL_ITER)&&(niter_Riemann<NMAX_ITER))
    {
        Pg_prev=Pg;
        if(Pg>L_p)
        {
            /* shock wave */
            a0 = GAMMA_G5 / L_rho;
            a1 = GAMMA_G6 * L_p;
            a2 = sqrt(a0 / (Pg+a1));
            W_L = (Pg-L_p) * a2;
            Z_L = a2 * (1.0 - 0.5*(Pg-L_p)/(a1+Pg));
        } else {
            /* rarefaction wave */
            pratio = Pg / L_p;
            W_L = GAMMA_G4 * cs_L * (pow(pratio, GAMMA_G1)-1);
            Z_L = 1 / (L_rho*cs_L) * pow(Pg/L_p, -GAMMA_G2);
        }
        if(Pg>R_p)
        {
            /* shock wave */
            a0 = GAMMA_G5 / R_rho;
            a1 = GAMMA_G6 * R_p;
            a2 = sqrt(a0 / (Pg+a1));
            W_R = (Pg-R_p) * a2;
            Z_R = a2 * (1.0 - 0.5*(Pg-R_p)/(a1+Pg));
        } else {
            /* rarefaction wave */
            pratio = Pg / R_p;
            W_R = GAMMA_G4 * cs_R * (pow(pratio, GAMMA_G1)-1);
            Z_R = 1 / (R_rho*cs_R) * pow(pratio, -GAMMA_G2);
        }
        if(niter_Riemann < NMAX_ITER / 2)
            Pg -= (W_L + W_R + dvel) / (Z_L + Z_R);
        else
            Pg -= 0.5 * (W_L + W_R + dvel) / (Z_L + Z_R);
        
        if(Pg < 0.1 * Pg_prev)
            Pg = 0.1 * Pg_prev;
        
        tol = 2.0 * fabs((Pg-Pg_prev)/(Pg+Pg_prev));
//        printf("tol %e \n", tol);
        niter_Riemann++;
    }
    if(niter_Riemann<NMAX_ITER)
    {
        *P_M = Pg;
        *S_M = 0.5*(v_line_L+v_line_R) + 0.5*(W_R-W_L);
//        printf("Convergence! P_M %e \t S_M %e \n", Riemann_out->P_M, Riemann_out->S_M);
        return niter_Riemann;
    } else {
        return 0;
    }
}




inline void sample_reimann_vaccum_internal(PKD pkd, double S, 
      double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
      double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f,
      double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    double S_L = v_line_L + GAMMA_G4 * cs_L;
    double S_R = v_line_R - GAMMA_G4 * cs_R;
    if(S <= S_L)
    {
        /* left fan */
       sample_reimann_vaccum_right(pkd,  S,
                R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                P_M,  S_M,  rho_f,  p_f,  v_f,
                n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
    }
    else if(S >= S_R)
    {
        /* right fan */
       sample_reimann_vaccum_left(pkd,  S,
                R_rho, R_p,  R_v,  L_rho, L_p,  L_v,
                P_M,  S_M,  rho_f,  p_f,  v_f,
                n_unit,  v_line_L,  v_line_R,  cs_L,  cs_R);
    }
    else
    {
        /* vacuum in between */
        *P_M = 0;
        *S_M = S;
#ifndef USE_MFM
        *rho_f = 0;
        *p_f = *P_M;
        int k;
        for(k=0;k<3;k++)
            v_f[k] = (L_v[k] + (R_v[k]-L_v[k]) * (S-S_L)/(S_R-S_L)) *
            (1-n_unit[k]) + S * n_unit[k];
#endif
    }
}

inline double guess_for_pressure(PKD pkd, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    double pmin, pmax;
    /* start with the usual lowest-order guess for the contact wave pressure */
    double pv = 0.5*(L_p+R_p) - 0.125*(v_line_R-v_line_L)*(L_p+R_p)*(cs_L+cs_R);
    pmin = DMIN(L_p,R_p);
    pmax = DMAX(L_p,R_p);
    
    /* if one side is vacuum, guess half the mean */
    if(pmin<=0)
        return 0.5*(pmin+pmax);

    /* if the two are sufficiently close, and pv is between both values, return it */
    double qrat = pmax / pmin;
    if(qrat <= 2.0 && (pmin <= pv && pv <= pmax))
        return pv;
    
    if(pv < pmin)
    {
        /* use two-rarefaction solution */
        double pnu = (cs_L+cs_R) - GAMMA_G7 * (v_line_R - v_line_L);
        double pde = cs_L / pow(L_p, GAMMA_G1) + cs_R / pow(R_p, GAMMA_G1);
        return pow(pnu / pde, GAMMA_G3);
    }
    else
    {
        /* two-shock approximation  */
        double gel = sqrt((GAMMA_G5 / L_rho) / (GAMMA_G6 * L_p + pv));
        double ger = sqrt((GAMMA_G5 / R_rho) / (GAMMA_G6 * R_p + pv));
        double x = (gel * L_p + ger * R_p - (v_line_R - v_line_L)) / (gel + ger);
        if(x < pmin || x > pmax)
            x = pmin;
        return x;
    }
}


/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/*  This is the "normal" Riemann fan, with no vacuum on L or R state! */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
inline void sample_reimann_standard(PKD pkd, double S, 
      double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
      double P_M, double S_M, double *rho_f, double *p_f, double *v_f, 
      double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
#ifndef HYDRO_MESHLESS_FINITE_VOLUME
    /* we don't actually need to evaluate the fluxes, and we already have P_M and S_M, which define the 
     contact discontinuity where the rho flux = 0; so can simply exit this routine */
    //return;
#endif
    int k; double C_eff,S_eff;
    if(S <= S_M)  /* sample point is left of contact discontinuity */
    {
        if(P_M <= L_p)	/* left fan (rarefaction) */
        {
            double S_check_L = v_line_L - cs_L;
            if(S <= S_check_L) /* left data state */
            {
                *p_f = L_p;
                *rho_f = L_rho;
                for(k=0;k<3;k++)
                    v_f[k] = L_v[k];
                return;
            }
            else
            {
                double C_eff_L = cs_L * pow(P_M / L_p, GAMMA_G1);
                double S_tmp_L = S_M - C_eff_L;
                
                if(S > S_tmp_L)	/* middle left state */
                {
                    *rho_f = L_rho * pow(P_M / L_p, GAMMA_G8);
                    *p_f = P_M;
                    for(k=0;k<3;k++)
                        v_f[k] = L_v[k] + (S_M-v_line_L)*n_unit[k];
                    return;
                }
                else		/* left state inside fan */
                {
                    S_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
                    C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
                    *rho_f = L_rho * pow(C_eff / cs_L, GAMMA_G4);
                    *p_f = L_p * pow(C_eff / cs_L, GAMMA_G3);
                    for(k=0;k<3;k++)
                        v_f[k] = L_v[k] + (S_eff-v_line_L)*n_unit[k];
                    return;
                }
            }
        }
        else			/* left shock */
        {
            if(L_p > 0)
            {
                double pml = P_M / L_p;
                double S_L = v_line_L - cs_L * sqrt(GAMMA_G2 * pml + GAMMA_G1);
                
                if(S <= S_L)	/* left data state */
                {
                    *p_f = L_p;
                    *rho_f = L_rho;
                    for(k=0;k<3;k++)
                        v_f[k] = L_v[k];
                    return;
                }
                else		/* middle left state behind shock */
                {
                    *rho_f = L_rho * (pml + GAMMA_G6) / (pml * GAMMA_G6 + 1.0);
                    *p_f = P_M;
                    for(k=0;k<3;k++)
                        v_f[k] = L_v[k] + (S_M-v_line_L)*n_unit[k];
                    return;
                }
            }
            else
            {
                *rho_f = L_rho / GAMMA_G6;
                *p_f = P_M;
                for(k=0;k<3;k++)
                    v_f[k] = L_v[k] + (S_M-v_line_L)*n_unit[k];
                return;
            }
        }
    }
    else    /* sample point is right of contact discontinuity */
    {
        if(P_M > R_p)	/* right shock */
        {
            if(R_p > 0)
            {
                double pmr = P_M / R_p;
                double S_R = v_line_R + cs_R * sqrt(GAMMA_G2 * pmr + GAMMA_G1);
                
                if(S >= S_R)	/* right data state */
                {
                    *p_f = R_p;
                    *rho_f = R_rho;
                    for(k=0;k<3;k++)
                        v_f[k] = R_v[k];
                    return;
                }
                else		/* middle right state behind shock */
                {
                    *rho_f = R_rho * (pmr + GAMMA_G6) / (pmr * GAMMA_G6 + 1.0);
                    *p_f = P_M;
                    for(k=0;k<3;k++)
                        v_f[k] = R_v[k] + (S_M-v_line_R)*n_unit[k];
                    return;
                }
            }
            else
            {
                *rho_f = R_rho / GAMMA_G6;
                *p_f = P_M;
                for(k=0;k<3;k++)
                    v_f[k] = R_v[k] + (S_M-v_line_R)*n_unit[k];
                return;
            }
        }
        else			/* right fan */
        {
            double S_check_R = v_line_R + cs_R;
            if(S >= S_check_R)		/* right data state */
            {
                *p_f = R_p;
                *rho_f = R_rho;
                for(k=0;k<3;k++)
                    v_f[k] = R_v[k];
                return;
            }
            else
            {
                double C_eff_R = cs_R * pow(P_M / R_p, GAMMA_G1);
                double S_tmp_R = S_M + C_eff_R;

                if(S <= S_tmp_R)	/* middle right state */
                {
                    *rho_f = R_rho * pow(P_M / R_p, GAMMA_G8);
                    *p_f = P_M;
                    for(k=0;k<3;k++)
                        v_f[k] = R_v[k] + (S_M-v_line_R)*n_unit[k];
                    return;
                }
                else		/* fan right state */
                {
                    S_eff = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
                    C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
                    *rho_f = R_rho * pow(C_eff / cs_R, GAMMA_G4);
                    *p_f = R_p * pow(C_eff / cs_R, GAMMA_G3);
                    for(k=0;k<3;k++)
                        v_f[k] = R_v[k] + (S_eff-v_line_R)*n_unit[k];
                    return;
                }
            }
        }
    }
}

/* -------------------------------------------------------------------------------------------------------------- */
/*  Part of exact Riemann solver: */
 /*    take the face state we have calculated from the exact Riemann solution and get the corresponding fluxes */
/*   (written by V. Springel for AREPO, with minor modifications) */
 /* -------------------------------------------------------------------------------------------------------------- */
void convert_face_to_flux(PKD pkd, 
      double *rho_f, double *p_f, double *v_f, 
      double n_unit[3])
{
    double rho, P, v[3], v_line=0, v_frame=0, h=0; int k;
    rho = *rho_f;
    P = *p_f;
    for(k=0;k<3;k++)
    {
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
    for(k=0;k<3;k++)
        v_f[k] = (*rho_f) * v[k] + P * n_unit[k];
    return;
}




/* --------------------------------------------------------------------------------- */
/* Part of exact Riemann solver: */
/* right state is a vacuum, but left state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_vaccum_right(PKD pkd, double S, 
                                double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
                                double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f, 
                                double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    //double S_L = v_line_L - GAMMA_G4 * cs_L;
    double S_L = v_line_L + GAMMA_G4 * cs_L; // above line was a sign error, caught by Bert Vandenbroucke
    //printf("Vaccuum right\n");
#ifdef USE_MFM
    /* in this code mode, we are -always- moving with the contact discontinuity so density flux = 0;
     this constrains where we reside in the solution fan */
    *P_M = 0;
    *S_M = S_L;
    return;
#endif
    
    if(S_L < S)
    {
        /* vacuum */
        *P_M = 0;
        *S_M = S_L;
        *rho_f = 0;
    } else {
        /* left fan */
        double S_L_check = v_line_L - cs_L;
        if(S_L_check < S)
        {
            /* rarefaction fan left state */
            double C_eff = GAMMA_G5 * (cs_L + GAMMA_G7 * (v_line_L - S));
            *P_M = L_p * pow(C_eff / cs_L, GAMMA_G3);
            *S_M = GAMMA_G5 * (cs_L + GAMMA_G7 * v_line_L + S);
            *rho_f = L_rho * pow(C_eff / cs_L, GAMMA_G4);
        } else {
            /* left data state */
            *P_M = L_p;
            *S_M = v_line_L;
            *rho_f = L_rho;
        }
    }
    *p_f = *P_M;
    int k;
    for(k=0;k<3;k++)
        v_f[k] = L_v[k] + (*S_M - v_line_L) * n_unit[k];
    return;
}



/* --------------------------------------------------------------------------------- */
/* part of exact Riemann solver: */
/* left state is a vacuum, but right state is not: sample the fan appropriately */
/*  (written by V. Springel for AREPO) */
/* --------------------------------------------------------------------------------- */
void sample_reimann_vaccum_left(PKD pkd, double S, 
                                double R_rho,double R_p, double R_v[3],double L_rho,double L_p, double L_v[3],
                                double *P_M, double *S_M, double *rho_f, double *p_f, double *v_f, 
                                double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    //printf("Vaccuum left\n");
    double S_R = v_line_R - GAMMA_G4 * cs_R;
#ifdef USE_MFM
    /* in this code mode, we are -always- moving with the contact discontinuity so density flux = 0; 
     this constrains where we reside in the solution fan */
    *P_M = 0;
    *S_M = S_R;
    return;
#endif

    if(S_R > S)
    {
        /* vacuum */
        *P_M = 0;
        *S_M = S_R;
        *rho_f = 0;
    } else {
        /* right fan */
        double S_R_check = v_line_R + cs_R;
        if(S_R_check > S)
        {
            /* rarefaction fan right state */
            double C_eff = GAMMA_G5 * (cs_R - GAMMA_G7 * (v_line_R - S));
            *P_M = R_p * pow(C_eff / cs_R, GAMMA_G3);
            *S_M = GAMMA_G5 * (-cs_R + GAMMA_G7 * v_line_R + S);
            *rho_f = R_rho * pow(C_eff / cs_R, GAMMA_G4);
        } else {
            /* right data state */
            *P_M = R_p;
            *S_M = v_line_R;
            *rho_f = R_rho;
        }
    }
    *p_f = *P_M;
    int k;
    for(k=0;k<3;k++)
        v_f[k] = R_v[k] + (*S_M - v_line_R) * n_unit[k];
    return;
}
