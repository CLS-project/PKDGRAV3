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

inline double DMAX(double a, double b) { return (a > b) ? a : b; }
inline double DMIN(double a, double b) { return (a < b) ? a : b; }
inline int Riemann_solver_exact(PKD pkd, 
      double R_rho, double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R);

inline int iterative_Riemann_solver(PKD pkd, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double v_line_L, double v_line_R, double cs_L, double cs_R);


inline double guess_for_pressure(PKD pkd, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double v_line_L, double v_line_R, double cs_L, double cs_R);

inline void sample_reimann_vaccum_internal(PKD pkd, double S, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
      double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R);



inline int Riemann_solver_exact(PKD pkd, 
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
            double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R, double h_L, double h_R)
{
    //printf("Going for the exact Riemann Solver \n");
    /* first, we need to check for all the special/exceptional cases that will cause things to go haywire */
   int niter = 0;
    if((L_p == 0 && L_p == 0) || (L_rho==0 && R_rho==0))
    {
        /* we're in a Vaccuum! */
        *P_M = *S_M = 0;
//#ifdef HYDRO_MESHLESS_FINITE_VOLUME
//        Riemann_out->Fluxes.rho = Riemann_out->Fluxes.p = Riemann_out->Fluxes.v[0] = Riemann_out->Fluxes.v[1] = Riemann_out->Fluxes.v[2] = 0;
//#endif
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

            // IA: For MFM, this do nothing
            //sample_reimann_standard(pkd, 0.0,R_rho, R_p, L_rho, L_p , *P_M, *S_M,n_unit,v_line_L,v_line_R,cs_L,cs_R);
            //
//           printf("Normal Riemann solution %e \n", Riemann_out->Fluxes.rho);
        }
        else
        {
            /* ICs lead to vacuum, need to sample vacuum solution */
            sample_reimann_vaccum_internal(pkd, 0.0,R_rho, R_p, L_rho, L_p,P_M, S_M,n_unit,v_line_L,v_line_R,cs_L,cs_R);
        }
    } else {
        /* one of the densities is zero or negative */
        if((L_rho<0)||(R_rho<0))
            assert(0);
        if(L_rho>0){
             double S_L = v_line_L + GAMMA_G4 * cs_L; // above line was a sign error, caught by Bert Vandenbroucke
             *P_M = 0;
             *S_M = S_L;
        }
        if(R_rho>0){
             double S_R = v_line_R - GAMMA_G4 * cs_R;
             *P_M = 0;
             *S_M = S_R;
        }
    }
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
      double R_rho,double R_p,double L_rho,double L_p,
      double *P_M, double *S_M, 
                                    double n_unit[3], double v_line_L, double v_line_R, double cs_L, double cs_R)
{
    double S_L = v_line_L + GAMMA_G4 * cs_L;
    double S_R = v_line_R - GAMMA_G4 * cs_R;
    if(S <= S_L)
    {
       *P_M = 0;
       *S_M = S_L;
    }
    else if(S >= S_R)
    {
       *P_M = 0;
       *S_M = S_R;
    }
    else
    {
        /* vacuum in between */
        *P_M = 0;
        *S_M = S;
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


