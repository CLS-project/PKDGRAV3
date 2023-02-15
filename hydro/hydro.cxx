/* File added by Isaac Alonso
 * the hydrodynamical part using a mesh-free
 * Hydrodynamical solver using the mesh free methods of
 * Lanson&Vila 2008, Gaburov&Nitadori 2011 and Hopkins 2015
 */
#include <algorithm>
#include "hydro.h"


#include <stdio.h>

using blitz::TinyVector;
using blitz::dot;

/* ----------------
 * HELPER FUNCTIONS
 * ----------------
 */
template<typename T>
int sign(T v) {return (v>0) - (v<0);}


/* We use a cubic spline kernel.
 * See, for example:
 * https://pysph.readthedocs.io/en/latest/reference/kernels.html
 */
double cubicSplineKernel(double r, double h) {
    double q;
    q = r/h;
    if (q<1.0) {
        return M_1_PI/(h*h*h)*( 1. - 1.5*q*q*(1.-0.5*q) );
    }
    else if (q<2.0) {
        return 0.25*M_1_PI/(h*h*h)*(2.-q)*(2.-q)*(2.-q);
    }
    else {
        return 0.0;
    }
}


void inverseMatrix(double *E, double *B) {
    double det;

    det = E[XX]*E[YY]*E[ZZ] + 2.0*E[XY]*E[YZ]*E[XZ] //+ E[XZ]*E[XY]*E[YZ]
          -E[XX]*E[YZ]*E[YZ] - E[ZZ]*E[XY]*E[XY] - E[YY]*E[XZ]*E[XZ];

    if (det==0) {
        printf("Singular matrix!\n");
        printf("XX %e \nXY %e \t YY %e \nXZ %e \t YZ %e \t ZZ %e \n",
               E[XX], E[XY], E[YY], E[XZ], E[YZ], E[ZZ]);
        abort();
        B[XX] = 0.;
        B[YY] = 0.;
        B[ZZ] = 0.;
        B[XY] = 0.;
        B[XZ] = 0.;
        B[YZ] = 0.;
    }

    det = 1./det;

    B[XX] = (E[YY]*E[ZZ] - E[YZ]*E[YZ])*det;
    B[YY] = (E[XX]*E[ZZ] - E[XZ]*E[XZ])*det;
    B[ZZ] = (E[YY]*E[XX] - E[XY]*E[XY])*det;

    B[XY] = -(E[XY]*E[ZZ] - E[YZ]*E[XZ])*det;
    B[XZ] = (E[XY]*E[YZ] - E[YY]*E[XZ])*det;
    B[YZ] = -(E[XX]*E[YZ] - E[XY]*E[XZ])*det;


}

double conditionNumber(double *E, double *B) {
    double modE = 0.;
    modE += E[XX]*E[XX];
    modE += E[YY]*E[YY];
    modE += E[ZZ]*E[ZZ];
    modE += 2.*E[XY]*E[XY];
    modE += 2.*E[XZ]*E[XZ];
    modE += 2.*E[YZ]*E[YZ];

    double modB = 0.;
    modB += B[XX]*B[XX];
    modB += B[YY]*B[YY];
    modB += B[ZZ]*B[ZZ];
    modB += 2.*B[XY]*B[XY];
    modB += 2.*B[XZ]*B[XZ];
    modB += 2.*B[YZ]*B[YZ];

    return sqrt(modB*modE)/3.;
}


void BarthJespersenLimiter(double *limVar, TinyVector<double,3> gradVar,
                           double var_max, double var_min,
                           TinyVector<double,3> dr) {
#ifdef DEBUG_FLUX_NOLIMITER
    *limVar = 1;
    return;
#endif
    double diff, lim;

    diff = dot(gradVar,dr);
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
void ConditionedBarthJespersenLimiter(double *limVar, TinyVector<myreal,3> gradVar,
                                      double var_max, double var_min,
                                      TinyVector<double,3> dr,
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


    diff = dot(gradVar,dr);
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

// Equation 10.39
inline void compute_Ustar(double rho_K, double S_K, double v_K,
                          double p_K, double h_K, double S_s,
                          double *rho_sK, double *rhov_sK, double *e_sK) {
    double fac = rho_K * (S_K - v_K)/(S_K-S_s);

    *rho_sK = fac;

    *rhov_sK = S_s * fac;

    double e_K = rho_K*h_K - p_K;
    *e_sK = fac * ( e_K/rho_K + (S_s - v_K)*(S_s + p_K/(rho_K*(S_K - v_K))) );
}






#define psi1 0.5
#define psi2 0.25
void genericPairwiseLimiter(double Lstate, double Rstate,
                            double *Lstate_face, double *Rstate_face) {
#ifdef DEBUG_FLUX_NOLIMITER
    return;
#endif
    double phi_max, phi_min, d1, d2, phi_mean, phi_p, phi_m;

    if (Lstate == Rstate) {
        *Lstate_face = Lstate;
        *Rstate_face = Rstate;
    }
    else {

        d1 = psi1*fabs(Lstate - Rstate);
        d2 = psi2*fabs(Lstate - Rstate);

        phi_mean = 0.5*(Lstate+Rstate);

        phi_min = std::min(Lstate, Rstate);
        phi_max = std::max(Lstate, Rstate);

        if (sign(phi_min - d1) == sign(phi_min) ) {
            phi_m = phi_min - d1;
        }
        else {
            phi_m = phi_min/(1. + d1/fabs(phi_min));
        }

        if (sign(phi_max + d1) == sign(phi_max) ) {
            phi_p = phi_max + d1;
        }
        else {
            phi_p = phi_max/(1. + d1/fabs(phi_max));
        }

        if (Lstate < Rstate) {
            *Lstate_face = std::max(phi_m, std::min(phi_mean+d2, *Lstate_face));
            *Rstate_face = std::min(phi_p, std::max(phi_mean-d2, *Rstate_face));
        }
        else {
            *Rstate_face = std::max(phi_m, std::min(phi_mean+d2, *Rstate_face));
            *Lstate_face = std::min(phi_p, std::max(phi_mean-d2, *Lstate_face));
        }

    }


}


void hydroSourceGravity(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                        double pDelta, double *pa, double dScaleFactor,
                        int bComove) {
    double gravE = 0.0;
    double gravE_dmdt = 0.0;
    double aFac_m2 = 1./(dScaleFactor*dScaleFactor);
    const auto &a = p.acceleration();

    if (bComove) {
        for (int j=0; j<3; j++) {
            // IA: One 1/a is from the definition of acceleration in pkdgrav3.
            //  The other comes from the shape of the source term, which is proportional to  1/a
            pa[j] = a[j]*aFac_m2; // TODO: Do 1/a2 only once
        }
    }
    else {
        for (int j=0; j<3; j++) {
            pa[j] = a[j];
        }
    }
    auto &v = p.velocity();
    for (int j=0; j<3; j++) {
        psph->mom[j] += 0.5*pDelta*(psph->lastMass*psph->lastAcc[j] + p.mass()*pa[j]);
        v[j] = psph->mom[j]/p.mass();
#ifndef USE_MFM
        // IA: Multiplying here by 'a' instead of doing it at hydro.c is simpler and more efficient.
        // However, it may hinder conservation properties. But doing the time average over two steps is not conservative anyway
        // In the Zeldovich case I have not found any relevant difference among both options
        gravE_dmdt +=  0.5*( psph->lastAcc[j]*psph->lastDrDotFrho[j] +  pa[j]*psph->drDotFrho[j]*dScaleFactor ) ;
#endif
        gravE += 0.5*pDelta*( psph->lastMom[j]*psph->lastAcc[j] + p.mass()*v[j]*pa[j]  );
    }
    if (pDelta==0.) gravE_dmdt = 0.;

    psph->E += gravE - 0.5*gravE_dmdt;
}




void hydroSourceExpansion(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                          double pDelta, double dScaleFactor, double dHubble,
                          int bComove, double dConstGamma) {
    //  E^{n+1} = E^{n} + dE_flux - dt*(H^{n} E^n + H^{n+1} E^{n+1})
    if (bComove) {
        psph->E = (psph->E - psph->lastHubble*pDelta*psph->lastE) /
                  (1.+pDelta*dHubble);
        psph->Uint = (psph->Uint -
                      psph->lastHubble*1.5*pDelta*psph->lastUint *
                      (dConstGamma - 1.)) /
                     (1.+1.5*pDelta*dHubble*(dConstGamma - 1.));
#ifdef ENTROPY_SWITCH
        psph->S = (psph->S -
                   psph->lastHubble*1.5*pDelta*psph->lastS *
                   (dConstGamma - 1.)) /
                  (1.+1.5*pDelta*dHubble*(dConstGamma - 1.));
#endif

        for (int j=0; j<3; j++) {
            psph->mom[j] = (psph->mom[j] -
                            0.5*psph->lastHubble*pDelta*psph->lastMom[j]) /
                           (1.+0.5*pDelta*dHubble);
        }
        psph->lastHubble = dHubble;
    }

}





void hydroSyncEnergies(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph, double pa[3], double dConstGamma) {
    double Ekin = 0.5*(psph->mom[0]*psph->mom[0] +
                       psph->mom[1]*psph->mom[1] +
                       psph->mom[2]*psph->mom[2])/p.mass();
    double Egrav = p.mass()*
                   sqrt(pa[0]*pa[0] + pa[1]*pa[1] + pa[2]*pa[2])*p.ball();

    if ( (Ekin+Egrav) > 100.*psph->Uint ) {
        // IA: The fluid is dominated by the kinetic energy,
        // so using Etot may be unreliable to compute the pressure
#ifdef ENTROPY_SWITCH
        if ((psph->S>0.) &&
                (psph->Uint<0.001*(psph->maxEkin+psph->Uint) || psph->Uint<0.001*Egrav )) {
            // The flow is smooth and/or dominated by gravity,
            // thus the entropy can be used to evolve the pressure
            psph->Uint = psph->S *
                         pow(pkdDensity(pkd,p), dConstGamma-1) /
                         (dConstGamma -1.);
        }
        else if (psph->Uint > 0.) {
            // The flow is NOT smooth, so the entropy can not be used
            psph->S = psph->Uint *
                      (dConstGamma -1.) *
                      pow(pkdDensity(pkd,p), -dConstGamma+1);
        }
        else {
            printf("WARNING %e \t S %e \t(%e) \t Uint %e ≤t(%e) \n",psph->P,
                   psph->S,
                   psph->S / pkdMass(pkd,p) * pow(pkdDensity(pkd,p), dConstGamma),
                   psph->Uint,
                   psph->Uint*psph->omega*(dConstGamma -1.) );
            psph->S=0.0;
            psph->Uint=0.0;
        }
#endif // ENTROPY_SWITCH
        psph->E = psph->Uint + Ekin;
    }
    else {
        // The fluid is dominated by pressure forces so the total energy
        // should be used
        psph->Uint = psph->E - Ekin;
#ifdef ENTROPY_SWITCH
        psph->S = psph->Uint *(dConstGamma -1.) *
                  pow(pkdDensity(pkd,p), -dConstGamma+1);
#endif
    }
}




void hydroSetPrimitives(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                        double dTuFac, double dConstGamma) {
    // Temperature minimum of T=0, but could be changed.
    // If cooling is used, the corresponding entropy floor
    // is applied inside cooling_cool_part
    double minUint = 0. * dTuFac * p.mass();
    if (psph->Uint < minUint) {
        double Ekin = 0.5*(psph->mom[0]*psph->mom[0] +
                           psph->mom[1]*psph->mom[1] +
                           psph->mom[2]*psph->mom[2])/p.mass();
        psph->Uint = minUint;
        psph->E = Ekin + minUint;
#ifdef ENTROPY_SWITCH
        psph->S = psph->Uint *(dConstGamma -1.) *
                  pow(p.density(), -dConstGamma+1);
#endif
    }
    psph->P = psph->Uint*psph->omega*(dConstGamma -1.);


    psph->c = sqrt(psph->P*dConstGamma/p.density());

    auto &v = p.velocity();
    for (int j=0; j<3; j++) {
        v[j] = psph->mom[j]/p.mass();

        /*IA: This is here for compatibility with hydro.c, as in there we use vPred.
         * I think that all could be changed to use only pkdVel instead.
         * But I am not sure if when adding gravity vPred would be needed,
         * thus I will *temporarly* keep it
         */
        psph->vPred[j] = v[j];
    }
}





void hydroSetLastVars(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph, double *pa,
                      double dScaleFactor, double dTime, double dDelta,
                      double dConstGamma) {
#ifndef USE_MFM
    for (int j=0; j<3; j++) {
        psph->lastDrDotFrho[j] = psph->drDotFrho[j]*dScaleFactor;
        psph->drDotFrho[j] = 0.;
    }
#endif
    for (int j=0; j<3; j++) {
        psph->lastAcc[j] = pa[j];
        psph->lastMom[j] = psph->mom[j];
    }
    psph->lastUpdateTime = dTime;
    psph->lastE = psph->E;
    psph->lastUint = psph->Uint;
    psph->lastMass = p.mass();
#ifdef ENTROPY_SWITCH
    if (dDelta <= 0) {
        // Initialize the entropy
        psph->S = p.mass() * psph->P *
                  pow(p.density(), -dConstGamma);
    }
    psph->lastS = psph->S;
#endif
}


