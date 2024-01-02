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

/* We use a cubic spline kernel.
 * See, for example:
 * https://pysph.readthedocs.io/en/latest/reference/kernels.html
 */

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

void hydroSourceGravity(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                        double pDelta, TinyVector<double,3> &pa, double dScaleFactor,
                        int bComove) {
    double aFac_m2 = 1./(dScaleFactor*dScaleFactor);
    const auto &a = p.acceleration();

    if (bComove) {
        // IA: One 1/a is from the definition of acceleration in pkdgrav3.
        //  The other comes from the shape of the source term, which is proportional to  1/a
        pa = a * aFac_m2; // TODO: Do 1/a2 only once
    }
    else {
        pa = a;
    }

    psph->mom += 0.5 * pDelta * (psph->lastMass*psph->lastAcc + p.mass()*pa);
    p.velocity() = psph->mom / p.mass();

    double gravE_dmdt = 0.;
#ifndef USE_MFM
    // IA: Multiplying here by 'a' instead of doing it at hydro.c is simpler and more efficient.
    // However, it may hinder conservation properties. But doing the time average over two steps is not conservative anyway
    // In the Zeldovich case I have not found any relevant difference among both options
    if (pDelta > 0.)
        gravE_dmdt = 0.5 * (dot(psph->lastAcc,psph->lastDrDotFrho) + dot(pa,psph->drDotFrho));
#endif

    const double gravE = 0.5 * pDelta * (dot(psph->lastMom,psph->lastAcc) +
                                         dot(psph->mom,pa));

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

        psph->mom = (psph->mom - 0.5*psph->lastHubble*pDelta*psph->lastMom) /
                    (1. + 0.5*pDelta*dHubble);
        psph->lastHubble = dHubble;
    }

}

void hydroSyncEnergies(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                       const TinyVector<double,3> &pa, double dConstGamma) {
    double Ekin = 0.5 * dot(psph->mom,psph->mom) / p.mass();
    double Egrav = p.mass() * sqrt(dot(pa,pa)) * p.ball();

    if ( (Ekin+Egrav) > 100.*psph->Uint ) {
        // IA: The fluid is dominated by the kinetic energy,
        // so using Etot may be unreliable to compute the pressure
#ifdef ENTROPY_SWITCH
        if ((psph->S>0.) &&
                (psph->Uint<0.001*(psph->maxEkin+psph->Uint) || psph->Uint<0.001*Egrav )) {
            // The flow is smooth and/or dominated by gravity,
            // thus the entropy can be used to evolve the pressure
            psph->Uint = psph->S *
                         pow(p.density(), dConstGamma-1) /
                         (dConstGamma -1.);
        }
        else if (psph->Uint > 0.) {
            // The flow is NOT smooth, so the entropy can not be used
            psph->S = psph->Uint *
                      (dConstGamma -1.) *
                      pow(p.density(), -dConstGamma+1);
        }
        else {
            printf("WARNING %e \t S %e \t(%e) \t Uint %e ≤t(%e) \n",psph->P,
                   psph->S,
                   psph->S / p.mass() * pow(p.density(), dConstGamma),
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
                  pow(p.density(), -dConstGamma+1);
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
        double Ekin = 0.5 * dot(psph->mom,psph->mom) / p.mass();
        psph->Uint = minUint;
        psph->E = Ekin + minUint;
#ifdef ENTROPY_SWITCH
        psph->S = psph->Uint *(dConstGamma -1.) *
                  pow(p.density(), -dConstGamma+1);
#endif
    }

    psph->P = psph->Uint*psph->omega*(dConstGamma -1.);
    psph->c = sqrt(psph->P*dConstGamma/p.density());
    p.velocity() = psph->mom / p.mass();
}

void hydroSetLastVars(PKD pkd, particleStore::ParticleReference &p, SPHFIELDS *psph,
                      const TinyVector<double,3> &pa, double dScaleFactor, double dTime,
                      double dDelta, double dConstGamma) {
#ifndef USE_MFM
    psph->lastDrDotFrho = psph->drDotFrho;
    psph->drDotFrho = 0.;
#endif
    psph->lastAcc = pa;
    psph->lastMom = psph->mom;
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

void hydroResetFluxes(SPHFIELDS *psph) {
    psph->Frho = 0.0;
    psph->Fene = 0.0;
    psph->Fmom = 0.0;
}
