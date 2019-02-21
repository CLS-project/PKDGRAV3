/* File added by Isaac Alonso for computing
 * the hydrodynamical part using a mesh-free
 * method, following the work of Hopkins 2015
 */

#include "pkd.h"
#include "smoothfcn.h"
#include "hydro.h" 
#include "riemann.h"
#include <stdio.h>


double cubicSplineKernel(double r, double h) {
   double q;
   q = r/h;
   if (q<1.0){
      return M_1_PI/(h*h*h)*( 1 - 1.5*q*q*(1.-0.5*q) );  
   }else if (q<2.0){
      return 0.25*M_1_PI/(h*h*h)*(2.-q)*(2.-q)*(2.-q);
   }else{
      return 0.0;
   }
}


#define XX 0
#define YY 3
#define ZZ 5
#define XY 1
#define XZ 2
#define YZ 4
void inverseMatrix(double* E, double* B){
   double det;

   det = E[XX]*E[YY]*E[ZZ] + 2.0*E[XY]*E[YZ]*E[XZ] //+ E[XZ]*E[XY]*E[YZ] 
        -E[XX]*E[YZ]*E[YZ] - E[ZZ]*E[XY]*E[XY] - E[YY]*E[XZ]*E[XZ];

   B[XX] = (E[YY]*E[ZZ] - 2.0*E[YZ])/det;
   B[YY] = (E[ZZ]*E[ZZ] - 2.0*E[XZ])/det;
   B[ZZ] = (E[YY]*E[XX] - 2.0*E[XY])/det;

   B[XY] = (E[XY]*E[ZZ] - E[YZ]*E[XZ])/det;
   B[XZ] = (E[XY]*E[YZ] - E[YY]*E[XZ])/det;
   B[YZ] = (E[XX]*E[YZ] - E[XY]*E[XZ])/det;

}



void initHydroForces(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    if (pkdIsActive(pkd,p)) {
	SPHFIELDS *psph = pkdSph(pkd,p);
	psph->uDot = 0;
	psph->fMetalsDot = 0;
	if (!pkd->param.bDoGravity) { /* Normally these are zeroed in Gravity */
	    p->uNewRung = 0;
	    pkdAccel(pkd,p)[0] = 0;
	    pkdAccel(pkd,p)[1] = 0;
	    pkdAccel(pkd,p)[2] = 0;
	    }
	}
    }

void initHydroForcesCached(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    if (pkdIsActive(pkd,p)) {
	SPHFIELDS *psph = pkdSph(pkd,p);
	psph->uDot = 0;
	psph->fMetalsDot = 0;
	p->uNewRung = 0;
	pkdAccel(pkd,p)[0] = 0;  
	pkdAccel(pkd,p)[1] = 0;  
	pkdAccel(pkd,p)[2] = 0; /* JW: Cached copies have zero accel! && rung */
	}
    }


void combHydroForces(void *vpkd, void *p1,void *p2) {
    PKD pkd = (PKD) vpkd;
    assert(!pkd->bNoParticleOrder);
    if (pkdIsActive(pkd,p1)) {
	SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);
	float *a1 = pkdAccel(pkd,p1), *a2 = pkdAccel(pkd,p2);
	psph1->uDot += psph2->uDot;
	psph1->fMetalsDot += psph2->fMetalsDot;
	a1[0] += a2[0];  
	a1[1] += a2[1];  
	a1[2] += a2[2]; 
	if (((PARTICLE *) p2)->uNewRung > ((PARTICLE *) p1)->uNewRung) 
	    ((PARTICLE *) p1)->uNewRung = ((PARTICLE *) p2)->uNewRung;
	}
    }


void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {

    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double E[6], B[6];  // IA: We are assumming 3D here!
    double ph, qh, rpq, hpq, dx,dy,dz,Wpq;
    double psi, psiTilde[3];
    double diff, pDensity;
    int i, j;

    /*
    assert(!pkd->bNoParticleOrder);

    aFac = (smf->a);        // comoving acceleration factor 
    vFac = (smf->bComove ? 1./(smf->a*smf->a) : 1.0); // converts v to xdot 
    dConstGammainv = 1/smf->dConstGamma;
    dtC = (1+0.6*smf->alpha)/(smf->a*smf->dEtaCourant);
    dtMu = (0.6*smf->beta)/(smf->a*smf->dEtaCourant);
    */

    /* Particle p data */
    psph = pkdSph(pkd,p);
    pDensity = pkdDensity(pkd,p); // TODO: check where is this computed
//    pMass = pkdMass(pkd,p); TODO: check where is computed.. is it useful for this implementation??? Probably not..
    ph = 0.5*fBall; /* IA: fBall seems to be the minimun distance to any neighbors ¿? Although it would be logical to be the maximum */
//    pActive = pkdIsActive(pkd,p);

    /* IA: I do not get this, so I will use my own cubic spline kernel */
//    ih2 = 1/(ph*ph);
//    fNorm = 0.5*M_1_PI*ih2/ph;
//    fNorm1 = fNorm*ih2;	/* converts to physical u */

    /* IA: Compute the \omega(x_i) normalization factor */
    psph->omega = 0.0;
    for (i=0; i<nSmooth; ++i){
       q = nnList[i].pPart;

       qh = 0.5*pkdBall(pkd,q); // TODO: Check this 

       rpq = sqrt(nnList[i].fDist2);
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

       psph->omega += cubicSplineKernel(rpq, hpq);
    }


    /* IA: We compute the density making use of Eq. 27 Hopkins 2015 */
    pkdSetDensity(pkd,p, pkdMass(pkd,p)*psph->omega);
    /* IA: TODO: having this information maybe some memory can be saved if we do not
     * store \omega */



    /* IA: Compute the E matrix (Hopkins 2015, eq 14) */
    for (i=0; i<6;++i){
       E[i] = 0.0; 
    }
    // TODO: Check initialization of E

    for (i=0; i<nSmooth;++i){
       q = nnList[i].pPart;

       dx = nnList[i].dx;
       dy = nnList[i].dy;
       dz = nnList[i].dz;

       qh = 0.5*pkdBall(pkd,q); // TODO: Check this 

       rpq = sqrt(nnList[i].fDist2);
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

       Wpq = cubicSplineKernel(rpq, hpq);

       E[0] += dx*dx*Wpq;
       E[3] += dy*dy*Wpq;
       E[5] += dz*dz*Wpq;

       E[1] += dy*dx*Wpq;
       E[2] += dz*dx*Wpq;
       E[4] += dy*dz*Wpq;

    }

    
    /* IA: Normalize the matrix */
    for (i=0; i<6;++i){
       E[i] /= psph->omega; 
    }

    /* IA: END of E matrix computation */

    /* IA: Now, we need to do the inverse */
    inverseMatrix(E, psph->B);
    
    /* IA: There is nothing more that we can do in this loop, as all the particle must have their densities
     * and B matrices computed */


    /* TODO: I should make sure that psph is not cleared after a smSmooth call!!! */

    }


/* IA: This routine will extrapolate the *already limited* gradient to the 'faces' and solve the 1D riemann problem. 
 * For now, the 1D riemann flux will be computed TWICE for a given face, one for each adjacent particles... This 
 * is highly sub-optimal, and I should think a way to 'mark' which riemann fluxes has been already computed.
 * IDEA: maybe a 64-bit with 0/1? But I think that I have no info about the neighbors of my neighbor..Or create a queue
 * of fluxes to be computed.. */
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double pv[3], qv[3], vFrame[3];
    double ph,qh, hpq, modApq, rpq, dx,dy,dz,Wpq,psi_p,psi_q,diff,pDensity;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3], dr[3]; 
    struct Input_vec_Riemann riemann_input;
    struct Riemann_outputs riemann_output; //TODO: do this outside the loop
    int i,j;

    psph = pkdSph(pkd, p);   

    for (i=0;i<nSmooth;++i){

	 q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);

	 dx = nnList[i].dx;
	 dy = nnList[i].dy;
	 dz = nnList[i].dz;

       qh = 0.5*pkdBall(pkd,q); // TODO: Check this 
       rpq = sqrt(nnList[i].fDist2);
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

       Wpq = cubicSplineKernel(rpq, hpq); 
       psi_p = Wpq/psph->omega;
       psi_q = Wpq/qsph->omega;

       /* Compute the gradients and extrapolate to the faces */
      psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi_p;
      psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi_p;
      psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi_p;

       psiTilde_q[0] = -(qsph->B[XX]*dx + qsph->B[XY]*dy + qsph->B[XZ]*dz)*psi_q;
       psiTilde_q[1] = -(qsph->B[XY]*dx + qsph->B[YY]*dy + qsph->B[YZ]*dz)*psi_q;
       psiTilde_q[2] = -(qsph->B[XZ]*dx + qsph->B[YZ]*dy + qsph->B[ZZ]*dz)*psi_q;


       // Velocity of the quadrature mid-point 
       for (j=0; j<3; j++){
          vFrame[j] = 0.5*(psph->vPred[j]+qsph->vPred[j]);
          // We boost to the reference of the p-q 'face'
          pv[j] = psph->vPred[j] - vFrame[j];
          qv[j] = qsph->vPred[j] - vFrame[j];
       }
 
       dr[0] = dx;
       dr[1] = dy;
       dr[2] = dz;
 
       // TODO: The next loops can be unified
       // TODO: Remove psph->*Grad variables
       // Density gradients
      diff = (pkdDensity(pkd,q) - pDensity);
      riemann_input.L.rho = pkdDensity(pkd,p);
      riemann_input.R.rho = pkdDensity(pkd,q);
      for (j=0; j<3; j++) {
         riemann_input.L.rho += 0.5*dr[j]*psiTilde_p[j]*diff;
         riemann_input.R.rho += 0.5*dr[j]*psiTilde_q[j]*diff; //Same sign (- from diff and another - from dr
      }

      // Momentum gradients 
      diff = qsph->vPred[0] - psph->vPred[0];// TODO :check this vPred!
      riemann_input.L.v[0] = pv[0];
      riemann_input.R.v[0] = qv[0];
      for (j=0; j<3; j++) {
//         psph->vxGrad[j] += psiTilde[j] * diff;
         riemann_input.L.v[0] += 0.5*dr[j]*psiTilde_p[j]*diff;
         riemann_input.R.v[0] += 0.5*dr[j]*psiTilde_q[j]*diff;
      }

      diff = qsph->vPred[1] - psph->vPred[1];// TODO :check this vPred!
      riemann_input.L.v[1] = pv[1];
      riemann_input.R.v[1] = qv[1];
      for (j=0; j<3; j++) {
//         psph->vyGrad[j] += psiTilde[j] * diff;
         riemann_input.L.v[1] += 0.5*dr[j]*psiTilde_p[j]*diff;
         riemann_input.R.v[1] += 0.5*dr[j]*psiTilde_q[j]*diff;
      }

      diff = qsph->vPred[2] - psph->vPred[2];// TODO :check this vPred!
      riemann_input.L.v[2] = pv[2];
      riemann_input.R.v[2] = qv[2];
      for (j=0; j<3; j++) {
//         psph->vzGrad[j] += psiTilde[j] * diff;
         riemann_input.L.v[2] += 0.5*dr[j]*psiTilde_p[j]*diff;
         riemann_input.R.v[2] += 0.5*dr[j]*psiTilde_q[j]*diff;
      }

      // Pressure gradients TODO: check when psph->cs is computed! is it predicted?
      riemann_input.L.p = psph->c*psph->c * pDensity / pkd->param.dConstGamma;
      riemann_input.R.p = qsph->c*qsph->c * pkdDensity(pkd,q) / pkd->param.dConstGamma;
      diff = riemann_input.R.p - riemann_input.L.p;
      for (j=0; j<3; j++){
//         psph->intEneGrad[j] += psiTilde[j] * diff;
         riemann_input.L.p += 0.5*dr[j]*psiTilde_p[j]*diff;
         riemann_input.R.p += 0.5*dr[j]*psiTilde_q[j]*diff;
      }

      // The cs and u values of the reconstructed states are computed inside the Riemann Solver!




       // Face where the riemann problem will be solved
       modApq = 0.0;
       for (j=0; j<3; j++){
          Apq[j] = psiTilde_q[j]/psph->omega - psiTilde_p[j]/qsph->omega;
          modApq += Apq[j]*Apq[j];
       }
       modApq = sqrt(modApq);

       for (j=0; j<3; j++){
          face_unit[j] = Apq[j]/modApq;
       }

       // It seems that in hydra_core_meshless.h they do not rotate any velocity, they just pass the face_unit
       // to the Riemann solver and let it do its stuff..

       
       // It uses as input parameters the density, pressure and velocities at both ends. They are assumed to be given in comoving coordinates, then they are
       // converted into physical units inside Riemann_Solver. 
       Riemann_solver(riemann_input, &riemann_output, face_unit, /*double press_tot_limiter TODO For now, just p>0: */ 0.0);

       psph->Frho = riemann_output.Fluxes.rho;
       psph->Fene = riemann_output.Fluxes.p;
       for (j=0;j<3;j++){ psph->Fmom[j] = riemann_output.Fluxes.v[j]; }

       // Now we de-boost the fluxes following Eq. A8 Hopkins 2015 (From hydra_core_meshless.h):
       /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame'
        which we do following Pakmor et al. 2011 */
       for(j=0;j<3;j++)
       {
           riemann_output.Fluxes.p += vFrame[j] * riemann_output.Fluxes.v[j];
           riemann_output.Fluxes.p += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
       }

       /* ok now we can actually apply this to the EOM */
       psph->Frho += riemann_output.Fluxes.rho*modApq; //IA TODO: Check if Face_Area_Norm is EQUAL to modAqp. It seems they are!
       psph->Fene += riemann_output.Fluxes.p*modApq; 
       for(j=0;j<3;j++) {psph->Fmom[j] += (riemann_output.Fluxes.v[j]+vFrame[j] * riemann_output.Fluxes.rho)*modApq;} // momentum flux (need to divide by mass) 

       // IA: we store this flux also for the q particle
       qsph->Frho -= psph->Frho;
       qsph->Fene -= psph->Fene;
       for(j=0;j<3;j++){ qsph->Fmom[j] -= psph->Fmom[j]; }



/*                
       // From hydra_evaluate.h //TODO: ALL
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                double dmass_holder = Fluxes.rho * dt_hydrostep, dmass_limiter;
                if(dmass_holder > 0) {dmass_limiter=P[j].Mass;} else {dmass_limiter=local.Mass;}
                dmass_limiter *= 0.1;
                if(fabs(dmass_holder) > dmass_limiter) {dmass_holder *= dmass_limiter / fabs(dmass_holder);}
                out.dMass += dmass_holder;
                out.DtMass += Fluxes.rho;
#ifndef BOX_SHEARING
                SphP[j].dMass -= dmass_holder;
#endif
                double gravwork[3]; gravwork[0]=Fluxes.rho*kernel.dp[0]; gravwork[1]=Fluxes.rho*kernel.dp[1]; gravwork[2]=Fluxes.rho*kernel.dp[2];
                for(k=0;k<3;k++) {out.GravWorkTerm[k] += gravwork[k];}
#endif
                for(k=0;k<3;k++) {out.Acc[k] += Fluxes.v[k];}
                out.DtInternalEnergy += Fluxes.p;

                // if this is particle j's active timestep, you should sent them the time-derivative information as well, for their subsequent drift operations 
                if(j_is_active_for_fluxes)
                {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                    SphP[j].DtMass -= Fluxes.rho;
                    for(k=0;k<3;k++) {SphP[j].GravWorkTerm[k] -= gravwork[k];}
#endif
                    for(k=0;k<3;k++) {SphP[j].HydroAccel[k] -= Fluxes.v[k];}
                    SphP[j].DtInternalEnergy -= Fluxes.p;


                //}

                /* if we have mass fluxes, we need to have metal fluxes if we're using them (or any other passive scalars) */
/*                
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                if(dmass_holder != 0)
                {
                }
#endif */

                /* --------------------------------------------------------------------------------- */
                /* don't forget to save the signal velocity for time-stepping! */
                /* --------------------------------------------------------------------------------- */
/*
                if(kernel.vsig > out.MaxSignalVel) out.MaxSignalVel = kernel.vsig;
                if(j_is_active_for_fluxes) {if(kernel.vsig > SphP[j].MaxSignalVel) SphP[j].MaxSignalVel = kernel.vsig;}
*/

//            } // for(n = 0; n < numngb; n++) //
    } // IA: End of loop over neighbors



}
