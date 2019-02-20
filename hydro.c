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
    
    /* IA: Now we are prepared for computing the gradients! */
    for (i=0;i<nSmooth;++i) {
	q = nnList[i].pPart;

	dx = nnList[i].dx;
	dy = nnList[i].dy;
	dz = nnList[i].dz;

      qh = 0.5*pkdBall(pkd,q); // TODO: Check this 

      rpq = sqrt(nnList[i].fDist2);
      hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

      Wpq = cubicSplineKernel(rpq, hpq); //TODO: Re-use this from the previous loop!!
      psi = Wpq/psph->omega;

      /* IA: By the way, we use this loop to compute the volume of each particle.
       * This may not be needed, depending on how the pkdDensity has been computed, so TODO: check! */
      qsph->V += psi; //Contribution of the 'volume' at position p, of the particle q. Note that in order
                      // to use this 'effective volume' value, a full loop around all particles must be done!!
                      // TODO: This must be initialized to zero outside the smooth loop!!!
                      // TODO: This can also be done using Eq. 27 Hopkins 2015
                      // Indeed, I am not sure my implementation is ok!!

      psiTilde[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
      psiTilde[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
      psiTilde[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

      // Density gradients
      diff = (pkdDensity(pkd,q) - pDensity);
      for (j=0; j<3; j++) {
         psph->densGrad[j] += psiTilde[j]*diff;
      }

      // Momentum gradients
      diff = qsph->vPred[0] - psph->vPred[0];// TODO :check this vPred!
      for (j=0; j<3; j++) {
         psph->vxGrad[j] += psiTilde[j] * diff;
      }

      diff = qsph->vPred[1] - psph->vPred[1];// TODO :check this vPred!
      for (j=0; j<3; j++) {
         psph->vyGrad[j] += psiTilde[j] * diff;
      }

      diff = qsph->vPred[2] - psph->vPred[2];// TODO :check this vPred!
      for (j=0; j<3; j++) {
         psph->vzGrad[j] += psiTilde[j] * diff;
      }

      // Internal energy gradients TODO: also, chech this uPred!
      diff = qsph->uPred + 0.5*pkdDensity(pkd,q)*(qsph->vPred[0]*qsph->vPred[0] + qsph->vPred[1]*qsph->vPred[1] + qsph->vPred[2]*qsph->vPred[2])
           - psph->uPred - 0.5*pDensity*(psph->vPred[0]*psph->vPred[0] + psph->vPred[1]*psph->vPred[1] + psph->vPred[2]*psph->vPred[2]);
      for (j=0; j<3; j++){
         psph->intEneGrad[j] += psiTilde[j] * diff;
      }


      // Dont know if needed, but I will also store the psiTildes. THERE IS ONE FOR EACH NEIGHBOR!!
//      psph->psiTildex = psiTildex;
//      psph->psiTildey = psiTildey;
//      psph->psiTildez = psiTildez;
    }

    /* TODO: I should make sure that psph is not cleared after a smSmooth call!!! */

    /* TODO: Here I should limit the gradients if needed */


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
    double ph,qh, hpq, modApq, rpq, dx,dy,dz,Wpq,psi_p,psi_q;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3]; 
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

       // Velocity of the quadrature mid-point 
       for (j=0; j<3; j++){
          vFrame[j] = 0.5*(psph->vPred[j]+qsph->vPred[j]);
          // We boost to the reference of the p-q 'face'
          pv[j] = psph->vPred[j] - vFrame[j];
          qv[j] = qsph->vPred[j] - vFrame[j];
       }


       // Extrapolation to the faces TODO: check dx/y/z sign!
       riemann_input.L.rho = pkdDensity(pkd,p) + 0.5*( psph->densGrad[0]*dx + psph->densGrad[1]*dy + psph->densGrad[2]*dz );
       riemann_input.L.v[0] = pv[0] + 0.5*( psph->vxGrad[0]*dx + psph->vxGrad[1]*dy + psph->vxGrad[2]*dz );
       riemann_input.L.v[1] = pv[1] + 0.5*( psph->vyGrad[0]*dx + psph->vyGrad[1]*dy + psph->vyGrad[2]*dz );
       riemann_input.L.v[2] = pv[2] + 0.5*( psph->vzGrad[0]*dx + psph->vzGrad[1]*dy + psph->vzGrad[2]*dz );
       riemann_input.L.u = psph->uPred + 0.5*( psph->intEneGrad[0]*dx + psph->intEneGrad[1]*dy + psph->intEneGrad[2]*dz );
       riemann_input.L.cs = sqrt(pkd->param.dConstGamma * riemann_input.L.p / riemann_input.L.rho) ;// This can be done inside Riemann solver
       riemann_input.L.p = riemann_input.L.u * (pkd->param.dConstGamma - 1.0) * riemann_input.L.rho; // Same

       riemann_input.R.rho = pkdDensity(pkd, q) - 0.5*( qsph->densGrad[0]*dx + qsph->densGrad[1]*dy + qsph->densGrad[2]*dz );
       riemann_input.R.v[0] = pv[0] - 0.5*( qsph->vxGrad[0]*dx + qsph->vxGrad[1]*dy + qsph->vxGrad[2]*dz );
       riemann_input.R.v[1] = pv[1] - 0.5*( qsph->vyGrad[0]*dx + qsph->vyGrad[1]*dy + qsph->vyGrad[2]*dz );
       riemann_input.R.v[2] = pv[2] - 0.5*( qsph->vzGrad[0]*dx + qsph->vzGrad[1]*dy + qsph->vzGrad[2]*dz );
       riemann_input.R.u = qsph->uPred - 0.5*( qsph->intEneGrad[0]*dx + qsph->intEneGrad[1]*dy + qsph->intEneGrad[2]*dz );
       riemann_input.R.cs = sqrt(pkd->param.dConstGamma * riemann_input.R.p / riemann_input.R.rho) ;// This can be done inside Riemann solver
       riemann_input.R.p = riemann_input.R.u * (pkd->param.dConstGamma - 1.0) * riemann_input.R.rho; // Same



       // Rotate the vectors to be aligned with the face. This means that the vx velocity will be aligned
       // with the Aij = Vi * psiTildej(xi) - Vj*psiTildei(xj)
       // [19/02/19] I need the psiTilde from BOTH sides!!!! Will I need to store such a big array (64 doubles/floats) per particle...
       // Other option would be to just store the B matrix, (6 doubles/floats) and reconstruct the psiTilde using the kernel function...
       // AS ALWAYS: MEMORY VS NB OF OPERATIONS!
       // [20/02/19] I think that I will go with saving the B matrix.. anyway this three different loop thing will have to be changed sooner or later
       rpq = sqrt(nnList[i].fDist2);
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

       Wpq = cubicSplineKernel(rpq, hpq); 
       psi_p = Wpq/psph->omega;
       psi_q = Wpq/qsph->omega;

       psiTilde_p[0] =  (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi_p;
       psiTilde_p[1] =  (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi_p;
       psiTilde_p[2] =  (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi_p;

       psiTilde_q[0] = -(qsph->B[XX]*dx + qsph->B[XY]*dy + qsph->B[XZ]*dz)*psi_q;
       psiTilde_q[1] = -(qsph->B[XY]*dx + qsph->B[YY]*dy + qsph->B[YZ]*dz)*psi_q;
       psiTilde_q[2] = -(qsph->B[XZ]*dx + qsph->B[YZ]*dy + qsph->B[ZZ]*dz)*psi_q;

       modApq = 0.0;
       for (j=0; j<3; j++){
          Apq[j] = psph->V*psiTilde_q[j] - qsph->V*psiTilde_p[j];
          modApq += Apq[j]*Apq[j];
       }
       modApq = sqrt(modApq);

       for (j=0; j<3; j++){
          face_unit[j] = Apq[j]/modApq;
       }

       // Now we rotate the velocity vectors in this new SoR TODO: Implement this using arrays and loops..
       // It seems that in hydra_core_meshless.h they do not rotate any velocity, they just pass the face_unit
       // to the Riemann solver and let it do its stuff..

       
       // Once we have the states in the correct SoR, we can compute the riemann flux
       // It uses as input parameters the density, pressure and velocities at both ends. They are assumed to be given in comoving coordinates, then they are
       // converted into physical units inside Riemann_Solver. 
 
       Riemann_solver(riemann_input, &riemann_output, face_unit, /*double press_tot_limiter TODO For now, just p>0: */ 0.0);

       psph->Frho = riemann_output.Fluxes.rho;
       psph->Fene = riemann_output.Fluxes.p;
       for (j=0;j<3;j++){ psph->Fmom[j] = riemann_output.Fluxes.v[j]; }

       // Now we de-boost the fluxes following Eq. A8 Hopkins 2015 (From hydra_core_meshless.h)
       /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame'
        which we do following Pakmor et al. 2011 */
       for(j=0;j<3;j++)
       {
           psph->Fene += vFrame[j] * riemann_output.Fluxes.v[j];
//#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
           /* Riemann_out->Fluxes.rho is un-modified */
           psph->Fene += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
           psph->Fmom[j] += vFrame[j] * riemann_output.Fluxes.rho; /* just boost by frame vel (as we would in non-moving frame) */
//#endif
       }

       /* ok now we can actually apply this to the EOM */
//#if defined(HYDRO_MESHLESS_FINITE_VOLUME)
       psph->Frho *= modApq; //Face_Area_Norm  ;//IA TODO: Check if Face_Area_Norm is EQUAL to modAqp
//#endif
       psph->Fene *= modApq; //Face_Area_Norm ; // this is really Dt of --total-- energy, need to subtract KE component for e */
       for(j=0;j<3;j++) {psph->Fmom[j] *= modApq;} //Face_Area_Norm;} // momentum flux (need to divide by mass) //




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
#ifdef MAGNETIC
#ifndef HYDRO_SPH
                for(k=0;k<3;k++) {out.Face_Area[k] += Face_Area_Vec[k];}
#endif
#ifndef FREEZE_HYDRO
                for(k=0;k<3;k++) {out.DtB[k]+=Fluxes.B[k];}
                out.divB += Fluxes.B_normal_corrected;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME) // mass-based phi-flux
                out.DtPhi += Fluxes.phi;
#endif
#ifdef HYDRO_SPH
                for(k=0;k<3;k++) {out.DtInternalEnergy+=magfluxv[k]*local.Vel[k]/All.cf_atime;}
                out.DtInternalEnergy += resistivity_heatflux;
#else
                double wt_face_sum = Face_Area_Norm * (-face_area_dot_vel+face_vel_i);
                out.DtInternalEnergy += 0.5 * kernel.b2_i*All.cf_a2inv*All.cf_a2inv * wt_face_sum;
#ifdef DIVBCLEANING_DEDNER
                for(k=0; k<3; k++)
                {
                    out.DtB_PhiCorr[k] += Riemann_out.phi_normal_db * Face_Area_Vec[k];
                    out.DtB[k] += Riemann_out.phi_normal_mean * Face_Area_Vec[k];
                    out.DtInternalEnergy += Riemann_out.phi_normal_mean * Face_Area_Vec[k] * local.BPred[k]*All.cf_a2inv;
                }
#endif
#ifdef MHD_NON_IDEAL
                for(k=0;k<3;k++) {out.DtInternalEnergy += local.BPred[k]*All.cf_a2inv*bflux_from_nonideal_effects[k];}
#endif
#endif
#endif
#endif  // magnetic //

                // if this is particle j's active timestep, you should sent them the time-derivative information as well, for their subsequent drift operations 
                if(j_is_active_for_fluxes)
                {
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                    SphP[j].DtMass -= Fluxes.rho;
                    for(k=0;k<3;k++) {SphP[j].GravWorkTerm[k] -= gravwork[k];}
#endif
                    for(k=0;k<3;k++) {SphP[j].HydroAccel[k] -= Fluxes.v[k];}
                    SphP[j].DtInternalEnergy -= Fluxes.p;

#ifdef MAGNETIC
#ifndef HYDRO_SPH
                    for(k=0;k<3;k++) {SphP[j].Face_Area[k] -= Face_Area_Vec[k];}
#endif
#ifndef FREEZE_HYDRO
                    for(k=0;k<3;k++) {SphP[j].DtB[k]-=Fluxes.B[k];}
                    SphP[j].divB -= Fluxes.B_normal_corrected;
#if defined(DIVBCLEANING_DEDNER) && defined(HYDRO_MESHLESS_FINITE_VOLUME) // mass-based phi-flux
                    SphP[j].DtPhi -= Fluxes.phi;
#endif
#ifdef HYDRO_SPH
                    for(k=0;k<3;k++) {SphP[j].DtInternalEnergy-=magfluxv[k]*VelPred_j[k]/All.cf_atime;}
                    SphP[j].DtInternalEnergy += resistivity_heatflux;
#else
                    double wt_face_sum = Face_Area_Norm * (-face_area_dot_vel+face_vel_j);
                    SphP[j].DtInternalEnergy -= 0.5 * kernel.b2_j*All.cf_a2inv*All.cf_a2inv * wt_face_sum;
#ifdef DIVBCLEANING_DEDNER
                    for(k=0; k<3; k++)
                    {
                        SphP[j].DtB_PhiCorr[k] -= Riemann_out.phi_normal_db * Face_Area_Vec[k];
                        SphP[j].DtB[k] -= Riemann_out.phi_normal_mean * Face_Area_Vec[k];
                        SphP[j].DtInternalEnergy -= Riemann_out.phi_normal_mean * Face_Area_Vec[k] * BPred_j[k]*All.cf_a2inv;
                    }
#endif
#ifdef MHD_NON_IDEAL
                    for(k=0;k<3;k++) {SphP[j].DtInternalEnergy -= BPred_j[k]*All.cf_a2inv*bflux_from_nonideal_effects[k];}
#endif
#endif
#endif
#endif // magnetic // */

                //}

                /* if we have mass fluxes, we need to have metal fluxes if we're using them (or any other passive scalars) */
/*                
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                if(dmass_holder != 0)
                {
#ifdef METALS
                    if(Fluxes.rho > 0)
                    {
                        // particle i gains mass from particle j 
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            out.Dyield[k] += (P[j].Metallicity[k] - local.Metallicity[k]) * dmass_holder;
                    } else {
                        // particle j gains mass from particle i 
                        dmass_holder /= -P[j].Mass;
                        for(k=0;k<NUM_METAL_SPECIES;k++)
                            P[j].Metallicity[k] += (local.Metallicity[k] - P[j].Metallicity[k]) * dmass_holder;
                    }
#endif 
                }
#endif */

                /* --------------------------------------------------------------------------------- */
                /* don't forget to save the signal velocity for time-stepping! */
                /* --------------------------------------------------------------------------------- */
/*
                if(kernel.vsig > out.MaxSignalVel) out.MaxSignalVel = kernel.vsig;
                if(j_is_active_for_fluxes) {if(kernel.vsig > SphP[j].MaxSignalVel) SphP[j].MaxSignalVel = kernel.vsig;}
#ifdef WAKEUP
                if(!(TimeBinActive[P[j].TimeBin]))
                {
                    if(kernel.vsig > WAKEUP*SphP[j].MaxSignalVel) PPPZ[j].wakeup = 1;
#if (SLOPE_LIMITER_TOLERANCE <= 0)
                    if(local.Timestep*WAKEUP < TimeStep_J) PPPZ[j].wakeup = 1;
#endif
                }
#endif
*/

//            } // for(n = 0; n < numngb; n++) //
    } // IA: End of loop over neighbors



}
