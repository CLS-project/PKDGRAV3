/* File added by Isaac Alonso for computing
 * the hydrodynamical part using a mesh-free
 * method, following the work of Hopkins 2015
 */


#include "hydro.h" 
#include <stdio.h>


void test(){
      printf("Test function inside hydro.c \n");
}


void cubicSplineKernel(double r, double h) {
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

   B[XX] = E[YY]*E[ZZ] - 2.0*E[YZ]
   B[YY] = E[ZZ]*E[ZZ] - 2.0*E[XZ]
   B[ZZ] = E[YY]*E[XX] - 2.0*E[YX]

   B[XY] = E[XY]*E[ZZ] - E[YZ]*E[XZ]
   B[XZ] = E[XY]*E[YZ] - E[YY]*E[XZ]
   B[YZ] = E[XX]*E[YZ] - E[XY]*E[XZ]

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


/* TODO: clean variables! */
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {

    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double ih2,r2,rs1,rq,rp;
    double dx,dy,dz,dvx,dvy,dvz,dvdotdr;
    double pPoverRho2,pPoverRho2f,pMass;
    double qPoverRho2,qPoverRho2f;
    double ph,pc,pDensity,visc,hav,absmu,Accp,Accq;
    double fNorm,fNorm1,aFac,vFac,gammainv,dtC,dtMu,dtEst,dt2;
    float *pa,*qa;
    double E[6], B[6];  // IA: We are assumming 3D here!
    int i,pActive,qActive;
    uint8_t uNewRung;

    assert(!pkd->bNoParticleOrder);

    aFac = (smf->a);        /* comoving acceleration factor */
    vFac = (smf->bComove ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */
    gammainv = 1/smf->gamma;
    dtC = (1+0.6*smf->alpha)/(smf->a*smf->dEtaCourant);
    dtMu = (0.6*smf->beta)/(smf->a*smf->dEtaCourant);

    /* Particle p data */
    psph = pkdSph(pkd,p);
    pc = psph->c;
    pDensity = pkdDensity(pkd,p); // TODO: check where is this computed
    pMass = pkdMass(pkd,p);
    pPoverRho2 = gammainv*psph->c*psph->c/pDensity;
    pPoverRho2f = pPoverRho2; //TODO: What is this for?!
    ph = 0.5*fBall; /* IA: fBall seems to be the minimun distance to any neighbors ¿? Although it would be logical to be the maximum */
    pActive = pkdIsActive(pkd,p);
    pa = pkdAccel(pkd,p);

    /* IA: I do not get this, so I will use my own cubic spline kernel */
//    ih2 = 1/(ph*ph);
//    fNorm = 0.5*M_1_PI*ih2/ph;
//    fNorm1 = fNorm*ih2;	/* converts to physical u */

    /* IA: Compute the \omega(x_i) normalization factor */
    omega = 0.0;
    for (i=0; i<nSmooth; ++i){
       q = nnList[i].pPart;

       qh = 0.5*pkdBall(pkd,q); // TODO: Check this 

       rpq = sqrt(nnList[i].dist2);
       hpq = 0.5*(qh+qp); // IA: We symmetrize the kernel size (probably needed, not sure)

       omega += cubicSplineKernel(rpq, hpq);
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

       rpq = sqrt(nnList[i].dist2);
       hpq = 0.5*(qh+qp); // IA: We symmetrize the kernel size (probably needed, not sure)

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
       E[i] /= omega; 
    }

    /* IA: END of E matrix computation */

    /* IA: Now, we need to do the inverse */
    inverseMatrix(&E, &B);
    
    /* IA: Now we are prepared for computing the gradients! */
    for (i=0;i<nSmooth;++i) {
	q = nnList[i].pPart;

	dx = nnList[i].dx;
	dy = nnList[i].dy;
	dz = nnList[i].dz;

      qh = 0.5*pkdBall(pkd,q); // TODO: Check this 

      rpq = sqrt(nnList[i].dist2);
      hpq = 0.5*(qh+qp); // IA: We symmetrize the kernel size (probably needed, not sure)

      Wpq = cubicSplineKernel(rpq, hpq); //TODO: Re-use this from the previous loop!!
      psi = Wpq/omega;

      psiTildex = (B[XX]*dx + B[XY]*dy + B[XZ]*dz)*psi;
      psiTildey = (B[XY]*dx + B[YY]*dy + B[YZ]*dz)*psi;
      psiTildez = (B[XZ]*dx + B[YZ]*dy + B[ZZ]*dz)*psi;

      // Density gradients
      psph->densGradx += psiTildex*(pkdDensity(pkd,q) - pDensity);
      psph->densGrady += psiTildey*(pkdDensity(pkd,q) - pDensity);
      psph->densGradz += psiTildez*(pkdDensity(pkd,q) - pDensity);

      // Momentum gradients
      pmom = pDensity*psph->vPred[0]; // TODO :check this vPred!
      dmom = pkdDensity(pkd,q)*qsph->vPred[0] - pmom;
      psph->momxGradx += psiTildex * dmom;
      psph->momxGrady += psiTildey * dmom;
      psph->momxGradz += psiTildez * dmom;

      pmom = pDensity*psph->vPred[1]; // TODO :check this vPred!
      dmom = pkdDensity(pkd,q)*qsph->vPred[1] - pmom;
      psph->momyGradx += psiTildex * dmom;
      psph->momyGrady += psiTildey * dmom;
      psph->momyGradz += psiTildez * dmom;

      pmom = pDensity*psph->vPred[2]; // TODO :check this vPred!
      dmom = pkdDensity(pkd,q)*qsph->vPred[2] - pmom;
      psph->momzGradx += psiTildex * dmom;
      psph->momzGrady += psiTildey * dmom;
      psph->momzGradz += psiTildez * dmom;

      // Internal energy gradients TODO: also, chech this uPred!
      dIntEne = qsph->uPred + 0.5*pkdDensity(pkd,q)*(qsph->vPred[0]*qsph->vPred[0] + qsph->vPred[1]*qsph->vPred[1] + qsph->vPred[2]*qsph->vPred[2])
              - psph->uPred - 0.5*pDensity*(psph->vPred[0]*psph->vPred[0] + psph->vPred[1]*psph->vPred[1] + psph->vPred[2]*psph->vPred[2]);
      psph->intEneGradx += psiTildex * dIntEne; 
      psph->intEneGrady += psiTildey * dIntEne;
      psph->intEneGradz += psiTildez * dIntEne;


      // Dont know if needed, but I will also store the psiTildes. THERE IS ONE FOR EACH NEIGHBOR!!
//      psph->psiTildex = psiTildex;
//      psph->psiTildey = psiTildey;
//      psph->psiTildez = psiTildez;
    }

    /* TODO: I should make sure that psph is not cleared after a smSmooth call!!! */

    /* TODO: Here I should limit the gradients if needed */


    }


/* IA: This routine will extrapolate the *already limited* gradient to the 'faces' and solve the 1D riemann problem 
 * For now, the 1D riemann flux will be computed TWICE for a given face, one for each adjacent particles... This 
 * is highly sub-optimal, and I should think a way to 'mark' which riemann fluxes has been already computed.
 * IDEA: maybe a 64-bit with 0/1? But I think that I have no info about the neighbors of my neighbor..*/
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    int i;

    psph = pkdSph(pkd, p);   

    for (i=0;i<nSmooth;++i){

	 q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);

	 dx = nnList[i].dx;
	 dy = nnList[i].dy;
	 dz = nnList[i].dz;

       qh = 0.5*pkdBall(pkd,q); // TODO: Check this 

       // Velocity of the quadrature mid-point 
       vFramex = 0.5*(psph->vPred[0]+qsph->vPred[0]);
       vFramey = 0.5*(psph->vPred[1]+qsph->vPred[1]);
       vFramez = 0.5*(psph->vPred[2]+qsph->vPred[2]);

       // We boost to the reference of the p-q 'face'
       pvx = psph->vPred[0] - vFramex;
       pvy = psph->vPred[1] - vFramey;
       pvz = psph->vPred[2] - vFramez;

       qvx = qsph->vPred[0] - vFramex;
       qvy = qsph->vPred[1] - vFramey;
       qvz = qsph->vPred[2] - vFramez;


       // Rotate the vectors to be aligned with the face. This means that the vx velocity will be aligned
       // with the Aij = Vi * psiTildej(xi) - Vj*psiTildei(xj)
       // I need the psiTilde from BOTH sides!!!! Will I need to store such a big array (64 doubles/floats) per particle...
       // Other option would be to just store the B matrix, (6 doubles/floats) and reconstruct the psiTilde using the kernel function...
       // AS ALWAYS: MEMORY VS NB OF OPERATIONS!
       
       // Once we have the fluxes in the correct SoR, we can compute the riemann flux
       // TODO: Use riemann.h from GIZMO!
 
    }



}
