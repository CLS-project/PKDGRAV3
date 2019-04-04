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

   B[XX] = (E[YY]*E[ZZ] - E[YZ]*E[YZ])/det;
   B[YY] = (E[XX]*E[ZZ] - E[XZ]*E[XZ])/det;
   B[ZZ] = (E[YY]*E[XX] - E[XY]*E[XY])/det;

   B[XY] = -(E[XY]*E[ZZ] - E[YZ]*E[XZ])/det;
   B[XZ] = (E[XY]*E[YZ] - E[YY]*E[XZ])/det;
   B[YZ] = -(E[XX]*E[YZ] - E[XY]*E[XZ])/det;

   if (det==0) {
      printf("Singular matrix!");
      abort();
   }

}


/* IA: We need to clear the SPHFIELD for all particles */
void initHydroLoop(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    int i;
    assert(!pkd->bNoParticleOrder);
//    if (pkdIsActive(pkd,p)) {
	SPHFIELDS *psph = pkdSph(pkd,p);
	psph->uDot = 0;
	psph->fMetalsDot = 0;
      for (i=0;i<6;i++) { psph->B[i] = 0.0; }
      psph->omega = 0.0;
      psph->Frho = 0.0;
      psph->Fene = 0.0;
      p->uNewRung = 0;
      for (i=0;i<3;i++) { 
         psph->Fmom[i] = 0.0;
	   pkdAccel(pkd,p)[i] = 0;
	   pkdAccel(pkd,p)[i] = 0;
	   pkdAccel(pkd,p)[i] = 0;
	}
	
    }

void initHydroLoopCached(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    int i;
//    if (pkdIsActive(pkd,p)) {
	SPHFIELDS *psph = pkdSph(pkd,p);
	psph->uDot = 0;
	psph->fMetalsDot = 0;
      for (i=0;i<6;i++) { psph->B[i] = 0.0; }
      psph->omega = 0.0;
      psph->Frho = 0.0;
      psph->Fene = 0.0;
      p->uNewRung = 0;
      for (i=0;i<3;i++) { 
         psph->Fmom[i] = 0.0;
	   pkdAccel(pkd,p)[i] = 0;
	   pkdAccel(pkd,p)[i] = 0;
	   pkdAccel(pkd,p)[i] = 0;
	}
    }

/* IA: I guess that this is called when we 'merge' information coming from different processors
 * about the same particle.. For the first hydro loop, all quantities are aditives so easy!
 */
void combFirstHydroLoop(void *vpkd, void *p1,void *p2) {
    PKD pkd = (PKD) vpkd;
    assert(!pkd->bNoParticleOrder);
    /*
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
	} */
    SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);
    int i;

    psph1->omega += psph2->omega;

    //IA: the fluxes are not added because they has not been yet computed!

    }


void combSecondHydroLoop(void *vpkd, void *p1,void *p2) {
    PKD pkd = (PKD) vpkd;
    assert(!pkd->bNoParticleOrder);
    SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);
    int i;

    // IA: Not sure about this.. what if one side has been limited?
    /*
    for (i=0;i<3;i++){ 
       psph1->gradRho[i] += psph2->gradRho[i];
       psph1->gradVx[i] += psph2->gradVx[i];
       psph1->gradVy[i] += psph2->gradVy[i];
       psph1->gradVz[i] += psph2->gradVz[i];
       psph1->gradP[i] += psph2->gradP[i];
       }
    for (i=0;i<6;i++){ psph1->B[i] += psph2->B[i]; }
    */
    }


void combThirdHydroLoop(void *vpkd, void *p1,void *p2) {
    PKD pkd = (PKD) vpkd;
    assert(!pkd->bNoParticleOrder);
    SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);
    int i;

    for (i=0;i<3;i++){ psph1->Fmom[i] += psph2->Fmom[i]; }
    psph1->Frho += psph2->Frho;
    psph1->Fene += psph2->Fene;
    if (((PARTICLE *) p2)->uNewRung > ((PARTICLE *) p1)->uNewRung) 
       ((PARTICLE *) p1)->uNewRung = ((PARTICLE *) p2)->uNewRung;
    }

void initHydroFluxes(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    int i;
//    if (pkdIsActive(pkd,p)) {
	SPHFIELDS *psph = pkdSph(pkd,p);
      psph->Frho = 0.0;
      psph->Fene = 0.0;
      psph->uNewRung = 0;
      for (i=0;i<3;i++) { 
         psph->Fmom[i] = 0.0;
	}
    }

void initHydroGradients(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    SPHFIELDS *psph = pkdSph(pkd,p);
    int j;
    for (j=0; j<3;j++){
      psph->gradRho[j] = 0.0;
      psph->gradVx[j] = 0.0;
      psph->gradVy[j] = 0.0;
      psph->gradVz[j] = 0.0;
      psph->gradP[j] = 0.0;
      }
    }



void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph;
    double ph, qh, rpq, hpq, dx,dy,dz,Wpq;
    int i;

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
    ph = fBall; /* IA: fBall seems to be the minimun distance to any neighbors ¿? Although it would be logical to be the maximum */


    /* IA: Compute the \omega(x_i) normalization factor */
    psph->omega = 0.0;
    for (i=0; i<nSmooth; ++i){
       q = nnList[i].pPart;

       qh = pkdBall(pkd,q); // TODO: Check this 

       rpq = sqrt(nnList[i].fDist2);
       hpq = ph;// 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

       psph->omega += cubicSplineKernel(rpq, hpq);
    }


    /* IA: We compute the density making use of Eq. 27 Hopkins 2015 */
    pkdSetDensity(pkd,p, pkdMass(pkd,p)*psph->omega);
/*    if(pkdDensity(pkd,p) > 3)*/  //  printf("mass %e \t omega %e \t density %e \n",pkdMass(pkd,p), psph->omega, pkdMass(pkd,p)*psph->omega);
}


void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {

    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double E[6];  // IA: We are assumming 3D here!
    double ph, qh, rpq, hpq,Wpq, dx,dy,dz, diff;
    double rho_max, rho_min, vx_max, vx_min, vy_max, vy_min, vz_max, vz_min, p_max, p_min;
    double limRho, limVx, limVy, limVz, limP;
    double psi, psiTilde_p[3];
    float  *pv, *qv;
    int i, j;
    psph = pkdSph(pkd,p);
    ph = fBall;

//    printf("(hydroGradients) BEGIN \n");

    /* IA: Compute the E matrix (Hopkins 2015, eq 14) */
    for (i=0; i<6;++i){
       E[i] = 0.0; 
    }

    for (i=0; i<nSmooth;++i){
       q = nnList[i].pPart;

       dx = -nnList[i].dx;
       dy = -nnList[i].dy;
       dz = -nnList[i].dz;

       qh = pkdBall(pkd,q);

       rpq = sqrt(nnList[i].fDist2);
       hpq = ph; //0.5*(qh+ph); 

       Wpq = cubicSplineKernel(rpq, hpq);

       E[XX] += dx*dx*Wpq;
       E[YY] += dy*dy*Wpq;
//       printf("dx %e dy %e dz %e Wpq %e \n", dx, dy, dz, Wpq);
       E[ZZ] += dz*dz*Wpq;

       E[XY] += dy*dx*Wpq;
       E[XZ] += dz*dx*Wpq;
       E[YZ] += dy*dz*Wpq;

    }

    /* IA: Normalize the matrix */
    for (i=0; i<6;++i){
       E[i] /= psph->omega; 
    }

    /* IA: END of E matrix computation */
//    printf("E_q [XX] %e \t [XY] %e \t [XZ] %e \n \t \t \t [YY] %e \t [YZ] %e \n \t\t\t \t \t \t [ZZ] %e \n", E[XX], E[XY], E[XZ], E[YY], E[YZ], E[ZZ]);

    /* IA: Now, we need to do the inverse */
    inverseMatrix(E, psph->B);
    
    /* IA: There is nothing more that we can do in this loop, as all the particle must have their densities
     * and B matrices computed */

    // IA: DEBUG: check if E^{-1} = B
    /*
    double *B = psph->B;
    double unityXX, unityYY, unityZZ, unityXY, unityXZ, unityYZ;
    unityXX = E[XX]*B[XX] + E[XY]*B[XY] + E[XZ]*B[XZ];
    unityYY = E[XY]*B[XY] + E[YY]*B[YY] + E[YZ]*B[YZ];
    unityZZ = E[XZ]*B[XZ] + E[YZ]*B[YZ] + E[ZZ]*B[ZZ];

    unityXY = E[XX]*B[XY] + E[XY]*B[YY] + E[XZ]*B[YZ]; 
    unityXZ = E[XX]*B[XZ] + E[XY]*B[YZ] + E[XZ]*B[ZZ]; 
    unityYZ = E[XY]*B[XZ] + E[YY]*B[YZ] + E[YZ]*B[ZZ]; 
    
    printf("XX %e \t YY %e \t ZZ %e \n", unityXX, unityYY, unityZZ);
    printf("XY %e \t XZ %e \t YZ %e \n", unityXY, unityXZ, unityYZ);
    */

    /* IA: Now we can compute the gradients 
     * This and the B matrix computation could be done in the first hydro loop (where omega is computed) but
     * at that time we do not know the density of all particles (because it depends on omega) */

//    printf("(hydroGradients) Begin GRADIENTS \n");
    for (i=0;i<nSmooth;++i){

	 q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);

	 dx = -nnList[i].dx;
	 dy = -nnList[i].dy;
	 dz = -nnList[i].dz;

       if (dx==0 && dy==0 && dz==0) continue;

       qh = pkdBall(pkd,q); 
       ph = fBall;
       rpq = sqrt(nnList[i].fDist2);
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)
//       hpq = ph;
//       hpq = qh > ph ? qh : ph; 

       Wpq = cubicSplineKernel(rpq, hpq); 
       psi = Wpq/psph->omega;
//       printf("omega_p %e omega_q %e \n", psph->omega, qsph->omega);
//       printf("psi_p %e \t psi_q %e \n", psi_p, psi_q);

       psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
       psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
       psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

       diff = pkdDensity(pkd, q) - pkdDensity(pkd, p);
       for (j=0; j<3; j++) {psph->gradRho[j] += diff*psiTilde_p[j]; }
       diff = qsph->vPred[0] - psph->vPred[0];
       for (j=0; j<3; j++) {psph->gradVx[j] += diff*psiTilde_p[j]; }
       diff = qsph->vPred[1] - psph->vPred[1];
       for (j=0; j<3; j++) {psph->gradVy[j] += diff*psiTilde_p[j]; }
       diff = qsph->vPred[2] - psph->vPred[2];
       for (j=0; j<3; j++) {psph->gradVz[j] += diff*psiTilde_p[j]; }
       diff = qsph->P - psph->P;
       for (j=0; j<3; j++) {psph->gradP[j] += diff*psiTilde_p[j]; }

       }

     /* IA: Now we can limit them */
//    printf("(hydroGradients) Begin LIMITER \n");
    
    // IA: First step, compute the maximum and minimum difference of each variable
    rho_min= vx_min= vy_min= vz_min= p_min =  HUGE_VAL;
    rho_max= vx_max= vy_max= vz_max= p_max = -HUGE_VAL;

    pv = pkdVel(pkd,p);
    for (i=0; i<nSmooth;++i){
//       if (nnList[i].dx==0 && nnList[i].dy==0 && nnList[i].dz==0) continue;
       q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);
       qv = pkdVel(pkd,q);

       if (pkdDensity(pkd,q) < rho_min) rho_min = pkdDensity(pkd,q);
       if (pkdDensity(pkd,q) > rho_max) rho_max = pkdDensity(pkd,q);

       if (qv[0] < vx_min) vx_min = qv[0];
       if (qv[0] > vx_max) vx_max = qv[0];

       if (qv[1] < vy_min) vy_min = qv[1];
       if (qv[1] > vy_max) vy_max = qv[1];

       if (qv[2] < vz_min) vz_min = qv[2];
       if (qv[2] > vz_max) vz_max = qv[2];

       if (qsph->P < p_min) p_min = qsph->P;
       if (qsph->P > p_max) p_max = qsph->P;
    }

    limRho= limVx= limVy= limVz= limP = 1.;
    double maxdx = 0.0; //IA TODO For debugging
    for (i=0; i<nSmooth;++i){
	 dx = -nnList[i].dx;
	 dy = -nnList[i].dy;
	 dz = -nnList[i].dz;
       if (dx==0 && dy==0 && dz==0) continue;
       q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);
       if (maxdx < fabs(dx)) maxdx = fabs(dx);

       // Las diferencias se pueden poner fuera, siempre son con p
       BarthJespersenLimiter(&limRho, psph->gradRho, rho_max-pkdDensity(pkd,p), rho_min-pkdDensity(pkd,p), dx, dy, dz);
       BarthJespersenLimiter(&limVx, psph->gradVx, vx_max-pv[0], vx_min-pv[0], dx, dy, dz);
       BarthJespersenLimiter(&limVy, psph->gradVy, vy_max-pv[1], vy_min-pv[1], dx, dy, dz);
       BarthJespersenLimiter(&limVz, psph->gradVz, vz_max-pv[2], vz_min-pv[2], dx, dy, dz);
       BarthJespersenLimiter(&limP, psph->gradP, p_max-psph->P, p_min-psph->P, dx, dy, dz);
    }

    for (j=0; j<3; j++){
       psph->gradRho[j] *= limRho;
       psph->gradVx[j] *= limVx;
       psph->fMetals = maxdx;
       psph->gradVy[j] *= limVy;
       psph->gradVz[j] *= limVz;
       psph->gradP[j] *= limP;
    }
    /* END OF LIMITER */
//    if (limRho<0.5) {
//      printf("Limiter: rho %e \t v %e %e %e \t p %e \n", limRho, limVx, limVy, limVz, limP);
//      printf("x %e y %e z %e \n", pkdPos(pkd,p,0), pkdPos(pkd,p,1), pkdPos(pkd,p,2));
//      printf("gradRho %e \n", psph->gradRho[0]);
//    }
    }


void BarthJespersenLimiter(double* limVar, double* gradVar, double var_max, double var_min, double dx, double dy, double dz){
    double diff, lim;

    diff = gradVar[0]*dx + gradVar[1]*dy + gradVar[2]*dz;
    if (var_min > 0) { var_min=0; } //IA: Can happen due to machine precision
    if (var_max < 0) { var_max=0; } //IA: Can happen due to machine precision
    if (diff > 0.) {
       lim = var_max/diff;
    }else if (diff < 0.){
       lim = var_min/diff;
    }else{
       lim = 1.;
    }
    if (lim > 1.) lim = 1.; // min(1,lim)
    if (lim < (*limVar) && lim>=0.) {/*printf("aa %e \n", lim);*/ *limVar = lim;} //IA: the second condition can happen due to machine precision errors when dealing with constant fields
    // FIXME IA: Option to avoid extrapolation or limiter
//    *limVar = 1.0;
//    *limVar = 0.0;
}





/* IA: This routine will extrapolate the primitives to the 'faces' and solve the 1D riemann problem. 
 * For now, the 1D riemann flux will be computed TWICE for a given face, one for each adjacent particles */ 
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double pv[3], qv[3], vFrame[3];
    double ph,qh, hpq, modApq, rpq, dx,dy,dz,Wpq,psi_p,psi_q,pDensity,pDeltaHalf, qDeltaHalf, dt2, dtEst, vsig_pq, dvDotd, dvDotdr;
    double vxdiff, vydiff, vzdiff, pdivv, qdivv, pdiffOverRhop, pdiffOverRhoq, psi;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3], dr[3]; 
    struct Input_vec_Riemann riemann_input;
    struct Riemann_outputs riemann_output; //TODO: do this outside the loop
    uint8_t uNewRung;
    int i,j;

    psph = pkdSph(pkd, p);   
    ph = fBall;

// printf("dDelta %e \n", smf->dDelta/(1<<p->uRung));
    pDeltaHalf = 0.0; //0.5*( smf->dDelta/(1<<p->uRung) );

    pDensity = pkdDensity(pkd,p);

    dtEst = HUGE_VAL;

    for (i=0;i<nSmooth;++i){

	 q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);
       qh = pkdBall(pkd,q); 

       qDeltaHalf = 0.0; //smf->dTime - qsph->lastUpdateTime + pDeltaHalf; 
//       printf("pDeltaHalf %e qDeltaHalf %e \n", pDeltaHalf, qDeltaHalf);

	 dx = nnList[i].dx;
	 dy = nnList[i].dy;
	 dz = nnList[i].dz;

       /* IA: in the nnList there is a 'copy' of the own particle, which we can omit as there is no fluxes
        * to be computed here */
       if (dx==0 && dy==0 && dz==0) continue;

       // Face where the riemann problem will be solved
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size 
//       hpq = qh > ph ? qh : ph; 
       rpq = sqrt(nnList[i].fDist2);

       Wpq = cubicSplineKernel(rpq, hpq); 
       if (Wpq==0.0){/*printf("hpq %e rpq %e \n", hpq, rpq);*/ continue; }

       psi = Wpq/psph->omega;
       psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
       psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
       psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

       psi = -Wpq/qsph->omega; // IA: minus because we are 'looking from' the other particle, thus -dr
       psiTilde_q[0] = (qsph->B[XX]*dx + qsph->B[XY]*dy + qsph->B[XZ]*dz)*psi;
       psiTilde_q[1] = (qsph->B[XY]*dx + qsph->B[YY]*dy + qsph->B[YZ]*dz)*psi;
       psiTilde_q[2] = (qsph->B[XZ]*dx + qsph->B[YZ]*dy + qsph->B[ZZ]*dz)*psi;

       modApq = 0.0;
       for (j=0; j<3; j++){
          Apq[j] = psiTilde_q[j]/psph->omega - psiTilde_p[j]/qsph->omega;
          modApq += Apq[j]*Apq[j];
       }
//       printf("modApq %e \n", modApq);
       modApq = sqrt(modApq);

       /* DEBUG
       if (modApq<=0.0) {
          printf("dx %e \t dy %e \t dz %e \n", dx, dy, dz);
          printf("rpq %e hpq %e ratio %e Wpq %e \n", rpq, hpq, rpq/hpq, Wpq);
       }
       assert(modApq>0.0); // Area should be positive!
       */


       for (j=0; j<3; j++){
          face_unit[j] = Apq[j]/modApq;
       }


       // Velocity of the quadrature mid-point 
       for (j=0; j<3; j++){
          vFrame[j] = 0.5*(psph->vPred[j]+qsph->vPred[j]);
          // We boost to the reference of the p-q 'face'
          pv[j] = psph->vPred[j] - vFrame[j];
          qv[j] = qsph->vPred[j] - vFrame[j];
       }
 
       // Mid-point rule 
       dr[0] = -0.5*dx;
       dr[1] = -0.5*dy;
       dr[2] = -0.5*dz;

      // Differences in the variables
      vxdiff = (qsph->vPred[0] - psph->vPred[0]);
      vydiff = (qsph->vPred[1] - psph->vPred[1]);
      vzdiff = (qsph->vPred[2] - psph->vPred[2]);


       // From Eqs 24,25 Hopkins 2015, to limit deltaT
       dvDotdr = (dx*vxdiff + dy*vydiff + dz*vzdiff);
       if (dvDotdr < 0) {
          vsig_pq = psph->c + qsph->c - dvDotdr/rpq;
       }else{
          vsig_pq = psph->c + qsph->c;
       }

       dt2 = 2.*smf->dEtaCourant * fBall /vsig_pq;	
	 if (dt2 < dtEst) dtEst=dt2;

      // Divergence of the velocity field for the forward in time prediction
      pdivv = psph->gradVx[0] + psph->gradVy[1] + psph->gradVz[2];
      qdivv = qsph->gradVx[0] + qsph->gradVy[1] + qsph->gradVz[2];

      pdivv *= pDeltaHalf;
      qdivv *= qDeltaHalf;

      riemann_input.L.rho = pkdDensity(pkd,p);
      riemann_input.R.rho = pkdDensity(pkd,q);
      riemann_input.L.v[0] = pv[0];
      riemann_input.R.v[0] = qv[0];
      riemann_input.L.v[1] = pv[1];
      riemann_input.R.v[1] = qv[1];
      riemann_input.L.v[2] = pv[2];
      riemann_input.R.v[2] = qv[2];
      riemann_input.L.p = psph->P;
      riemann_input.R.p = qsph->P;

//      printf("1) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("1) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

      // We add the gradients terms (from extrapolation and forward prediction)
      for (j=0; j<3; j++) {
//         if (fabs(riemann_input.L.rho - riemann_input.R.rho) > 0.1 && psph->gradRho[j]>1.0){
//           printf("L.rho %e L.rho_face %e R.rho %e R.rho_face %e %e\n", riemann_input.L.rho, riemann_input.L.rho+dr[j]*psph->gradRho[j], riemann_input.R.rho, riemann_input.R.rho-dr[j]*qsph->gradRho[j], psph->gradRho[j]);
//         }
         riemann_input.L.rho += ( dr[j] - pDeltaHalf*pv[j])*psph->gradRho[j];
         riemann_input.R.rho += (-dr[j] - qDeltaHalf*qv[j])*qsph->gradRho[j];

         riemann_input.L.v[0] += ( dr[j]*psph->gradVx[j]);
         riemann_input.R.v[0] += (-dr[j]*qsph->gradVx[j]); 
                                                   
         riemann_input.L.v[1] += ( dr[j]*psph->gradVy[j]);
         riemann_input.R.v[1] += (-dr[j]*qsph->gradVy[j]);
                                                   
         riemann_input.L.v[2] += ( dr[j]*psph->gradVz[j]);
         riemann_input.R.v[2] += (-dr[j]*qsph->gradVz[j]);

         riemann_input.L.p += ( dr[j] - pDeltaHalf*pv[j])*psph->gradP[j];
         riemann_input.R.p += (-dr[j] - qDeltaHalf*qv[j])*qsph->gradP[j];
      }
//      printf("2) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("2) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

      for (j=0; j<3; j++){ // Forward extrapolation of velocity
         riemann_input.L.v[j] -= pv[j]*pdivv + psph->gradP[j]/pDensity*pDeltaHalf;
         riemann_input.R.v[j] -= qv[j]*qdivv + qsph->gradP[j]/pkdDensity(pkd,q)*qDeltaHalf;
      }

      riemann_input.L.rho -= pDensity*pdivv;
      riemann_input.R.rho -= pkdDensity(pkd,q)*qdivv;
      riemann_input.L.p -= pkd->param.dConstGamma*psph->P*pdivv;
      riemann_input.R.p -= pkd->param.dConstGamma*qsph->P*qdivv;
//      printf("3) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("3) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);
//      printf("p - gradRho %e %e %e \n", psph->gradRho[0], psph->gradRho[1], psph->gradRho[2]);
//      printf("q - gradRho %e %e %e \n", qsph->gradRho[0], qsph->gradRho[1], qsph->gradRho[2]);


       // It seems that in hydra_core_meshless.h they do not rotate any velocity, they just pass the face_unit
       // to the Riemann solver and let it do its stuff..

       
       // It uses as input parameters the density, pressure and velocities at both ends. They are assumed to be given in comoving coordinates, then they are
       // converted into physical units inside Riemann_Solver. 
       
       // IA: DEBUG: Tests for the riemann solver extracted from Toro (10.1007/b79761)
       // Test 1
//       riemann_input.L.rho = 1.0; riemann_input.L.p = 1.0; riemann_input.L.v[0] = 0.0;
//       riemann_input.L.rho = 0.125; riemann_input.L.p = 0.1; riemann_input.L.v[0] = 0.0;

       Riemann_solver(riemann_input, &riemann_output, face_unit, /*double press_tot_limiter TODO For now, just p>0: */ 0.0);

       // IA: DEBUG
//       printf("Riemann_output: rho %e \tp %e \tv %e %e %e \n", riemann_output.Fluxes.rho, riemann_output.Fluxes.p, riemann_output.Fluxes.v[0], riemann_output.Fluxes.v[1], riemann_output.Fluxes.v[2]);
//       abort();

       // Check for NAN fluxes
       if (riemann_output.Fluxes.rho!=riemann_output.Fluxes.rho) abort();
       if (riemann_output.Fluxes.p!=riemann_output.Fluxes.p) abort(); 


       // Now we de-boost the fluxes following Eq. A8 Hopkins 2015 (From hydra_core_meshless.h):
       /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame'
        which we do following Pakmor et al. 2011 */
       for(j=0;j<3;j++)
       {
           riemann_output.Fluxes.p += vFrame[j] * riemann_output.Fluxes.v[j];
           riemann_output.Fluxes.p += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
       }

       // IA: Now we just multiply by the face area
       riemann_output.Fluxes.p *= modApq;
       riemann_output.Fluxes.rho *= modApq;
       for (j=0;j<3;j++) {riemann_output.Fluxes.v[j] *= modApq; riemann_output.Fluxes.v[j] += vFrame[j]*riemann_output.Fluxes.rho;  } // De-boost (modApq included in Fluxes.rho)

       // IA: Own contribution is always added
       psph->Frho += riemann_output.Fluxes.rho; 
       psph->Fene += riemann_output.Fluxes.p; 
       for(j=0;j<3;j++) {psph->Fmom[j] += riemann_output.Fluxes.v[j];} 

       /* IA: If the other particle is not active, we just add the proportional part of the flux in the given timestep. For example, if p is in rung 2 and q in rung 1,
        * then we only add half of the fluxes in this step, as the we need two timesteps of particle p to be syncronized with q. */
       double dtFracDueToRungDiff; 
       dtFracDueToRungDiff = 1./(1<<(p->uRung - q->uRung));
       

       if (!pkdIsActive(pkd,q)){
          qsph->Frho -= dtFracDueToRungDiff*riemann_output.Fluxes.rho;
          qsph->Fene -= dtFracDueToRungDiff*riemann_output.Fluxes.p;
          for(j=0;j<3;j++){ qsph->Fmom[j] -= dtFracDueToRungDiff*riemann_output.Fluxes.v[j]; }
       }else{
          if (2.*qh < rpq) {  // q is active but p is not in its neighbors list
             qsph->Frho -= riemann_output.Fluxes.rho;
             qsph->Fene -= riemann_output.Fluxes.p;
             for(j=0;j<3;j++){ qsph->Fmom[j] -= riemann_output.Fluxes.v[j]; }
          }
       }




                
       // From hydra_evaluate.h //TODO: Understand this mass_holder thingy... Apart from that, is only updating the Dt structs, which I do not need to do
/*
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
*/
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

    } // IA: End of loop over neighbors


    // IA: Mass change limiter
    
    double dmass_limiter = 0.1*pkdMass(pkd,p), dmass_holder = psph->Frho * ( smf->dDelta/(1<<p->uRung) );

    if (fabs(dmass_holder) > dmass_limiter) {
       printf("Limiting! \n");
       psph->Frho *= dmass_limiter / fabs(dmass_holder);
    }
    


    // IA: Timestep criteria based on the hydro accelerations
    double a[3], acc;
    double cfl = 0.001, dtAcc;
    
    for (j=0;j<3;j++) { a[j] = (pkdVel(pkd,p)[j]*psph->Frho + psph->Fmom[j])/pkdMass(pkd,p); }
    acc = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    dtAcc = sqrt(2*cfl*cfl*fBall/acc);
    if (dtAcc < dtEst) dtEst = dtAcc;
    

//    dtEst = 1.e-7; //IA FIXME forced same rungs for all particles
    psph->uNewRung = pkdDtToRung(dtEst,smf->dDelta,MAX_RUNG);


    //IA: TODO FIXME Now this is done temporarly on pkdHydroStep
    //if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung; 

}
