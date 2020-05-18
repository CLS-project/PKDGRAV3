/* File added by Isaac Alonso for computing
 * the hydrodynamical part using a mesh-free
 * method, following the work of Hopkins 2015
 */

#include "pkd.h"
#include "smoothfcn.h"
#include "hydro.h" 
#include "riemann.h"
#include <stdio.h>

//IA: Ref https://pysph.readthedocs.io/en/latest/reference/kernels.html
double cubicSplineKernel(double r, double h) {
   double q;
   q = r/h;
   if (q<1.0){
      return M_1_PI/(h*h*h)*( 1. - 1.5*q*q*(1.-0.5*q) );  
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
      printf("Singular matrix!\n");
      printf("XX %e \nXY %e \t YY %e \nXZ %e \t YZ %e \t ZZ %e \n", E[XX], E[XY], E[YY], E[XZ], E[YZ], E[ZZ]);
      abort();
      B[XX] = 0.;
      B[YY] = 0.;
      B[ZZ] = 0.;
      B[XY] = 0.;
      B[XZ] = 0.;
      B[YZ] = 0.;
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
//      for (i=0;i<6;i++) { psph->B[i] = 0.0; }
//      p->uNewRung = 0;
//      for (i=0;i<3;i++) { 
//         psph->Fmom[i] = 0.0;
//	   pkdAccel(pkd,p)[i] = 0;
//	   pkdAccel(pkd,p)[i] = 0;
//	   pkdAccel(pkd,p)[i] = 0;
//	}
	

//      psph->Frho = 0.0;
//      psph->Fene = 0.0;
//      psph->uNewRung = 0;
//      for (i=0;i<3;i++) { 
//         psph->Fmom[i] = 0.0;
//	}
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
//      for (i=0;i<6;i++) { psph->B[i] = 0.0; }
//      p->uNewRung = 0;
//      for (i=0;i<3;i++) { 
//         psph->Fmom[i] = 0.0;
//	   pkdAccel(pkd,p)[i] = 0;
//	   pkdAccel(pkd,p)[i] = 0;
//	   pkdAccel(pkd,p)[i] = 0;
//	}
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
   
    
    float *p1mass = pkdField(p1,pkd->oMass);
    float *p2mass = pkdField(p2,pkd->oMass);
    *p1mass += *p2mass;

    psph1->mom[0] += psph2->mom[0];
    psph1->mom[1] += psph2->mom[1];
    psph1->mom[2] += psph2->mom[2];
    psph1->E += psph2->E;
//    if (((PARTICLE *) p2)->uNewRung > ((PARTICLE *) p1)->uNewRung) 
//       ((PARTICLE *) p1)->uNewRung = ((PARTICLE *) p2)->uNewRung;
    }

void initHydroFluxes(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    int i;
//    if (pkdIsActive(pkd,p)) {
	SPHFIELDS *psph = pkdSph(pkd,p);
//      psph->Frho = 0.0;
//      psph->Fene = 0.0;
//      psph->uNewRung = 0;
//      for (i=0;i<3;i++) { 
//         psph->Fmom[i] = 0.0;
//	}
    }

void initHydroStep(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;

}

void combHydroStep(void *vpkd, void *p1,void *p2) {
    PKD pkd = (PKD) vpkd;
    assert(!pkd->bNoParticleOrder);
    SPHFIELDS *psph1 = pkdSph(pkd,p1), *psph2 = pkdSph(pkd,p2);

       if (((PARTICLE *) p2)->uNewRung > ((PARTICLE *) p1)->uNewRung)
           ((PARTICLE *) p1)->uNewRung = ((PARTICLE *) p2)->uNewRung;
}


/* IA: If this works as I think it works, I should put to zero all the 
 * conserved quantities, which will be updated during the hydro loop.
 * Then those will be merged with the actual particle information inside
 * combThirdHydroLoop
 * However, the comb does not seem to be called... But everything seems
 * to work just fine.. (chuckles) I'm in danger */
void initHydroFluxesCached(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    SPHFIELDS *psph = pkdSph(pkd,p);
    int i;

    float *pmass = pkdField(p,pkd->oMass);
    *pmass = 0.0;
    psph->mom[0] = 0.;
    psph->mom[1] = 0.;
    psph->mom[2] = 0.;
    psph->E = 0.;

//    if (pkdIsActive(pkd,p)) {
//      psph->Frho = 0.0;
//      psph->Fene = 0.0;
//      psph->uNewRung = 0;
//      for (i=0;i<3;i++) { 
//         psph->Fmom[i] = 0.0;
//	}
    }

void initHydroGradients(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    SPHFIELDS *psph = pkdSph(pkd,p);
    int j;
//    for (j=0; j<3;j++){
//      psph->gradRho[j] = 0.0;
//      psph->gradVx[j] = 0.0;
//      psph->gradVy[j] = 0.0;
//      psph->gradVz[j] = 0.0;
//      psph->gradP[j] = 0.0;
//      }
    }



void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph;
    double ph, qh, rpq, hpq, dx,dy,dz,Wpq, c;
    int i;

    /* Particle p data */
    psph = pkdSph(pkd,p);
    ph = fBall;

#ifndef FIXED_NSMOOTH_STRICT
    /* IA: Compute the \omega(x_i) normalization factor */
    psph->omega = 0.0;
    ph = pkdBall(pkd,p);
    for (i=0; i<nSmooth; ++i){
       q = nnList[i].pPart;

       rpq = sqrt(nnList[i].fDist2);
       hpq = ph; 

       psph->omega += cubicSplineKernel(rpq, hpq);
    }
    
    /* IA: If we are using a iterative procedure for computing the smoothing length, then:
     *    - Particles marked are those which still needs iterations
     *    - Particles not marked are those with a correct smoothing length
     */
    if (pkd->param.bIterativeSmoothingLength && p->bMarked){
#ifndef FIXED_NSMOOTH_RELAXED
       c = 4.*M_PI/3. * psph->omega *fBall*fBall*fBall*8.; 
#else
       c = nSmooth;
#endif
       if (fabs(c-pkd->param.nSmooth) < pkd->param.dNeighborsStd0){
          p->bMarked = 0;
       }else{
          float newBall;
          newBall = fBall * pow(  pkd->param.nSmooth/c  ,0.3333333333);
       //   if (nSmooth <= 1) newBall *= 2.*fBall;
       
          pkdSetBall(pkd,p, 0.5*(newBall+fBall));
          psph->fLastBall = fBall;
          
#else // FIXED_NSMOOTH_STRICT
    double minR2, lastMin;
    if (pkd->param.bIterativeSmoothingLength && p->bMarked){
       c = nSmooth;
       // Check if we have enough neighbors, otherwise increse fBall
       if (c <  pkd->param.nSmooth){
          pkdSetBall(pkd,p,fBall * pow(  1.2*pkd->param.nSmooth/c  ,0.3333333333));
       }else if (c >= pkd->param.nSmooth){
          // Now we look for the distance to the nSmooth-th neighbor
          minR2 = HUGE_VAL;
          lastMin = 0.;

          for (int n; n<pkd->param.nSmooth-2; n++){
             minR2 = HUGE_VAL;
             for (i=0; i<nSmooth; i++){
                if (nnList[i].fDist2 < minR2)
                   if (nnList[i].fDist2 > lastMin)
                      minR2 = nnList[i].fDist2;
             }
             lastMin = minR2;

          }

          pkdSetBall(pkd,p, 0.5*sqrt(lastMin));

          p->bMarked=0;
#endif
       }
    }

#ifdef FIXED_NSMOOTH_STRICT
    psph->omega = 0.0;
    ph = pkdBall(pkd,p);
    for (i=0; i<nSmooth; ++i){
       q = nnList[i].pPart;

       rpq = sqrt(nnList[i].fDist2);
       hpq = ph; 

       psph->omega += cubicSplineKernel(rpq, hpq);
    }
#endif

    /* IA: We compute the density making use of Eq. 27 Hopkins 2015 */
    pkdSetDensity(pkd,p, pkdMass(pkd,p)*psph->omega);
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

    psph->nLastNeighs = nSmooth;
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

      //printf("rpq %e hpq %e omega %e Wpq %e \n", rpq, hpq, psph->omega, Wpq);
    }

    /* IA: Normalize the matrix */
    for (i=0; i<6;++i){
       E[i] /= psph->omega; 
    }


    

    /* IA: END of E matrix computation */
    //printf("nSmooth %d fBall %e E_q [XX] %e \t [XY] %e \t [XZ] %e \n \t \t \t [YY] %e \t [YZ] %e \n \t\t\t \t \t \t [ZZ] %e \n", nSmooth, fBall, E[XX], E[XY], E[XZ], E[YY], E[YZ], E[ZZ]);

    /* IA: Now, we need to do the inverse */
    inverseMatrix(E, psph->B);
    

    /* IA: Computation of the condition number */
    double modE = 0.;
    modE += E[XX]*E[XX]; 
    modE += E[YY]*E[YY]; 
    modE += E[ZZ]*E[ZZ]; 
    modE += 2.*E[XY]*E[XY]; 
    modE += 2.*E[XZ]*E[XZ]; 
    modE += 2.*E[YZ]*E[YZ]; 

    double modB = 0.;
    modB += psph->B[XX]*psph->B[XX]; 
    modB += psph->B[YY]*psph->B[YY]; 
    modB += psph->B[ZZ]*psph->B[ZZ]; 
    modB += 2.*psph->B[XY]*psph->B[XY]; 
    modB += 2.*psph->B[XZ]*psph->B[XZ]; 
    modB += 2.*psph->B[YZ]*psph->B[YZ]; 

    psph->Ncond = sqrt(modB*modE)/3.;



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
    for (j=0; j<3;j++){
       psph->gradRho[j] = 0.0;
       psph->gradVx[j] = 0.0;
       psph->gradVy[j] = 0.0;
       psph->gradVz[j] = 0.0;
       psph->gradP[j] = 0.0;
#if defined(MAKE_GLASS) || defined(REGULARIZE_MESH)
       psph->cellCM[j] = 0.0;
#endif
    }
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
       hpq = ph; 

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
#if defined(MAKE_GLASS) || defined(REGULARIZE_MESH)
       for (j=0; j<3; j++) {psph->cellCM[j] += psi*psiTilde_p[j]; }
#endif
    }
#if defined(MAKE_GLASS) || defined(REGULARIZE_MESH)
    double CMfactor = - fBall*fBall*M_PI*fBall*fBall*fBall * psph->omega / 3.;
    for (j=0; j<3; j++) { psph->cellCM[j] *= CMfactor; }
#endif

     /* IA: Now we can limit them */
    
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
#if defined(LIMITER_BARTH) || defined(LIMITER_CONDBARTH)
    for (i=0; i<nSmooth;++i){
	 dx = -nnList[i].dx; //Vector from p to q
	 dy = -nnList[i].dy;
	 dz = -nnList[i].dz;
       if (dx==0 && dy==0 && dz==0) continue;
       q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);

       // TODO: The differences could be computed outside of this loop
#ifdef LIMITER_BARTH
       BarthJespersenLimiter(&limRho, psph->gradRho, rho_max-pkdDensity(pkd,p), rho_min-pkdDensity(pkd,p), dx, dy, dz);
       BarthJespersenLimiter(&limVx, psph->gradVx, vx_max-pv[0], vx_min-pv[0], dx, dy, dz);
       BarthJespersenLimiter(&limVy, psph->gradVy, vy_max-pv[1], vy_min-pv[1], dx, dy, dz);
       BarthJespersenLimiter(&limVz, psph->gradVz, vz_max-pv[2], vz_min-pv[2], dx, dy, dz);
       BarthJespersenLimiter(&limP, psph->gradP, p_max-psph->P, p_min-psph->P, dx, dy, dz);
#endif

#ifdef LIMITER_CONDBARTH
       ConditionedBarthJespersenLimiter(&limRho, psph->gradRho, rho_max-pkdDensity(pkd,p), rho_min-pkdDensity(pkd,p), dx, dy, dz, 100., psph->Ncond);
       ConditionedBarthJespersenLimiter(&limVx, psph->gradVx, vx_max-pv[0], vx_min-pv[0], dx, dy, dz, 100., psph->Ncond);
       ConditionedBarthJespersenLimiter(&limVy, psph->gradVy, vy_max-pv[1], vy_min-pv[1], dx, dy, dz, 100., psph->Ncond);
       ConditionedBarthJespersenLimiter(&limVz, psph->gradVz, vz_max-pv[2], vz_min-pv[2], dx, dy, dz, 100., psph->Ncond);
       ConditionedBarthJespersenLimiter(&limP, psph->gradP, p_max-psph->P, p_min-psph->P, dx, dy, dz, 100., psph->Ncond);
#endif
    }
#endif

    for (j=0; j<3; j++){
       psph->gradRho[j] *= limRho;
       psph->gradVx[j] *= limVx;
       psph->gradVy[j] *= limVy;
       psph->gradVz[j] *= limVz;
       psph->gradP[j] *= limP;
    }
    /* END OF LIMITER */
    /*
    if (p->iOrder==55) {
      printf("Limiter: rho %e \t v %e %e %e \t p %e \n", limRho, limVx, limVy, limVz, limP);
      printf("x %e y %e z %e \n", pkdPos(pkd,p,0), pkdPos(pkd,p,1), pkdPos(pkd,p,2));
      printf("gradRho %e \n", psph->gradRho[0]);
      
      for (i=0; i<nSmooth;++i){
       q = nnList[i].pPart;
         printf("\t p %" PRId64 " \t q %" PRId64 " \n", p->iOrder, q->iOrder);
      }
    }
    */
    }


void BarthJespersenLimiter(double* limVar, double* gradVar, double var_max, double var_min, double dx, double dy, double dz){
    double diff, lim;

    diff = (gradVar[0]*dx + gradVar[1]*dy + gradVar[2]*dz);
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
    if (lim < 0.) lim = 0.;
    if (lim < (*limVar)) { *limVar = lim;} 
    // FIXME IA: Option to avoid extrapolation or limiter
//    *limVar = 1.0;
//    *limVar = 0.0;
}

// IA: In this version we take into account the condition number, which give us an idea about how 'well aligned' are the particles
void ConditionedBarthJespersenLimiter(double* limVar, double* gradVar, double var_max, double var_min, double dx, double dy, double dz, double Ncrit, double Ncond){
    double diff, lim, beta;

    diff = Ncrit/Ncond;
    diff = diff < 1. ? diff : 1.;
    diff *= 2.;
    beta = (1. < diff) ? diff : 1.; 


    diff = (gradVar[0]*dx + gradVar[1]*dy + gradVar[2]*dz);
    if (var_min > 0) { var_min=0; } //IA: Can happen due to machine precision
    if (var_max < 0) { var_max=0; } //IA: Can happen due to machine precision
    if (diff > 0.) {
       lim = var_max/diff;
    }else if (diff < 0.){
       lim = var_min/diff;
    }else{
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



/* IA: This routine will extrapolate the primitives to the 'faces' and solve the 1D riemann problem. 
 * For now, the 1D riemann flux will be computed TWICE for a given face, one for each adjacent particles */ 
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
   //IA TODO Clean unused variables!
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    float *pmass, *qmass; 
    double minDt;
    double pv[3], qv[3], vFrame[3];
    double ph,qh, hpq, modApq, rpq, dx,dy,dz,Wpq,psi_p,psi_q,pDensity,pDeltaHalf, qDeltaHalf, dt2, dtEst, vsig_pq, dvDotd, dvDotdr;
    double vxdiff, vydiff, vzdiff, pdivv, qdivv, pdiffOverRhop, pdiffOverRhoq, psi;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3], dr[3]; 
    struct Input_vec_Riemann riemann_input;
    struct Riemann_outputs riemann_output; 
    uint8_t uNewRung;
    int i,j;

    psph = pkdSph(pkd, p);   
    ph = pkdBall(pkd, p); 

    pDensity = pkdDensity(pkd,p);

    dtEst = HUGE_VAL;

    for (i=0;i<nSmooth;++i){

	 q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);
       qh = pkdBall(pkd,q); 


       /* IA: We update the conservatives variables taking the minimum timestep between the particles, as in AREPO */
       if (!pkdIsActive(pkd,q)) { // If q is not active we now that p has the smallest dt 
          minDt = smf->dDelta/(1<<p->uRung) ;
       } else { // Otherwise we need to explicitly check
          if (p->uRung > q->uRung) {
             minDt = smf->dDelta/(1<<p->uRung) ;
          }else{
             minDt = smf->dDelta/(1<<q->uRung) ;
          }
       }

       if (smf->dDelta > 0) {
          pDeltaHalf = smf->dTime - psph->lastUpdateTime + 0.5*smf->dDelta/(1<<p->uRung);
          qDeltaHalf = smf->dTime - qsph->lastUpdateTime + 0.5*smf->dDelta/(1<<q->uRung);
       }else{
          qDeltaHalf = 0.0; // For the initialization step we do not extrapolate because we dont have a reliable dDelta
          pDeltaHalf = 0.0;
       }
       if(pkd->csm->val.bComove)
       {
          qDeltaHalf /= smf->a;
          pDeltaHalf /= smf->a;
       }


//       printf("pDeltaHalf %e qDeltaHalf %e \n", pDeltaHalf, qDeltaHalf);

	 dx = nnList[i].dx;
	 dy = nnList[i].dy;
	 dz = nnList[i].dz;
#ifdef FORCE_1D
       if (dz!=0) continue;
       if (dy!=0) continue;
#endif
#ifdef FORCE_2D
       if (dz!=0) continue;
#endif

       /* IA: in the nnList there is a 'copy' of the own particle, which we can omit as there are no fluxes
        * to be computed here */
       if (dx==0 && dy==0 && dz==0) continue;

       // Face where the riemann problem will be solved
       hpq = ph;
       rpq = sqrt(nnList[i].fDist2);
       if (qh/0.50 < rpq) continue;

       Wpq = cubicSplineKernel(rpq, hpq); 
       if (Wpq==0.0){/*printf("hpq %e rpq %e \t %e \n", hpq, rpq, rpq/hpq); */continue; }

       // \tilde{\psi}_j (x_i)
       psi = -cubicSplineKernel(rpq, ph)/psph->omega;   
       psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
       psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
       psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

       // \tilde{\psi}_i (x_j)
       psi = cubicSplineKernel(rpq, qh)/qsph->omega; // IA: minus because we are 'looking from' the other particle, thus -dr
       psiTilde_q[0] = (qsph->B[XX]*dx + qsph->B[XY]*dy + qsph->B[XZ]*dz)*psi;
       psiTilde_q[1] = (qsph->B[XY]*dx + qsph->B[YY]*dy + qsph->B[YZ]*dz)*psi;
       psiTilde_q[2] = (qsph->B[XZ]*dx + qsph->B[YZ]*dy + qsph->B[ZZ]*dz)*psi;

       modApq = 0.0;
       for (j=0; j<3; j++){
          Apq[j] = psiTilde_p[j]/psph->omega - psiTilde_q[j]/qsph->omega;
          modApq += Apq[j]*Apq[j];
       }
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
      
      // IA: Placing this here solved the convergence problem for the comoving soundwaves.
      //   This problem may be caused because we do not use the time extrapolated cell-centered states in 
      //   this limiter
      //
      //   TODO: One possible solution is extrapolate the variables at PrimVar (or equivalent) and use them
      //   directly here. This also would reduce the computational cost, as we are doing the exact same
      //   temporal extrapolation ~2 nSmooth times instead of just once.
      /*
      genericPairwiseLimiter(pkdDensity(pkd,p), pkdDensity(pkd,q), &riemann_input.L.rho, &riemann_input.R.rho);
      genericPairwiseLimiter(psph->P, qsph->P, &riemann_input.L.p, &riemann_input.R.p);
      genericPairwiseLimiter(pv[0], qv[0], &riemann_input.L.v[0], &riemann_input.R.v[0]);
      genericPairwiseLimiter(pv[1], qv[1], &riemann_input.L.v[1], &riemann_input.R.v[1]);
      genericPairwiseLimiter(pv[2], qv[2], &riemann_input.L.v[2], &riemann_input.R.v[2]);
      */

      
      double temp;
      if(pkd->csm->val.bComove){
         
         for (j=0;j<3;j++){
            temp = smf->H * pDeltaHalf * smf->a * pv[j];
            riemann_input.L.v[j] -= temp;
            vFrame[j] -= 0.5*temp;

            temp = smf->H * qDeltaHalf * smf->a * qv[j];
            riemann_input.R.v[j] -= temp;
            vFrame[j] -= 0.5*temp;
         }
         
         riemann_input.L.p -= 3. * smf->H * pDeltaHalf * smf->a * (pkd->param.dConstGamma - 1.) * psph->P;
         riemann_input.R.p -= 3. * smf->H * qDeltaHalf * smf->a * (pkd->param.dConstGamma - 1.) * qsph->P;
         
      }




      for (j=0; j<3; j++){ // Forward extrapolation of velocity
         temp = pv[j]*pdivv + psph->gradP[j]/pDensity*pDeltaHalf;
         riemann_input.L.v[j] -= temp;
         vFrame[j] -= 0.5*temp;

         temp = qv[j]*qdivv + qsph->gradP[j]/pkdDensity(pkd,q)*qDeltaHalf;
         riemann_input.R.v[j] -= temp;
         vFrame[j] -= 0.5*temp;
      }

      for (j=0; j<3; j++){
         temp = psph->lastAcc[j]*pDeltaHalf;
         riemann_input.L.v[j] += temp;
         vFrame[j] += 0.5*temp;

         temp = qsph->lastAcc[j]*qDeltaHalf;
         riemann_input.R.v[j] += temp;
         vFrame[j] += 0.5*temp;
      }

      riemann_input.L.rho -= pDensity*pdivv;
      riemann_input.R.rho -= pkdDensity(pkd,q)*qdivv;
      riemann_input.L.p -= pkd->param.dConstGamma*psph->P*pdivv;
      riemann_input.R.p -= pkd->param.dConstGamma*qsph->P*qdivv;

      genericPairwiseLimiter(pkdDensity(pkd,p), pkdDensity(pkd,q), &riemann_input.L.rho, &riemann_input.R.rho);
      genericPairwiseLimiter(psph->P, qsph->P, &riemann_input.L.p, &riemann_input.R.p);
      genericPairwiseLimiter(pv[0], qv[0], &riemann_input.L.v[0], &riemann_input.R.v[0]);
      genericPairwiseLimiter(pv[1], qv[1], &riemann_input.L.v[1], &riemann_input.R.v[1]);
      genericPairwiseLimiter(pv[2], qv[2], &riemann_input.L.v[2], &riemann_input.R.v[2]);
       
       // IA: DEBUG: Tests for the riemann solver extracted from Toro (10.1007/b79761)
       // Test 1
//       riemann_input.L.rho = 1.0; riemann_input.L.p = 1.0; riemann_input.L.v[0] = 0.0;
//       riemann_input.L.rho = 0.125; riemann_input.L.p = 0.1; riemann_input.L.v[0] = 0.0;

       if (riemann_input.L.rho < 0) {riemann_input.L.rho = pkdDensity(pkd,p); /* printf("WARNING, L.rho < 0 : using first-order scheme \n");*/ }
       if (riemann_input.R.rho < 0) {riemann_input.R.rho = pkdDensity(pkd,q); /* printf("WARNING, R.rho < 0 : using first-order scheme \n");*/ }
       if (riemann_input.L.p < 0) {riemann_input.L.p = psph->P;  /*  printf("WARNING, L.p < 0 : using first-order scheme \n");*/ }
       if (riemann_input.R.p < 0) {riemann_input.R.p = qsph->P;  /*  printf("WARNING, R.p < 0 : using first-order scheme \n");*/ }


       Riemann_solver(pkd, riemann_input, &riemann_output, face_unit, /*double press_tot_limiter TODO For now, just p>0: */ 0.0);
      
       // IA: MFM
#ifdef USE_MFM
       if (riemann_output.Fluxes.rho != 0){
          printf("Frho %e \n",riemann_output.Fluxes.rho);
          abort();
       }
        riemann_output.Fluxes.rho = 0.;
        riemann_output.Fluxes.p = riemann_output.P_M * riemann_output.S_M;
        for(j=0;j<3;j++)
            riemann_output.Fluxes.v[j] = riemann_output.P_M * face_unit[j];
#endif       
       // IA: End MFM

       // Force 2D
#ifdef FORCE_1D
       riemann_output.Fluxes.v[2] = 0.;
       riemann_output.Fluxes.v[1] = 0.;
#endif
#ifdef FORCE_2D
       riemann_output.Fluxes.v[2] = 0.;
#endif
      


       // IA: DEBUG
//       abort();

       // Check for NAN fluxes
       if (riemann_output.Fluxes.rho!=riemann_output.Fluxes.rho) riemann_output.Fluxes.rho = 0.;//abort();
       if (riemann_output.Fluxes.p!=riemann_output.Fluxes.p) riemann_output.Fluxes.p = 0.;//abort(); 


       if(pkd->csm->val.bComove){ 
          minDt /= smf->a; // 1/a term before \nabla
          /*
          modApq *= smf->a*smf->a; 
          for (j=0;j<3;j++){
             vFrame[j] /= smf->a;
          }
          */
       }


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


       if (smf->dDelta > 0){
       pmass = pkdField(p,pkd->oMass);
       qmass = pkdField(q,pkd->oMass);

       /*
       if(pkd->csm->val.bComove){ 
          for (j=0;j<3;j++) riemann_output.Fluxes.v[j] *= smf->a;
          riemann_output.Fluxes.p *= smf->a*smf->a;
       }
       */

       
       if ( (2.*qh < rpq) | !pkdIsActive(pkd,q)){  
            *qmass += minDt * riemann_output.Fluxes.rho ;

            qsph->mom[0] += minDt * riemann_output.Fluxes.v[0] ;
            qsph->mom[1] += minDt * riemann_output.Fluxes.v[1] ;
            qsph->mom[2] += minDt * riemann_output.Fluxes.v[2] ;

            qsph->E += minDt * riemann_output.Fluxes.p;

            qsph->Uint += minDt * ( riemann_output.Fluxes.p - riemann_output.Fluxes.v[0]*qsph->vPred[0] 
                                                            - riemann_output.Fluxes.v[1]*qsph->vPred[1]
                                                            - riemann_output.Fluxes.v[2]*qsph->vPred[2]
                                  + 0.5*(qsph->vPred[0]*qsph->vPred[0] + qsph->vPred[1]*qsph->vPred[1] + qsph->vPred[2]*qsph->vPred[2]) * riemann_output.Fluxes.rho );

            qsph->drDotFrho[0] += minDt * riemann_output.Fluxes.rho * dx;
            qsph->drDotFrho[1] += minDt * riemann_output.Fluxes.rho * dy;
            qsph->drDotFrho[2] += minDt * riemann_output.Fluxes.rho * dz;
       } 
            *pmass -= minDt * riemann_output.Fluxes.rho ;

            psph->mom[0] -= minDt * riemann_output.Fluxes.v[0] ;
            psph->mom[1] -= minDt * riemann_output.Fluxes.v[1] ;
            psph->mom[2] -= minDt * riemann_output.Fluxes.v[2] ;

            psph->E -= minDt * riemann_output.Fluxes.p;

            psph->Uint -= minDt * ( riemann_output.Fluxes.p - riemann_output.Fluxes.v[0]*psph->vPred[0] 
                                                            - riemann_output.Fluxes.v[1]*psph->vPred[1]
                                                            - riemann_output.Fluxes.v[2]*psph->vPred[2]
                                  + 0.5*(psph->vPred[0]*psph->vPred[0] + psph->vPred[1]*psph->vPred[1] + psph->vPred[2]*psph->vPred[2]) * riemann_output.Fluxes.rho );

            psph->drDotFrho[0] += minDt * riemann_output.Fluxes.rho * dx;
            psph->drDotFrho[1] += minDt * riemann_output.Fluxes.rho * dy;
            psph->drDotFrho[2] += minDt * riemann_output.Fluxes.rho * dz;

       }
/*IA:  Old fluxes update (see 15/04/19 ) 
 * TODO: This is not needed for the update of the conserved variables. Instead,
 * it is now only used for the acceleleration criteria. Memory-wise, this can be
 * substantially improved
 */
       // IA: Own contribution is always added
       psph->Frho += riemann_output.Fluxes.rho; 
       psph->Fene += riemann_output.Fluxes.p; 
       for(j=0;j<3;j++) {psph->Fmom[j] += riemann_output.Fluxes.v[j];} 

       if (!pkdIsActive(pkd,q)){
          qsph->Frho -= riemann_output.Fluxes.rho;
          qsph->Fene -= riemann_output.Fluxes.p;
          for(j=0;j<3;j++){ qsph->Fmom[j] -= riemann_output.Fluxes.v[j]; }
       }
    } // IA: End of loop over neighbors





}



/* IA: Compute the hydrodynamical time step of this particle, based on two criterias: acceleration and signal velocity */
void hydroStep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double dt2, dtEst, vsig_pq, dvDotdr, dx, dy, dz;
    uint8_t uNewRung;
    int i,j;

    psph = pkdSph(pkd, p);   

    dtEst = HUGE_VAL;

    /* IA: Signal velocity criteria */
    for (i=0;i<nSmooth;++i){

	 dx = nnList[i].dx;
	 dy = nnList[i].dy;
	 dz = nnList[i].dz;

       if (dx==0 && dy==0 && dz==0) continue;

	 q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);

       // From Eqs 24,25 Hopkins 2015, to limit deltaT
       dvDotdr = (dx*(qsph->vPred[0] - psph->vPred[0]) + 
                  dy*(qsph->vPred[1] - psph->vPred[1]) +
                  dz*(qsph->vPred[2] - psph->vPred[2]));

       if (dvDotdr < 0) {
          vsig_pq = psph->c + qsph->c - dvDotdr/sqrt(nnList[i].fDist2);
       }else{
          vsig_pq = psph->c + qsph->c;
       }

       dt2 = 2.*smf->dEtaCourant * fBall /vsig_pq;	
	 if (dt2 < dtEst) dtEst=dt2;

    }

#ifdef HERNQUIST_POTENTIAL
    // IA: Timestep criteria based on the Hernsquist potential
    const double const_reduced_hubble_cgs = 3.2407789e-18;
    //const double H0 = 0.704 * const_reduced_hubble_cgs * pkd->param.dSecUnit;
    const double H0 = 70.4/ pkd->param.dKmPerSecUnit * ( pkd->param.dKpcUnit / 1e3);

    const double concentration = 9.0;
    const double M200 = 135.28423603962767; //137.0 ; // / pkd->param.dMsolUnit;
    const double V200 = cbrt(M200*H0);
    //const double R200 = V200/(H0);
    const double R200 = cbrt(M200/(100.*H0*H0));
    const double RS = R200 / concentration;

    const double al = RS * sqrt(2. * (log(1. + concentration) -
                                    concentration / (1. + concentration)));

    const double mass = M200;(1.-0.041);

  /* Calculate the relative potential with respect to the centre of the
   * potential */
  dx = pkdPos(pkd,p,0); //- potential->x[0];
  dy = pkdPos(pkd,p,1); //- potential->x[1];
  dz = pkdPos(pkd,p,2); //- potential->x[2];

  /* calculate the radius  */
    const double epsilon =  0.2/pkd->param.dKpcUnit;
    const double epsilon2 = epsilon*epsilon;
  const float r = sqrtf(dx * dx + dy * dy + dz * dz + epsilon2);
  const float sqrtgm_inv = 1.f / sqrtf(mass);

  /* Calculate the circular orbital period */
  const float period = 2.f * M_PI * sqrtf(r) * al *
                       (1 + r / al) * sqrtgm_inv;

  /* Time-step as a fraction of the circular orbital time */
  double time_step = 0.01 * period;

  if (time_step < dtEst) dtEst = time_step;
#endif





    // IA: Timestep criteria based on the hydro+grav accelerations

    double a[3], acc;
    float* pa = pkdAccel(pkd,p);
    double cfl = smf->dCFLacc, dtAcc;
    
    double pdt = smf->dDelta/(1<<p->uRung) ;
    for (j=0;j<3;j++) { a[j] = pa[j] + (-pkdVel(pkd,p)[j]*psph->Frho + psph->Fmom[j])/pkdMass(pkd,p); }
    acc = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    float h = pkd->param.bDoGravity ? (fBall < pkdSoft(pkd,p) ? fBall : pkdSoft(pkd,p) ) : fBall;
    dtAcc = cfl*sqrt(2*h/acc);
    

    if (dtAcc < dtEst) dtEst = dtAcc;
    uNewRung = pkdDtToRung(dtEst,smf->dDelta,MAX_RUNG);
    if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung; 

    if ( p->uNewRung < pkd->param.dMinDt ) p->uNewRung = pkd->param.dMinDt;
    // IA: Timestep limiter that imposes that I must have a dt which is at most, 
    // four times (i.e., 2 rungs) the smallest dt of my neighbours

    if (smf->dDelta >= 0){
    for (i=0; i<nSmooth; ++i){
        q = nnList[i].pPart;

        //if ( (q->uNewRung - p->uNewRung) > 2) uNewRung = q->uNewRung-1;
        //if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung; 


        if ( (p->uNewRung - q->uNewRung) > 2) q->uNewRung = p->uNewRung-1;

    }
    
    }


}

#define psi1 0.5
#define psi2 0.25
#define SIGN(x) (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0) )
#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y)  ((X) > (Y) ? (X) : (Y))
void genericPairwiseLimiter(double Lstate, double Rstate, double *Lstate_face, double *Rstate_face){
   double phi_max, phi_min, d1, d2, phi_mean, phi_p, phi_m;

   if (Lstate == Rstate){
      *Lstate_face = Lstate;
      *Rstate_face = Rstate;
   }else{

      d1 = psi1*fabs(Lstate - Rstate);
      d2 = psi2*fabs(Lstate - Rstate);

      phi_mean = 0.5*(Lstate+Rstate);

      phi_min = MIN(Lstate, Rstate);
      phi_max = MAX(Lstate, Rstate);

      if (SIGN(phi_min - d1) == SIGN(phi_min) ){
         phi_m = phi_min - d1;
      }else{
         phi_m = phi_min/(1. + d1/fabs(phi_min));
      }

      if (SIGN(phi_max + d1) == SIGN(phi_max) ){
         phi_p = phi_max + d1;
      }else{
         phi_p = phi_max/(1. + d1/fabs(phi_max));
      }

      if (Lstate < Rstate){
         *Lstate_face = MAX(phi_m, MIN(phi_mean+d2, *Lstate_face));
         *Rstate_face = MIN(phi_p, MAX(phi_mean-d2, *Rstate_face));
      }else{
         *Rstate_face = MAX(phi_m, MIN(phi_mean+d2, *Rstate_face));
         *Lstate_face = MIN(phi_p, MAX(phi_mean-d2, *Lstate_face));
      }

   }


}
