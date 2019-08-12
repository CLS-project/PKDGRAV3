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
      printf("Singular matrix!\n");
      printf("XX %e \nXY %e \t YY %e \nXZ %e \t YZ %e \t ZZ %e \n", E[XX], E[XY], E[YY], E[XZ], E[YZ], E[ZZ]);
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
//    *pmass = 0.0;
//    psph->mom[0] = 0.;
//    psph->mom[1] = 0.;
//    psph->mom[2] = 0.;
//    psph->E = 0.;

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

    /* IA: Compute the \omega(x_i) normalization factor */
    psph->omega = 0.0;
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
       c = 4.*M_PI/3. * psph->omega ; 
       //if (fabs(nSmooth-pkd->param.nSmooth) < pkd->param.iNeighborsStd){
       if (fabs(c*fBall*fBall*fBall*8.-pkd->param.nSmooth) < pkd->param.iNeighborsStd){
          p->bMarked = 0;
       }else{
          if (psph->fLastBall == 0.0) { // IA: We have just read the input file. So we can only do one Newton iteration
             float newBall = pow( pkd->param.nSmooth/c, 1./3.)/2.;      
//             printf("p %" PRId64 " omega %e Nngb %e nSmooth %d fBall %e newBall %e NewNngb %e \n", p->iOrder, psph->omega, c*fBall*fBall*fBall*8., nSmooth, fBall, newBall, c*newBall*newBall*newBall*8.);
             pkdSetBall(pkd,p, newBall); 
             psph->fLastBall = fBall;
             psph->nLastNeighs = nSmooth;
          }else{ //IA: We have the last two points, thus we can apply the bisection rule
             if (psph->fLastBall == fBall) {p->bMarked=0; return; } // We got stuck

             int minNeighs = psph->nLastNeighs < nSmooth ? psph->nLastNeighs : nSmooth;
             float minBall = psph->fLastBall < fBall ? psph->fLastBall : fBall;


             int maxNeighs = psph->nLastNeighs > nSmooth ? psph->nLastNeighs : nSmooth;
             float maxBall = psph->fLastBall > fBall ? psph->fLastBall : fBall;

             if (minNeighs > pkd->param.nSmooth){ 
                float newBall = pow( pkd->param.nSmooth/c, 1./3.)/2.;      
                pkdSetBall(pkd,p, newBall); 
                psph->fLastBall = fBall;
                psph->nLastNeighs = nSmooth;
//                pkdSetBall(pkd, p, 0.92*minBall); 
//                psph->fLastBall = minBall;
//                psph->nLastNeighs = minNeighs;
             }else if (maxNeighs < pkd->param.nSmooth) {
                float newBall = pow( pkd->param.nSmooth/c, 1./3.)/2.;      
                pkdSetBall(pkd,p, newBall); 
                psph->fLastBall = fBall;
                psph->nLastNeighs = nSmooth;
//                pkdSetBall(pkd, p, 1.23*maxBall);
//                psph->fLastBall = maxBall;
//                psph->nLastNeighs = maxNeighs;
             }else{
                pkdSetBall(pkd, p, 0.5*(minBall + maxBall));
                psph->fLastBall = 0.0;
                psph->nLastNeighs = minNeighs;
             }
/*
             if ( psph->nLastNeighs > nSmooth ){ // We had more neighbors before
                if (nSmooth > pkd->param.nSmooth) { // But still more than the expected
                   // We keep iterating using a Newton method
                   float newBall = 0.9*fBall; //pow( pkd->param.nSmooth/c, 1./3.) ;      
                   pkdSetBall(pkd,p, newBall ); 
                   psph->fLastBall = fBall;
                   psph->nLastNeighs = nSmooth;
                   
                }else{ // And the target is in between. Apply middle point rule
                   pkdSetBall(pkd,p, (psph->fLastBall + fBall)*0.5 );
                   psph->fLastBall = fBall;
                   psph->nLastNeighs = nSmooth;
                }

             }else{  // We had less neighbors

                if (pkd->param.nSmooth < psph->nLastNeighs) {  // And we still have too much compared to the target
                   // So we can only do a newton iteration
                   float newBall = 0.9*pow( pkd->param.nSmooth/c, 1./3.) ;      
                   pkdSetBall(pkd,p, newBall ); 
                   //psph->fLastBall = fBall;
                   //psph->nLastNeighs = nSmooth;
                }else{ // The target is in between.
                   pkdSetBall(pkd,p, (psph->fLastBall + fBall)*0.5 );
                   psph->fLastBall = fBall;
                   psph->nLastNeighs = nSmooth;
                }

             }
                */
             
          }
       }
    }


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
    //printf("E_q [XX] %e \t [XY] %e \t [XZ] %e \n \t \t \t [YY] %e \t [YZ] %e \n \t\t\t \t \t \t [ZZ] %e \n", E[XX], E[XY], E[XZ], E[YY], E[YZ], E[ZZ]);

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
    for (j=0; j<3;j++){
       psph->gradRho[j] = 0.0;
       psph->gradVx[j] = 0.0;
       psph->gradVy[j] = 0.0;
       psph->gradVz[j] = 0.0;
       psph->gradP[j] = 0.0;
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

    //IA: Test for Gaburov limiter
   // double irho_min, irho_max, ip_max, ip_min;
   // irho_min=  ip_min =  HUGE_VAL;
   // irho_max=  ip_max = -HUGE_VAL;

    pv = pkdVel(pkd,p);
    for (i=0; i<nSmooth;++i){
//       if (nnList[i].dx==0 && nnList[i].dy==0 && nnList[i].dz==0) continue;
       q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);
       qv = pkdVel(pkd,q);

    //IA: Test for Gaburov limiter
      // double irho = pkdDensity(pkd,p) + (psph->gradRho[0]*dx + psph->gradRho[1]*dy + psph->gradRho[2]*dz);
      // if (irho < irho_min) irho_min = irho;
      // if (irho > irho_max) irho_max = irho;
       if (pkdDensity(pkd,q) < rho_min) rho_min = pkdDensity(pkd,q);
       if (pkdDensity(pkd,q) > rho_max) rho_max = pkdDensity(pkd,q);

       if (qv[0] < vx_min) vx_min = qv[0];
       if (qv[0] > vx_max) vx_max = qv[0];

       if (qv[1] < vy_min) vy_min = qv[1];
       if (qv[1] > vy_max) vy_max = qv[1];

       if (qv[2] < vz_min) vz_min = qv[2];
       if (qv[2] > vz_max) vz_max = qv[2];

    //IA: Test for Gaburov limiter
      // double ip = psph->P + (psph->gradP[0]*dx + psph->gradP[1]*dy + psph->gradP[2]*dz);
      // if (ip < ip_min) ip_min = ip;
      // if (ip > ip_max) ip_max = ip;
       if (qsph->P < p_min) p_min = qsph->P;
       if (qsph->P > p_max) p_max = qsph->P;
    }
 
    /*
    //IA: Test for Gaburov limiter
    double a = ( rho_max-pkdDensity(pkd,p) )/ ( irho_max-pkdDensity(pkd,p) );
    double b = ( pkdDensity(pkd,p)-rho_min )/ ( pkdDensity(pkd,p)-irho_min );
    limRho = a<b ? a : b;
    limRho = 1<limRho ? 1 : limRho;

    a = ( p_max-psph->P )/ ( ip_max-psph->P );
    b = ( psph->P-p_min )/ ( psph->P-ip_min );
    limP = a<b ? a : b;
    limP = 1<limP ? 1 : limP;
    */

    limRho= limVx= limVy= limVz= limP = 1.;
    for (i=0; i<nSmooth;++i){
	 dx = -nnList[i].dx; //Vector from p to q
	 dy = -nnList[i].dy;
	 dz = -nnList[i].dz;
       if (dx==0 && dy==0 && dz==0) continue;
       q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);

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

       if (smf->dTime > 0) {
          pDeltaHalf = smf->dTime - psph->lastUpdateTime + minDt*0.5;
          qDeltaHalf = smf->dTime - qsph->lastUpdateTime + minDt*0.5; //smf->dTime - qsph->lastUpdateTime + pDeltaHalf; 
       }else{
          qDeltaHalf = 0.0; // For the initialization step we do not extrapolate because we dont have a reliable dDelta
          pDeltaHalf = 0.0;
       }
//       printf("pDeltaHalf %e qDeltaHalf %e \n", pDeltaHalf, qDeltaHalf);

	 dx = nnList[i].dx;
	 dy = nnList[i].dy;
	 dz = nnList[i].dz;

       /* IA: in the nnList there is a 'copy' of the own particle, which we can omit as there are no fluxes
        * to be computed here */
       if (dx==0 && dy==0 && dz==0) continue;

       // Face where the riemann problem will be solved
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size 
//       hpq = qh > ph ? qh : ph; 
       rpq = sqrt(nnList[i].fDist2);
       if (qh/0.50 < rpq) continue;

       Wpq = cubicSplineKernel(rpq, hpq); 
//       Wpq = 0.5*( cubicSplineKernel(rpq, ph) + cubicSplineKernel(rpq, qh) ); 
       if (Wpq==0.0){/*printf("hpq %e rpq %e \t %e \n", hpq, rpq, rpq/hpq); */continue; }

       // \tilde{\psi}_j (x_i)
       psi = Wpq/psph->omega;   
       psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi;
       psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi;
       psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi;

       // \tilde{\psi}_i (x_j)
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
          psph->vPred[j] = pkdVel(pkd,p)[j];
          psph->vPred[j] += pkdAccel(pkd,p)[j]*pDeltaHalf;

          qsph->vPred[j] = pkdVel(pkd,q)[j];
          qsph->vPred[j] += pkdAccel(pkd,q)[j]*qDeltaHalf;

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
/*
      uint64_t A, B;
      A = 4;
      B = 55;

if (p->iOrder == A || p->iOrder==B){
      printf("%d \t p %" PRId64 "; %e \t q %" PRId64 "; %e \t %e \n", i, p->iOrder, 2.*pkdBall(pkd,p), q->iOrder, 2.*pkdBall(pkd,q), rpq);
}
*/

      riemann_input.L.rho -= pDensity*pdivv;
      riemann_input.R.rho -= pkdDensity(pkd,q)*qdivv;
      riemann_input.L.p -= pkd->param.dConstGamma*psph->P*pdivv;
      riemann_input.R.p -= pkd->param.dConstGamma*qsph->P*qdivv;
//if ((p->iOrder == B && q->iOrder==A) || (p->iOrder == A && q->iOrder==B)){     
//      printf("3) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("3) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);
//      printf("p %" PRId64 " - gradRho %e %e %e \n", p->iOrder, psph->gradRho[0], psph->gradRho[1], psph->gradRho[2]);
//      printf("q %" PRId64 " - gradRho %e %e %e \n", q->iOrder, qsph->gradRho[0], qsph->gradRho[1], qsph->gradRho[2]);
//}

       // It seems that in hydra_core_meshless.h they do not rotate any velocity, they just pass the face_unit
       // to the Riemann solver and let it do its stuff..

       
       // It uses as input parameters the density, pressure and velocities at both ends. They are assumed to be given in comoving coordinates, then they are
       // converted into physical units inside Riemann_Solver. 
       
       // IA: DEBUG: Tests for the riemann solver extracted from Toro (10.1007/b79761)
       // Test 1
//       riemann_input.L.rho = 1.0; riemann_input.L.p = 1.0; riemann_input.L.v[0] = 0.0;
//       riemann_input.L.rho = 0.125; riemann_input.L.p = 0.1; riemann_input.L.v[0] = 0.0;

       if (riemann_input.L.rho <= 0) {riemann_input.L.rho = pkdDensity(pkd,p); printf("WARNING, L.rho < 0 : using first-order scheme \n"); }
       if (riemann_input.R.rho <= 0) {riemann_input.R.rho = pkdDensity(pkd,q); printf("WARNING, R.rho < 0 : using first-order scheme \n"); }
       if (riemann_input.L.p <= 0) {riemann_input.L.p = psph->P;    printf("WARNING, L.p < 0 : using first-order scheme \n"); }
       if (riemann_input.R.p <= 0) {riemann_input.R.p = qsph->P;    printf("WARNING, R.p < 0 : using first-order scheme \n");}
       Riemann_solver(pkd, riemann_input, &riemann_output, face_unit, /*double press_tot_limiter TODO For now, just p>0: */ 0.0);
      
       // IA: MFM
#ifdef USE_MFM
       if (riemann_output.Fluxes.rho != 0){
          printf("Frho %e \n",riemann_output.Fluxes.rho);
          abort();
       }
       for (j=0;j<3;j++){
          vFrame[j] += riemann_output.S_M*face_unit[j];
       }
#endif       
       // IA: End MFM

      


       // IA: DEBUG
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

//if ((p->iOrder == A && q->iOrder==B) || (p->iOrder == B && q->iOrder==A)){     
//         printf("Riemann_output: rho %e \tp %e \tv %e %e %e \n", riemann_output.Fluxes.rho, riemann_output.Fluxes.p, riemann_output.Fluxes.v[0], riemann_output.Fluxes.v[1], riemann_output.Fluxes.v[2]);
//}

       if (smf->dTime > 0){
       pmass = pkdField(p,pkd->oMass);
       qmass = pkdField(q,pkd->oMass);


       if (!pkdIsActive(pkd,q)){  
            *qmass += minDt * riemann_output.Fluxes.rho ;

            qsph->mom[0] += minDt * riemann_output.Fluxes.v[0] ;
            qsph->mom[1] += minDt * riemann_output.Fluxes.v[1] ;
            qsph->mom[2] += minDt * riemann_output.Fluxes.v[2] ;

            qsph->E += minDt * riemann_output.Fluxes.p;

            qsph->Uint += minDt * ( riemann_output.Fluxes.p - riemann_output.Fluxes.v[0]*qsph->vPred[0] 
                                                            - riemann_output.Fluxes.v[1]*qsph->vPred[1]
                                                            - riemann_output.Fluxes.v[2]*qsph->vPred[2]
                                  + 0.5*(qsph->vPred[0]*qsph->vPred[0] + qsph->vPred[1]*qsph->vPred[1] + qsph->vPred[2]*qsph->vPred[2]) * riemann_output.Fluxes.rho );
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


    // IA: Mass change limiter. For now, disabled as it does not make a lot of sense. To keep the conservative
    // properties of the code, this could be added as a time limiter
   /* 
    double dmass_limiter = 0.1*pkdMass(pkd,p), dmass_holder = psph->Frho * ( smf->dDelta/(1<<p->uRung) );

    if (fabs(dmass_holder) > dmass_limiter) {
    //   printf("Limiting! \n");
       psph->Frho *= dmass_limiter / fabs(dmass_holder);
    }
    */



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





    // IA: Timestep criteria based on the hydro+grav accelerations

    double a[3], acc;
    float* pa = pkdAccel(pkd,p);
    double cfl = smf->dCFLacc, dtAcc;
    
    for (j=0;j<3;j++) { a[j] = pa[j] + (pkdVel(pkd,p)[j]*psph->Frho + psph->Fmom[j])/pkdMass(pkd,p); }
    acc = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    dtAcc = cfl*sqrt(2*fBall/acc);
    
    // IA: New dtAcc in the case for circular orbits (keplerian disk!)
    //double r = sqrt(pkdPos(pkd,p,0)*pkdPos(pkd,p,0) + pkdPos(pkd,p,1)*pkdPos(pkd,p,1));
    //dtAcc = 2.*3.14159*sqrt(r*r*r) * cfl;
    //printf("r %f \t dtAcc %e \n", sqrt(pkdPos(pkd,p,0)*pkdPos(pkd,p,0) + pkdPos(pkd,p,1)*pkdPos(pkd,p,1)), dtAcc);

    if (dtAcc < dtEst) dtEst = dtAcc;
    uNewRung = pkdDtToRung(dtEst,smf->dDelta,MAX_RUNG);
    if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung; 

    // IA: Timestep limiter that imposes that I must have a dt which is at most, 
    // four times (i.e., 2 rungs) the smallest dt of my neighbours

    if (smf->dTime >= 0){
    for (i=0; i<nSmooth; ++i){
        q = nnList[i].pPart;

        if ( (q->uNewRung - p->uNewRung) > 2) uNewRung = q->uNewRung-1;

        if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung; 
    }
    
    }


}
