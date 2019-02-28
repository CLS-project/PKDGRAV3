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

    for (i=0;i<6;i++){ psph1->B[i] += psph2->B[i]; }
    psph1->omega += psph2->omega;

    //IA: the fluxes are not added because they has not been yet computed!

    }

void combSecondHydroLoop(void *vpkd, void *p1,void *p2) {
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
      p->uNewRung = 0;
      for (i=0;i<3;i++) { 
         psph->Fmom[i] = 0.0;
	}
    }

/* We build the geometrical information of the cell, this is:
 *  - omega
 *  - B matrix
 */
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {

    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double E[6];  // IA: We are assumming 3D here!
    double ph, qh, rpq, hpq, dx,dy,dz,Wpq;
    double psi, psiTilde[3];
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
    ph = fBall; /* IA: fBall seems to be the minimun distance to any neighbors ¿? Although it would be logical to be the maximum */


    /* IA: Compute the \omega(x_i) normalization factor */
    psph->omega = 0.0;
    for (i=0; i<nSmooth; ++i){
       q = nnList[i].pPart;

       qh = pkdBall(pkd,q); // TODO: Check this 

       rpq = sqrt(nnList[i].fDist2);
       hpq = ph;// 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)
       //printf("rpq %e ph %e qh %e Wpq %e \n", rpq, ph, qh, cubicSplineKernel(rpq, hpq));

       psph->omega += cubicSplineKernel(rpq, hpq);
    }


    /* IA: We compute the density making use of Eq. 27 Hopkins 2015 */
    pkdSetDensity(pkd,p, pkdMass(pkd,p)*psph->omega);
/*    if(pkdDensity(pkd,p) > 3)*/  //  printf("mass %e \t omega %e \t density %e \n",pkdMass(pkd,p), psph->omega, pkdMass(pkd,p)*psph->omega);



    /* IA: Compute the E matrix (Hopkins 2015, eq 14) */
    for (i=0; i<6;++i){
       E[i] = 0.0; 
    }

    for (i=0; i<nSmooth;++i){
       q = nnList[i].pPart;

       dx = nnList[i].dx;
       dy = nnList[i].dy;
       dz = nnList[i].dz;

       qh = pkdBall(pkd,q); // TODO: Check this 

       rpq = sqrt(nnList[i].fDist2);
       hpq = ph; //0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

       Wpq = cubicSplineKernel(rpq, hpq);

       E[XX] += dx*dx*Wpq;
       E[YY] += dy*dy*Wpq;
    //   printf("dx %e dy %e dz %e Wpq %e \n", dx, dy, dz, Wpq);
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
//    printf("x %e y %e z %e \n", pkdPos(pkd,p,0), pkdPos(pkd,p,1), pkdPos(pkd,p,2));
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

    }


/* IA: This routine will extrapolate the primitives to the 'faces' and solve the 1D riemann problem. 
 * For now, the 1D riemann flux will be computed TWICE for a given face, one for each adjacent particles... This 
 * is highly sub-optimal, and I should think a way to 'mark' which riemann fluxes has been already computed.
 * IDEA: maybe a 64-bit with 0/1? But I think that I have no info about the neighbors of my neighbor..Or create a queue
 * of fluxes to be computed.. */
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double pv[3], qv[3], vFrame[3];
    double ph,qh, hpq, modApq, rpq, dx,dy,dz,Wpq,psi_p,psi_q,pDensity,pDeltaHalf, qDeltaHalf, dt2, dtEst, vsig_pq, dvDotd, dvDotdr;
    double rhodiff, vxdiff, vydiff, vzdiff, pdiff, pdivv, qdivv, pdiffOverRhop, pdiffOverRhoq;
    double psiTilde_p[3], psiTilde_q[3], Apq[3], face_unit[3], dr[3]; 
    struct Input_vec_Riemann riemann_input;
    struct Riemann_outputs riemann_output; //TODO: do this outside the loop
    uint8_t uNewRung;
    int i,j;

    psph = pkdSph(pkd, p);   

    //TODO: I think that this assummes that both particles are syncronized, this is, that both had the last primitive update
    // at the same time. If we have 'p' active and 'q' is not, then the forward prediction for the set of primitives variables
    // will be different for each particle!!!
    //dtp = smf->pkd->param.dDelta/(1<<p->uRung);/* size of step just completed by star */
// printf("dDelta %e \n", smf->dDelta/(1<<p->uRung));
    pDeltaHalf = 0.5*( smf->dDelta/(1<<p->uRung) );

    pDensity = pkdDensity(pkd,p);

     

    // The conversion to primitive variables for the p-particle has only to be done once
    if (pkdIsActive(pkd,p)){ //FIXME TODO there is no need for doing this inside the q-particles loop!!
       // Rho has already been computed from the first hydro loop
       riemann_input.L.p = psph->E;
       for (j=0; j<3; j++){
          psph->vPred[j] = psph->mom[j]/pkdMass(pkd,p);
          riemann_input.R.p -= 0.5*psph->mom[j]*psph->mom[j]/pkdMass(pkd,p);
       }
       riemann_input.L.p *= psph->omega*(pkd->param.dConstGamma - 1.);

       // We set new particle velocities 
       for (j=0;j<3; j++){
          pkdVel(pkd,p)[j] = psph->vPred[j];
       }
    }
    psph->c = sqrt(riemann_input.L.p*pkd->param.dConstGamma/pDensity);
    dtEst = HUGE_VAL;

    for (i=0;i<nSmooth;++i){

	 q = nnList[i].pPart;
       qsph = pkdSph(pkd, q);

       qDeltaHalf = pDeltaHalf; //0.5*( smf->dDelta/(1<<q->uRung) );
//       printf("pDeltaHalf %e qDeltaHalf %e \n", pDeltaHalf, qDeltaHalf);

	 dx = nnList[i].dx;
	 dy = nnList[i].dy;
	 dz = nnList[i].dz;

       /* IA: in the nnList there is a 'copy' of the own particle, which we can omit as there is no fluxes
        * to be computed here */
       if (dx==0 && dy==0 && dz==0) continue;

       qh = pkdBall(pkd,q); // TODO: Check this, it may or not be updated 
       ph = fBall;
       rpq = sqrt(nnList[i].fDist2);
       hpq = 0.5*(qh+ph); // IA: We symmetrize the kernel size (probably needed, not sure)

       Wpq = cubicSplineKernel(rpq, hpq); 
       psi_p = Wpq/psph->omega;
       psi_q = Wpq/qsph->omega;
//       printf("omega_p %e omega_q %e \n", psph->omega, qsph->omega);
//       printf("psi_p %e \t psi_q %e \n", psi_p, psi_q);

       /* Compute the gradients and extrapolate to the faces */
      psiTilde_p[0] = (psph->B[XX]*dx + psph->B[XY]*dy + psph->B[XZ]*dz)*psi_p;
      psiTilde_p[1] = (psph->B[XY]*dx + psph->B[YY]*dy + psph->B[YZ]*dz)*psi_p;
      psiTilde_p[2] = (psph->B[XZ]*dx + psph->B[YZ]*dy + psph->B[ZZ]*dz)*psi_p;

//      printf("dx %e dy %e dz %e \n", dx, dy, dz);
       psiTilde_q[0] = -(qsph->B[XX]*dx + qsph->B[XY]*dy + qsph->B[XZ]*dz)*psi_q;
       psiTilde_q[1] = -(qsph->B[XY]*dx + qsph->B[YY]*dy + qsph->B[YZ]*dz)*psi_q;
       psiTilde_q[2] = -(qsph->B[XZ]*dx + qsph->B[YZ]*dy + qsph->B[ZZ]*dz)*psi_q;

//       printf("B_p [XX] %e \t [XY] %e \t [XZ] %e \n \t \t \t [YY] %e \t [YZ] %e \n \t\t\t \t \t \t [ZZ] %e \n", psph->B[XX], psph->B[XY], psph->B[XZ], psph->B[YY], psph->B[YZ], psph->B[ZZ]);
//       printf("B_q [XX] %e \t [XY] %e \t [XZ] %e \n \t \t \t [YY] %e \t [YZ] %e \n \t\t\t \t \t \t [ZZ] %e \n", qsph->B[XX], qsph->B[XY], qsph->B[XZ], qsph->B[YY], qsph->B[YZ], qsph->B[ZZ]);

//       printf("dx %e \t dy %e \t dz %e \n", dx, dy, dz);
//       printf("psiTilde_p %e \t %e \t %e \n", psiTilde_p[0], psiTilde_p[1], psiTilde_p[2]);
//       printf("psiTilde_q %e \t %e \t %e \n", psiTilde_q[0], psiTilde_q[1], psiTilde_q[2]);

       // Compute primitives variables from conserved ones:
       // m, mv, E -> rho, v, p
       // From AREPO, pag 820, last paragraph, they say that the primitives variables are only updated
       // for the active particles, but i do not see the advantage of that.
//       if (pkdIsActive(pkd,q)){ //FIXME This is repeated too many times!! This should be done in another particle loop? W/o neighs of course
          // Rho has already been computed from the first hydro loop
          riemann_input.R.p = qsph->E;
//          printf("qsph->E %e \n", qsph->E);
          for (j=0; j<3; j++){
             qsph->vPred[j] = qsph->mom[j]/pkdMass(pkd,q);
             riemann_input.R.p -= 0.5*qsph->mom[j]*qsph->mom[j]/pkdMass(pkd,q);
          }
          riemann_input.R.p *= qsph->omega*(pkd->param.dConstGamma - 1.);

          // We set new particle velocities TODO: is this the right place to do this?
          for (j=0;j<3; j++){
             pkdVel(pkd,q)[j] = psph->vPred[j];
          }
//       }


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

      // Differences in the variables
      rhodiff = (pkdDensity(pkd,q) - pDensity);
      vxdiff = (qsph->vPred[0] - psph->vPred[0]);
      vydiff = (qsph->vPred[1] - psph->vPred[1]);
      vzdiff = (qsph->vPred[2] - psph->vPred[2]);
      pdiff = (riemann_input.R.p - riemann_input.L.p);

 // FIXME This cs computation is also done inside riemann_solver
       qsph->c = sqrt(riemann_input.R.p*pkd->param.dConstGamma/pkdDensity(pkd,q));

       // From Eqs 24,25 Hopkins 2015, to limit deltaT

       dvDotdr = (dx*vxdiff + dy*vydiff + dz*vzdiff);
       if (dvDotdr < 0) {
          vsig_pq = psph->c + qsph->c - dvDotdr/sqrt(nnList[i].fDist2);
       }else{
          vsig_pq = psph->c + qsph->c;
       }

       dt2 = 2.*smf->dEtaCourant * fBall /vsig_pq;	
	 if (dt2 < dtEst) dtEst=dt2;

      // Divergence of the velocity field for the forward in time prediction
      pdivv = 0.0;
      qdivv = 0.0;

      pdivv += vxdiff*psiTilde_p[0];
      pdivv += vydiff*psiTilde_p[1];
      pdivv += vzdiff*psiTilde_p[2];

      qdivv -= vxdiff*psiTilde_q[0];
      qdivv -= vydiff*psiTilde_q[1];
      qdivv -= vzdiff*psiTilde_q[2];

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

//      printf("1) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("1) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

      // This notation may be confusing... sorry
      pdiffOverRhop = -pdiff/pDensity * pDeltaHalf;  // We include here the minus comming from eq. A4
      pdiffOverRhoq = pdiff/pkdDensity(pkd,q) * qDeltaHalf;

      // TODO: Has this to be done in physical units? I guess..
      // We to add the terms on the divergence of v from the forward prediction TODO: I have commented them for simplicity
      riemann_input.L.rho -= pkdDensity(pkd,p)*pdivv;
      riemann_input.R.rho -= pkdDensity(pkd,q)*qdivv;
//      printf("2) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
      for (j=0; j<3; j++) { 
         riemann_input.L.v[j] -= pv[j]*pdivv;
         riemann_input.R.v[j] -= qv[j]*qdivv;
      }
      riemann_input.L.p -= pkd->param.dConstGamma*riemann_input.L.p*pdivv;
      riemann_input.R.p -= pkd->param.dConstGamma*riemann_input.R.p*qdivv;
//      printf("2) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);

      // We add the gradients terms (from extrapolation and forward prediction)
      for (j=0; j<3; j++) {
         riemann_input.L.rho += ( 0.5*dr[j] - pDeltaHalf*pv[j])*psiTilde_p[j]*rhodiff;
         riemann_input.R.rho += (-0.5*dr[j] - qDeltaHalf*qv[j])*psiTilde_q[j]*(-rhodiff);

         riemann_input.L.v[0] += (0.5*dr[j]*vxdiff + pdiffOverRhop)*psiTilde_p[j];
         riemann_input.R.v[0] += (0.5*dr[j]*vxdiff + pdiffOverRhoq)*psiTilde_q[j]; // The two minus at the first term compensate each other
                                                   
         riemann_input.L.v[1] += (0.5*dr[j]*vydiff + pdiffOverRhop)*psiTilde_p[j];
         riemann_input.R.v[1] += (0.5*dr[j]*vydiff + pdiffOverRhoq)*psiTilde_q[j];
                                                   
         riemann_input.L.v[2] += (0.5*dr[j]*vzdiff + pdiffOverRhop)*psiTilde_p[j];
         riemann_input.R.v[2] += (0.5*dr[j]*vzdiff + pdiffOverRhoq)*psiTilde_q[j];

         riemann_input.L.p += ( 0.5*dr[j] - pDeltaHalf*pv[j])*psiTilde_p[j]*pdiff;
         riemann_input.R.p += (-0.5*dr[j] - qDeltaHalf*qv[j])*psiTilde_q[j]*(-pdiff);
      }
//      printf("3) L.rho %e \t R.rho %e \n", riemann_input.L.rho, riemann_input.R.rho);
//      printf("3) L.p %e \t R.p %e \n", riemann_input.L.p, riemann_input.R.p);


      // The cs and u values of the reconstructed states are computed inside the Riemann Solver!




       // Face where the riemann problem will be solved
       modApq = 0.0;
       for (j=0; j<3; j++){
          Apq[j] = psiTilde_q[j]/psph->omega - psiTilde_p[j]/qsph->omega;
//          printf("Apq %d %e \n", j, psph->omega);
          modApq += Apq[j]*Apq[j];
       }
//       printf("modApq %e \n", modApq);
       modApq = sqrt(modApq);

       for (j=0; j<3; j++){
          face_unit[j] = Apq[j]/modApq;
       }

       // It seems that in hydra_core_meshless.h they do not rotate any velocity, they just pass the face_unit
       // to the Riemann solver and let it do its stuff..

       
       // It uses as input parameters the density, pressure and velocities at both ends. They are assumed to be given in comoving coordinates, then they are
       // converted into physical units inside Riemann_Solver. 
//       printf("unit_face %e %e %e \n", face_unit[0], face_unit[1], face_unit[2]);
       riemann_input.L.p = riemann_input.R.p;     
       riemann_input.L.rho = riemann_input.R.rho;     
       Riemann_solver(riemann_input, &riemann_output, face_unit, /*double press_tot_limiter TODO For now, just p>0: */ 0.0);
//       printf("riemann_output.Fluxes.v %e %e %e \n", riemann_output.Fluxes.v[0], riemann_output.Fluxes.v[1], riemann_output.Fluxes.v[2]);

       // Check for NAN fluxes
       if (riemann_output.Fluxes.rho!=riemann_output.Fluxes.rho) printf("riemann_output.Fluxes.rho %e \n", riemann_output.Fluxes.rho);
       if (riemann_output.Fluxes.p!=riemann_output.Fluxes.p) abort();

//       psph->Frho = riemann_output.Fluxes.rho;
//       psph->Fene = riemann_output.Fluxes.p;
//       for (j=0;j<3;j++){ psph->Fmom[j] = riemann_output.Fluxes.v[j]; }

       // Now we de-boost the fluxes following Eq. A8 Hopkins 2015 (From hydra_core_meshless.h):
       /* the fluxes have been calculated in the rest frame of the interface: we need to de-boost to the 'simulation frame'
        which we do following Pakmor et al. 2011 */
       for(j=0;j<3;j++)
       {
           riemann_output.Fluxes.p += vFrame[j] * riemann_output.Fluxes.v[j];
           riemann_output.Fluxes.p += (0.5*vFrame[j]*vFrame[j])*riemann_output.Fluxes.rho;
       }

       // IA: Now we just multiply by the face area
       psph->Frho += riemann_output.Fluxes.rho*modApq; 
       psph->Fene += riemann_output.Fluxes.p*modApq; 
       for(j=0;j<3;j++) {psph->Fmom[j] += (riemann_output.Fluxes.v[j]+vFrame[j] * riemann_output.Fluxes.rho)*modApq;} // momentum flux (need to divide by mass) 

       /* IA: we will compute the same flux from p->q and q->p so we now halve the fluxes */
       psph->Frho *= 0.5;
       psph->Fene *= 0.5;
       psph->Fmom[0] *= 0.5;
       psph->Fmom[1] *= 0.5;
       psph->Fmom[2] *= 0.5;


       // IA: we store this flux also for the q particle
       qsph->Frho -= psph->Frho;
       qsph->Fene -= psph->Fene;
       for(j=0;j<3;j++){ qsph->Fmom[j] -= psph->Fmom[j]; }



                
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


    uNewRung = pkdDtToRung(dtEst,smf->dDelta,MAX_RUNG);
    if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung; 
//    printf("dtEst %e uNewRung %d \n", dtEst, p->uNewRung);

}
