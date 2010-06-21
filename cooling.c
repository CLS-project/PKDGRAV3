#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include "cooling.h"
#include "param.h"

#define SET_MODEL FC_MOD_FUNC_(cooling_module,COOLING_MODULE,set_model,SET_MODEL)
void SET_MODEL(int *Nmodel, double *J0in_in, double *J0min_in, 
	       double *alpha_in, double *normfacJ0_in, double *zreioniz_in,
	       int *correct_cooling, double *realistic_ne, 
	       double *h, double *omegab, double *omega0, double *omegaL, 
	       double *astart_sim, double *pT2_sim);
void RAMSES_set_model(int Nmodel, double J0in_in, double J0min_in, 
		      double alpha_in, double normfacJ0_in, double zreioniz_in,
		      int correct_cooling, double realistic_ne, 
		      double h, double omegab, double omega0, double omegaL, 
		      double astart_sim, double *pT2_sim)  {
    SET_MODEL(&Nmodel,&J0in_in,&J0min_in,&alpha_in,
	      &normfacJ0_in
	      ,&zreioniz_in, &correct_cooling,
	      &realistic_ne, 
	      &h,&omegab,&omega0,&omegaL,&astart_sim,
	      pT2_sim);
    }

#define SET_TABLE FC_MOD_FUNC_(cooling_module,COOLING_MODULE,set_table,SET_TABLE)
void SET_TABLE(double *aexp);
void RAMSES_set_table(double aexp) {
    SET_TABLE(&aexp);
    }

#define INTERPOLATE_TABLE FC_MOD_FUNC_(cooling_module,COOLING_MODULE,interpolate_table,INTERPOLATE_TABLE)
void INTERPOLATE_TABLE(double *nHin,double *T2in,double *n_spec,
		       double *pT,double *pmu);
void RAMSES_interpolate_table(double nHin,double T2in,double *n_spec,
			      double *pT,double *pmu) {
    INTERPOLATE_TABLE(&nHin,&T2in,n_spec,pT,pmu);
    }

#define SOLVE_COOLING FC_MOD_FUNC_(cooling_module,COOLING_MODULE,solve_cooling,SOLVE_COOLING)
void SOLVE_COOLING(double *nH,double *T2,double *ZSolar, double *dt,
		   double *deltaT2, int *ncell);
void RAMSES_solve_cooling(double *nH,double *T2,double *ZSolar, double dt,
			  double *deltaT2, int ncell) {
    SOLVE_COOLING(nH,T2,ZSolar,&dt,deltaT2,&ncell);
    }

COOL *CoolInit() 
    {
    COOL *cl;
    cl = (COOL *) malloc(sizeof(COOL));
    assert(cl!=NULL);
    return cl;
    }

void CoolSetup( COOL *cl, double dGmPerCcUnit,double dComovingGmPerCcUnit,
		double dErgPerGmUnit,double dSecUnit,double dKpcUnit,
		double dOmega0,double dHubble0,double dLambda,
		double dOmegab,double dOmegaRad,
		double a,double z,double dTime, COOLPARAM CoolParam)
    {
    double z_reion=8.5, J21=0.0, a_spec=1.0, h0=75;
    double omega_b=.045, omega_m= 0.25, omega_l = 0.75,aexp_ini = 0.01;
    double T2_sim;
    int Nmodel;
    int haardt_madau=1; // default is false, true for clus run

    cl->CoolParam = CoolParam;

    cl->dGmPerCcUnit = dGmPerCcUnit;
    cl->dErgPerGmUnit = dErgPerGmUnit;
    cl->dSecUnit = dSecUnit;
    cl->dErgPerGmPerSecUnit = cl->dErgPerGmUnit / cl->dSecUnit;
    cl->diErgPerGmUnit = 1./dErgPerGmUnit;
    cl->dKpcUnit = dKpcUnit;
    cl->dComovingGmPerCcUnit = cl->dGmPerCcUnit*pow(1+z,3.);

    cl->dnHUnit = cl->dComovingGmPerCcUnit*(0.76/1.66e-24); /* Hardwire to Ramses X=0.76 */
    cl->dT2Unit = cl->dErgPerGmUnit/(1.5*1.38062e-16/1.66e-24);

    aexp_ini = a;
    h0=dHubble0/cl->dSecUnit*3.0857e19; /* Convert to km/s per Mpc */
    /* nCOM = del_star*omega_b*rhoc/aexp**3*X/mH -- nISM is max of nCOM and nstar */
    cl->nCOMFac = 200.*dOmegab*1.88e-29*pow(h0/100.,2.0)*0.76/1.66e-24;
    omega_l = dLambda;
    omega_b = dOmegab;
    omega_m = dOmega0;
    haardt_madau = CoolParam.bUV;

    Nmodel=-1;
    if(!haardt_madau)Nmodel=2;
    RAMSES_set_model(Nmodel,J21*1e-21,-1.0,a_spec,-1.0,z_reion, 
		     -1,2, h0/100.,omega_b,omega_m,omega_l, 
		     aexp_ini,&T2_sim);

    cl->dT2Init = T2_sim;

    printf("CoolSetup %d %d %g %g %g %g %g %g %g %g\n",haardt_madau,Nmodel,aexp_ini,omega_l,h0,omega_b,omega_m,omega_l,z_reion,T2_sim);
    }

void CoolSetTime( COOL *cl, double dTime, double z, int bUpdateTable ) {

    cl->z = z;
    cl->dTime = dTime;
    cl->dComovingGmPerCcUnit = cl->dGmPerCcUnit*pow(1+z,3.);
    cl->dnHUnit = cl->dComovingGmPerCcUnit*(0.76/1.66e-24); /* Hardwire to Ramses */
    
    cl->nISM = cl->CoolParam.dTFloor_nstar;
    {
    double nCOM = cl->nCOMFac*pow(1+z,3.0);
    printf("Cool Set Time: %g %g %g %d\n",z,nCOM,cl->nISM,bUpdateTable);
    if (nCOM > cl->nISM) cl->nISM = nCOM;
    }

    if (bUpdateTable) {
	double aexp = 1/(1+z);
	RAMSES_set_table(aexp); /* Generate new cooling tables */
	printf("Cool: Called RAMSES_set_table\n");
	}
    }

void clApplyT2Floor( COOL *cl, double nH, double *pT2 ) {
    if (cl->CoolParam.bTFloor) {
	double T2Floor = cl->CoolParam.dTFloor_T2star*pow(nH/cl->nISM,
							  cl->CoolParam.dTFloor_gstar-1)
	    +cl->CoolParam.dTFloor_T2min;  /* How Ramses does floor Tf1+Tf2 */
	if (*pT2 < T2Floor) *pT2 = T2Floor;
	}
    }
    

void CoolIntegrateEnergyCode(COOL *cl, COOLPARTICLE *cp, double *ECode, 
			     double ExternalHeatingCode, double dRhoCode, 
			     double ZMetal, double *posCode, double tStep ) {
    double dt;
    double nH,T2,deltaT2,ZSolar;

    /* How Ramses does PdV - operator split -- 
       better to do predictor corrector or include cts. PdV heating in cooling */
    *ECode += ExternalHeatingCode*tStep; 

    T2 = cl->dT2Unit*(*ECode); /* cgs for Ramses */
    nH = cl->dnHUnit*dRhoCode;
    ZSolar = ZMetal/0.02;
    dt = cl->dSecUnit*tStep; 
    RAMSES_solve_cooling(&nH,&T2,&ZSolar,dt,&deltaT2,1); /* Radiative Cooling/Heating */
    T2 += deltaT2;
    clApplyT2Floor( cl, nH, &T2 );
    *ECode = T2/cl->dT2Unit;
    }

   
void CoolAddParams( COOLPARAM *CoolParam, void *prm ) {
    CoolParam->bUV = 1;
    prmAddParam(prm,"bUV",0,&CoolParam->bUV,sizeof(int),"UV",
		"use UV (Haard Madau) = +UV");
    CoolParam->bTFloor = 1;
    prmAddParam(prm,"bTFloor",0,&CoolParam->bTFloor,sizeof(int),"tfloor",
		"use TFloor = +tfloor");
    CoolParam->dTFloor_T2min = 0.01;
    prmAddParam(prm,"dTFloor_T2min",2,&CoolParam->dTFloor_T2min,sizeof(double),"t2min",
		"dTFloor_T2min = 0.01" );
    CoolParam->dTFloor_T2star = 1e4;
    prmAddParam(prm,"dTFloor_T2star ",2,&CoolParam->dTFloor_T2star,sizeof(double),"t2fstar",
		"dTFloor_star = 1e4" );
    CoolParam->dTFloor_nstar = 0.1;
    prmAddParam(prm,"dTFloor_nstar ",2,&CoolParam->dTFloor_nstar ,sizeof(double),"t2fnstar",
		"dTFloor_nstar = 0.1" );
    CoolParam->dTFloor_gstar = 5/3.;
    prmAddParam(prm,"dTFloor",2,&CoolParam->dTFloor_gstar,sizeof(double),"t2floorgstar",
		"dTFloor_gstar = 5/3." );
    }
	
void CoolLogParams( COOLPARAM *CoolParam, FILE *fp ) {
    fprintf(fp,"\n# Cooling: bUV: %d",CoolParam->bUV);
    fprintf(fp," bTFloor: %d",CoolParam->bTFloor);
    fprintf(fp," dTFloor: %g",CoolParam->dTFloor_T2min);
    fprintf(fp," dTFloor_T2star: %g",CoolParam->dTFloor_T2star);
    fprintf(fp," dTFloor_nstar: %g",CoolParam->dTFloor_nstar);
    fprintf(fp," dTFloor_gstar: %g",CoolParam->dTFloor_gstar);
    }


void CoolTempFromEnergyCode( COOL *cl, COOLPARTICLE *cp, double *pE, double *pT,double dRhoCode, double fMetal )
    {
    double T2, nH, nspec[6], mu;
    
    T2 = cl->dT2Unit*(*pE);
    nH = cl->dnHUnit*dRhoCode;
    clApplyT2Floor( cl, nH, &T2 );
    RAMSES_interpolate_table(nH,T2,nspec,pT,&mu);
    }

void CoolEnergyCodeFromTemp( COOL *cl, COOLPARTICLE *cp, double *pE, double *pT,double dRhoCode, double fMetal )
    {
    double T, nH, nspec[6], mu;
    double T2Left, T2Right, T2Mid, TLeft, TRight, TMid;

    nH = cl->dnHUnit*dRhoCode;
    T = *pT;
    /* mu = 1.21961; X=.76, Y=0.24, zero metals, neutral */
    /* mu = 0.588235; X=.76, Y=0.24, zero metals, ionized */
    T2Left = T/1.22; 
    T2Right = T/0.58; 
    RAMSES_interpolate_table(nH,T2Left,nspec,&TLeft,&mu);
    assert(mu < 1.22);
    RAMSES_interpolate_table(nH,T2Right,nspec,&TRight,&mu);
    assert(mu > 0.58);
    for (;;) {
	T2Mid = 0.5*(T2Left+T2Right);
	if (T2Right-T2Left < 1e-4*T2Mid) break;
	RAMSES_interpolate_table(nH,T2Mid,nspec,&TMid,&mu);
	if (T > TMid) {
	    T2Left = T2Mid; 
	    TLeft = TMid;
	    }
	else {
	    T2Right = T2Mid;
	    TRight = TMid;
	    }
	}
    clApplyT2Floor( cl, nH, &T2Mid );
    *pE = T2Mid/cl->dT2Unit;
    }

void CoolInitEnergyCode( COOL *cl, COOLPARTICLE *cp, double *pE, double *pT, double dRhoCode, double fMetal )
    {
    double T2 = cl->dT2Init;
    clApplyT2Floor( cl, cl->dnHUnit*dRhoCode, &T2 );
    *pE = T2/cl->dT2Unit;

/*    fprintf(stderr,"SPECIAL INIT: %g  %g  %g %g  %g %g\n",cl->dnHUnit*dRhoCode, *pT,*pT*1.71101e-08,*pE,cl->dT2Init,T2);
      assert(0);*/
    }
