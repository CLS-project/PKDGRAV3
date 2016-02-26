#ifndef COOLING_HINCLUDED
#define COOLING_HINCLUDED

typedef struct CoolingParametersStruct {
    int    bUV;
    int    bTFloor;
    double dTFloor_T2min;
    double dTFloor_T2star;
    double dTFloor_nstar;
    double dTFloor_gstar;
    } COOLPARAM;

typedef struct CoolingContextStruct { 
    double     z; /* Redshift */
    double     dTime;
    
    COOLPARAM  CoolParam;
    
    double     dGmPerCcUnit;
    double     dComovingGmPerCcUnit;
    double     dErgPerGmUnit;
    double     dSecUnit;
    double     dErgPerGmPerSecUnit;
    double     diErgPerGmUnit;
    double     dKpcUnit;

    double     dnHUnit;
    double     nISM;
    double     nCOMFac;
    double     dT2Unit;    
    double     dT2Init;
    
    double     Y_H; /* not used right now */
    double     Y_He;
    double     Y_eMAX;
    
/* Diagnostic */
    int        its;
#if defined(COOLDEBUG) 
    MDL        mdl; /* For diag/debug outputs */
    struct particle *p; /* particle pointer needed for SN feedback */
#endif
    } COOL;

typedef struct CoolingParticleStruct {
    double Y_HI,Y_HeI,Y_HeII;	/* Abundance of ions */
    } COOLPARTICLE;

COOL *CoolInit();
void CoolSetup( COOL *cl, double dGmPerCcUnit,double dComovingGmPerCcUnit,
		double dErgPerGmUnit,double dSecUnit,double dKpcUnit,
		double dOmega0,double dHubble0,double dLambda,
		double dOmegab,double dOmegaRad,
		double a,double z,double dTime, COOLPARAM CoolParam);
void CoolSetTime( COOL *cl, double dTime, double z, int bUpdateTable );
void CoolIntegrateEnergyCode(COOL *cl, COOLPARTICLE *cp, double *ECode, 
			     double ExternalHeatingCode, double rhoCode, 
			     double ZMetal, double *posCode, double tStep );
void CoolAddParams( COOLPARAM *CoolParam, void *prm );
void CoolLogParams( COOLPARAM *CoolParam, FILE *fp );
void CoolTempFromEnergyCode( COOL *cl, COOLPARTICLE *cp, double *pE, double *pT,double dDensity, double fMetal );
void CoolEnergyCodeFromTemp( COOL *cl, COOLPARTICLE *cp, double *pE, double *pT,double dDensity, double fMetal );
void CoolInitEnergyCode( COOL *cl, COOLPARTICLE *cp, double *pE, double *pT, double dDensity, double fMetal );

#endif
