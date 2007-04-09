#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"
#include "floattype.h"

#define MSR_INIT_E		1
#define MSR_STEP_E		0

typedef struct msrContext {
    PRM prm;
    PST pst;
    MDL mdl;
    LCL lcl;
    FLOAT fCenter[3];
    /*
    ** Parameters.
    */
    struct parameters param;	   
    /*
    ** Other stuff...
    */
    int nThreads;
    int N;
    int nDark;
    int nGas;
    int nStar;
    int nMaxOrder;		/* Order number of last particle */
    int nMaxOrderGas;
    int nMaxOrderDark;
    int iCurrMaxRung;
    double dCrit;
    /*
    ** Comoving coordinate variables.
    */
    double dEcosmo;
    double dUOld;
    double dTimeOld;
    /*
    ** Redshift output points.
    */
    int nMaxOuts;
    int nOuts;
    double *pdOutTime;
    int iOut;
    /*
    ** Processor mapping for one-node-output functions.
    */
    int *pMap;
    /*
     * Domain Decomposition Done
     */
    int *nRung;
    int bDoneDomainDecomp;
    int iLastRungDD;
    int iLastRungSD;
    int nActive;
    int nSmoothActive;
    int nGroups;
    int nBins;
    int bAntiGrav;
    } * MSR;

void msrInitialize(MSR *,MDL,int,char **);
void msrLogParams(MSR msr, FILE *fp);
int msrGetLock(MSR msr);
int msrCheckForStop(MSR msr);
void msrFinish(MSR);
double msrReadTipsy(MSR);
void msrWriteTipsy(MSR,char *,double);
void msrSetSoft(MSR msr,double);
void msrDomainDecomp(MSR,int iRung,int bGreater,int bSplitVA);
void msrBuildTree(MSR msr,double dMass,double dTime);
void msrBuildTreeExcludeVeryActive(MSR msr,double dMass,double dTime);

#ifdef GASOLINE
void msrCalcBallBound(MSR,double fBallFactor);
#endif

void msrDomainColor(MSR);
void msrReorder(MSR);
void msrOutArray(MSR,char *,int);
void msrOutVector(MSR,char *,int);
void msrSmooth(MSR,double,int,int,int);
void msrReSmooth(MSR,double,int,int,int);
void msrUpdateSoft(MSR,double);
void msrGravity(MSR,double,double,int *,double *,double *,double *,int *);
void msrCalcEandL(MSR,int,double,double *,double *,double *,double *,double *);
void msrDrift(MSR,double,double);
void msrKick(MSR,double,double);
double msrReadCheck(MSR,int *);
void msrWriteCheck(MSR,double,int);
int msrOutTime(MSR,double);
void msrReadOuts(MSR,double);
void msrTopStepKDK(MSR msr,
		   double dStep,	/* Current step */
		   double dTime,	/* Current time */
		   double dDelta,	/* Time step */
		   int iRung,		/* Rung level */
		   int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		   int iRungVeryActive, /* rung *below which* very active particles are */
		   int iAdjust,		/* Do an adjust? */
		   double *pdActiveSum,
		   double *pdWMax,
		   double *pdIMax,
		   double *pdEMax,
		   int *piSec);
void msrStepVeryActiveKDK(MSR msr, double dStep, double dTime, double dDelta,
			  int iRung);

void msrBallMax(MSR msr, int iRung, int bGreater);
/*------------------*/
/* Active Functions */
/*------------------*/
void msrActiveRung(MSR msr, int iRung, int bGreater);
void msrActiveOrder(MSR msr);

/* Replacement functions */
void msrActiveExactType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask);
void msrActiveType(MSR msr, unsigned int iTestMask, unsigned int iSetMask);
void msrSetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask);
void msrResetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask);
int  msrCountType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask);
void msrActiveMaskRung(MSR msr, unsigned int iSetMask, int iRung, int bGreater);
void msrActiveTypeRung(MSR msr, unsigned int iTestMask, unsigned int iSetMask, int iRung, int bGreater);
void msrActiveTypeOrder(MSR msr, unsigned int iTestMask );
/*------------------*/
/* Active Functions */
/*------------------*/

void msrVelocityRung(MSR msr,int iRung,double dDelta,double dTime,int bAll);
void msrCalcWriteStart(MSR);
void msrAddDelParticles(MSR msr);
void msrGravStep(MSR msr, double dTime);
void msrAccelStep(MSR msr, double dTime);
void msrDensityStep(MSR msr, double dTime);
void msrInitDt(MSR msr);
int msrDtToRung(MSR msr, int iRung, double dDelta, int bAll);

/*
** Interface functions.
*/
int msrSteps(MSR);
char *msrOutName(MSR);
double msrDelta(MSR);
int msrLogInterval(MSR);
int msrCheckInterval(MSR);
int msrOutInterval(MSR);
int msrRestart(MSR);
int msrComove(MSR);
double msrSoft(MSR);
int msrDoDensity(MSR);
int msrDoGravity(MSR msr);
void msrInitStep(MSR msr);
void msrSetRung(MSR msr, int iRung);
void msrInitAccel(MSR msr);
void msrSwitchTheta(MSR msr,double);
int msrMaxOrder(MSR msr);

void msrInitTimeSteps(MSR,double,double);

void msrFof(MSR msr,int nFOFsDone,int iSmoothType,int bSymmetric);
void msrGroupMerge(MSR msr);
void msrGroupProfiles(MSR msr,int nFOFsDone,int iSmoothType,int bSymmetric);
void msrOutGroups(MSR msr,char *,int,double dTime);
void msrDeleteGroups(MSR msr);
#ifdef RELAXATION
void msrInitRelaxation(MSR msr);
void msrRelaxation(MSR msr,double dTime,double deltaT,int iSmoothType,int bSymmetric);
#endif /* RELAXATION  */

#endif
