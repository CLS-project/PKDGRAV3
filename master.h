#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED

#include <stdint.h>

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#define MSR_INIT_E		1
#define MSR_STEP_E		0

typedef struct msrContext {
    PRM prm;
    PST pst;
    MDL mdl;
    LCL lcl;
    FLOAT fCenter[6];
    /*
    ** Parameters.
    */
    struct parameters param;
    /*
    ** Other stuff...
    */
    int nThreads;
    uint64_t N;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    uint64_t nMaxOrder;		/* Order number of last particle */
    int iCurrMaxRung;
    double dThetaMin;

    /*
    ** Tree moments (for Ewald)
    */
    MOMC momTreeRoot[IRUNGMAX+1];
    double momTreeCom[IRUNGMAX+1][3];

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

    uint8_t iRungVeryActive;    /* NOTE: The first very active particle is at iRungVeryActive + 1 */

    /*
     * Domain Decomposition Done
     */
    uint64_t *nRung;
    int iRungDD;
    int iLastRungRT,iLastRungDD;
    uint64_t nActive;
    int nGroups;
    int nBins;
    int bAntiGrav;

    int bSavePending;

    /* Values for a restore from checkpoint */
    double dCheckpointTime;
    int iCheckpointStep;
    int nCheckpointThreads;
    char achCheckpointName[PST_FILENAME_SIZE];
    int nCheckpointClasses;
    PARTCLASS aCheckpointClasses[PKD_MAX_CLASSES];
    } * MSR;

int msrInitialize(MSR *,MDL,int,char **);
void msrLogParams(MSR msr, FILE *fp);
void msrprintf(MSR msr, const char *Format, ... );
int msrGetLock(MSR msr);
int msrCheckForStop(MSR msr);
void msrFinish(MSR);
void msrInitializePStore(MSR msr, uint64_t *nSpecies);
double msrGenerateIC(MSR);
double msrRead(MSR msr,const char *achInFile);
void msrWrite(MSR,const char *,double, int bCheckpoint );
void msrSetSoft(MSR msr,double);
void msrDomainDecomp(MSR,int iRung,int bOthers,int bSplitVA);
void msrBuildTree(MSR msr,double dTime,int bNeedEwald);
void msrBuildTreeVeryActive(MSR msr,double dTime,int bNeedEwald,uint8_t uRungDD);
void msrBuildTreeByRung(MSR msr,double dTime,int bNeedEwald,int iRung);
void msrBuildTreeExcludeVeryActive(MSR msr,double dTime);
void msrBuildTreeMarked(MSR msr,double dTime);
void msrCalcBound(MSR msr,BND *pbnd);
void msrCalcVBound(MSR msr,BND *pbnd);
void msrDomainColor(MSR);
void msrReorder(MSR);
void msrOutArray(MSR,const char *,int);
void msrOutVector(MSR,const char *,int);
void msrSmoothSetSMF(MSR msr, SMF *smf, double dTime);
void msrSmooth(MSR,double,int,int,int);
int msrDoGas(MSR msr);
void msrFastGasPhase1(MSR,double,int);
void msrFastGasPhase2(MSR,double,int);
void msrReSmooth(MSR,double,int,int);
void msrUpdateSoft(MSR,double);
uint8_t msrGravity(MSR msr,uint8_t uRungLo, uint8_t uRungHi,int iRoot1,int iRoot2,
    double dTime,double dStep,int bKickClose,int bKickOpen,int bEwald,int nGroup,int *piSec,uint64_t *pnActive);
void msrCalcEandL(MSR msr,int bFirst,double dTime,double *E,double *T,double *U,double *Eth,double *L,double *F,double *W);
void msrDrift(MSR,double dTime,double dDelta,int iRoot);
void msrScaleVel(MSR msr,double dvFac);
double msrAdjustTime(MSR msr, double aOld, double aNew);
void msrKick(MSR,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
double msrReadCheck(MSR,int *);
void msrWriteCheck(MSR,double,int);
int msrOutTime(MSR,double);
void msrReadOuts(MSR,double);
void msrNewTopStepKDK(MSR msr,
    uint8_t uRung,		/* Rung level */
    double *pdStep,	/* Current step */
    double *pdTime,	/* Current time */
    uint8_t *puRungMax,
    int *piSec);
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
		   int *piSec);
void msrTopStepHSDKD(MSR msr,
		   double dStep,	/* Current step */
		   double dTime,	/* Current time */
		   double dDelta,	/* Time step */
		   int iRung,		/* Rung level */
		   int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		   int iRungVeryActive, /* rung *below which* very active particles are */
		   int iAdjust,		/* Do an adjust? */
		   double *pdActiveSum,
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
void msrActiveMaskRung(MSR msr, unsigned int iSetMask, int iRung, int bGreater);
/*------------------*/
/* Active Functions */
/*------------------*/

void msrVelocityRung(MSR msr,int iRung,double dDelta,double dTime,int bAll);
uint64_t msrCalcWriteStart(MSR);
void msrGetNParts(MSR msr);
void msrAddDelParticles(MSR msr);
void msrGravStep(MSR msr, double dTime);
void msrAccelStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime);
void msrDensityStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime);
int msrUpdateRung(MSR msr, uint8_t uRung);
int msrCountRungs(MSR msr, uint64_t *nRungs);

/*
** Interface functions.
*/
int msrSteps(MSR);
void msrOutputPk(MSR msr,int iStep,double dTime);
void msrCheckpoint(MSR msr, int iStep, double dTime);
double msrRestore(MSR msr);
void msrOutput(MSR msr, int iStep, double dTime, int bCheckpoint);
char *msrOutName(MSR);
char *msrBuildName(MSR msr,char *achFile,int iStep);
char *msrBuildIoName(MSR msr,char *achFile,int iStep);
double msrDelta(MSR);
int msrLogInterval(MSR);
int msrCheckInterval(MSR);
const char *msrCheckTypes(MSR msr);
int msrOutInterval(MSR);
const char *msrOutTypes(MSR msr);
int msrRestart(MSR);
int msrComove(MSR);
double msrSoft(MSR);
int msrDoDensity(MSR);
#ifdef USE_PNG
int msrPNGResolution(MSR msr);
#endif
int msrDoGravity(MSR msr);
void msrInitStep(MSR msr);
void msrSetRung(MSR msr, uint8_t uRungLo, uint8_t uRungHi, int uRung);
void msrZeroNewRung(MSR msr, uint8_t uRungLo, uint8_t uRungHi, int uRung);
int msrMaxRung(MSR msr);
void msrSwitchTheta(MSR msr,double);
uint64_t msrMaxOrder(MSR msr);

void msrFof(MSR msr, double exp);
void msrHop(MSR msr, double exp);
void msrHopWrite(MSR msr, const char *fname);
void msrGroupMerge(MSR msr, double exp);
void msrGroupProfiles(MSR msr, double exp);
void msrOutGroups(MSR msr,const char *,int,double dTime);
void msrDeleteGroups(MSR msr);
void msrInitRelaxation(MSR msr);
void msrRelaxation(MSR msr,double dTime,double deltaT,int iSmoothType,int bSymmetric);
/* Gas routines */
void msrInitSph(MSR,double);
void msrSph(MSR msr, double dTime, double dStep);
void msrSphStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime);
void msrCoolSetup(MSR msr, double);
void msrCooling(MSR msr,double dTime,double dStep,int bUpdateState, int bUpdateTable,int bInterateDt);
void msrStarForm( MSR, double, int);
/* END Gas routines */

void msrHostname(MSR msr);
void msrMemStatus(MSR msr);


void msrSelSrcAll(MSR msr);
void msrSelDstAll(MSR msr);
void msrSelSrcGas(MSR msr);
void msrSelDstGas(MSR msr);
void msrSelSrcStar(MSR msr);
void msrSelDstStar(MSR msr, int, double);
void msrSelSrcDeleted(MSR msr);
void msrSelDstDeleted(MSR msr);
uint64_t msrSelSrcMass(MSR msr,double dMinMass,double dMaxMass,int setIfTrue,int ClearIfFalse);
uint64_t msrSelDstMass(MSR msr,double dMinMass,double dMaxMass,int setIfTrue,int ClearIfFalse);
uint64_t msrSelSrcById(MSR msr,uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstById(MSR msr,uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcPhaseDensity(MSR msr,double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstPhaseDensity(MSR msr,double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcBox(MSR msr,double *dCenter, double *dSize,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstBox(MSR msr,double *dCenter, double *dSize,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcSphere(MSR msr,double *r, double dRadius,int setIfTrue,int clearIfFalse);
uint64_t msrSelDstSphere(MSR msr,double *r, double dRadius,int setIfTrue,int clearIfFalse);
uint64_t msrSelSrcCylinder(MSR msr,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse );
uint64_t msrSelDstCylinder(MSR msr,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse );

void msrDeepestPot(MSR msr,double *r, float *fPot);
double msrTotalMass(MSR msr);
void msrProfile(
    MSR msr, const PROFILEBIN **pBins, int *pnBins, double *r,
    double dMinRadius, double dLogRadius, double dMaxRadius,
    int nPerBin, int nBins, int nAccuracy );
void msrDeleteProfile(MSR msr);

void msrCalcCOM(MSR msr,const double *dCenter, double dRadius,
		double *com, double *vcm, double *L, double *M);
void msrInitGrid(MSR msr,int x,int y,int z);
void msrGridProject(MSR msr,double x,double y,double z);
#ifdef MDL_FFTW
void msrMeasurePk(MSR msr,double *dCenter,double dRadius,int nGrid,int nBins,float *fK,float *fPk);
#endif
void msrPSGroupFinder(MSR msr);
void msrOutPsGroups(MSR msr,const char *pszFile,int iOutType, double dTime);
void msrUnbind(MSR msr);
void msrSetPSGroupIds(MSR msr);

#endif
