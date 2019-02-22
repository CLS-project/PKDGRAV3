/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED

#include <stdint.h>
#include <signal.h>
#include <time.h>

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"

#define MSR_INIT_E		1
#define MSR_STEP_E		0

static time_t timeGlobalSignalTime = 0;
static int bGlobalOutput = 0;

#ifndef _MSC_VER
static inline void USR1_handler(int signo) {
    signal(SIGUSR1,USR1_handler);
    timeGlobalSignalTime = time(0);
    }

static inline void USR2_handler(int signo) {
    signal(SIGUSR2,USR2_handler);
    bGlobalOutput = 1;
    }
#endif

typedef struct msrContext {
    PRM prm;
    PST pst;
    MDL mdl;
    LCL lcl;
    double fCenter[6];
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
    int iRungDD, iRungDT;
    int iLastRungRT,iLastRungDD;
    uint64_t nActive;
    int nGroups;
    int nBins;
    int bAntiGrav;

    int bSavePending;

    long lStart; /* starting time of job */


    /* Values for a restore from checkpoint */
    double dCheckpointTime;
    int iCheckpointStep;
    int nCheckpointThreads;
    char achCheckpointName[PST_FILENAME_SIZE];
    int nCheckpointClasses;
    PARTCLASS aCheckpointClasses[PKD_MAX_CLASSES];
    } * MSR;
#ifdef __cplusplus
extern "C" {
#endif
double msrTime();
int msrInitialize(MSR *,MDL,int,char **);
void msrLogParams(MSR msr, FILE *fp);
void msrprintf(MSR msr, const char *Format, ... );
int msrGetLock(MSR msr);
int msrCheckForStop(MSR msr, const char *achStopFile);
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
int msrMeshlessHydro(MSR msr);
int msrFirstHydroLoop(MSR msr);
void msrSetFirstHydroLoop(MSR msr, int value);
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
void msrCheckForOutput(MSR msr,int iStep,double dTime,int *pbDoCheckpoint,int *pbDoOutput);
int msrNewTopStepKDK(MSR msr,
    int bDualTree,      /* Should be zero at rung 0! */
    uint8_t uRung,	/* Rung level */
    double *pdStep,	/* Current step */
    double *pdTime,	/* Current time */
    uint8_t *puRungMax,int *piSec,int *pbDoCheckpoint,int *pbDoOutput);
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

void msrLightConeOpen(MSR msr, int iStep);
void msrLightConeClose(MSR msr, int iStep);
void msrLightConeVel(MSR msr);
void msrInflate(MSR msr,int iStep);

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
void msrOutputLinPk(MSR msr, int iStep, double dTime);
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

void msrNewFof(MSR msr, double exp);
void msrGroupStats(MSR msr);
void msrHop(MSR msr, double exp);
void msrHopWrite(MSR msr, const char *fname);
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


void msrSelAll(MSR msr);
void msrSelGas(MSR msr);
void msrSelStar(MSR msr);
void msrSelDeleted(MSR msr);
uint64_t msrSrcMass(MSR msr,double dMinMass,double dMaxMass,int setIfTrue,int ClearIfFalse);
uint64_t msrSrcById(MSR msr,uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse);
uint64_t msrSrcPhaseDensity(MSR msr,double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse);
uint64_t msrSrcBox(MSR msr,double *dCenter, double *dSize,int setIfTrue,int clearIfFalse);
uint64_t msrSrcSphere(MSR msr,double *r, double dRadius,int setIfTrue,int clearIfFalse);
uint64_t msrSrcCylinder(MSR msr,double *dP1, double *dP2, double dRadius, int setIfTrue, int clearIfFalse );

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
void msrGridCreateFFT(MSR msr, int nGrid);
void msrGridDeleteFFT(MSR msr);
void msrAssignMass(MSR msr,int iAssignment,int nGrid);
void msrMeasurePk(MSR msr,int iAssignment,int bInterlace,int nGrid,int nBins,uint64_t *nPk,float *fK,float *fPk);
void msrMeasureLinPk(MSR msr,int nGridLin,double a,double dBoxSize,
                uint64_t *nPk,float *fK,float *fPk);
void msrSetLinGrid(MSR msr,double dTime, int nGrid, int bKickClose, int bKickOpen);
void msrLinearKick(MSR msr, double dTime, int bKickClose, int bKickOpen);
#endif
void msrPSGroupFinder(MSR msr);
void msrOutPsGroups(MSR msr,const char *pszFile,int iOutType, double dTime);
void msrUnbind(MSR msr);
void msrSetPSGroupIds(MSR msr);
int msrGetParticles(MSR msr, int nIn, uint64_t *ID, struct outGetParticles *out);
void msrOutputOrbits(MSR msr,int iStep,double dTime);
#ifdef __cplusplus
}
#endif
#endif
