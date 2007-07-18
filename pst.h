#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "parameters.h"
#include "pkd.h"
#include "mdl.h"
#include "smoothfcn.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "moments.h"

#include "parameters.h"
#include "cosmo.h"
#ifdef PLANETS
#include "collision.h"
#endif


typedef struct lclBlock {
    char *pszDataPath;
    PKD	pkd;
    int nPstLvl;
    int iWtFrom;
    int iWtTo;
    int iPart;
    int iOrdSplit;
    FLOAT fSplit;
    FLOAT fWtLow;
    FLOAT fWtHigh;
    FLOAT fLow;
    FLOAT fHigh;
    int nSplit;
    int nWriteStart;
    } LCL;

typedef struct pstContext {
    struct pstContext *pstLower;
    MDL mdl;
    LCL *plcl;
    int idSelf;
    int idUpper;
    int nLeaves;
    int nLower;
    int nUpper;
    int iLvl;
    BND bnd;
    int iSplitDim;
    int iOrdSplit;
    FLOAT fSplit;
    FLOAT fSplitInactive;
    int nTotal;
    int iVASplitSide;
    int nLowTot;    /* total number of particles in the lower subset of processors */
    int nHighTot;   /* total number of particles in the upper subset of processors */
    int ioIndex;
    } * PST;


#define PST_SERVICES		100
#define PST_FILENAME_SIZE	512

enum pst_service {
    PST_SRV_STOP,
    PST_SETADD,
    PST_LEVELIZE,
    PST_READTIPSY,
    PST_DOMAINDECOMP,
    PST_CALCBOUND,
    PST_RUNGDDWEIGHT,
    PST_WEIGHT,
    PST_COUNTVA,
    PST_WEIGHTWRAP,
    PST_FREESTORE,
    PST_COLREJECTS,
    PST_SWAPREJECTS,
    PST_COLORDREJECTS,
    PST_DOMAINORDER,
    PST_LOCALORDER,
    PST_OUTARRAY,
    PST_OUTVECTOR,
    PST_WRITETIPSY,
    PST_BUILDTREE,
    PST_DISTRIBCELLS,
    PST_CALCBOUNDBALL,
    PST_DISTRIBBOUNDBALL,
    PST_DISTRIBROOT,
    PST_CALCROOT,
    PST_SMOOTH,
    PST_GRAVITY,
    PST_GRAVEXTERNAL,
    PST_CALCEANDL,
    PST_CALCEANDLEXT,
    PST_DRIFT,
    PST_DRIFTINACTIVE,
    PST_CACHEBARRIER,
    PST_ROPARTICLECACHE,
    PST_PARTICLECACHEFINISH,
    PST_STEPVERYACTIVE,
    PST_STEPVERYACTIVEH,
    PST_COPY0,
    PST_PREDICTOR,
    PST_CORRECTOR,
    PST_SUNCORRECTOR,
    PST_PREDICTORINACTIVE,
    PST_KICK,
    PST_SETSOFT,
    PST_PHYSICALSOFT,
    PST_PREVARIABLESOFT,
    PST_POSTVARIABLESOFT,
    PST_SETTOTAL,
    PST_ONENODEREADINIT,
    PST_SWAPALL,
    PST_MASSCHECK,
    PST_ACTIVEORDER,
    PST_INITSTEP,
    PST_SETRUNG,
    PST_ACTIVERUNG,
    PST_CURRRUNG,
    PST_GRAVSTEP,
    PST_ACCELSTEP,
    PST_DENSITYSTEP,
    PST_GETMAP,
    PST_SETRUNGVERYACTIVE,
    PST_SETPARTICLETYPES,
    PST_MARKSMOOTH,
    PST_RESMOOTH,
    PST_INITACCEL,
    PST_DTTORUNG,
    PST_INITDT,
    PST_ORDWEIGHT,
    PST_SETWRITESTART,
    PST_COLNPARTS,
    PST_NEWORDER,
    PST_SETNPARTS,
    PST_DENSCHECK,
    PST_FOF,
    PST_GROUPMERGE, 
    PST_GROUPPROFILES,
    PST_INITRELAXATION,
    PST_CLEARTIMER,
    PST_FINDIOS,
    PST_STARTIO,
#ifdef PLANETS
    PST_READSS,
    PST_WRITESS,
    PST_SUNINDIRECT,
    PST_GRAVSUN,
    PST_NEXTCOLLISION,
    PST_GETCOLLIDERINFO,
    PST_DOCOLLISION,
    PST_GETVARIABLEVERYACTIVE,
    PST_CHECKHELIODIST,
#ifdef SYMBA
    PST_STEPVERYACTIVES,
    PST_DRMINTORUNG,
    PST_MOMSUN,
    PST_DRIFTSUN,
    PST_KEPLERDRIFT,
#endif /* SYMBA */
#endif /* PLANETS */
#ifdef HERMITE
    PST_AARSETHSTEP,
    PST_FIRSTDT,
#endif
    };

void pstAddServices(PST,MDL);
void pstInitialize(PST *,MDL,LCL *);
void pstFinish(PST);

/* PST_SETADD */

struct inSetAdd {
    int id;
    };
void pstSetAdd(PST,void *,int,void *,int *);

/* PST_LEVELIZE */
struct inLevelize {
    int iLvl;
    };
void pstLevelize(PST,void *,int,void *,int *);

/* PST_READTIPSY */
struct inReadTipsy {
    int nFileStart;
    int nFileEnd;
    int nDark;	
    int nGas;
    int nStar;
    float fExtraStore;
    FLOAT fPeriod[3];
    int bStandard;
    double dvFac;
    double dTuFac;
    int bDoublePos;
    char achInFile[PST_FILENAME_SIZE];
    char achOutName[PST_FILENAME_SIZE];
    };
void pstReadTipsy(PST,void *,int,void *,int *);

/* PST_DOMAINDECOMP */
struct inDomainDecomp {
    BND bnd;
    int nBndWrap[3];
    int bDoRootFind;
    int bDoSplitDimFind;
    int bSplitVA;
    int nActive;
    int nTotal;
    };
void pstDomainDecomp(PST,void *,int,void *,int *);

/* PST_CALCBOUND */
struct outCalcBound {
    BND bnd;
    };
void pstCalcBound(PST,void *,int,void *,int *);

struct inRungDDWeight {
    int iMaxRung;
    double dWeight; 
    };

void pstRungDDWeight(PST,void *,int,void *,int *);

/* PST_WEIGHT */
struct inWeight {
    int iSplitDim;
    FLOAT fSplit;
    int iSplitSide;
    int ittr;
    int pFlag;
    };
struct outWeight {
    int nLow;
    int nHigh;
    FLOAT fLow;
    FLOAT fHigh;
    };
void pstWeight(PST,void *,int,void *,int *);

/* PST_COUNTVA */
struct inCountVA {
    int iSplitDim;
    FLOAT fSplit;
    };
struct outCountVA {
    int nLow;
    int nHigh;
    };
void pstCountVA(PST,void *,int,void *,int *);

/* PST_WEIGHTWRAP */
struct inWeightWrap {
    int iSplitDim;
    FLOAT fSplit;
    FLOAT fSplit2;
    int iSplitSide;
    int ittr;
    int iVASplitSide;
    };
struct outWeightWrap {
    int nLow;
    int nHigh;
    };
void pstWeightWrap(PST,void *,int,void *,int *);

/* PST_FREESTORE */
struct outFreeStore {
    int nFreeStore;
    };
void pstFreeStore(PST,void *,int,void *,int *);

/*
** This structure is used by reject collectors and SwapRejects
*/
typedef struct outReject {
    int id;
    int nRejects;
    int nSpace;
    int nLocal;
    } OREJ;

/* PST_COLREJECTS */
void pstColRejects(PST,void *,int,void *,int *);

/* PST_SWAPREJECTS */
void pstSwapRejects(PST,void *,int,void *,int *);

/* PST_COLORDREJECTS */
struct inColOrdRejects {
    int iOrdSplit;
    int iSplitSide;
    };
void pstColOrdRejects(PST,void *,int,void *,int *);

/* PST_DOMAINORDER */
struct inDomainOrder {
    int iMaxOrder;
    };
void pstDomainOrder(PST,void *,int,void *,int *);

/* PST_LOCALORDER */
void pstLocalOrder(PST,void *,int,void *,int *);

/* PST_OUTARRAY */
struct inOutArray {
    char achOutFile[PST_FILENAME_SIZE];
    int iType;
    };
void pstOutArray(PST,void *,int,void *,int *);

/* PST_OUTVECTOR */
struct inOutVector {
    char achOutFile[PST_FILENAME_SIZE];
    int iDim;
    int iType;
    };
void pstOutVector(PST,void *,int,void *,int *);

/* PST_WRITETIPSY */
struct inWriteTipsy {
    double dTime;
    double dvFac;
    double duTFac;
    int bDoublePos;
    int bStandard;
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstWriteTipsy(PST,void *,int,void *,int *);

#ifdef USE_MDL_IO

/* PST_FINDIOS */
struct inFindIOS {
    int nLower; /* Number of particles to the left */
    int N;
    };
struct outFindIOS {
    int nCount[MDL_MAX_IO_PROCS];
    };
void pstFindIOS(PST,void *,int,void *,int *);

/* PST_STARTIO */
struct inStartIO {
#ifdef IO_SPLIT
    int nCount[MDL_MAX_IO_PROCS];
#endif
    double dvFac;
    double duTFac;
    double dTime;
    int bDoublePos;
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstStartIO(PST,void *,int,void *,int *);
#endif

/* PST_BUILDTREE */
struct inBuildTree {
    double diCrit2;
    double dTimeStamp;
    int nBucket;
    int iCell;
    int nCell;
    int bTreeSqueeze;
    int bExcludeVeryActive;
    };
void pstBuildTree(PST,void *,int,void *,int *);

/* PST_DISTRIBCELLS */
void pstDistribCells(PST,void *,int,void *,int *);

#ifdef GASOLINE
/* PST_CALCBOUNDBALL */
struct inCalcBoundBall {
    double fBallFactor;
    int iCell;
    int nCell;
    };
void pstCalcBoundBall(PST,void *,int,void *,int *);

/* PST_DISTRIBBOUNDBALL */
void pstDistribBoundBall(PST,void *,int,void *,int *);
#endif /* of GASOLINE */

/* PST_CALCROOT */
struct ioCalcRoot {
    MOMC momc;
    };
void pstCalcRoot(PST,void *,int,void *,int *);

/* PST_DISTRIBROOT */
void pstDistribRoot(PST,void *,int,void *,int *);

/* PST_SMOOTH */
struct inSmooth {
    int nSmooth;
    int bGasOnly;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    int eParticleTypes; /* Smooth over which particle types */
    double dfBall2OverSoft2;
    SMF smf;
    };
void pstSmooth(PST,void *,int,void *,int *);

/* PST_GRAVITY */
struct inGravity {
  double dTime;
    int nReps;
    int bPeriodic;
    int bEwald;
    int iEwOrder;
    double dEwCut;
    double dEwhCut;
    };
struct outGravity {
    int nActive;
    int nTreeActive;
    double aSun[3];
    double dPartSum;
    double dCellSum;
    double dFlop;
    /*	
    ** Collected CPU time stats.
    */
    double dWSum;
    double dWMax;
    double dWMin;
    double dISum;
    double dIMax;
    double dIMin;
    double dESum;
    double dEMax;
    double dEMin;
    /*
    ** Cache Statistics.
    */
    double dpASum;
    double dpMSum;
    double dpCSum;
    double dpTSum;
    double dcASum;
    double dcMSum;
    double dcCSum;
    double dcTSum;
    };
void pstGravity(PST,void *,int,void *,int *);

/* PST_CALCEANDL */
struct outCalcEandL {
    double T;
    double U;
    double Eth;
    double L[3];
    };
void pstCalcEandL(PST,void *,int,void *,int *);

/* PST_DRIFT */
struct inDrift {
    double dTime;
    double dDelta;
    FLOAT fCenter[3];
    int bPeriodic;
    int bFandG;
    FLOAT fCentMass;
    };
void pstDrift(PST,void *,int,void *,int *);
void pstDriftInactive(PST,void *,int,void *,int *);

/* PST_ROPARTICLECACHE */

void pstROParticleCache(PST, void *, int, void *, int *);

/* PST_PARTICLECACHEFINISH */

void pstParticleCacheFinish(PST, void *, int, void *, int *);

/* PST_CACHEBARRIER */
void pstCacheBarrier(PST, void *, int, void *, int *);

/* PST_STEPVERYACTIVE */
struct inStepVeryActive 
    {
    double dStep;
    double dTime;
    double dDelta;
    int iRung;
    int nMaxRung;
    double diCrit2;
    double aSunInact[3];
    double adSunInact[3];
    double dSunMass;
    };
struct outStepVeryActive
    {
    int nMaxRung;
    };
void pstStepVeryActiveKDK(PST,void *,int,void *,int *);

#ifdef HERMITE
/* PST_STEPVERYACTIVEH */
struct inStepVeryActiveH 
    {
    double dStep;
    double dTime;
    double dDelta;
    int iRung;
    int nMaxRung;
    double diCrit2;
    double aSunInact[3];
    double adSunInact[3];
    double dSunMass;
    };
struct outStepVeryActiveH
    {
    int nMaxRung;
    };
void pstStepVeryActiveHermite(PST,void *,int,void *,int *);

/* PST_COPY0 */
struct inCopy0 {
    double dTime;   
    };
void pstCopy0(PST,void *,int,void *,int *);

/* PST_PREDICTOR */
struct inPredictor {
    double dTime;
    };
void pstPredictor(PST,void *,int,void *,int *);

/* PST_CORRECTOR */
struct inCorrector {
  double dTime;
    };
void pstCorrector(PST,void *,int,void *,int *);

/* PST_SUNCORRECTOR */
struct inSunCorrector {
    double dTime;
    double dSunMass;
    };
void pstSunCorrector(PST,void *,int,void *,int *);

/* PST_PREDICTORINACTIVE */
struct inPredictorInactive {
    double dTime;
    };
void pstPredictorInactive(PST,void *,int,void *,int *);

/* PST_AARSETHSTEP */
struct inAarsethStep {
    double dEta;
    };
void pstAarsethStep(PST,void *,int,void *,int *);

/* PST_FIRSTDT */
void pstFirstDt(PST,void *,int,void *,int *);

#endif /* Hermite*/

/* PST_UPDATEUDOT */
struct inUpdateuDot {
    double duDelta;
    double z;
    int iGasModel;
    int bUpdateY;
    };
struct outUpdateuDot {
    double Time;
    double MaxTime;
    double SumTime;
    int nSum;
    };

void pstUpdateuDot(PST,void *,int,void *,int *);

/* PST_KICK */
struct inKick {
    double dvFacOne;
    double dvFacTwo;
    double dvPredFacOne;
    double dvPredFacTwo;
    double duDelta;
    double duPredDelta;
    double duDotLimit;
    int iGasModel;
    double z;
    };
struct outKick {
    double Time;
    double MaxTime;
    double SumTime;
    int nSum;
    };

void pstKick(PST,void *,int,void *,int *);

/* PST_SETSOFT */
struct inSetSoft {
    double dSoft;
    };
void pstSetSoft(PST,void *,int,void *,int *);

#ifdef CHANGESOFT
/* PST_PHYSICALSOFT */
struct inPhysicalSoft {
    double dSoftMax;
    double dFac;
    int bSoftMaxMul;
    };
void pstPhysicalSoft(PST,void *,int,void *,int *);

/* PST_PREVARIABLESOFT */
void pstPreVariableSoft(PST,void *,int,void *,int *);

/* PST_POSTVARIABLESOFT */
struct inPostVariableSoft {
    double dSoftMax;
    int bSoftMaxMul;
    };
void pstPostVariableSoft(PST,void *,int,void *,int *);
#endif

/* PST_SETTOTAL */
struct outSetTotal {
    int nTotal;
    };
void pstSetTotal(PST,void *,int,void *,int *);

/* PST_ONENODEREADINIT */
void pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SWAPALL */
void pstSwapAll(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_ACTIVEORDER */
void pstActiveOrder(PST,void *,int,void *,int *);

/* PST_SETRUNG */
struct inSetRung {
    int iRung;
    };
void pstSetRung(PST,void *,int,void *,int *);

/* PST_INITSTEP */
struct inInitStep {
    struct parameters param;
    struct csmContext csm;
    };
void pstInitStep(PST,void *,int,void *,int *);

struct inDensCheck {
    int iRung;
    int bGreater;
    int iMeasure;
    };
struct outDensCheck {
    double dMaxDensError;
    double dAvgDensError;
    int nError;
    int nTotal;
    };

void pstDensCheck(PST,void *,int,void *,int *);

/* PST_ACTIVERUNG */
struct inActiveRung {
    int iRung;
    int bGreater;
    };

void pstActiveRung(PST,void *,int,void *,int *);

/* PST_CURRRUNG */
struct inCurrRung {
    int iRung;
    };
struct outCurrRung {
    int iCurrent;
    };
void pstCurrRung(PST,void *,int,void *,int *);

/* PST_GRAVSTEP */
struct inGravStep {
    double dEta;
    double dRhoFac;
    };
void pstGravStep(PST,void *,int,void *,int *);

/* PST_ACCELSTEP */
struct inAccelStep {
    double dEta;
    double dVelFac;
    double dAccFac;
    int    bDoGravity;
    int    bEpsAcc;
    int    bSqrtPhi;
    double dhMinOverSoft;
    };
void pstAccelStep(PST,void *,int,void *,int *);

/* PST_DENSITYSTEP */
struct inDensityStep {
    double dEta;
    double dRhoFac;
    };
void pstDensityStep(PST,void *,int,void *,int *);

/* PST_GETMAP */
struct inGetMap {
    int nStart;
    };
void pstGetMap(PST,void *,int,void *,int *);

void pstSetRungVeryActive(PST,void *,int,void *,int *);

struct inSetParticleTypes {
    int iDummy;
    };

void pstSetParticleTypes(PST,void *,int,void *,int *);

/* PST_RESMOOTH */
struct inReSmooth {
    int nSmooth;
    int bGasOnly;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    int eParticleTypes; /* Smooth over which particle types */
    double dfBall2OverSoft2;
    SMF smf;
    };
void pstReSmooth(PST,void *,int,void *,int *);

/* PST_DTTORUNG */
struct inDtToRung {
    int iRung;
    double dDelta;
    int iMaxRung;
    int bAll;
    };
struct outDtToRung {
    int nRungCount[256];
    };
void pstDtToRung(PST,void *,int,void *,int *);

/* PST_INITDT */
struct inInitDt {
    double dDelta;
    };
void pstInitDt(PST,void *,int,void *,int *);

/* PST_ORDWEIGHT */
struct inOrdWeight {
    int iOrdSplit;
    int iSplitSide;
    int ittr;
    };
struct outOrdWeight {
    int nLow;
    int nHigh;
    };
void pstOrdWeight(PST,void *,int,void *,int *);

/* PST_SETWRITESTART */
struct inSetWriteStart {
    int nWriteStart;
    };
void pstSetWriteStart(PST,void *,int,void *,int *);

/* PST_COLNPARTS */
struct outColNParts {
    int nNew;
    int nDeltaGas;
    int nDeltaDark;
    int nDeltaStar;
    };
void pstColNParts(PST, void *, int, void *, int *);

/* PST_NEWORDER */
void pstNewOrder(PST, void *, int, void *, int *);

/* PST_SETNPARTS */
struct inSetNParts {
    int nGas;
    int nDark;
    int nStar;
    int nMaxOrderGas;
    int nMaxOrderDark;
    double dSunMass;
    };
void pstSetNParts(PST, void *, int, void *, int *);

/* PST_MARKSMOOTH */
struct inMarkSmooth {
    int nSmooth;
    int bGasOnly;
    int bPeriodic;
    int bSymmetric;
    int iMarkType;
    double dfBall2OverSoft2;
    SMF smf;
    };
void pstMarkSmooth(PST,void *,int,void *,int *);


/* PST_CLEARTIMER */
struct inClearTimer 
    {
    int iTimer;
    };

void pstClearTimer(PST,void *,int,void *,int *);

struct inFof {
    int nFOFsDone;
    int nSmooth;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    int eParticleTypes; /* Smooth over which particle types */
    SMF smf;
    };
struct inGroupMerge{
    int bPeriodic;
    SMF smf;
    };
struct inGroupProfiles{
    int nSmooth;
    int nFOFsDone;
    int bPeriodic;
    int nTotalGroups;
    int bLogBins;
    int bSymmetric;
    int iSmoothType;
    int eParticleTypes; /* Smooth over which particle types */
    SMF smf;
    };
void pstFof(PST,void *,int,void *,int *);
void pstGroupMerge(PST,void *,int,void *,int *);
void pstGroupProfiles(PST,void *,int,void *,int *);
#ifdef RELAXATION
void pstInitRelaxation(PST,void *,int,void *,int *);
#endif


#ifdef PLANETS
/* PLANETS begin */

/* PST_WRITESS */
struct inWriteSS {
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteSS(PST,void *,int,void *,int *);

/* PST_READSS */
struct inReadSS {
	int nFileStart;
	int nFileEnd;
	int nDark;
	int nGas;			/* always zero */
	int nStar;			/* always zero */
	int iOrder;
	float fExtraStore;
	FLOAT fPeriod[3];	/* for compatability */
	char achInFile[PST_FILENAME_SIZE];
        double dSunMass;
	};
void pstReadSS(PST,void *,int,void *,int *);

/* PST_SUNINDIRECT */
struct inSunIndirect{
      int iFlag;
     }; 
struct outSunIndirect{
      double aSun[3];
      double adSun[3];
     }; 
void pstSunIndirect(PST,void *,int,void *,int *);

/* PST_GRAVSUN */
struct inGravSun{
      double aSun[3];
      double adSun[3];
      double dSunMass;
      };
void pstGravSun(PST,void *,int,void *,int *);

/* PST_NEXTCOLLISION */
struct outNextCollision {
	double dt;
	int iOrder1,iOrder2;
	};
void pstNextCollision(PST,void *,int,void *,int *);

/* PST_GETCOLLIDERINFO */
struct inGetColliderInfo {
	int iOrder;
	};
struct outGetColliderInfo {
	COLLIDER Collider;
	};
void pstGetColliderInfo(PST,void *,int,void *,int *);

/* PST_DOCOLLISION */
struct inDoCollision {
	double dt;
	COLLIDER Collider1,Collider2;
        int bPeriodic;
	COLLISION_PARAMS CP;
	};
struct outDoCollision {
	COLLIDER Out[MAX_NUM_FRAG];
	double dT;
        int iOutcome,nOut;
        };
void pstDoCollision(PST,void *,int,void *,int *);

/* PST_GETVARIABLEVERYACTIVE */
struct outGetVariableVeryActive {
        double dDeltaEcoll;
        };
void pstGetVariableVeryActive(PST,void *,int,void *,int *);

/* PST_CHECKHELIODIST */
struct outCheckHelioDist {
        double dT;
        double dSM;
        };
void pstCheckHelioDist(PST,void *,int,void *,int *);

#ifdef SYMBA

/* PST_STEPVERYACTIVESYMBA */
struct inStepVeryActiveS 
    {
    double dStep;
    double dTime;
    double dDelta;
    int iRung;
    int nMaxRung;
    double diCrit2;
    double dSunMass;
    };
struct outStepVeryActiveS
    {
    int nMaxRung;
    };

void pstStepVeryActiveSymba(PST,void *,int,void *,int *);

/* PST_DRMINTORUNG */
struct inDrminToRung {
    int iRung;
    int iMaxRung;
    };
struct outDrminToRung {
    int nRungCount[256];
    };
void pstDrminToRung(PST,void *,int,void *,int *);

/* PST_MOMSUN */
struct outMomSun{
    double momSun[3];
};

void pstMomSun(PST,void *,int,void *,int *); 

/* PST_DRIFTSUN */
struct inDriftSun{
    double vSun[3];
    double dDelta;
};
void pstDriftSun(PST,void *,int,void *,int *); 

/* PST_KEPLERDRIFT */
struct inKeplerDrift{
    double dDelta;
    double dSunMass;
};

void pstKeplerDrift(PST,void *,int,void *,int *); 

#endif /* SYMBA */
#endif /* PLANETS */

#endif 
