#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "parameters.h"
#include "pkd.h"
#include "mdl.h"
#include "smoothfcn.h"
#include "floattype.h"
#include "moments.h"

#include "parameters.h"
#include "cosmo.h"

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
    PST_KICK,
    PST_SETSOFT,
    PST_PHYSICALSOFT,
    PST_PREVARIABLESOFT,
    PST_POSTVARIABLESOFT,
    PST_SETTOTAL,
    PST_ONENODEREADINIT,
    PST_SWAPALL,
    PST_MASSCHECK,
    PST_ACTIVETYPEORDER,
    PST_ACTIVEORDER,
    PST_INITSTEP,
    PST_SETRUNG,
    PST_ACTIVERUNG,
    PST_CURRRUNG,
    PST_GRAVSTEP,
    PST_ACCELSTEP,
    PST_DENSITYSTEP,
    PST_GETMAP,
    PST_ACTIVEEXACTTYPE,
    PST_ACTIVETYPE,
    PST_SETTYPE,
    PST_RESETTYPE,
    PST_COUNTTYPE,
    PST_ACTIVEMASKRUNG,
    PST_ACTIVETYPERUNG,
    PST_SETPARTICLETYPES,
    PST_MARKSMOOTH,
    PST_RESMOOTH,
    PST_RESMOOTHWALK,
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
    int bStandard;
    double dvFac;
    double duTFac;
    int bDoublePos;
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstWriteTipsy(PST,void *,int,void *,int *);

/* PST_BUILDTREE */
struct inBuildTree {
    int nBucket;
    double diCrit2;
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
    double dfBall2OverSoft2;
    SMF smf;
    };
void pstSmooth(PST,void *,int,void *,int *);

/* PST_GRAVITY */
struct inGravity {
    int nReps;
    int bPeriodic;
    int bEwald;
    int iEwOrder;
    int bDoSun;
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
    struct parameters param;
    struct csmContext csm;
    };
struct outStepVeryActive
    {
    int nMaxRung;
    };
void pstStepVeryActiveKDK(PST,void *,int,void *,int *);

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

struct inActiveTypeOrder {
    unsigned int iTestMask;
    };

void pstActiveTypeOrder(PST,void *,int,void *,int *);

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

struct inActiveType {
    unsigned int iFilterMask;
    unsigned int iTestMask;
    unsigned int iSetMask;
    int iRung;
    int bGreater;
    };

void pstActiveExactType(PST,void *,int,void *,int *);
void pstActiveType(PST,void *,int,void *,int *);
void pstSetType(PST,void *,int,void *,int *);
void pstResetType(PST,void *,int,void *,int *);
void pstCountType(PST,void *,int,void *,int *);
void pstActiveMaskRung(PST,void *,int,void *,int *);
void pstActiveTypeRung(PST,void *,int,void *,int *);

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
    double dfBall2OverSoft2;
    SMF smf;
    };
void pstReSmooth(PST,void *,int,void *,int *);

#ifdef GASOLINE
/* PST_RESMOOTHWALK */
void pstReSmoothWalk(PST,void *,int,void *,int *);
#endif

/* PST_INITACCEL */
void pstInitAccel(PST,void *,int,void *,int *);

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
    SMF smf;
    };
void pstFof(PST,void *,int,void *,int *);
void pstGroupMerge(PST,void *,int,void *,int *);
void pstGroupProfiles(PST,void *,int,void *,int *);
#ifdef RELAXATION
void pstInitRelaxation(PST,void *,int,void *,int *);
#endif

#endif
