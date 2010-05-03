#ifndef PST_HINCLUDED
#define PST_HINCLUDED
#define PST_H_MODULE_ID "$Id$"

#include "pkd.h"
#include "mdl.h"
#include "smoothfcn.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "moments.h"
#include "outtype.h"

#include "parameters.h"
#include "cosmo.h"
#ifdef PLANETS
#include "collision.h"
#endif


typedef struct lclBlock {
    char *pszDataPath;
    PKD	pkd;
    int iWtFrom;
    int iWtTo;
    int iPart;
    uint64_t iOrdSplit;
    FLOAT fSplit;
    FLOAT fWtLow;
    FLOAT fWtHigh;
    FLOAT fLow;
    FLOAT fHigh;
    int nSplit;
    uint64_t nWriteStart;
    PKDOUT pkdout;
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
    uint64_t iOrdSplit;
    FLOAT fSplit;
    FLOAT fSplitInactive;
    uint64_t nTotal;
    int iVASplitSide;
    uint64_t nLowTot;    /* total number of particles in the lower subset of processors */
    uint64_t nHighTot;   /* total number of particles in the upper subset of processors */
    int ioIndex;
    } * PST;


#define PST_SERVICES		100
#define PST_FILENAME_SIZE	512
#define PST_MAX_FILES           800

enum pst_service {
    PST_SRV_STOP,
    PST_SPLITIO,
    PST_SETADD,
    PST_READFILE,
    PST_PEANOHILBERTCOUNT,
    PST_DOMAINDECOMP,
    PST_CALCBOUND,
    PST_COMBINEBOUND,
    PST_WEIGHT,
    PST_COUNTVA,
    PST_WEIGHTWRAP,
    PST_FREESTORE,
    PST_COLREJECTS,
    PST_SWAPREJECTS,
    PST_COLORDREJECTS,
    PST_DOMAINORDER,
    PST_LOCALORDER,
    PST_COMPRESSASCII,
    PST_WRITEASCII,
    PST_WRITETIPSY,
    PST_BUILDTREE,
    PST_DISTRIBCELLS,
    PST_CALCROOT,
    PST_DISTRIBROOT,
    PST_ENFORCEPERIODIC,
    PST_TREENUMSRCACTIVE,
    PST_BOUNDSWALK,
    PST_SMOOTH,
    PST_FASTGASPHASE1,
    PST_FASTGASPHASE2,
    PST_FASTGASCLEANUP,
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
    PST_SETTOTAL,
    PST_ONENODEREADINIT,
    PST_SWAPALL,
    PST_MASSCHECK,
    PST_ACTIVEORDER,
    PST_INITSTEP,
    PST_SETRUNG,
    PST_ZERONEWRUNG,
    PST_ACTIVERUNG,
    PST_CURRRUNG,
    PST_GRAVSTEP,
    PST_ACCELSTEP,
    PST_SPHSTEP,
    PST_STARFORM,
    PST_DENSITYSTEP,
    PST_COOLSETUP,
    PST_COOLING,
    PST_CORRECTENERGY,
    PST_GETMAP,
    PST_SETRUNGVERYACTIVE,
    PST_MARKSMOOTH,
    PST_RESMOOTH,
    PST_INITACCEL,
    PST_UPDATERUNG,
    PST_INITDT,
    PST_ORDWEIGHT,
    PST_SETWRITESTART,
    PST_ADDWRITESTART,
    PST_COLNPARTS,
    PST_NEWORDER,
    PST_GETNPARTS,
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
    PST_HANDSUNMASS,
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
#ifdef USE_GRAFIC
    PST_GENERATEIC,
#endif
    PST_HOSTNAME,
    PST_MEMSTATUS,
    PST_GETCLASSES,
    PST_SETCLASSES,
    PST_SWAPCLASSES,

    PST_SELSRCALL,
    PST_SELDSTALL,
    PST_SELSRCGAS,
    PST_SELDSTGAS,
    PST_SELSRCSTAR,
    PST_SELDSTSTAR,
    PST_SELSRCDELETED,
    PST_SELDSTDELETED,
    PST_SELSRCMASS,
    PST_SELDSTMASS,
    PST_SELSRCBYID,
    PST_SELDSTBYID,
    PST_SELSRCPHASEDENSITY,
    PST_SELDSTPHASEDENSITY,

    PST_SELSRCBOX,
    PST_SELDSTBOX,
    PST_SELSRCSPHERE,
    PST_SELDSTSPHERE,
    PST_SELSRCCYLINDER,
    PST_SELDSTCYLINDER,

    PST_DEEPESTPOT,
    PST_PROFILE,
    PST_CALCDISTANCE,
    PST_CALCCOM,
    PST_COUNTDISTANCE,
    PST_PEAKVC,
    PST_INITGRID,
    PST_GRIDPROJECT,
#ifdef MDL_FFTW
    PST_MEASUREPK,
#endif
    PST_TOTALMASS,
    PST_BUILDPSDTREE,
    PST_PSD,
    PST_PSFOF,
    };

void pstAddServices(PST,MDL);
void pstInitialize(PST *,MDL,LCL *);
void pstFinish(PST);

/* PST_SETADD */

struct inSetAdd {
    int idLower;
    int idUpper;
    };
void pstSetAdd(PST,void *,int,void *,int *);

/* PST_READFILE */
struct inReadFile {
    uint64_t nNodeStart; /* First particle to read (of total) */
    uint64_t nNodeEnd;   /* Last particle to read (of total) */
    uint64_t nSpecies[FIO_SPECIES_LAST];
    uint64_t mMemoryModel;
    double fPeriod[3];
    double dvFac;
    double dTuFac;
    double dOmega0;
    double dOmegab;
    float fExtraStore;
    int nTreeBitsLo;
    int nTreeBitsHi;
    int nBucket;
    //int nFiles;
    int iCacheSize;
    int nProcessors;
    //uint8_t bStandard;
    //uint8_t bDoublePos;
    //short   eFileType;
    char achFilename[PST_FILENAME_SIZE];
    };
void pstReadFile(PST,void *,int,void *,int *);


struct inPeanoHilbertCount {
    };
struct outPeanoHilbertCount {
    };
void pstPeanoHilbertCount(PST,void *,int,void *,int *);

/* PST_DOMAINDECOMP */
struct inDomainDecomp {
    BND bnd;
    int nBndWrap[3];
    int bDoRootFind;
    int bDoSplitDimFind;
    int bSplitVA;
    uint64_t nActive;
    uint64_t nTotal;
    };
void pstDomainDecomp(PST,void *,int,void *,int *);

/* PST_CALCBOUND */
void pstCalcBound(PST,void *,int,void *,int *);

/* PST_COMBINEBOUND */
void pstCombineBound(PST,void *,int,void *,int *);

/* PST_WEIGHT */
struct inWeight {
    int iSplitDim;
    FLOAT fSplit;
    int iSplitSide;
    int ittr;
    int pFlag;
    };
struct outWeight {
    uint64_t nLow;
    uint64_t nHigh;
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
    uint64_t nLow;
    uint64_t nHigh;
    };
void pstWeightWrap(PST,void *,int,void *,int *);

/* PST_FREESTORE */
struct outFreeStore {
    uint64_t nFreeStore;
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
    uint64_t iOrdSplit;
    int iSplitSide;
    };
void pstColOrdRejects(PST,void *,int,void *,int *);

/* PST_DOMAINORDER */
struct inDomainOrder {
    uint64_t iMaxOrder;
    };
void pstDomainOrder(PST,void *,int,void *,int *);

/* PST_LOCALORDER */
void pstLocalOrder(PST,void *,int,void *,int *);

/* PST_COMPRESSASCII */
struct inCompressASCII {
    uint64_t nTotal;
    int iFile;
    int iType;
    int iDim;
    };
struct outCompressASCII {
    uint64_t nBytes;
    };
void pstCompressASCII(PST,void *,int,void *,int *);

/* PST_WRITEASCII */
struct inWriteASCII {
    uint64_t nFileOffset;
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstWriteASCII(PST,void *,int,void *,int *);

/* PST_WRITETIPSY */
struct inWriteTipsy {
    double dTime;
    double dvFac;
    int bDoublePos;
    int bStandard;
    int nProcessors;
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
    uint64_t N;
    double dvFac;
    double dTime;
    double dEcosmo;
    double dTimeOld;
    double dUOld;
    int bDoublePos;
    char achOutName[PST_FILENAME_SIZE];
    };
void pstStartIO(PST,void *,int,void *,int *);

/* PST_IO_LOAD */
struct inIOLoad {
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    double dvFac;
    double dTuFac;
    int nBucket;
    float fExtraStore;
    FLOAT fPeriod[3];
    };

void pstIOLoad(PST,void *,int,void *,int *);

#endif

/* PST_BUILDTREE */
struct inBuildTree {
    double diCrit2;
    int nBucket;
    int iCell;
    int nCell;
    int bExcludeVeryActive;
    };
void pstBuildTree(PST,void *,int,void *,int *);

/* PST_DISTRIBCELLS */
void pstDistribCells(PST,void *,int,void *,int *);

/* PST_CALCROOT */
struct ioCalcRoot {
    MOMC momc;
    };
void pstCalcRoot(PST,void *,int,void *,int *);

/* PST_DISTRIBROOT */
void pstDistribRoot(PST,void *,int,void *,int *);

/* PST_ENFORCEPERIODIC */
void pstEnforcePeriodic(PST,void *,int,void *,int *);

/* PST_TREENUMSRCACTIVE */
struct inTreeNumSrcActive {
    uint8_t uRungLo;
    uint8_t uRungHi;
};
void pstTreeNumSrcActive(PST,void *,int,void *,int *);

/* PST_BOUNDSWALK */
struct inBoundsWalk {
    BND bnd;
    uint8_t uRungLo;
    uint8_t uRungHi;
};
struct outBoundsWalk {
    uint64_t nActive;
    uint64_t nContained;
};
void pstBoundsWalk(PST,void *,int,void *,int *);

/* PST_SMOOTH */
struct inSmooth {
    int nSmooth;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    SMF smf;
    };
void pstSmooth(PST,void *,int,void *,int *);

/* PST_RESMOOTH */
void pstReSmooth(PST,void *,int,void *,int *);

/* PST_FASTGASPHASE1 */
void pstFastGasPhase1(PST,void *,int,void *,int *);

/* PST_FASTGASPHASE2 */
void pstFastGasPhase2(PST,void *,int,void *,int *);

/* PST_FASTGASCLEANUP */
void pstFastGasCleanup(PST,void *,int,void *,int *);


/* PST_GRAVITY */
struct inGravity {
    double dTime;
    int nReps;
    int bPeriodic;
    int bEwald;
    double dEwCut;
    double dEwhCut;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
struct outGravity {
    int nActive;
    double dPartSum;
    double dCellSum;
    double dFlop;
    /*
    ** Collected CPU time.
    */
    double dWalkTime;
    /*
    ** Cache Statistics.
    */
    CASTAT cs;
#ifdef INSTRUMENT
    double dComputing;
    double dSynchronizing;
    double dWaiting;
#endif
    };
void pstGravity(PST,void *,int,void *,int *);

/* PST_CALCEANDL */
struct outCalcEandL {
    double T;
    double U;
    double Eth;
    double L[3];
    double F[3];
    double W;
    };
void pstCalcEandL(PST,void *,int,void *,int *);

/* PST_DRIFT */
struct inDrift {
    double dTime;
    double dDelta;
    double dDeltaVPred;
    double dDeltaUPred;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstDrift(PST,void *,int,void *,int *);

/* PST_ROPARTICLECACHE */

void pstROParticleCache(PST, void *, int, void *, int *);

/* PST_PARTICLECACHEFINISH */

void pstParticleCacheFinish(PST, void *, int, void *, int *);

/* PST_CACHEBARRIER */
void pstCacheBarrier(PST, void *, int, void *, int *);

/* PST_STEPVERYACTIVE */
struct inStepVeryActive {
    double dStep;
    double dTime;
    double dDelta;
    int iRung;
    int nMaxRung;
    double diCrit2;
    double aSunInact[3];
    double adSunInact[3];
    double dSunMass;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
struct outStepVeryActive {
    int nMaxRung;
    };
void pstStepVeryActiveKDK(PST,void *,int,void *,int *);

#ifdef HERMITE
/* PST_STEPVERYACTIVEH */
struct inStepVeryActiveH {
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
struct outStepVeryActiveH {
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

/* PST_KICK */
struct inKick {
    double dTime;
    double dDelta;
    double dDeltaVPred;
    double dDeltaU;
    double dDeltaUPred;
    uint8_t uRungLo;
    uint8_t uRungHi;
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

/* PST_PHYSICALSOFT */
struct inPhysicalSoft {
    double dSoftMax;
    double dFac;
    int bSoftMaxMul;
    };
void pstPhysicalSoft(PST,void *,int,void *,int *);

/* PST_SETTOTAL */
struct outSetTotal {
    uint64_t nTotal;
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
    uint8_t uRung;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstSetRung(PST,void *,int,void *,int *);

/* PST_ZERONEWRUNG */
struct inZeroNewRung {
    uint8_t uRung;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstZeroNewRung(PST,void *,int,void *,int *);

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

/* PST_ACCELSTEP */
struct inAccelStep {
    double dEta;
    double dVelFac;
    double dAccFac;
    int    bDoGravity;
    int    bEpsAcc;
    int    bSqrtPhi;
    double dhMinOverSoft;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstAccelStep(PST,void *,int,void *,int *);

/* PST_SPHSTEP */
struct inSphStep {
    double dAccFac;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstSphStep(PST,void *,int,void *,int *);

/* PST_STARFORM */
struct inStarForm
    {
    double dRateCoeff;
    double dTMax;
    double dDenMin;
    double dDelta;
    
    double dTime;
    double dInitStarMass;
    double dESNPerStarMass;
    double dtCoolingShutoff;

    double dtFeedbackDelay;
    double dMassLossPerStarMass;
    double dZMassPerStarMass;
    double dMinGasMass;
    int bdivv;
    };

struct outStarForm 
    {
    int nFormed;
    int nDeleted;
    double dMassFormed;
    };

void pstStarForm(PST,void *,int,void *,int *);

/* PST_DENSITYSTEP */
struct inDensityStep {
    double dEta;
    double dRhoFac;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstDensityStep(PST,void *,int,void *,int *);

/* PST_COOLSETUP */
struct inCoolSetup {
    double dGmPerCcUnit;
    double dComovingGmPerCcUnit;
    double dErgPerGmUnit;
    double dSecUnit;
    double dKpcUnit;

    double dOmega0;
    double dHubble0;
    double dLambda;
    double dOmegab;
    double dOmegaRad;

    double a;
    double z;
    double dTime;
    COOLPARAM CoolParam;
    };

void pstCoolSetup(PST,void *,int,void *,int *);

/* PST_COOLING */
struct inCooling {
    double dTime;	
    double z;
    int bUpdateState;
    int bUpdateTable;
    int bIterateDt;
    int bIsothermal;
    };
struct outCooling {
    double Time;
    double MaxTime;
    double SumTime;
    int nSum;
    };
void pstCooling(PST,void *,int,void *,int *);

/* PST_CORRECTENERGY */
struct inCorrectEnergy {
    double dTuFac;
    double z;
    double dTime;
    double iDirection;
    };
void pstCorrectEnergy(PST, void *,int,void *,int *);

/* PST_GETMAP */
struct inGetMap {
    int nStart;
    };
void pstGetMap(PST,void *,int,void *,int *);

void pstSetRungVeryActive(PST,void *,int,void *,int *);

/* PST_UPDATERUNG */
struct inUpdateRung {
    uint8_t uRungLo;  /* Minimum Rung to modify */
    uint8_t uRungHi;  /* Maximum Rung to modify */
    uint8_t uMinRung; /* Minimum it can be set to */
    uint8_t uMaxRung; /* Maximum it can be set to */
    };
struct outUpdateRung {
    uint64_t nRungCount[256];
    };
void pstUpdateRung(PST,void *,int,void *,int *);

/* PST_ORDWEIGHT */
struct inOrdWeight {
    uint64_t iOrdSplit;
    int iSplitSide;
    int ittr;
    };
struct outOrdWeight {
    uint64_t nLow;
    uint64_t nHigh;
    };
void pstOrdWeight(PST,void *,int,void *,int *);

/* PST_SETWRITESTART */
struct inSetWriteStart {
    uint64_t nWriteStart;
    };
void pstSetWriteStart(PST,void *,int,void *,int *);

/* PST_ADDWRITESTART */
struct inAddWriteStart {
    uint64_t nWriteStart;
    };
void pstAddWriteStart(PST,void *,int,void *,int *);

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

/* PST_GETNPARTS */
/* see pkd.h
 struct outGetNParts { 
	int n;
    int nGas;
    int nDark;
    int nStar;
    int iMaxOrderGas;
    int iMaxOrderDark;
    int iMaxOrderStar;
    };
*/
void pstGetNParts(PST, void *, int, void *, int *);

/* PST_SETNPARTS */
struct inSetNParts {
    uint64_t nGas;
    uint64_t nDark;
    uint64_t nStar;
    uint64_t nMaxOrderGas;
    uint64_t nMaxOrderDark;
    uint64_t nMaxOrder;
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
struct inClearTimer {
    int iTimer;
    };

void pstClearTimer(PST,void *,int,void *,int *);

struct inFof {
    int nSmooth;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    int iCenterType;
    SMF smf;
    };
struct inGroupMerge {
    int bPeriodic;
    int iCenterType;
    SMF smf;
  };
struct inGroupProfiles {
    int nSmooth;
    int bPeriodic;
    int nTotalGroups;
    int bSymmetric;
    int iSmoothType;
    SMF smf;
    };
void pstFof(PST,void *,int,void *,int *);
void pstGroupMerge(PST,void *,int,void *,int *);
void pstGroupProfiles(PST,void *,int,void *,int *);
void pstInitRelaxation(PST,void *,int,void *,int *);

#ifdef PLANETS
/* PLANETS begin */

/* PST_WRITESS */
struct inWriteSS {
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstWriteSS(PST,void *,int,void *,int *);

/* PST_READSS */
struct inReadSS {
    uint64_t nFileStart;
    uint64_t nFileEnd;
    uint64_t nDark;
    uint64_t nGas;			/* always zero */
    uint64_t nStar;			/* always zero */
    uint64_t iOrder;
    int nBucket;
    float fExtraStore;
    FLOAT fPeriod[3];	/* for compatability */
    char achInFile[PST_FILENAME_SIZE];
    double dSunMass;
    };
void pstReadSS(PST,void *,int,void *,int *);

/* PST_SUNINDIRECT */
struct inSunIndirect {
    int iFlag;
    };
struct outSunIndirect {
    double aSun[3];
    double adSun[3];
    };
void pstSunIndirect(PST,void *,int,void *,int *);

/* PST_GRAVSUN */
struct inGravSun {
    double aSun[3];
    double adSun[3];
    double dSunMass;
    };
void pstGravSun(PST,void *,int,void *,int *);

/* PST_HANDSUNMASS */
struct inHandSunMass {
    double dSunMass;
    };
void pstHandSunMass(PST,void *,int,void *,int *);


/* PST_NEXTCOLLISION */
struct outNextCollision {
    double dt;
    uint64_t iOrder1,iOrder2;
    };
void pstNextCollision(PST,void *,int,void *,int *);

/* PST_GETCOLLIDERINFO */
struct inGetColliderInfo {
    uint64_t iOrder;
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
struct inStepVeryActiveS {
    double dStep;
    double dTime;
    double dDelta;
    int iRung;
    int nMaxRung;
    double diCrit2;
    double dSunMass;
    };
struct outStepVeryActiveS {
    int nMaxRung;
    };

void pstStepVeryActiveSymba(PST,void *,int,void *,int *);

/* PST_DRMINTORUNG */
struct inDrminToRung {
    int iRung;
    int iMaxRung;
    };
struct outDrminToRung {
    uint64_t nRungCount[256];
    };
void pstDrminToRung(PST,void *,int,void *,int *);

/* PST_MOMSUN */
struct outMomSun {
    double momSun[3];
    };

void pstMomSun(PST,void *,int,void *,int *);

/* PST_DRIFTSUN */
struct inDriftSun {
    double vSun[3];
    double dDelta;
    };
void pstDriftSun(PST,void *,int,void *,int *);

/* PST_KEPLERDRIFT */
struct inKeplerDrift {
    double dDelta;
    double dSunMass;
    };

void pstKeplerDrift(PST,void *,int,void *,int *);

#endif /* SYMBA */
#endif /* PLANETS */

#ifdef USE_GRAFIC
/* PST_GENERATEIC */
struct inGenerateIC {
    double h;
    double dBoxSize;
    double omegac;
    double omegab;
    double omegav;
    int iSeed;
    int nGrid;
    int nBucket;
    int bComove;
    float fExtraStore;
    FLOAT fPeriod[3];
    };
struct outGenerateIC {
    double dExpansion;
    };

void pstGenerateIC(PST,void *,int,void *,int *);
#endif

/* PST_HOSTNAME */
struct outHostname {
    int  iMpiID;
    char szHostname[20];
    };
void pstHostname(PST,void *,int,void *,int *);

/* PST_MEMSTATUS */
struct outMemStatus {
#ifdef __linux__
    uint64_t minflt;
    uint64_t majflt;
    uint64_t vsize;
    uint64_t rss;
#endif
    uint64_t nCheck;
    };
void pstMemStatus(PST,void *,int,void *,int *);

/* PST_GETCLASSES - Output PARTCLASS[] */
void pstGetClasses(PST,void *,int,void *,int *);

/* PST_SETCLASSES - Input PARTCLASS[] */
void pstSetClasses(PST,void *,int,void *,int *);

/* PST_SWAPCLASSES - Input PARTCLASS[] - Output PARTCLASS[] */
void pstSwapClasses(PST,void *,int,void *,int *);


/* PST_SELSRCALL */
void pstSelSrcAll(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstAll(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCGAS */
void pstSelSrcGas(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstGas(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCSTAR */
struct inSelDstStar {
    int bFB;
    double dTimeFB;
    };
void pstSelSrcStar(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstStar(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCDELETED */
void pstSelSrcDeleted(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstDeleted(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCMASS */
struct inSelMass {
    double dMinMass;
    double dMaxMass;
    int setIfTrue;
    int clearIfFalse;
    };
struct outSelMass {
    uint64_t nSelected;
    };
void pstSelSrcMass(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstMass(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCBYID */
struct inSelById {
    double idStart;
    double idEnd;
    int setIfTrue;
    int clearIfFalse;
    };
struct outSelById {
    uint64_t nSelected;
    };
void pstSelSrcById(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstById(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCPHASEDENSITY */
struct inSelPhaseDensity {
    double dMinDensity;
    double dMaxDensity;
    int setIfTrue;
    int clearIfFalse;
    };
struct outSelPhaseDensity {
    uint64_t nSelected;
    };
void pstSelSrcPhaseDensity(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstPhaseDensity(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCBOX */
struct inSelBox {
    double dCenter[3];
    double dSize[3];
    int setIfTrue;
    int clearIfFalse;
    };
struct outSelBox {
    uint64_t nSelected;
    };
void pstSelSrcBox(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstBox(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCSPHERE */
struct inSelSphere {
    double r[3];
    double dRadius;
    int setIfTrue;
    int clearIfFalse;
    };
struct outSelSphere {
    uint64_t nSelected;
    };
void pstSelSrcSphere(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstSphere(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SELSRCCYLINDER */
struct inSelCylinder {
    double dP1[3];
    double dP2[3];
    double dRadius;
    int setIfTrue;
    int clearIfFalse;
    };
struct outSelCylinder {
    uint64_t nSelected;
    };
void pstSelSrcCylinder(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstCylinder(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_DEEPESTPOT - Input inDeepestPot - Output outDeepestPot */
struct inDeepestPot {
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
struct outDeepestPot {
    double   r[3];
    uint64_t nChecked;
    float    fPot;
    };
void pstDeepestPot(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_PROFILE */
#define PST_MAX_PROFILE_BINS 1000000
struct inProfile {
    double dCenter[3];
    double com[3];
    double vcm[3];
    double L[3];
    uint32_t nBins;
    uint8_t uRungLo;
    uint8_t uRungHi;
    double dRadii[PST_MAX_PROFILE_BINS];
    };
void pstProfile(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_CALCDISTANCE */
struct inCalcDistance {
    double dCenter[3];
    double dRadius;
    };
void pstCalcDistance(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_CALCCOM */
struct inCalcCOM {
    double dCenter[3];
    double dRadius;
    };
struct outCalcCOM {
    double com[3];
    double vcm[3];
    double L[3];
    double M;
    uint64_t N;
    };
void pstCalcCOM(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_COUNTDISTANCE */
struct inCountDistance {
    double dRadius2Inner;
    double dRadius2Outer;
    };
struct outCountDistance {
    uint64_t nCount;
    };
void pstCountDistance(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_PEAKVC */
#define PST_MAX_PEAKVC 1000000
struct inPeakVc {
    double dCenter[3];
    int iGroup;
    short iProcessor;
    };
struct outPeakVc {
    double dPeakVc;
    };
void pstPeakVc(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_INITGRID */
struct inInitGrid {
    int n1, n2, n3, a1;
    int s, n;
    };
void pstInitGrid(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_GRIDPROJECT */
struct inGridProject {
    double r[3];
    };
void pstGridProject(PST pst,void *vin,int nIn,void *vout,int *pnOut);

#ifdef MDL_FFTW
#define PST_MAX_K 4096
/* PST_MEASUREPK */
struct inMeasurePk {
    double dCenter[3];
    double dRadius;
    int nGrid;
    };
struct outMeasurePk {
    float fPower[PST_MAX_K];
    int   nPower[PST_MAX_K];
    };
void pstMeasurePk(PST pst,void *vin,int nIn,void *vout,int *pnOut);
#endif

/* PST_TOTALMASS */
struct outTotalMass {
    double dMass;
    };
void pstTotalMass(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_PSD */
struct inPSD {
    int nSmooth;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    SMF smf;
    };
struct outPSD {
    int dummy;
    };
void pstBuildPsdTree(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSD(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_PSFOF */
void pstPsFof(PST pst,void *vin,int nIn,void *vout,int *pnOut);

#endif
