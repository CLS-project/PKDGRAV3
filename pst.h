#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "pkd.h"
#include "mdl.h"
#include "psd.h"
#include "smoothfcn.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "moments.h"
#include "outtype.h"

#include "parameters.h"
#include "cosmo.h"
#include "ic.h"

typedef struct lclBlock {
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
    uint64_t nLowerStore;
    uint64_t nUpperStore;
    uint64_t nLowTot;    /* total number of particles in the lower subset of processors */
    uint64_t nHighTot;   /* total number of particles in the upper subset of processors */
    int ioIndex;
    } * PST;


#define PST_SERVICES		100
#define PST_FILENAME_SIZE	512
#define PST_MAX_FILES           16384

enum pst_service {
    PST_SRV_STOP,
    PST_SPLITIO,
    PST_SETADD,
    PST_READFILE,
#ifdef MPI_VERSION
    PST_ORB_BEGIN,
    PST_ORB_SELECT_RUNG,
    PST_ORB_UPDATE_RUNG,
    PST_ORB_DECOMP,
    PST_ORB_FINISH,
    PST_ORB_ROOT_FIND,
    PST_ORB_SPLIT,
    PST_PEANOHILBERTDECOMP,
    PST_RUNGORDER,
#endif
    PST_DOMAINDECOMP,
    PST_CALCBOUND,
    PST_CALCVBOUND,
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
    PST_WRITE,
    PST_BUILDTREE,
    PST_DUMPTREES,
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
    PST_SCALEVEL,
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
    PST_KICKTREE,
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
    PST_SETRUNGVERYACTIVE,
    PST_MARKSMOOTH,
    PST_RESMOOTH,
    PST_INITACCEL,
    PST_UPDATERUNG,
    PST_UPDATERUNGBYTREE,
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
    PST_HOP_LINK,
    PST_HOP_JOIN,
    PST_HOP_FINISH_UP,
    PST_HOP_TREE_BUILD,
    PST_HOP_GRAVITY,
    PST_HOP_UNBIND,
    PST_HOP_ASSIGN_GID,
    PST_HOP_SEND_STATS,
    PST_GROUPMERGE,
    PST_GROUPPROFILES,
    PST_INITRELAXATION,
    PST_CLEARTIMER,
    PST_INITIALIZEPSTORE,
    PST_GETFFTMAXSIZES,
    PST_GENERATEIC,
    PST_CONSTRUCTIC,
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
    PST_SELSRCGROUP,
    PST_SELDSTGROUP,

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
    PST_INITGRID,
    PST_GRIDPROJECT,
#ifdef MDL_FFTW
    PST_MEASUREPK,
#endif
    PST_TOTALMASS,
    PST_BUILDPSDTREE,
    PST_PSD,
    PST_PSDLINK,
    PST_PSD_INIT,
    PST_PSD_JOINBRIDGES,
    PST_PSD_CLG,
    PST_PSD_ASSIGN_GLOBAL_IDS,
    PST_PSD_MERGENOISYGROUPS,
    PST_PSD_JOINGROUPBRIDGES,
    PST_PSD_SETGLOBALID,
    PST_PSD_FINISH,
    PST_WRITE_PSGROUPS,
    PST_UNBIND,
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
    int nGroup;
    int nDomainRungs;
    int iCacheSize;
    int iWorkQueueSize;
    int iCUDAQueueSize;
    int nProcessors;
    /*char achFilename[PST_FILENAME_SIZE];*/
    };
typedef char inReadFileFilename[PST_FILENAME_SIZE];
void pstReadFile(PST,void *,int,void *,int *);

#ifdef MPI_VERSION

/* PST_PEANOHILBERTDECOMP */
struct inPeanoHilbertDecomp {
    int     nRungs;
    int     iMethod;
    };
struct outPeanoHilbertDecomp {
    int x;
    };
void pstPeanoHilbertDecomp(PST,void *,int,void *,int *);

/* PST_ORB_BEGIN,PST_ORB_FINISH */
struct inOrbBegin {
    int nRungs;
    };
void pstOrbBegin(PST,void *,int,void *,int *);
void pstOrbFinish(PST,void *,int,void *,int *);

/* PST_SELECT_RUNG */
struct inOrbSelectRung {
    int iRung;
    };
struct outOrbSelectRung {
    uint64_t nActive;
    };
void pstOrbSelectRung(PST,void *,int,void *,int *);

/* PST_UPDATE_RUNG */
void pstOrbUpdateRung(PST,void *,int,void *,int *);

/* PST_ORB_DECOMP */
struct inOrbDecomp {
    BND bnd;
    };
void pstOrbDecomp(PST,void *,int,void *,int *);


/* PST_ORB_ROOT_FIND */
struct inOrbRootFind {
    BND     bnd;
    double  dFraction;
    double  dReserveFraction; /* [0,1) */
    uint64_t nLower, nUpper;
    };
struct outOrbRootFind {
    double dSplit;
    int    iDim;
    int    nDomains;
    };
void pstOrbRootFind(PST,void *,int,void *,int *);

/* PST_ORB_SPLIT */
void pstOrbSplit(PST,void *,int,void *,int *);

/* PST_RUNGORDER */
struct inRungOrder {
    int     iRung;
    };
struct outRungOrder {
    total_t nMoved;
    BND     bnd;
    };
void pstRungOrder(PST,void *,int,void *,int *);

#endif

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

/* PST_CALCVBOUND */
void pstCalcVBound(PST,void *,int,void *,int *);

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

/* PST_WRITE */
struct inWrite {
    BND bnd;
    double dTime;
    double dEcosmo;
    double dTimeOld;
    double dUOld;
    double dvFac;
    double dBoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    uint64_t nSph;
    uint64_t nDark;
    uint64_t nStar;
    int bStandard;
    int iIndex;
    int nProcessors;
    int iLower, iUpper;
    int bHDF5;
    int mFlags;
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstWrite(PST,void *,int,void *,int *);

/* PST_BUILDTREE */
struct inBuildTree {
    int nBucket;
    int nTrees;
    }; /* followed by an array of TREESPEC */
void pstBuildTree(PST,void *,int,void *,int *);

/* PST_DUMPTREES */
void pstDumpTrees(PST,void *,int,void *,int *);

/* PST_CALCROOT */
struct ioCalcRoot {
    MOMC momc;
    };
void pstCalcRoot(PST,void *,int,void *,int *);

/* PST_DISTRIBROOT */
struct ioDistribRoot {
    double r[3];
    MOMC momc;
    };
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

/* PST_HOP_LINK */
struct inHopLink {
    double dHopTau;
    int nSmooth;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    SMF smf;
    };
void pstHopLink(PST,void *,int,void *,int *);

/* PST_HOP_JOIN */
struct outHopJoin {
    uint64_t nGroups;
    int bDone;
    };
void pstHopJoin(PST,void *,int,void *,int *);

/* PST_HOP_FINISH_UP */
struct inHopFinishUp{
    double fPeriod[3];
    int bPeriodic;
    int nMinGroupSize;
    };
void pstHopFinishUp(PST,void *,int,void *,int *);

/* PST_HOP_TREE_BUILD */
/* PST_BUILDTREE */
struct inHopTreeBuild {
    int nBucket;
    };
void pstHopTreeBuild(PST,void *,int,void *,int *);

/* PST_HOP_GRAVITY */
struct inHopGravity {
    double dTime;
    double dEwCut;
    double dEwhCut;
    double dThetaMin;
    int bPeriodic;
    int nGroup;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstHopGravity(PST,void *,int,void *,int *);

/* PST_HOP_UNBIND */
struct inHopUnbind {
    double dTime;
    double fPeriod[3];
    int bPeriodic;
    int nMinGroupSize;
    int iIteration;
    };
struct outHopUnbind {
    uint64_t nEvaporated;
    uint64_t nGroups;
    };
void pstHopUnbind(PST,void *,int,void *,int *);

/* PST_HOP_ASSIGN_GID */
void pstHopAssignGID(PST,void *,int,void *,int *);

/* PST_HOP_SEND_STATS */
void pstHopSendStats(PST,void *,int,void *,int *);

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
    double dEwCut;
    double dEwhCut;
    double dThetaMin;
    int nReps;
    int bPeriodic;
    int bEwald;
    int nGroup;
    int iRoot1;
    int iRoot2;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
struct outGravity {
    int nActive;
    int nLocal;
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
    double dDelta;
    double dDeltaVPred;
    double dDeltaUPred;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstDrift(PST,void *,int,void *,int *);

/* PST_DRIFT */
struct inScaleVel {
    double dvFac;
    };
void pstScaleVel(PST,void *,int,void *,int *);

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
    double dThetaMin;
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

/* PST_KICKTREE */
struct inKickTree {
    double dTime;
    double dDelta;
    double dDeltaVPred;
    double dDeltaU;
    double dDeltaUPred;
    int iRoot;
    };
struct outKickTree {
    double Time;
    double MaxTime;
    double SumTime;
    int nSum;
    };
void pstKickTree(PST,void *,int,void *,int *);

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

#ifdef COOLING
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
#endif

/* PST_CORRECTENERGY */
struct inCorrectEnergy {
    double dTuFac;
    double z;
    double dTime;
    int    iDirection;
    };
void pstCorrectEnergy(PST, void *,int,void *,int *);


void pstSetRungVeryActive(PST,void *,int,void *,int *);

/* PST_UPDATERUNG */
struct inUpdateRung {
    uint8_t uRungLo;  /* Minimum Rung to modify */
    uint8_t uRungHi;  /* Maximum Rung to modify */
    uint8_t uMinRung; /* Minimum it can be set to */
    uint8_t uMaxRung; /* Maximum it can be set to */
    };
struct outUpdateRung {
    uint64_t nRungCount[MAX_RUNG];
    };
void pstUpdateRung(PST,void *,int,void *,int *);

/* PST_UPDATE_RUNGBYTREE */
struct inUpdateRungByTree {
    int iRoot;
    uint8_t uMinRung;
    uint8_t uMaxRung;
    };
void pstUpdateRungByTree(PST,void *,int,void *,int *);

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
    };
*/
void pstGetNParts(PST, void *, int, void *, int *);

/* PST_SETNPARTS */
struct inSetNParts {
    uint64_t nGas;
    uint64_t nDark;
    uint64_t nStar;
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

/* PST_INITIALIZEPSTORE */
struct inInitializePStore {
    uint64_t mMemoryModel;
    uint64_t nStore;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    double fPeriod[3];
    int nBucket;
    int nGroup;
    int nTreeBitsLo;
    int nTreeBitsHi;
    int nDomainRungs;
    int iCacheSize;
    int iWorkQueueSize;
    int iCUDAQueueSize;
    };
void pstInitializePStore(PST,void *,int,void *,int *);

#ifdef MDL_FFTW
/* PST_GETFFTMAXSIZES */
struct inGetFFTMaxSizes {
    int nx,ny,nz;
    };
struct outGetFFTMaxSizes {
    uint64_t nMaxLocal;
    int nMaxZ;
    int nMaxY;
    };
void pstGetFFTMaxSizes(PST,void *,int,void *,int *);

/* PST_GENERATEIC */
#define MAX_TF 4096
struct inGenerateIC {
    struct inInitializePStore ps;
    uint64_t nPerNode;
    double h;
    double dBoxSize;
    double omegac;
    double omegab;
    double omegav;
    double sigma8;
    double spectral;
    double dExpansion;
    float fExtraStore;
    int iSeed;
    int nGrid;
    int bComove;
    int nTf;
    double k[MAX_TF];
    double tf[MAX_TF];
    };
struct outGenerateIC {
    double dExpansion;
    };
void pstGenerateIC(PST,void *,int,void *,int *);

/* PST_CONSTRUCTIC */
struct inConstructIC {
    MDLFFT fft;
    gridptr dic[6];
    gridpos *pos;
    /* Threads get parts of each slab */
    uint64_t iBegYr, iEndYr;
    uint64_t iBegZk, iEndZk;
    int iSeed;
    double dBoxSize;
    double omegam;
    double omegav;
    double a;
    };
struct outConstructIC {
    double dExpansion;
    };
void pltConstructIC(PST,void *,int,void *,int *);
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
    uint64_t nBytesTree;
    uint64_t nBytesCl;
    uint64_t nBytesIlp;
    uint64_t nBytesIlc;
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
    uint64_t idStart;
    uint64_t idEnd;
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

/* PST_SECSRCGROUP */
void pstSelSrcGroup(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstSelDstGroup(PST pst,void *vin,int nIn,void *vout,int *pnOut);

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
    int nBucket;
    int iCell;
    int nCell;
    int nSmooth;
    int bPeriodic;
    };
struct outPSD {
    int dummy;
    };

struct inAssignGlobalIds
{
    int offs;
    int count;
};

struct inWritePsGroups
{
    char achOutFile[PST_FILENAME_SIZE];
    int iType;
};

void pstBuildPsdTree(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSD(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDLink(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDInit(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDFinish(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDJoinBridges(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDCountLocalGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDAssignGlobalIds(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDUpdateRemoteGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDUpdateParticleGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDMergeNoisyGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDJoinGroupBridges(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstPSDSetGlobalId(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstWritePsGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut);
void pstUnbind(PST pst,void *vin,int nIn,void *vout,int *pnOut);

#endif
