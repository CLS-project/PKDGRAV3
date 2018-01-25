#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "pkd.h"
#include "mdl.h"
#include "smoothfcn.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "moments.h"
#include "outtype.h"
#include "output.h"

#include "parameters.h"
#include "cosmo.h"
#include "ic.h"

typedef struct lclBlock {
    PKD	pkd;
    int iWtFrom;
    int iWtTo;
    int iPart;
    uint64_t iOrdSplit;
    double fSplit;
    double fWtLow;
    double fWtHigh;
    double fLow;
    double fHigh;
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
    double fSplit;
    double fSplitInactive;
    uint64_t nTotal;
    int iVASplitSide;
    uint64_t nLowerStore;
    uint64_t nUpperStore;
    uint64_t nLowTot;    /* total number of particles in the lower subset of processors */
    uint64_t nHighTot;   /* total number of particles in the upper subset of processors */

    /*
    ** Analysis information. Sticking this stuff in the PST is perhaps a bad idea,
    ** but in avoids order(nThread) storage on the master.
    */
    uint64_t nGroupsLower;

    } * PST;


#define PST_SERVICES		100
#define PST_FILENAME_SIZE	512
#define PST_MAX_FILES           16384

enum pst_service {
    PST_SRV_STOP,
    PST_SPLITIO,
    PST_SETADD,
    PST_READFILE,
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
    PST_OUTPUT,
    PST_OUTPUT_SEND,
    PST_CHECKPOINT,
    PST_RESTORE,
    PST_BUILDTREE,
    PST_DISTRIBTOPTREE,
    PST_DUMPTREES,
    PST_TREEINITMARKED,
    PST_CALCROOT,
    PST_DISTRIBROOT,
    PST_ENFORCEPERIODIC,
    PST_SMOOTH,
    PST_FASTGASPHASE1,
    PST_FASTGASPHASE2,
    PST_FASTGASCLEANUP,
    PST_GRAVITY,
    PST_GRAVEXTERNAL,
    PST_LIGHTCONE,
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
    PST_COUNTRUNGS,
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
    PST_NEW_FOF,
    PST_FOF_PHASES,
    PST_FOF_FINISH_UP,
    PST_FOFMERGE,
    PST_HOP_LINK,
    PST_HOP_JOIN,
    PST_HOP_FINISH_UP,
    PST_HOP_TREE_BUILD,
    PST_HOP_GRAVITY,
    PST_HOP_UNBIND,
    PST_GROUP_RELOCATE,
    PST_GROUP_COUNT_GID,
    PST_GROUP_ASSIGN_GID,
    PST_GROUP_STATS,
    PST_HOP_SEND_STATS,
    PST_GROUPMERGE,
    PST_GROUPPROFILES,
    PST_INITRELAXATION,
    PST_CLEARTIMER,
    PST_INITIALIZEPSTORE,
    PST_GETFFTMAXSIZES,
    PST_GENERATEIC,
    PLT_GENERATEIC,
    PST_MOVEIC,
    PLT_MOVEIC,
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
    PST_LIGHTCONE_OPEN,
    PST_LIGHTCONE_CLOSE,
    PST_INFLATE,
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

/* PST_INITIALIZEPSTORE */
struct inInitializePStore {
    uint64_t mMemoryModel;
    uint64_t nStore;
    uint64_t nSpecies[FIO_SPECIES_LAST];
    double fPeriod[3];
    uint64_t nMinEphemeral;
    uint64_t nMinTotalStore;
    int nEphemeralBytes;
    int nTreeBitsLo;
    int nTreeBitsHi;
    int iCacheSize;
    int iWorkQueueSize;
    int iCUDAQueueSize;
    int bLightCone;
    int bLightConeParticles;
    };
void pstInitializePStore(PST,void *,int,void *,int *);


/* PST_READFILE */
struct inReadFile {
    uint64_t nNodeStart; /* First particle to read (of total) */
    uint64_t nNodeEnd;   /* Last particle to read (of total) */
    double dvFac;
    double dTuFac;
    double dOmega0;
    double dOmegab;
    int nProcessors;
    /*char achFilename[PST_FILENAME_SIZE];*/
    };
typedef char inReadFileFilename[PST_FILENAME_SIZE];
void pstReadFile(PST,void *,int,void *,int *);

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
    double fSplit;
    int iSplitSide;
    int ittr;
    int pFlag;
    };
struct outWeight {
    uint64_t nLow;
    uint64_t nHigh;
    double fLow;
    double fHigh;
    };
void pstWeight(PST,void *,int,void *,int *);

/* PST_COUNTVA */
struct inCountVA {
    int iSplitDim;
    double fSplit;
    };
struct outCountVA {
    int nLow;
    int nHigh;
    };
void pstCountVA(PST,void *,int,void *,int *);

/* PST_WEIGHTWRAP */
struct inWeightWrap {
    int iSplitDim;
    double fSplit;
    double fSplit2;
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
    uint64_t iMinOrder;
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

/* PST_RESTORE */
struct inRestore {
    int nProcessors;
    char achInFile[PST_FILENAME_SIZE];
    };
void pstRestore(PST,void *,int,void *,int *);

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

/* PST_CHECKPOINT */
void pstCheckpoint(PST,void *,int,void *,int *);

struct inOutput {
    int iProcessor;   /* Output number: 0 to nParaWrite */
    int nProcessor;   /* Number of processors left in parallel */
    int iPartner;     /* Who to send the data to */
    int nPartner;     /* How many partners there are */
    outType eOutputType;  /* What kind of output */
    char achOutFile[PST_FILENAME_SIZE];
    };

/* PST_OUTPUT */
void pstOutput(PST,void *,int,void *,int *);

/* PST_OUTPUT_SEND */
void pstOutputSend(PST,void *,int,void *,int *);

/* PST_BUILDTREE */
struct inBuildTree {
    int nBucket;      /* Bucket Size */
    int nGroup;       /* Group Size */
    uint32_t uRoot;   /* Which root node to use */
    uint32_t utRoot;  /* Template tree */
    };
void pstBuildTree(PST,void *,int,void *,int *);

/* PST_DISTRIBTOPTREE */
struct inDistribTopTree {
    uint32_t uRoot; /* Which root node to use */
    uint32_t nTop;
    };
void pstDistribTopTree(PST,void *,int,void *,int *);

/* PST_DUMPTREES */
struct inDumpTrees {
    int bOnlyVA;
    uint8_t uRungDD; /* Domain DD was done on this rung */
    };
void pstDumpTrees(PST,void *,int,void *,int *);

/* PST_TREEINITMARKED */
void pstTreeInitMarked(PST,void *,int,void *,int *);

/* PST_CALCROOT */
struct inCalcRoot {
    double com[3];
    uint32_t uRoot;
    };
struct outCalcRoot {
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
    int nGroup;
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

/* PST_GROUP_RELOCATE */
void pstGroupRelocate(PST,void *,int,void *,int *);

/* PST_GROUP_COUNT_GID */
struct outGroupCountGID {
    uint64_t nGroups;
    };
void pstGroupCountGID(PST,void *,int,void *,int *);

/* PST_GROUP_ASSIGN_GID */
struct inGroupAssignGID {
    uint64_t iStartGID;
    };
void pstGroupAssignGID(PST,void *,int,void *,int *);

/* PST_GROUP_STATS */
struct inGroupStats {
    int bPeriodic;
    double dPeriod[3];
    double rEnvironment[2];
    };
void pstGroupStats(PST,void *,int,void *,int *);


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
    int bKickClose;
    int bKickOpen;
    double dAccFac;
    vel_t dtClose[IRUNGMAX+1];
    vel_t dtOpen[IRUNGMAX+1];
    double dtLCDrift[IRUNGMAX+1];
    double dtLCKick[IRUNGMAX+1];
    double dLookbackFac;
    double dLookbackFacLCP;
    };


typedef struct StatsCollector {
    double dSum;
    double dSum2;
    double dMax;
    int idMax;
    int n;
    } STAT;

/*
** The outGravityReduct structure is at the beginning of the output message, 
** followed by number of threads times the outGravityPerProc structure.
*/
struct outGravityReduct {
    STAT sLocal;
    STAT sActive;
    STAT sPart;
    STAT sPartNumAccess;
    STAT sPartMissRatio;
    STAT sCell;
    STAT sCellNumAccess;
    STAT sCellMissRatio;
    STAT sFlop;
#ifdef INSTRUMENT
    STAT sComputing;
    STAT sWaiting;
    STAT sSynchronizing;
#endif
#ifdef __linux__
    STAT sRSS;
    STAT sFreeMemory;
#endif
    double dFlopSingleCPU;
    double dFlopDoubleCPU;
    double dFlopSingleGPU;
    double dFlopDoubleGPU;
    uint64_t nActive;
    uint64_t nRung[IRUNGMAX+1];
    };
struct outGravityPerProc {
    /*
    ** Collected CPU time.
    */
    double dWalkTime;
    };
void pstGravity(PST,void *,int,void *,int *);

/* PST_LIGHTCONE */
struct inLightCone {
    double dtLCDrift[IRUNGMAX+1];
    double dtLCKick[IRUNGMAX+1];
    double dLookbackFac;
    double dLookbackFacLCP;
    uint8_t uRungLo;
    uint8_t uRungHi;
    };
void pstLightCone(PST,void *,int,void *,int *);

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
    int iRoot;
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
    struct csmVariables cosmo;
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

/* PST_COUNTRUNGS */
struct outCountRungs {
    uint64_t nRungs[MAX_RUNG+1];
    };
void pstCountRungs(PST,void *,int,void *,int *);

/* PST_ACCELSTEP */
struct inAccelStep {
    double dEta;
    double dVelFac;
    double dAccFac;
    int    bDoGravity;
    int    bEpsAcc;
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
void pstFof(PST,void *,int,void *,int *);

/* PST_NEW_FOF */
struct inNewFof {
    double dTau2;
    int nMinMembers;
    };
void pstNewFof(PST,void *,int,void *,int *);

/* PST_FOF_PHASES */
struct outFofPhases {
    int bMadeProgress;
    };   
void pstFofPhases(PST,void *,int,void *,int *);

/* PST_FOF_FINISH_UP */
struct inFofFinishUp{
    int nMinGroupSize;
    };
void pstFofFinishUp(PST,void *,int,void *,int *);


struct inGroupMerge {
    int bPeriodic;
    int iCenterType;
    SMF smf;
  };
void pstGroupMerge(PST,void *,int,void *,int *);

struct inGroupProfiles {
    int nSmooth;
    int bPeriodic;
    int nTotalGroups;
    int bSymmetric;
    int iSmoothType;
    SMF smf;
    };
void pstGroupProfiles(PST,void *,int,void *,int *);
void pstInitRelaxation(PST,void *,int,void *,int *);

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
    uint64_t nPerNode;
    double dBoxSize;
//    double omegam;
//    double omegab;
//    double omegav;
//    double sigma8;
//    double normalization;
//    double spectral;
    double dExpansion;
    int iSeed;
    int bFixed;
    float fPhase;
    int nGrid;
    int b2LPT;
    int bComove;
    int nTf;
    int nInflateFactor;
    struct csmVariables cosmo;
    double k[MAX_TF];
    double tf[MAX_TF];
    };
struct outGenerateIC {
    double dExpansion;
    uint64_t N;
    double noiseMean;
    double noiseCSQ;
    };
void pstGenerateIC(PST,void *,int,void *,int *);

struct inGenerateICthread {
    struct inGenerateIC *ic;
    MDLFFT fft;
    };
void pltGenerateIC(PST,void *,int,void *,int *);

/* PLT_MOVEIC */
struct inMoveIC {
    overlayedParticle *pBase;
    uint64_t iStart;
    uint64_t nMove;
    float fMass;
    float fSoft;
    int nGrid;
    int nInflateFactor;
    };
void pltMoveIC(PST,void *,int,void *,int *);
void pstMoveIC(PST,void *,int,void *,int *);
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
    uint64_t rss;
    uint64_t freeMemory;
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
#define PST_MAX_K_BINS 2500
/* PST_MEASUREPK */
struct inMeasurePk {
    double dTotalMass;
    int nGrid;
    int nBins;
    };
struct outMeasurePk {
    double fK[PST_MAX_K_BINS];
    double fPower[PST_MAX_K_BINS];
    uint64_t nPower[PST_MAX_K_BINS];
    };
void pstMeasurePk(PST pst,void *vin,int nIn,void *vout,int *pnOut);
#endif

/* PST_TOTALMASS */
struct outTotalMass {
    double dMass;
    };
void pstTotalMass(PST pst,void *vin,int nIn,void *vout,int *pnOut);

struct inLightConeOpen {
    int nSideHealpix;
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstLightConeOpen(PST pst,void *vin,int nIn,void *vout,int *pnOut);
struct inLightConeClose {
    char achOutFile[PST_FILENAME_SIZE];
    };
void pstLightConeClose(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_INFLATE */
struct inInflate {
    int nInflateReps;
    };
void pstInflate(PST pst,void *vin,int nIn,void *vout,int *pnOut);


#endif
