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

#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "pkd.h"
#include "mdl.h"
#include "smooth/smoothfcn.h"
#include "gravity/moments.h"
#include "io/outtype.h"

#include "cosmo.h"
#include "ic/ic.h"

struct LCL {
    PKD pkd;
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
};

struct pstNode {
    pstNode() = delete;
    explicit pstNode(mdl::mdlClass* mdl, LCL* lcl)
        : pstLower_(nullptr), mdl_(mdl), lcl_(lcl), idUpper_(-1), nLeaves_(1), nLower_(0), nUpper_(0), iLvl_(0), iSplitDim(-1), iOrdSplit(0), fSplit(0.0),
          fSplitInactive(0.0), nTotal(0), nLowerGroups(0) {}

  private:
    pstNode* pstLower_;  // Left child node: found locally
    mdl::mdlClass* mdl_; // MDL class
    LCL* lcl_;           // Variables that are identical on all PST nodes
    int idUpper_;        // ID of the thread with the right child node
    int nLeaves_;        // Number of leaves in this subtree (nLower_ + nUpper_, or 1 for a leaf)
    int nLower_;         // Number of leaves in the lower subtree (0 if leaf)
    int nUpper_;         // Number of leaves in the upper subtree (0 if leaf)
    int iLvl_;           // Depth in the PST tree

  public:
    Bound bnd;
    int iSplitDim;
    uint64_t iOrdSplit;
    double fSplit;
    double fSplitInactive;
    uint64_t nTotal;
    uint64_t nLowerGroups; /* count of number of groups in the lower sub-tree, used for making global group ids. */

    [[nodiscard]] auto OffNode() const { return nLeaves_ > mdl_->Cores(); }
    [[nodiscard]] auto OnNode() const { return nLeaves_ <= mdl_->Cores(); }
    [[nodiscard]] auto AmNode() const { return nLeaves_ == mdl_->Cores(); }
    [[nodiscard]] auto NotNode() const { return nLeaves_ != mdl_->Cores(); }
    [[nodiscard]] auto AmCore() const { return nLeaves_ == 1; }
    [[nodiscard]] auto NotCore() const { return nLeaves_ > 1; }

    [[nodiscard]] auto Self() const { return mdl_->Self(); }
    [[nodiscard]] auto n_lower() const { return nLower_; }
    [[nodiscard]] auto n_upper() const { return nUpper_; }
    [[nodiscard]] auto n_total() const { return nLeaves_; }
    [[nodiscard]] auto lower() const { return pstLower_; }
    [[nodiscard]] auto upper() const { return idUpper_; }
    [[nodiscard]] auto depth() const { return iLvl_; }

    [[nodiscard]] auto mdl() const { return mdl_; }
    [[nodiscard]] auto lcl() const { return lcl_; }

    auto create_lower(int idMiddle, int idUp) {
        assert(pstLower_ == nullptr);
        idUpper_ = idMiddle;
        nLower_ = idMiddle - Self();
        nUpper_ = idUp - idMiddle;
        nLeaves_ = nLower_ + nUpper_;
        pstLower_ = new pstNode(mdl_, lcl_);
        pstLower_->iLvl_ = iLvl_ + 1;
        return pstLower_;
    }
};
using PST = pstNode*;

#define PST_SERVICES 100
#define PST_MAX_FILES 16384

namespace Legacy {
inline constexpr std::string_view ReadFile = "ReadFile";
inline constexpr std::string_view CompressASCII = "CompressASCII";
inline constexpr std::string_view SendParticles = "SendParticles";
inline constexpr std::string_view SendArray = "SendArray";
inline constexpr std::string_view WriteASCII = "WriteASCII";
inline constexpr std::string_view Write = "Write";
inline constexpr std::string_view Output = "Output";
inline constexpr std::string_view OutputSend = "OutputSend";
inline constexpr std::string_view Checkpoint = "Checkpoint";
inline constexpr std::string_view Smooth = "Smooth";
inline constexpr std::string_view Gravity = "Gravity";
inline constexpr std::string_view LightCone = "LightCone";
inline constexpr std::string_view CalcEandL = "CalcEandL";
inline constexpr std::string_view Drift = "Drift";
inline constexpr std::string_view ComputePrimVars = "ComputePrimVars";
inline constexpr std::string_view WakeParticles = "WakeParticles";
inline constexpr std::string_view FluxStats = "FluxStats";
inline constexpr std::string_view GrackleInit = "GrackleInit";
inline constexpr std::string_view CoolingUpdate = "CoolingUpdate";
inline constexpr std::string_view CoolingUpdateZ = "CoolingUpdateZ";
inline constexpr std::string_view CoolingInit = "CoolingInit";
inline constexpr std::string_view CoolingHydReion = "CoolingHydReion";
inline constexpr std::string_view ChemCompInit = "ChemCompInit";
inline constexpr std::string_view BhPlaceSeed = "BhPlaceSeed";
inline constexpr std::string_view BhReposition = "BhReposition";
inline constexpr std::string_view BhInit = "BhInit";
inline constexpr std::string_view BhAccretion = "BhAccretion";
inline constexpr std::string_view ROParticleCache = "ROParticleCache";
inline constexpr std::string_view ParticleCacheFinish = "ParticleCacheFinish";
inline constexpr std::string_view Kick = "Kick";
inline constexpr std::string_view KickTree = "KickTree";
inline constexpr std::string_view PhysicalSoft = "PhysicalSoft";
inline constexpr std::string_view SetTotal = "SetTotal";
inline constexpr std::string_view OneNodeReadInit = "OneNodeReadInit";
inline constexpr std::string_view ActiveOrder = "ActiveOrder";
inline constexpr std::string_view AccelStep = "AccelStep";
inline constexpr std::string_view StarForm = "StarForm";
inline constexpr std::string_view StarFormInit = "StarFormInit";
inline constexpr std::string_view DensityStep = "DensityStep";
inline constexpr std::string_view CorrectEnergy = "CorrectEnergy";
inline constexpr std::string_view ReSmooth = "ReSmooth";
inline constexpr std::string_view ReSmoothNode = "ReSmoothNode";
inline constexpr std::string_view ReorderInNodes = "ReorderInNodes";
inline constexpr std::string_view StellarEvolutionInit = "StellarEvolutionInit";
inline constexpr std::string_view SetWriteStart = "SetWriteStart";
inline constexpr std::string_view AddWriteStart = "AddWriteStart";
inline constexpr std::string_view ColNParts = "ColNParts";
inline constexpr std::string_view NewOrder = "NewOrder";
inline constexpr std::string_view SetNParts = "SetNParts";
inline constexpr std::string_view HopLink = "HopLink";
inline constexpr std::string_view HopJoin = "HopJoin";
inline constexpr std::string_view HopFinishUp = "HopFinishUp";
inline constexpr std::string_view HopTreeBuild = "HopTreeBuild";
inline constexpr std::string_view HopGravity = "HopGravity";
inline constexpr std::string_view HopUnbind = "HopUnbind";
inline constexpr std::string_view GroupRelocate = "GroupRelocate";
inline constexpr std::string_view GroupStats = "GroupStats";
inline constexpr std::string_view InitializePStore = "InitializePStore";
inline constexpr std::string_view GenerateIC = "GenerateIC";
inline constexpr std::string_view LocalGenerateIC = "LocalGenerateIC";
inline constexpr std::string_view MoveIC = "MoveIC";
inline constexpr std::string_view LocalMoveIC = "LocalMoveIC";
inline constexpr std::string_view MemStatus = "MemStatus";
inline constexpr std::string_view GetClasses = "GetClasses";
inline constexpr std::string_view SetClasses = "SetClasses";
inline constexpr std::string_view SwapClasses = "SwapClasses";
inline constexpr std::string_view Profile = "Profile";
inline constexpr std::string_view CalcDistance = "CalcDistance";
inline constexpr std::string_view CalcCOM = "CalcCOM";
inline constexpr std::string_view SetSPHoptions = "SetSPHoptions";
inline constexpr std::string_view InitializeEOS = "InitializeEOS";
inline constexpr std::string_view DensityContrast = "DensityContrast";
inline constexpr std::string_view Interlace = "Interlace";
inline constexpr std::string_view GridBinK = "GridBinK";
inline constexpr std::string_view BispectrumSelect = "BispectrumSelect";
inline constexpr std::string_view BispectrumCalculate = "BispectrumCalculate";
inline constexpr std::string_view GetMinDt = "GetMinDt";
inline constexpr std::string_view SetGlobalDt = "SetGlobalDt";
inline constexpr std::string_view PlaceBHSeed = "PlaceBHSeed";
inline constexpr std::string_view BHInit = "BHInit";
inline constexpr std::string_view BHReposition = "BHReposition";
inline constexpr std::string_view BHAccretion = "BHAccretion";
} // namespace Legacy

// Legacy service IDs. These will be removed in the future.
enum pst_service {
    PST_SRV_STOP = 0, /* service 0 is always STOP and handled by MDL */
    PST_READFILE,
    PST_COMPRESSASCII,
    PST_SENDPARTICLES,
    PST_SENDARRAY,
    PST_WRITEASCII,
    PST_WRITE,
    PST_OUTPUT,
    PST_OUTPUT_SEND,
    PST_CHECKPOINT,
    PST_SMOOTH,
    PST_GRAVITY,
    PST_CALCEANDL,
    PST_DRIFT,
#ifdef GRACKLE
    PST_GRACKLEINIT,
#endif
    PST_CHEMCOMPINIT,
    PST_BH_PLACESEED,
    PST_BH_REPOSITION,
    PST_BH_INIT,
    PST_BH_ACCRETION,
    PST_KICK,
    PST_KICKTREE,
    PST_PHYSICALSOFT,
    PST_SETTOTAL,
    PST_ONENODEREADINIT,
    PST_ACTIVEORDER,
    PST_ACCELSTEP,
    PST_STARFORM,
    PST_STARFORMINIT,
    PST_DENSITYSTEP,
    PST_CORRECTENERGY,
    PST_RESMOOTH,
#ifdef OPTIM_SMOOTH_NODE
    PST_RESMOOTHNODE,
#endif
#ifdef OPTIM_REORDER_IN_NODES
    PST_REORDERINNODES,
#endif
#ifdef STELLAR_EVOLUTION
    PST_STELLAREVOLUTIONINIT,
#endif
    PST_SETWRITESTART,
    PST_ADDWRITESTART,
    PST_COLNPARTS,
    PST_NEWORDER,
    PST_SETNPARTS,
    PST_HOP_LINK,
    PST_HOP_JOIN,
    PST_HOP_FINISH_UP,
    PST_HOP_TREE_BUILD,
    PST_HOP_GRAVITY,
    PST_HOP_UNBIND,
    PST_GROUP_RELOCATE,
    PST_GROUP_STATS,
    PST_INITIALIZEPSTORE,
    PST_GENERATEIC,
    PLT_GENERATEIC,
    PST_MOVEIC,
    PLT_MOVEIC,
    PST_MEMSTATUS,
    PST_GETCLASSES,
    PST_SETCLASSES,
    PST_SWAPCLASSES,
    PST_PROFILE,
    PST_CALCDISTANCE,
    PST_CALCCOM,
    PST_SETSPHOPTIONS,
    PST_INITIALIZEEOS,
    PST_DENSITY_CONTRAST,
    PST_INTERLACE,
    PST_GRID_BIN_K,
    PST_BISPECTRUM_SELECT,
    PST_BISPECTRUM_CALCULATE,
    PST_GETMINDT,
    PST_SETGLOBALDT,
};

void pstAddServices(PST, mdl::mdlClass*);

/* PST_INITIALIZEPSTORE */
struct inInitializePStore {
    uint64_t mMemoryModel;
    uint64_t nStore;
    uint64_t nDecoratorStore;
    fioSpeciesList nSpecies;
    blitz::TinyVector<double, 3> fPeriod;
    uint64_t nMinEphemeral;
    uint64_t nMinTotalStore;
    uint32_t nIntegerFactor;
    uint32_t nBucket;
    int nEphemeralBytes;
    int nTreeBitsLo;
    int nTreeBitsHi;
    int iCacheSize;
    int iCacheMaxInflight;
    int iWorkQueueSize;
    mdl::MPI_RECV_MODE mpi_recv_mode;
    int mpi_recv_count;
    bool mpi_RDMA;
};
struct outInitializePStore {
    int nSizeParticle;
    int nSizeDecorator;
    int nSizeNode;
};
int pstInitializePStore(PST, void*, int, void*, int);

/* PST_READFILE */
struct inReadFile {
    fioSpeciesList nSpeciesStart;
    fioSpeciesList nSpeciesEnd;
    double dvFac;
    double dTuFac;
    double dOmega0;
    double dOmegab;
    int nProcessors;
    /*char achFilename[PATH_MAX];*/
};
typedef char inReadFileFilename[PATH_MAX];
int pstReadFile(PST, void*, int, void*, int);

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
int pstCompressASCII(PST, void*, int, void*, int);

/* PST_WRITEASCII */
struct inWriteASCII {
    uint64_t nFileOffset;
    char achOutFile[PATH_MAX];
};
int pstWriteASCII(PST, void*, int, void*, int);

/* PST_WRITE */
struct inWrite {
    Bound bnd;
    UNITS units;
    double dTime;
    double dExp;
    double dEcosmo;
    double dTimeOld;
    double dUOld;
    double dvFac;
    double dTuFac;
    double dBoxSize;
    double Omega0;
    double OmegaLambda;
    double HubbleParam;
    uint64_t nGas;
    uint64_t nDark;
    uint64_t nStar;
    uint64_t nBH;
    int bStandard;
    int iIndex;
    int nProcessors;
    int nProcessorsTotal;
    int iLower, iUpper;
    int bHDF5;
    int mFlags;
    char achOutFile[PATH_MAX];
};
int pstWrite(PST, void*, int, void*, int);

/* PST_SENDPARTICLES */
int pstSendParticles(PST, void*, int, void*, int);

/* PST_SENDARRAY */
struct inSendArray {
    double dvFac;
    PKD_FIELD field;
    int iTo;
    int iUnitSize;
    int bMarked;
};
int pstSendArray(PST, void*, int, void*, int);

/* PST_CHECKPOINT */
int pstCheckpoint(PST, void*, int, void*, int);

/* PST_OUTPUT_SEND */
int pstOutputSend(PST, void*, int, void*, int);

/* PST_HOP_LINK */
struct inHopLink {
    double dHopTau;
    int nSmooth;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    SMF smf;
};
int pstHopLink(PST, void*, int, void*, int);

/* PST_HOP_JOIN */
struct outHopJoin {
    uint64_t nGroups;
    int bDone;
};
int pstHopJoin(PST, void*, int, void*, int);

/* PST_HOP_FINISH_UP */
struct inHopFinishUp {
    blitz::TinyVector<double, 3> fPeriod;
    int bPeriodic;
    int nMinGroupSize;
};
int pstHopFinishUp(PST, void*, int, void*, int);

/* PST_HOP_TREE_BUILD */
/* PST_BUILDTREE */
struct inHopTreeBuild {
    int nBucket;
    int nGroup;
};
int pstHopTreeBuild(PST, void*, int, void*, int);

/* PST_HOP_GRAVITY */
struct inHopGravity {
    double dTime;
    double dEwCut;
    double dEwhCut;
    double dTheta;
    int bPeriodic;
    uint8_t uRungLo;
    uint8_t uRungHi;
};
int pstHopGravity(PST, void*, int, void*, int);

/* PST_HOP_UNBIND */
struct inHopUnbind {
    double dTime;
    blitz::TinyVector<double, 3> fPeriod;
    int bPeriodic;
    int nMinGroupSize;
    int iIteration;
};
struct outHopUnbind {
    uint64_t nEvaporated;
    uint64_t nGroups;
};
int pstHopUnbind(PST, void*, int, void*, int);

/* PST_GROUP_RELOCATE */
int pstGroupRelocate(PST, void*, int, void*, int);

/* PST_GROUP_STATS */
struct inGroupStats {
    int bPeriodic;
    blitz::TinyVector<double, 3> dPeriod;
    uint64_t iGlobalStart; /* initially set to 1 at the master level */
    double rEnvironment[2];
};
int pstGroupStats(PST, void*, int, void*, int);

/* PST_SMOOTH */
struct inSmooth {
    int nSmooth;
    int bPeriodic;
    int bSymmetric;
    int iSmoothType;
    SMF smf;
};
struct outSmooth {
    uint64_t nSmoothed;
};
int pstSmooth(PST, void*, int, void*, int);

#ifdef GRACKLE
struct inGrackleInit {
    char achCoolingTable[PATH_MAX];
    UNITS units;
    int bComove;
    double dScaleFactor;
};
int pstGrackleInit(PST, void*, int, void*, int);
#endif
int pstChemCompInit(PST, void*, int, void*, int);
struct inPlaceBHSeed {
    double dTime;
    double dScaleFactor;
    double dDenMin;
    double dBHMhaloMin;
    double dTau;
    double dBHSeedMass;
    uint8_t uRungMax;
};
struct outPlaceBHSeed {
    int nBHs;
};
int pstPlaceBHSeed(PST, void*, int, void*, int);
int pstBHInit(PST, void*, int, void*, int);
int pstBHReposition(PST, void*, int, void*, int);
struct inBHAccretion {
    double dScaleFactor;
};
int pstBHAccretion(PST, void*, int, void*, int);

/* PST_RESMOOTH */
int pstReSmooth(PST, void*, int, void*, int);
#ifdef OPTIM_SMOOTH_NODE
int pstReSmoothNode(PST, void*, int, void*, int);
#endif
#ifdef OPTIM_REORDER_IN_NODES
int pstReorderWithinNodes(PST, void*, int, void*, int);
#endif

/* PST_GRAVITY */
struct inGravity {
    double dTime;
    double dEwCut;
    double dEwhCut;
    double dTheta;
    int nReps;
    int bPeriodic;
    int bEwald;
    int bGPU;
    int iRoot1;
    int iRoot2;
    struct pkdKickParameters kick;
    struct pkdLightconeParameters lc;
    struct pkdTimestepParameters ts;
    SPHOptions SPHoptions;
};

using STAT = struct StatsCollector {
    double dSum;
    double dSum2;
    double dMax;
    int idMax;
    int n;
};

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
    uint64_t nRung[IRUNGMAX + 1];
    uint64_t nTilesTotal;
    uint64_t nTilesCPU;
};
int pstGravity(PST, void*, int, void*, int);

/* PST_CALCEANDL */
struct outCalcEandL {
    double T;
    double U;
    double Eth;
    blitz::TinyVector<double, 3> L;
    blitz::TinyVector<double, 3> F;
    double W;
};
int pstCalcEandL(PST, void*, int, void*, int);

/* PST_DRIFT */
struct inDrift {
    double dTime;
    double dDelta;
    double dDeltaVPred;
    double dDeltaUPred;
    int iRoot;
    int bDoGas;
    int bGasEvolveDensity;
};
int pstDrift(PST, void*, int, void*, int);
int pstEndTimestepIntegration(PST, void*, int, void*, int);
int pstWakeParticles(PST, void*, int, void*, int);
struct outFluxStats {
    int nAvoided;
    int nComputed;
};

/* PST_KICK */
struct inKick {
    double dTime;
    double dDelta;
    double dDeltaVPred;
    double dDeltaU;
    double dDeltaUPred;
    int bDoGas;
    uint8_t uRungLo;
    uint8_t uRungHi;
};
int pstKick(PST, void*, int, void*, int);

/* PST_KICKTREE */
struct inKickTree {
    double dTime;
    double dDelta;
    double dDeltaVPred;
    double dDeltaU;
    double dDeltaUPred;
    int iRoot;
};
int pstKickTree(PST, void*, int, void*, int);

/* PST_PHYSICALSOFT */
struct inPhysicalSoft {
    double dSoftMax;
    double dFac;
    int bSoftMaxMul;
};
int pstPhysicalSoft(PST, void*, int, void*, int);

/* PST_SETTOTAL */
struct outSetTotal {
    uint64_t nTotal;
};
int pstSetTotal(PST, void*, int, void*, int);

/* PST_ONENODEREADINIT */
int pstOneNodeReadInit(PST pst, void* vin, int nIn, void* vout, int nOut);

/* PST_ACTIVEORDER */
int pstActiveOrder(PST, void*, int, void*, int);

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

void pstDensCheck(PST, void*, int, void*, int*);

/* PST_ACCELSTEP */
struct inAccelStep {
    double dEta;
    double dVelFac;
    double dAccFac;
    double dDelta;
    int iMaxRung;
    int bDoGravity;
    int bEpsAcc;
    double dhMinOverSoft;
    uint8_t uRungLo;
    uint8_t uRungHi;
};
int pstAccelStep(PST, void*, int, void*, int);

/* PST_DENSITYSTEP */
struct inDensityStep {
    double dDelta;
    double dEta;
    double dRhoFac;
    int iMaxRung;
    uint8_t uRungLo;
    uint8_t uRungHi;
};
int pstDensityStep(PST, void*, int, void*, int);

/* PST_CORRECTENERGY */
struct inCorrectEnergy {
    double dTuFac;
    double z;
    double dTime;
    int iDirection;
};
int pstCorrectEnergy(PST, void*, int, void*, int);

/* PST_SETWRITESTART */
struct inSetWriteStart {
    uint64_t nWriteStart;
};
int pstSetWriteStart(PST, void*, int, void*, int);

/* PST_ADDWRITESTART */
struct inAddWriteStart {
    uint64_t nWriteStart;
};
int pstAddWriteStart(PST, void*, int, void*, int);

/* PST_COLNPARTS */
struct outColNParts {
    int nNew;
    int nDeltaGas;
    int nDeltaDark;
    int nDeltaStar;
    int nDeltaBH;
};
int pstColNParts(PST, void*, int, void*, int);

/* PST_NEWORDER */
int pstNewOrder(PST, void*, int, void*, int);

/* PST_SETNPARTS */
struct inSetNParts {
    uint64_t nGas;
    uint64_t nDark;
    uint64_t nStar;
    uint64_t nBH;
};
int pstSetNParts(PST, void*, int, void*, int);

#ifdef MDL_FFTW
    /* PST_GENERATEIC */
    #define MAX_TF 4096
struct inGenerateIC {
    uint64_t nPerNode;
    double dBoxSize;
    double dBoxMass;
    double dBaryonFraction;
    double dExpansion;
    int iSeed;
    int bFixed;
    float fPhase;
    int nGrid;
    int iLPT;
    int bICgas;
    int nBucket;
    double dInitialT;
    double dInitialH;
    double dInitialHe;
    double dInitialC;
    double dInitialN;
    double dInitialO;
    double dInitialNe;
    double dInitialMg;
    double dInitialSi;
    double dInitialFe;
    double dInitialMetallicity;
    double dTuFac;
    int bComove;
    int nTf;
    double k[MAX_TF];
    double tf[MAX_TF];
};
struct outGenerateIC {
    double dExpansion;
    uint64_t N;
    double noiseMean;
    double noiseCSQ;
};
int pstGenerateIC(PST, void*, int, void*, int);

struct inGenerateICthread {
    struct inGenerateIC* ic;
    MDLFFT fft;
};
int pltGenerateIC(PST, void*, int, void*, int);

/* PLT_MOVEIC */
struct inMoveIC {
    overlayedParticle* pBase;
    uint64_t iStart;
    uint64_t nMove;
    float fMass;
    float fSoft;
    int nGrid;
    int bICgas;
    int nBucket;
    double dInitialT;
    double dInitialH;
    double dInitialHe;
    double dInitialC;
    double dInitialN;
    double dInitialO;
    double dInitialNe;
    double dInitialMg;
    double dInitialSi;
    double dInitialFe;
    double dInitialMetallicity;
    double dExpansion;
    double dBaryonFraction;
    double dTuFac;
};
int pltMoveIC(PST, void*, int, void*, int);
int pstMoveIC(PST, void*, int, void*, int);
#endif

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
int pstMemStatus(PST, void*, int, void*, int);

/* PST_GETCLASSES - Output PARTCLASS[] */
int pstGetClasses(PST, void*, int, void*, int);

/* PST_SETCLASSES - Input PARTCLASS[] */
int pstSetClasses(PST, void*, int, void*, int);

/* PST_SWAPCLASSES - Input PARTCLASS[] - Output PARTCLASS[] */
int pstSwapClasses(PST, void*, int, void*, int);

/* PST_PROFILE */
#define PST_MAX_PROFILE_BINS 1000000
struct inProfile {
    blitz::TinyVector<double, 3> dCenter;
    blitz::TinyVector<double, 3> com;
    blitz::TinyVector<double, 3> vcm;
    blitz::TinyVector<double, 3> L;
    uint32_t nBins;
    uint8_t uRungLo;
    uint8_t uRungHi;
    double dRadii[PST_MAX_PROFILE_BINS];
};
int pstProfile(PST pst, void* vin, int nIn, void* vout, int nOut);

/* PST_CALCDISTANCE */
struct inCalcDistance {
    blitz::TinyVector<double, 3> dCenter;
    double dRadius;
    int bPeriodic;
};
int pstCalcDistance(PST pst, void* vin, int nIn, void* vout, int nOut);

/* PST_CALCCOM */
struct inCalcCOM {
    blitz::TinyVector<double, 3> dCenter;
    double dRadius;
    int bPeriodic;
};
struct outCalcCOM {
    blitz::TinyVector<double, 3> com;
    blitz::TinyVector<double, 3> vcm;
    blitz::TinyVector<double, 3> L;
    double M;
    uint64_t N;
};
int pstCalcCOM(PST pst, void* vin, int nIn, void* vout, int nOut);

/* PST_SETSPHOPTIONS */
struct inSetSPHoptions {
    SPHOptions SPHoptions;
};
int pstSetSPHoptions(PST pst, void* vin, int nIn, void* vout, int nOut);
/* PST_INITIALIZEEOS */
int pstInitializeEOS(PST pst, void* vin, int nIn, void* vout, int nOut);

#ifdef MDL_FFTW
    #define PST_MAX_K_BINS 2500
/* PST_GRID_BIN_K */
struct inGridBinK {
    int nBins;
    int iGrid;
};
struct outGridBinK {
    double fK[PST_MAX_K_BINS];
    double fPower[PST_MAX_K_BINS];
    uint64_t nPower[PST_MAX_K_BINS];
};
int pstGridBinK(PST pst, void* vin, int nIn, void* vout, int nOut);
/* PST_DENSITY_CONTRAST */
struct inDensityContrast {
    int iGrid;
    int k;
    double dTotalMass;
};
int pstDensityContrast(PST pst, void* vin, int nIn, void* vout, int nOut);
/* PST_INTERLACE */
struct inInterlace {
    int iGridTarget;
    int iGridSource;
};
int pstInterlace(PST pst, void* vin, int nIn, void* vout, int nOut);
/* PST_BISPECTRUM_SELECT */
struct inBispectrumSelect {
    int iGridTarget;
    int iGridSource;
    double kmin, kmax;
};
int pstBispectrumSelect(PST pst, void* vin, int nIn, void* vout, int nOut);
/* PST_BISPECTRUM_CALCULATE */
struct inBispectrumCalculate {
    int iGrid1, iGrid2, iGrid3;
};
int pstBispectrumCalculate(PST pst, void* vin, int nIn, void* vout, int nOut);
#endif

/* PST_GETMINDT */
struct outGetMinDt {
    uint8_t uMinDt;
};
int pstGetMinDt(PST pst, void* vin, int nIn, void* vout, int nOut);
int pstSetGlobalDt(PST pst, void* vin, int nIn, void* vout, int nOut);

#endif
