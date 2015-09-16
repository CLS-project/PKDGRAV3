#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include <assert.h>
#include <stdint.h>
#include <string.h>
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#else
#define PRIu64 "llu"
#endif
#ifdef __linux__
#include <unistd.h>
#endif
#ifdef USE_ITT
#include "ittnotify.h"
#endif
#include "mdl.h"
#include "pst.h"
#include "pkd.h"
#include "smooth.h"
#include "hop.h"
#ifdef USE_PSD
#include "psdtree.h"
#include "unbind.h"
#endif

#define pstOffNode(pst) ((pst)->nLeaves > mdlCores((pst)->mdl))
#define pstOnNode(pst) ((pst)->nLeaves <= mdlCores((pst)->mdl))
#define pstAmNode(pst) ((pst)->nLeaves == mdlCores((pst)->mdl))
#define pstNotNode(pst) ((pst)->nLeaves != mdlCores((pst)->mdl))
#define pstAmCore(pst) ((pst)->nLeaves == 1)
#define pstNotCore(pst) ((pst)->nLeaves > 1)

/*
** Order:
**  nLocal  - depends on all local particles, regardless of rung
**  nActive - depends only on active particles
**  1       - constant time
**
** Input:
**  Bcast   - The same data is sent to all processors
**  Scatter - A unique value is sent to each processor
**  Fan Out
**
** Output:
**  Reduce  - The return value is reduced (sum, etc.)
**  Gather  - A unique value is received from all processors
**
**
**  Service               Order   Input   Output |
**  -------               -----   -----   ------ |
**  SetAdd                        yes     -      | fan out
**  DomainDecomp                  yes     -      |
**  CalcBound                     -       Reduce | Custom reduce: BND_COMBINE
**  CombineBound                  -       Reduce | Custom reduce: BND_COMBINE
**  Weight                        Bcast   Reduce | Sum several fields
**  CountVA                       Bcast   Reduce | Sum two fields
**  WeightWrap                    Bcast   Reduce | Sum two fields
**  OrdWeight                     Bcast   Reduce | Sum two fields
**  FreeStore                     -       Reduce | Sum a field
**  ColRejects                    -       Gather |
**  SwapRejects                   Scatter Gather |
**  ColOrdRejects                 Yes     Gather |
**  DomainOrder                   Yes     -      |
**  LocalOrder                    -       -      | Bcast
**  WriteTipsy                    Yes     -      |
**  BuildTree                     Yes     Many   | Multiple cells
**  CalcRoot              nLocal  -       Yes    |
**  DistribRoot           1       Bcast   -      |
**  EnforcePeriodic               Bcast   -      |
**  Smooth                        Yes     -      |
**  FastGasPhase1                 Yes     -      |
**  FastGasPhase2                 Yes     -      |
**  FastGasCleanup                Yes     -      |
**  Gravity                       Yes     Gather |
**  CalcEandL                     -       Reduce |
**  Drift                         Bcast   -      |
**  CacheBarrier                  -       -      |
**  StepVeryActiveKDK             Yes     Yes    |
**  Copy0                         Yes     -      |
**  Predictor                     Yes     -      |
**  Corrector                     Yes     -      |
**  SunCorrector                  Yes     -      |
**  PredictorInactive             Yes     -      |
**  AarsethStep                   Yes     -      |
**  FirstDt                       -       -      |
**  ROParticleCache               -       -      |
**  ParticleCacheFinish           -       -      |
**  Kick                  nLocal  Yes     Yes    |
**  KickTree              nActive Yes     Yes    |
**  SetSoft                       Yes     -      |
**  PhysicalSoft                  Yes     -      |
**  SetTotal                      -       Yes    |
**  SetWriteStart                 Yes     -      |
**  OneNodeReadInit               Yes     Gather |
**  SwapAll                       Yes     -      |
**  ActiveOrder                   -       Yes    |
**  InitStep                      Yes     -      |
**  SetRung                       Yes     -      |
**  ZeroNewRung                   Yes     -      |
**  ActiveRung                    Yes     -      |
**  CurrRung                      Yes     Yes    |
**  DensityStep                   Yes     -      |
**  CoolSetup                     Yes     -      |
**  Cooling                       Yes     Yes    |
**  CorrectEnergy                 Yes     -      |
**  AccelStep                     Yes     -      |
**  SphStep                       Yes     -      |
**  SetRungVeryActive             Yes     -      |
**  ReSmooth                      Yes     -      |
**  UpdateRung                    Yes     Yes    |
**  ColNParts                     -       Gather |
**  NewOrder                      Scatter -      |
**  GetNParts                     -       Gather |
**  SetNParts                     Yes     -      |
**  ClearTimer                    Yes     -      |
**  Fof                           Yes     -      |
**  GroupMerge                    Yes     Yes    |
**  GroupProfiles                 Yes     Yes    |
**  InitRelaxation                -       -      |
**  FindIOS                       Yes     Yes    |
**  StartIO                       Yes     -      |
**  ReadSS                        Yes     -      |
**  WriteSS                       Yes     -      |
**  SunIndirect                   Yes     Yes    |
**  GravSun                       Yes     -      |
**  HandSunMass                   Yes     -      |
**  NextCollision                 -       Yes    |
**  GetColliderInfo               Yes     Yes    |
**  DoCollision                   Yes     Yes    |
**  GetVariableVeryActive         -       Yes    |
**  CheckHelioDist                -       Yes    |
**  StepVeryActiveSymba           Yes     Yes    |
**  DrminToRung                   Yes     Yes    |
**  MomSun                        -       Yes    |
**  DriftSun                      Yes     -      |
**  KeplerDrift                   Yes     -      |
**  GenerateIC                    Yes     Yes    |
**  Hostname                      -       Gather |
**  MemStatus                     -       Gather |
*/

void pstAddServices(PST pst,MDL mdl) {
    int nThreads;

    nThreads = mdlThreads(mdl);

    mdlAddService(mdl,PST_SETADD,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetAdd,
		  sizeof(struct inSetAdd),0);
    mdlAddService(mdl,PST_READFILE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstReadFile,
	          sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE),0);
#ifdef MPI_VERSION
    mdlAddService(mdl,PST_ORB_BEGIN,pst,
	(void (*)(void *,void *,int,void *,int *)) pstOrbBegin,
	sizeof(struct inOrbBegin),0);
    mdlAddService(mdl,PST_ORB_SELECT_RUNG,pst,
	(void (*)(void *,void *,int,void *,int *)) pstOrbSelectRung,
	sizeof(struct inOrbSelectRung),sizeof(struct outOrbSelectRung));
    mdlAddService(mdl,PST_ORB_UPDATE_RUNG,pst,
	(void (*)(void *,void *,int,void *,int *)) pstOrbUpdateRung,0,0);
    mdlAddService(mdl,PST_ORB_FINISH,pst,
	(void (*)(void *,void *,int,void *,int *)) pstOrbFinish,0,0);
    mdlAddService(mdl,PST_ORB_DECOMP,pst,
                  (void (*)(void *,void *,int,void *,int *)) pstOrbDecomp,
                  sizeof(struct inOrbDecomp),0);
    mdlAddService(mdl,PST_ORB_ROOT_FIND,pst,
                  (void (*)(void *,void *,int,void *,int *)) pstOrbRootFind,
                  sizeof(struct inOrbRootFind),sizeof(struct outOrbRootFind));
    mdlAddService(mdl,PST_ORB_SPLIT,pst,
                  (void (*)(void *,void *,int,void *,int *)) pstOrbSplit,
                  sizeof(int),0);
    mdlAddService(mdl,PST_PEANOHILBERTDECOMP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPeanoHilbertDecomp,
		  sizeof(struct inPeanoHilbertDecomp),sizeof(struct outPeanoHilbertDecomp));
    mdlAddService(mdl,PST_RUNGORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstRungOrder,
	          sizeof(struct inRungOrder),sizeof(struct outRungOrder));
#endif
    mdlAddService(mdl,PST_DOMAINDECOMP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDomainDecomp,
		  sizeof(struct inDomainDecomp),0);
    mdlAddService(mdl,PST_CALCBOUND,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcBound,
		  0,sizeof(BND));
    mdlAddService(mdl,PST_CALCVBOUND,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcVBound,
		  0,sizeof(BND));
    mdlAddService(mdl,PST_COMBINEBOUND,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCombineBound,
		  0,sizeof(BND));
    mdlAddService(mdl,PST_WEIGHT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWeight,
		  sizeof(struct inWeight),sizeof(struct outWeight));
    mdlAddService(mdl,PST_COUNTVA,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCountVA,
		  sizeof(struct inCountVA),sizeof(struct outCountVA));
    mdlAddService(mdl,PST_WEIGHTWRAP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWeightWrap,
		  sizeof(struct inWeightWrap),sizeof(struct outWeightWrap));
    mdlAddService(mdl,PST_ORDWEIGHT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOrdWeight,
		  sizeof(struct inOrdWeight),sizeof(struct outOrdWeight));
    mdlAddService(mdl,PST_FREESTORE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFreeStore,
		  0,sizeof(struct outFreeStore));
    mdlAddService(mdl,PST_COLREJECTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstColRejects,
		  0,nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_SWAPREJECTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSwapRejects,
		  nThreads*sizeof(int),nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_COLORDREJECTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstColOrdRejects,
		  sizeof(struct inColOrdRejects),nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_DOMAINORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDomainOrder,
		  sizeof(struct inDomainOrder),0);
    mdlAddService(mdl,PST_LOCALORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstLocalOrder,
		  sizeof(struct inDomainOrder),0);
    mdlAddService(mdl,PST_COMPRESSASCII,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCompressASCII,
		  sizeof(struct inCompressASCII),sizeof(struct outCompressASCII));
    mdlAddService(mdl,PST_WRITEASCII,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWriteASCII,
		  sizeof(struct inWriteASCII),0);
    mdlAddService(mdl,PST_WRITE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWrite,
		  sizeof(struct inWrite),0);
    /*
    ** Calculate the number of levels in the top tree and use it to
    ** define the size of the messages.
    */
    mdlAddService(mdl,PST_BUILDTREE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstBuildTree,
	sizeof(struct inBuildTree)+MAX_RUNG*sizeof(TREESPEC),2*pkdMaxNodeSize());
    mdlAddService(mdl,PST_DUMPTREES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDumpTrees,
		  0,0);
    mdlAddService(mdl,PST_CALCROOT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcRoot,
	          3*sizeof(double),sizeof(struct ioCalcRoot));
    mdlAddService(mdl,PST_DISTRIBROOT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDistribRoot,
		  sizeof(struct ioDistribRoot),0);
    mdlAddService(mdl,PST_ENFORCEPERIODIC,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstEnforcePeriodic,
		  sizeof(BND),0);
    mdlAddService(mdl,PST_HOP_LINK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopLink,
	          sizeof(struct inHopLink),sizeof(uint64_t));
    mdlAddService(mdl,PST_HOP_JOIN,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopJoin,
		  sizeof(struct inHopLink),sizeof(struct outHopJoin));
    mdlAddService(mdl,PST_HOP_FINISH_UP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopFinishUp,
	          sizeof(struct inHopFinishUp),sizeof(uint64_t));
    mdlAddService(mdl,PST_HOP_TREE_BUILD,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopTreeBuild,
	          sizeof(struct inHopTreeBuild),0);
    mdlAddService(mdl,PST_HOP_GRAVITY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopGravity,
	          sizeof(struct inHopGravity),0);
    mdlAddService(mdl,PST_HOP_UNBIND,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopUnbind,
	          sizeof(struct inHopUnbind),sizeof(struct outHopUnbind));
    mdlAddService(mdl,PST_HOP_ASSIGN_GID,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopAssignGID,
		  0,0);
    mdlAddService(mdl,PST_HOP_SEND_STATS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopSendStats,
		  0,0);
    mdlAddService(mdl,PST_SMOOTH,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSmooth,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_FASTGASPHASE1,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFastGasPhase1,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_FASTGASPHASE2,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFastGasPhase2,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_FASTGASCLEANUP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFastGasCleanup,
		  0,0);
    mdlAddService(mdl,PST_GRAVITY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGravity,
		  sizeof(struct inGravity),nThreads*sizeof(struct outGravity));
    mdlAddService(mdl,PST_CALCEANDL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcEandL,
		  0,sizeof(struct outCalcEandL));
    mdlAddService(mdl,PST_DRIFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDrift,
		  sizeof(struct inDrift),0);
    mdlAddService(mdl,PST_SCALEVEL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstScaleVel,
		  sizeof(struct inScaleVel),0);
    mdlAddService(mdl,PST_CACHEBARRIER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCacheBarrier,
		  0,0);
    mdlAddService(mdl,PST_STEPVERYACTIVE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstStepVeryActiveKDK,
		  sizeof(struct inStepVeryActive),
		  sizeof(struct outStepVeryActive));
    mdlAddService(mdl,PST_ROPARTICLECACHE,pst,
		  (void (*)(void *,void *,int,void *,int *))pstROParticleCache,
		  0,0);
    mdlAddService(mdl,PST_PARTICLECACHEFINISH,pst,
		  (void (*)(void *,void *,int,void *,int *))pstParticleCacheFinish,
		  0,0);
    mdlAddService(mdl,PST_KICK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstKick,
		  sizeof(struct inKick),sizeof(struct outKick));
    mdlAddService(mdl,PST_KICKTREE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstKickTree,
		  sizeof(struct inKickTree),sizeof(struct outKickTree));
    mdlAddService(mdl,PST_SETSOFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetSoft,
		  sizeof(struct inSetSoft),0);
    mdlAddService(mdl,PST_PHYSICALSOFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPhysicalSoft,
		  sizeof(struct inPhysicalSoft),0);
    mdlAddService(mdl,PST_SETTOTAL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetTotal,
		  0,sizeof(struct outSetTotal));
    mdlAddService(mdl,PST_SETWRITESTART,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetWriteStart,
		  sizeof(struct inSetWriteStart),0);
    mdlAddService(mdl,PST_ADDWRITESTART,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstAddWriteStart,
		  sizeof(struct inAddWriteStart),0);
    mdlAddService(mdl,PST_ONENODEREADINIT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOneNodeReadInit,
	          sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE),
	          nThreads*sizeof(int));
    mdlAddService(mdl,PST_SWAPALL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSwapAll,
		  sizeof(int),0);
    mdlAddService(mdl,PST_ACTIVEORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveOrder,
		  0,sizeof(uint64_t));
    mdlAddService(mdl,PST_INITSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitStep,
		  sizeof(struct inInitStep),0);
    mdlAddService(mdl,PST_SETRUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetRung,
		  sizeof(struct inSetRung),0);
    mdlAddService(mdl,PST_ZERONEWRUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstZeroNewRung,
		  sizeof(struct inZeroNewRung),0);
    mdlAddService(mdl,PST_ACTIVERUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveRung,
		  sizeof(struct inActiveRung),0);
    mdlAddService(mdl,PST_CURRRUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCurrRung,
		  sizeof(struct inCurrRung),sizeof(struct outCurrRung));
    mdlAddService(mdl,PST_DENSITYSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDensityStep,
		  sizeof(struct inDensityStep),0);
#ifdef COOLING
    mdlAddService(mdl,PST_COOLSETUP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCoolSetup,
		  sizeof(struct inCoolSetup),0);
    mdlAddService(mdl,PST_COOLING,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCooling,
		  sizeof(struct inCooling),sizeof(struct outCooling));
#endif
    mdlAddService(mdl,PST_CORRECTENERGY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCorrectEnergy,
		  sizeof(struct inCorrectEnergy),0);
    mdlAddService(mdl,PST_ACCELSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstAccelStep,
		  sizeof(struct inAccelStep), 0);
    mdlAddService(mdl,PST_SPHSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSphStep,
		  sizeof(struct inSphStep), 0);
    mdlAddService(mdl,PST_STARFORM,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstStarForm,
		  sizeof(struct inStarForm),sizeof(struct outStarForm));
    mdlAddService(mdl,PST_SETRUNGVERYACTIVE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetRungVeryActive,
		  sizeof(struct inSetRung),0);
    mdlAddService(mdl,PST_RESMOOTH,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstReSmooth,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_UPDATERUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstUpdateRung,
		  sizeof(struct inUpdateRung),sizeof(struct outUpdateRung));
    mdlAddService(mdl,PST_UPDATERUNGBYTREE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstUpdateRungByTree,
		  sizeof(struct inUpdateRungByTree),sizeof(struct outUpdateRung));
    mdlAddService(mdl,PST_COLNPARTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstColNParts,
		  0,nThreads*sizeof(struct outColNParts));
    mdlAddService(mdl,PST_NEWORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstNewOrder,
		  nThreads*sizeof(uint64_t),0);
    mdlAddService(mdl,PST_GETNPARTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetNParts,
		  0,sizeof(struct outGetNParts));
    mdlAddService(mdl,PST_SETNPARTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetNParts,
		  sizeof(struct inSetNParts),0);
    mdlAddService(mdl,PST_CLEARTIMER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstClearTimer,
		  sizeof(struct inClearTimer),0);
    mdlAddService(mdl,PST_FOF,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFof,
		  sizeof(struct inFof),0);
    mdlAddService(mdl,PST_GROUPMERGE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGroupMerge,
		  sizeof(struct inGroupMerge),sizeof(int));
    mdlAddService(mdl,PST_GROUPPROFILES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGroupProfiles,
		  sizeof(struct inGroupProfiles),sizeof(int));
    mdlAddService(mdl,PST_INITRELAXATION,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitRelaxation,0,0);
#ifdef MDL_FFTW
    mdlAddService(mdl,PST_INITIALIZEPSTORE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitializePStore,
		  sizeof(struct inInitializePStore),0);
    mdlAddService(mdl,PST_GETFFTMAXSIZES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetFFTMaxSizes,
		  sizeof(struct inGetFFTMaxSizes),sizeof(struct outGetFFTMaxSizes));
    mdlAddService(mdl,PST_GENERATEIC,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGenerateIC,
		  sizeof(struct inGenerateIC),sizeof(struct outGenerateIC));
    mdlAddService(mdl,PST_CONSTRUCTIC,pst,
		  (void (*)(void *,void *,int,void *,int *)) pltConstructIC,
		  sizeof(struct inConstructIC),sizeof(struct outConstructIC));
#endif
    mdlAddService(mdl,PST_HOSTNAME,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHostname,
		  0,nThreads*sizeof(struct outHostname));
    mdlAddService(mdl,PST_MEMSTATUS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstMemStatus,
		  0,nThreads*sizeof(struct outMemStatus));
    mdlAddService(mdl,PST_GETCLASSES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetClasses,
		  0, PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlAddService(mdl,PST_SETCLASSES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetClasses,
		  PKD_MAX_CLASSES*sizeof(PARTCLASS), 0);
    mdlAddService(mdl,PST_SWAPCLASSES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSwapClasses,
		  PKD_MAX_CLASSES*sizeof(PARTCLASS),
		  PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlAddService(mdl,PST_SELSRCALL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcAll,
		  0, 0 );
    mdlAddService(mdl,PST_SELDSTALL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstAll,
		  0, 0 );
    mdlAddService(mdl,PST_SELSRCGAS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcGas,
		  0, 0 );
    mdlAddService(mdl,PST_SELDSTGAS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstGas,
		  0, 0 );
    mdlAddService(mdl,PST_SELSRCSTAR,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcStar,
		  0, 0 );
    mdlAddService(mdl,PST_SELDSTSTAR,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstStar,
		  sizeof(struct inSelDstStar), 0 );
    mdlAddService(mdl,PST_SELSRCDELETED,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcDeleted,
		  0, 0 );
    mdlAddService(mdl,PST_SELDSTDELETED,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstDeleted,
		  0, 0 );
    mdlAddService(mdl,PST_SELSRCBYID,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcById,
		  sizeof(struct inSelById), sizeof(struct outSelById));
    mdlAddService(mdl,PST_SELDSTBYID,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstById,
		  sizeof(struct inSelById), sizeof(struct outSelById));
    mdlAddService(mdl,PST_SELSRCMASS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcMass,
		  sizeof(struct inSelMass), sizeof(struct outSelMass));
    mdlAddService(mdl,PST_SELDSTMASS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstMass,
		  sizeof(struct inSelMass), sizeof(struct outSelMass));
    mdlAddService(mdl,PST_SELSRCPHASEDENSITY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcPhaseDensity,
		  sizeof(struct inSelPhaseDensity), sizeof(struct outSelPhaseDensity));
    mdlAddService(mdl,PST_SELDSTPHASEDENSITY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstPhaseDensity,
		  sizeof(struct inSelPhaseDensity), sizeof(struct outSelPhaseDensity));
    mdlAddService(mdl,PST_SELSRCBOX,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcBox,
		  sizeof(struct inSelBox), sizeof(struct outSelBox));
    mdlAddService(mdl,PST_SELDSTBOX,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstBox,
		  sizeof(struct inSelBox), sizeof(struct outSelBox));
    mdlAddService(mdl,PST_SELSRCSPHERE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcSphere,
		  sizeof(struct inSelSphere), sizeof(struct outSelSphere));
    mdlAddService(mdl,PST_SELDSTSPHERE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstSphere,
		  sizeof(struct inSelSphere), sizeof(struct outSelSphere));
    mdlAddService(mdl,PST_SELSRCCYLINDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcCylinder,
		  sizeof(struct inSelCylinder), sizeof(struct outSelCylinder));
    mdlAddService(mdl,PST_SELDSTCYLINDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstCylinder,
		  sizeof(struct inSelCylinder), sizeof(struct outSelCylinder));
    mdlAddService(mdl,PST_SELSRCGROUP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSrcGroup,
		  sizeof(int), 0);
    mdlAddService(mdl,PST_SELDSTGROUP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDstGroup,
		  sizeof(int), 0);
    mdlAddService(mdl,PST_DEEPESTPOT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDeepestPot,
		  sizeof(struct inDeepestPot), sizeof(struct outDeepestPot));
    mdlAddService(mdl,PST_PROFILE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstProfile,
		  sizeof(struct inProfile), 0); 
    mdlAddService(mdl,PST_CALCDISTANCE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcDistance,
		  sizeof(struct inCalcDistance), 0);
    mdlAddService(mdl,PST_CALCCOM,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcCOM,
		  sizeof(struct inCalcCOM), sizeof(struct outCalcCOM));
    mdlAddService(mdl,PST_COUNTDISTANCE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCountDistance,
		  sizeof(struct inCountDistance), sizeof(struct outCountDistance));
    mdlAddService(mdl,PST_INITGRID,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitGrid,
		  sizeof(struct inInitGrid), 0);
    mdlAddService(mdl,PST_GRIDPROJECT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGridProject,
		  sizeof(struct inGridProject), 0);
#ifdef MDL_FFTW
    mdlAddService(mdl,PST_MEASUREPK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstMeasurePk,
		  sizeof(struct inMeasurePk), sizeof(struct outMeasurePk));
#endif
    mdlAddService(mdl,PST_TOTALMASS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstTotalMass,
		  0, sizeof(struct outTotalMass));

#ifdef USE_PSD
    mdlAddService(mdl,PST_BUILDPSDTREE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstBuildPsdTree,
		  sizeof(struct inPSD),pkdMaxNodeSize());
    mdlAddService(mdl,PST_PSD,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSD,
		  sizeof(struct inPSD), 0);
    mdlAddService(mdl,PST_PSDLINK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDLink,
		  sizeof(struct inPSD), 0);
    mdlAddService(mdl,PST_PSD_INIT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDInit,
		  sizeof(struct inPSD), 0);
    mdlAddService(mdl,PST_PSD_JOINBRIDGES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDJoinBridges,
		  0, sizeof(int));

    mdlAddService(mdl,PST_PSD_CLG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDCountLocalGroups,
		  sizeof(int), nThreads*sizeof(struct inAssignGlobalIds));
    mdlAddService(mdl,PST_PSD_ASSIGN_GLOBAL_IDS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDAssignGlobalIds,
		  nThreads*sizeof(struct inAssignGlobalIds), 0);
    mdlAddService(mdl,PST_PSD_MERGENOISYGROUPS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDMergeNoisyGroups,
		  0, 0);
    mdlAddService(mdl,PST_PSD_JOINGROUPBRIDGES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDJoinGroupBridges,
		  0, sizeof(int));
    mdlAddService(mdl,PST_PSD_SETGLOBALID,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDSetGlobalId,
		  0, 0);
    mdlAddService(mdl,PST_PSD_FINISH,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPSDFinish,
		  sizeof(struct inPSD), 0);
    mdlAddService(mdl,PST_WRITE_PSGROUPS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWritePsGroups,
		  sizeof(struct inWritePsGroups), 0);
    mdlAddService(mdl,PST_UNBIND,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstUnbind,
		  sizeof(int), 0);
#endif
    mdlCommitServices(mdl);
   }

void pstInitialize(PST *ppst,MDL mdl,LCL *plcl) {
    PST pst;

    pst = (PST)malloc(sizeof(struct pstContext));
    mdlassert(mdl,pst != NULL);
    *ppst = pst;
    pst->plcl = plcl;
    pst->mdl = mdl;
    pst->idSelf = mdlSelf(mdl);
    pst->pstLower = NULL;
    pst->idUpper = -1;	/* invalidate upper 'id' */
    pst->nLeaves = 1;
    pst->nLower = 0;
    pst->nUpper = 0;
    pst->iSplitDim = -1;
    pst->iVASplitSide = 0;
    pst->nLowTot = 0;
    pst->nHighTot = 0;
    }


void pstFinish(PST pst) {
    PST pstKill;

    while (pst) {
	pstKill = pst;
	if (pst->nLeaves == 1 && pst->plcl->pkd)
	    pkdFinish(pst->plcl->pkd);
	pst = pst->pstLower;
	free(pstKill);
	}
    }

void pstSetAdd(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    PST pstNew;
    struct inSetAdd *in = vin;
    int n, idMiddle,iProcLower,iProcUpper;
    mdlassert(pst->mdl,nIn == sizeof(struct inSetAdd));
    mdlassert(pst->mdl,pst->nLeaves==1);
    mdlassert(pst->mdl,in->idLower==mdlSelf(pst->mdl));
    n = in->idUpper - in->idLower;
    idMiddle = (in->idUpper + in->idLower) / 2;
    if ( n > 1 ) {
	int rID;
	/* Make sure that the pst lands on core zero */
	iProcLower = mdlThreadToProc(pst->mdl,in->idLower);
	iProcUpper = mdlThreadToProc(pst->mdl,in->idUpper-1);
	if (iProcLower!=iProcUpper) {
	    idMiddle = mdlProcToThread(pst->mdl,mdlThreadToProc(pst->mdl,idMiddle));
	    }
	pst->nLeaves += n - 1;
	pst->nLower = idMiddle - in->idLower;
	pst->nUpper = in->idUpper - idMiddle;

	in->idLower = idMiddle;
	pst->idUpper = in->idLower;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETADD,in,nIn);
	in->idLower = mdlSelf(pst->mdl);
	in->idUpper = idMiddle;
	pstInitialize(&pstNew,pst->mdl,pst->plcl);
	pst->pstLower = pstNew;
	pstSetAdd(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    if (pnOut) *pnOut = 0;
    }

void pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inReadFile *in = vin;
    int *pout = vout;
    uint64_t nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(struct inReadFile));
    nFileStart = in->nNodeStart;
    nFileEnd = in->nNodeEnd;
    nFileTotal = nFileEnd - nFileStart + 1;
    
    if (pst->nLeaves > 1) {
	int rID;
	nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
	in->nNodeStart = nFileSplit;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_ONENODEREADINIT,in,nIn);
	in->nNodeStart = nFileStart;
	in->nNodeEnd = nFileSplit - 1;
	pstOneNodeReadInit(pst->pstLower,in,nIn,pout,pnOut);
	in->nNodeEnd = nFileEnd;
	mdlGetReply(pst->mdl,rID,pout+pst->nLower,pnOut);
	}
    else {
	/*
	** Determine the size of the local particle store.
	*/
	nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
	if (plcl->pkd) pkdFinish(plcl->pkd);
	pkdInitialize(
	    &plcl->pkd,pst->mdl,nStore,in->nBucket,in->nGroup,
	    in->nTreeBitsLo,in->nTreeBitsHi,
	    in->iCacheSize,in->iWorkQueueSize,in->iCUDAQueueSize,in->fPeriod,
	    in->nSpecies[FIO_SPECIES_DARK],
	    in->nSpecies[FIO_SPECIES_SPH],
	    in->nSpecies[FIO_SPECIES_STAR],
	    in->mMemoryModel, in->nDomainRungs);
	*pout = nFileTotal; /* Truncated: okay */
	}
    if (pnOut) *pnOut = sizeof(*pout) * pst->nLeaves;
    }

static void _SwapClasses(PKD pkd, int id) {
    PARTCLASS *pClass;
    int n;
    int rID;

    pClass = malloc(PKD_MAX_CLASSES*sizeof(PARTCLASS));
    assert(pClass!=NULL);

    n = pkdGetClasses( pkd, PKD_MAX_CLASSES, pClass );
    rID = mdlReqService(pkd->mdl,id,PST_SWAPCLASSES,pClass,n*sizeof(PARTCLASS));
    mdlGetReply(pkd->mdl,rID,pClass,&n);
    n = n / sizeof(PARTCLASS);
    pkdSetClasses( pkd, n, pClass, 0 );
    free(pClass);
    }

void pstReadFile(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inReadFile *in = vin;
    FIO fio;
    uint64_t nNodeStart,nNodeEnd,nNodeTotal,nNodeSplit,nStore;
    int rID;

    mdlassert(pst->mdl,nIn >= sizeof(struct inReadFile));

    nNodeStart = in->nNodeStart;
    nNodeEnd = in->nNodeEnd;
    nNodeTotal = nNodeEnd - nNodeStart + 1;
    if (pst->nLeaves > 1 && in->nProcessors > 1) {
	int nProcessors = in->nProcessors;
	int nProcUpper = pst->nUpper * nProcessors / pst->nLeaves;
	int nProcLower = nProcessors - nProcUpper;
	nNodeSplit = nNodeStart + pst->nLower*(nNodeTotal/pst->nLeaves);
	if ( nProcessors > 1 ) {
	    in->nProcessors = nProcUpper;
	    in->nNodeStart = nNodeSplit;
	    in->nNodeEnd = nNodeEnd;
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_READFILE,in,nIn);
	    }

	in->nProcessors = nProcLower;
	in->nNodeStart = nNodeStart;
	in->nNodeEnd = nNodeSplit - 1;
	pstReadFile(pst->pstLower,in,nIn,NULL,NULL);

	if ( nProcessors <= 1 ) {
	    in->nProcessors = nProcUpper;
	    in->nNodeStart = nNodeSplit;
	    in->nNodeEnd = nNodeEnd;
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_READFILE,in,nIn);
	    }

	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	int *nParts = malloc(pst->nLeaves * sizeof(nParts));
	int nid, i;
	uint64_t nStart;
	PKD pkd;
	MDL mdl;
	PST pst0;

	assert(nParts!=NULL);
	pstOneNodeReadInit(pst,in,sizeof(*in),nParts,&nid);
	pkd = plcl->pkd;
	mdl = pkd->mdl;

	fio = fioLoad(in+1,in->dOmega0,in->dOmegab);
	assert(fio!=NULL);

	nStart = nNodeStart + nParts[0];
	for(i=1; i<pst->nLeaves; ++i) {
	    int id = mdlSelf(mdl) + i;
	    int inswap;
	    /*
	     * Read particles into the local storage.
	     */
	    assert(pkd->nStore >= nParts[i]);
	    pkdReadFIO(pkd, fio, nStart, nParts[i], in->dvFac,in->dTuFac);
	    nStart += nParts[i];
	    /*
	     * Now shove them over to the remote processor.
	     */
	    _SwapClasses(pkd,id);
	    inswap = mdlSelf(mdl);
	    rID = mdlReqService(mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	    pkdSwapAll(pkd, id);
	    mdlGetReply(mdl,rID,NULL,NULL);
	    }
	pkdReadFIO(pkd, fio, nNodeStart, nParts[0], in->dvFac,in->dTuFac);
	free(nParts);
	fioClose(fio);
	}
    if (pnOut) *pnOut = 0;
    }

int _pstRejMatch(PST pst,int n1,OREJ *p1,int n2,OREJ *p2,int *pidSwap) {
    int id,i,i1=-1,i2=-1,id1,id2;
    local_t nLarge;
    total_t s1,s2,r1,r2;

    /*
    ** Check to see if there is enough space...
    */
    s1 = 0;
    r1 = 0;
    for (i=0;i<n1;++i) {
	s1 += p1[i].nSpace;
	r1 += p1[i].nRejects;
	}
    s2 = 0;
    r2 = 0;
    for (i=0;i<n2;++i) {
	s2 += p2[i].nSpace;
	r2 += p2[i].nRejects;
	}
    mdlassert(pst->mdl,r1 <= s2);
    mdlassert(pst->mdl,r2 <= s1);
    /*
    ** First invalidate the pidSwap array.
    */
    for (id=0;id<mdlThreads(pst->mdl);++id) pidSwap[id] = -1;
    /*
    ** Now map largest nReject of p1 to largest nSpace of p2.
    */
    while (1) {
	nLarge = 0;
	for (i=0;i<n1;++i) {
	    if (p1[i].nRejects > nLarge) {
		nLarge = p1[i].nRejects;
		i1 = i;
		}
	    }
	if (nLarge == 0) break;
	nLarge = 0;
	for (i=0;i<n2;++i) {
	    if (p2[i].nSpace > nLarge) {
		nLarge = p2[i].nSpace;
		i2 = i;
		}
	    }
	if (nLarge == 0) break;
	p1[i1].nRejects = 0;
	p1[i1].nSpace = 0;
	p2[i2].nRejects = 0;
	p2[i2].nSpace = 0;
	id1 = p1[i1].id;
	id2 = p2[i2].id;
	pidSwap[id1] = id2;
	pidSwap[id2] = id1;
	}
    /*
    ** Now map largest nReject of p2 to largest nSpace of p1.
    ** However, already mapped stuff is ignored, by the above!
    */
    while (1) {
	nLarge = 0;
	for (i=0;i<n2;++i) {
	    if (p2[i].nRejects > nLarge) {
		nLarge = p2[i].nRejects;
		i2 = i;
		}
	    }
	if (nLarge == 0) break;
	nLarge = 0;
	for (i=0;i<n1;++i) {
	    if (p1[i].nSpace > nLarge) {
		nLarge = p1[i].nSpace;
		i1 = i;
		}
	    }
	if (nLarge == 0) break;
	p1[i1].nRejects = 0;
	p1[i1].nSpace = 0;
	p2[i2].nRejects = 0;
	p2[i2].nSpace = 0;
	id1 = p1[i1].id;
	id2 = p2[i2].id;
	pidSwap[id1] = id2;
	pidSwap[id2] = id1;
	}
    for (i=0;i<mdlThreads(pst->mdl);++i)
	if (pidSwap[i] != -1) return(1);
    return(0);
    }


#define MAX_ITTR	64

void _pstRootSplit(PST pst,int iSplitDim,int bDoRootFind,int bDoSplitDimFind,
		   int bSplitVA) {
    int NUM_SAFETY = 4;			/* slop space when filling up memory */
    uint64_t nSafeTot;			/* total slop space we have to play with */
    uint64_t margin;			/* more slop */
    int d,ittr,nOut;
    /*
    ** Why are these initialized to -1 here???
    int nLow=-1,nHigh=-1;
    */
    uint64_t nLow,nHigh;
    uint64_t nLowerStore,nUpperStore;
    uint64_t nLowTot,nHighTot;
    uint64_t nLast;					/* number of particles at the last split iteration */
    uint64_t nTotalActive;
    int nDiff=0;				/* Difference between one iteration and the next, seems to only be used to warn. */
    FLOAT fLow,fHigh;
    FLOAT fl,fu,fm=-1,fmm;
    struct outFreeStore outFree;
    struct inWeight inWt;
    struct inWeightWrap inWtWrap;
    struct outWeight outWtLow;
    struct outWeight outWtHigh;
    struct outWeightWrap outWtWrLow;
    struct outWeightWrap outWtWrHigh;
    struct inCountVA inCtVA;
    struct outCountVA outCtVA;
    OREJ *pLowerRej,*pUpperRej;
    int *pidSwap,iRet;
    char ach[256];				/* Debug */
    mdlTimer t;
    int pFlag;					/* 0 => we are splitting all particles at once. 1 => we first split active, and then inactive. */
    int dBnd;
    int rID;

    mdlZeroTimer(pst->mdl,&t);
    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    rID = mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
    pstFreeStore(pst->pstLower,NULL,0,&outFree,NULL);
    nLowerStore = outFree.nFreeStore;
    mdlGetReply(pst->mdl,rID,&outFree,NULL);
    nUpperStore = outFree.nFreeStore;

    mdlprintf(pst->mdl,"_pstRootSplit: id %d Level %d\n",pst->idSelf,pst->iLvl);
    mdlPrintTimer(pst->mdl,"TIME START _pstRootSplit ",&t);
    mdlprintf(pst->mdl,"_pstRootSplit: fA0 %f fIA0 %f RS? %d DC? %d\n",
	      pst->fSplit,pst->fSplitInactive,bDoRootFind, bDoSplitDimFind);
    /* Debug */
    /*
      sprintf(ach,"id: %d _pstRootSplit\n", pst->idSelf );
      mdlDiag(pst->mdl,ach);
    */

    if (bDoSplitDimFind || pst->iSplitDim == -1) {
	pst->iSplitDim = iSplitDim;
	}

    d = dBnd = pst->iSplitDim;
    fl = pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd];
    fu = pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd];
    fm = pst->fSplit;
    ittr = -1;

    if (bDoRootFind || fm<fl || fm>fu || (bSplitVA && (pst->iVASplitSide == 0))) {
	/*
	** First order the particles into active/inactive order...
	*/
	uint64_t nActiveOrder;
	pstActiveOrder(pst, NULL, 0, &nActiveOrder, NULL); /* SOON NO MORE ACTIVE ORDER */

	fmm = (fl + fu)/2;
	/*
	 * First find total number of active particles.
	 */
	inWt.iSplitDim = d;
	inWt.fSplit = fmm;
	inWt.ittr = 0;
	inWt.iSplitSide = 1;
	inWt.pFlag = 1;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	mdlGetReply(pst->mdl,rID,&outWtHigh,NULL);
	nTotalActive = outWtLow.nLow + outWtHigh.nLow
		       + outWtLow.nHigh + outWtHigh.nHigh;
	mdlassert(pst->mdl,nActiveOrder == nTotalActive);
	pFlag = 1;
	if (nTotalActive <=1) {
	    pFlag = 0;			/* Divide them all */
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	    mdlGetReply(pst->mdl,rID,&outWtHigh,NULL);
	    nTotalActive = outWtLow.nLow + outWtHigh.nLow
			   + outWtLow.nHigh + outWtHigh.nHigh;
	    }

	/*
	** Now start the ROOT finder based on balancing active weight ALONE!
	** (unless pFlag == 0)
	*/
	ittr = 0;
	while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
	    fm = fmm;
	    inWt.iSplitDim = d;
	    inWt.fSplit = fm;
	    inWt.ittr = ittr;
	    inWt.iSplitSide = 1;
	    inWt.pFlag = pFlag;
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	    mdlGetReply(pst->mdl,rID,&outWtHigh,NULL);
	    /*
	    ** Add lower and Upper subsets weights and numbers
	    */
	    nLow = outWtLow.nLow + outWtHigh.nLow;
	    nHigh = outWtLow.nHigh + outWtHigh.nHigh;
	    fLow = outWtLow.fLow + outWtHigh.fLow;
	    fHigh = outWtLow.fHigh + outWtHigh.fHigh;
	    /*
	      printf("ittr:%d l:%d u:%d lw:%f uw:%f\n",ittr,nLow,nHigh,fLow,fHigh);
	    */
	    if (nLow == 1 && nHigh == 1) /* break on trivial case */
		break;
	    if (pFlag) {			/* split on work */
		if (fLow/pst->nLower > fHigh/pst->nUpper) fu = fm;
		else if (fLow/pst->nLower < fHigh/pst->nUpper) fl = fm;
		else break;
		}
	    else {				/* split on number */
		if (nLow/(double)pst->nLower >
			nHigh/(double)pst->nUpper) fu = fm;
		else if (nLow/(double)pst->nLower <
			 nHigh/(double)pst->nUpper) fl = fm;
		else break;
		}
	    fmm = (fl + fu)/2;
	    ++ittr;
	    }
	/*
	** We have now found a new split. We need to decide which side of this split
	** we want to put the very active particles.
	*/
	inCtVA.iSplitDim = d;
	inCtVA.fSplit = fm;
	pstCountVA(pst,&inCtVA,sizeof(inCtVA),&outCtVA,NULL);
	if (outCtVA.nLow == 0 && outCtVA.nHigh == 0) {
	    pst->iVASplitSide = 0;
	    }
	else if (outCtVA.nHigh > outCtVA.nLow) {
	    pst->iVASplitSide = -1;
	    }
	else {
	    pst->iVASplitSide = 1;
	    }
	mdlPrintTimer(pst->mdl,"TIME active split _pstRootSplit ",&t);
	}

    pst->fSplit = fm;

    mdlprintf(pst->mdl, "id: %d (%d) Chose split: %f (%f,%f) %d %d\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd], pst->nLower, pst->nUpper);
    if (ittr != -1)
	mdlprintf(pst->mdl, "  Low %"PRIu64" %f,  High %"PRIu64" %f, ittr=%d\n",
		  nLow,outWtLow.fLow + outWtHigh.fLow, nHigh,
		  outWtLow.fHigh + outWtHigh.fHigh,ittr);
    nLow = 0;
    nHigh = 0;
    fLow = 0.0;
    fHigh = 0.0;

    /*
    ** Now we see if the TOTAL number of particles in the lower and upper
    ** subsets exceeds the local particle stores. If so then we need to
    ** find a new boundary to distribute the INACTIVE particles so that
    ** everything fits.
    */
    inWtWrap.iSplitDim = d;
    fl = pst->fSplit + 1e-6*pst->bnd.fMax[dBnd];
    fu = pst->fSplit - 1e-6*pst->bnd.fMax[dBnd];

    if (!bDoSplitDimFind) fm = pst->fSplitInactive;
    else {
	fm = 0.5*(fl+fu);
	if (fm < pst->bnd.fCenter[dBnd]) fm = pst->bnd.fCenter[dBnd] + 1.000001*pst->bnd.fMax[dBnd];
	else fm = pst->bnd.fCenter[dBnd] - 1.000001*pst->bnd.fMax[dBnd];
	}
    mdlprintf(pst->mdl, "id: %d (%d) Zeroeth guess reverse split: %f (%f,%f)\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd]);
    inWtWrap.fSplit = fm;
    inWtWrap.fSplit2 = pst->fSplit;
    inWtWrap.ittr = 0;
    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
    inWtWrap.iSplitSide = 1;
    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
    inWtWrap.iSplitSide = 0;
    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
    mdlGetReply(pst->mdl,rID,&outWtWrHigh,NULL);
    /*
    ** Add lower and Upper subsets numbers of particles
    */
    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;

    nSafeTot = nLowerStore + nUpperStore - (nLowTot + nHighTot);
    if (nSafeTot/pst->nLeaves < NUM_SAFETY) {
	NUM_SAFETY = nSafeTot/pst->nLeaves;
	sprintf(ach,"id: %d tripped inactive NUM_SAFETY %d  Low %"PRIu64"/%"PRIu64"  High %"PRIu64"/%"PRIu64"\n",
		pst->idSelf, NUM_SAFETY, nLowTot, nLowerStore, nHighTot, nUpperStore);
	mdlDiag(pst->mdl,ach);
	mdlprintf(pst->mdl,"id: %d tripped inactive NUM_SAFETY %d  Low %%"PRIu64"/%"PRIu64"  High %"PRIu64"/%"PRIu64"\n",
		  pst->idSelf, NUM_SAFETY, nLowTot, nLowerStore, nHighTot, nUpperStore);
	}

    margin = nSafeTot/pst->nLeaves/20;
    if (margin < NUM_SAFETY/2) margin = NUM_SAFETY/2;

    mdlprintf(pst->mdl,"id: %d  %d Low %"PRIu64"/%"PRIu64"   %d High %"PRIu64"/%"PRIu64"  NUM_SAFETY %d margin %d\n",
	      pst->idSelf, pst->nLower,nLowTot, nLowerStore, pst->nUpper,nHighTot, nUpperStore,NUM_SAFETY,margin);


    if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) {
	sprintf(ach,"id: %d: nLowTot > nLowerStore-NUM_SAFETY*pst->nLower %"PRIu64" %"PRIu64" %d %d\n",
		pst->idSelf, nLowTot, nLowerStore, NUM_SAFETY, pst->nLower);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fl = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1;
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,rID,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*
	      mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
	      pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if (nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    /*
	      if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) fl = fm;
	      else if (nLowTot < nLowerStore-2*NUM_SAFETY*pst->nLower) fu = fm;
	    */
	    /*
	      if (nLowTot/pst->nLower > nHighTot/pst->nUpper) fl = fm;
	      else if (nLowTot/pst->nLower < nHighTot/pst->nUpper-NUM_SAFETY) fu = fm;
	    */
	    if (nLowTot > nLowerStore-margin*pst->nLower) fl = fm;
	    else if (nLowTot < nLowerStore-2*margin*pst->nLower) fu = fm;
	    else {
		fl = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %"PRIu64" High %"PRIu64"\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
	    if (abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nLowTot <= nLowerStore);
	mdlPrintTimer(pst->mdl,"TIME fix lower II _pstRootSplit ",&t);
	}
    else if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) {
	sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %"PRIu64" %"PRIu64" %d %d\n",
		pst->idSelf, nHighTot, nUpperStore, NUM_SAFETY, pst->nUpper);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fu = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1;
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,rID,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*
			mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
			pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
			pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if (nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    /*
	      if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) fu = fm;
	      else if (nHighTot < nUpperStore-2*NUM_SAFETY*pst->nUpper) fl = fm;
	    */
	    /*
	      if (nHighTot/pst->nUpper > nLowTot/pst->nLower) fu = fm;
	      else if (nHighTot/pst->nUpper < nLowTot/pst->nLower-NUM_SAFETY) fl = fm;
	    */
	    if (nHighTot > nUpperStore-margin*pst->nUpper) fu = fm;
	    else if (nHighTot < nUpperStore-2*margin*pst->nUpper) fl = fm;
	    else {
		fu = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %"PRIu64" High %"PRIu64"\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
	    if (abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFETY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nHighTot <= nUpperStore);
	mdlPrintTimer(pst->mdl,"TIME fix upper II _pstRootSplit ",&t);
	}

    if (nLowTot < NUM_SAFETY*pst->nLower) {
	sprintf(ach,"id: %d: nLowTot < NUM_SAFETY*pst->nLower %"PRIu64" %"PRIu64" %d %d\n",
		pst->idSelf, nLowTot, nLowerStore, NUM_SAFETY, pst->nLower);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fu = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1;
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,rID,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*
	      mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
	      pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if (nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    if (nLowTot > margin*pst->nLower) fl = fm;
	    else if (nLowTot < NUM_SAFETY*pst->nLower) fu = fm;
	    else {
		fl = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix too few Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %"PRIu64" High %"PRIu64"\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
	    if (abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nLowTot <= nLowerStore);
	mdlPrintTimer(pst->mdl,"TIME fix lower II _pstRootSplit ",&t);
	}
    if (nHighTot < NUM_SAFETY*pst->nUpper) {
	sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %"PRIu64" %"PRIu64" %d %d\n",
		pst->idSelf, nHighTot, nUpperStore, NUM_SAFETY, pst->nUpper);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fl = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1;
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,rID,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*
			mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
			pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
			pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if (nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    if (nHighTot > margin*pst->nUpper) fu = fm;
	    else if (nHighTot < NUM_SAFETY*pst->nUpper) fl = fm;
	    else {
		fu = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] &&
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix Too few High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %"PRIu64" High %"PRIu64"\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
	    if (abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFETY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nHighTot <= nUpperStore);
	mdlPrintTimer(pst->mdl,"TIME fix upper II _pstRootSplit ",&t);
	}

    mdlassert(pst->mdl, nLowTot >= pst->nLower);
    mdlassert(pst->mdl, nHighTot >= pst->nUpper);
    mdlassert(pst->mdl, nLowTot <= nLowerStore);
    mdlassert(pst->mdl, nHighTot <= nUpperStore);

    pst->nLowTot = nLowTot;
    pst->nHighTot = nHighTot;

    mdlPrintTimer(pst->mdl,"TIME Total Split _pstRootSplit ",&t);

    mdlprintf(pst->mdl, "id: %d (%d) Chose reverse split: %f (%f,%f)\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
    pst->fSplitInactive = fm;

    /*
    ** First Collect rejects.
    **
    ** Careful, SERVICE PST_COLREJECTS does NOT conform strictly to
    ** the proper use of MDL. This should be fixed in the future.
    ** FIXED -- MDL modified.
    */
    pLowerRej = malloc(pst->nLower*sizeof(OREJ));
    mdlassert(pst->mdl,pLowerRej != NULL);
    pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
    mdlassert(pst->mdl,pUpperRej != NULL);
    pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
    mdlassert(pst->mdl,pidSwap != NULL);

    rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,NULL,0);
    pstColRejects(pst->pstLower,NULL,0,pLowerRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
    mdlGetReply(pst->mdl,rID,pUpperRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);

    mdlPrintTimer(pst->mdl,"TIME Collected Rejects _pstRootSplit ",&t);


    ittr = 0;
    while (1) {
	iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,
			    pst->nUpper,pUpperRej,pidSwap);
	if (!iRet) break;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
		      mdlThreads(pst->mdl)*sizeof(int));
	pstSwapRejects(pst->pstLower,pidSwap,
		       mdlThreads(pst->mdl)*sizeof(int),pLowerRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,rID,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);

	++ittr;
	}


    free(pLowerRej);
    free(pUpperRej);
    free(pidSwap);

    mdlPrintTimer(pst->mdl,"TIME (FINISH) Swapped Rejects _pstRootSplit ",&t);

    }

#ifdef MPI_VERSION
void pstOrbBegin(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inOrbBegin *in = (struct inOrbBegin *)vin;
    int rID;
    if (pst->nLeaves > 1) {
        rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_BEGIN,vin,nIn);
        pstOrbBegin(pst->pstLower,vin,nIn,vout,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        pkdOrbBegin(plcl->pkd,in->nRungs);
        }
    }

void pstOrbSelectRung(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    uint64_t nActive;
    struct inOrbSelectRung *in = (struct inOrbSelectRung *)vin;
    struct outOrbSelectRung *out = (struct outOrbSelectRung *)vout;
    int rID;
    if (pst->nLeaves > 1) {
	struct outOrbSelectRung outLower, outUpper;
        rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_SELECT_RUNG,vin,nIn);
        pstOrbSelectRung(pst->pstLower,vin,nIn,&outLower,NULL);
        mdlGetReply(pst->mdl,rID,&outUpper,NULL);
	nActive = outLower.nActive + outUpper.nActive;
        }
    else {
	nActive = pkdOrbSelectRung(plcl->pkd,in->iRung);
        }
    if (out) out->nActive = nActive;
    if (pnOut) *pnOut = sizeof(struct outOrbSelectRung);
    }

void pstOrbUpdateRung(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_UPDATE_RUNG,vin,nIn);
        pstOrbUpdateRung(pst->pstLower,vin,nIn,vout,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        pkdOrbUpdateRung(plcl->pkd);
        }
    }

void pstOrbFinish(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_FINISH,vin,nIn);
        pstOrbFinish(pst->pstLower,vin,nIn,vout,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        pkdOrbFinish(plcl->pkd);
        }
    }

void pstOrbSplit(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    int *in = (int*)vin;

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_SPLIT,NULL,0);
        pstOrbSplit(pst->pstLower,vin,nIn,vout,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        pkdOrbSplit(plcl->pkd,nIn==0?-1:*in);
        }
    }

void pstOrbRootFind(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inOrbRootFind *in = (struct inOrbRootFind *)vin;
    struct outOrbRootFind *out = (struct outOrbRootFind *)vout;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inOrbRootFind));
    if (pst->nLeaves > 1) {
	struct inOrbRootFind no;
	struct outOrbRootFind outno;
	no.bnd = in->bnd;
	no.dFraction = 0.0; /* Do not split anyone but us (and maybe not even us) */
        rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_ROOT_FIND,&no,sizeof(no));
        pstOrbRootFind(pst->pstLower,vin,nIn,vout,pnOut);
        mdlGetReply(pst->mdl,rID,&outno,NULL);
        }
    else {
        out->nDomains = pkdOrbRootFind(plcl->pkd,in->dFraction,
	    in->nLower,in->nUpper,in->dReserveFraction,
	    &in->bnd,&out->dSplit,&out->iDim);
        }
    }

/* Our domain must be split -- bounds given as input */
void pstOrbDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inOrbDecomp *in = (struct inOrbDecomp *)vin;
    struct inOrbDecomp in1, in2;
    struct inOrbRootFind irf;
    struct outOrbRootFind out;
    int iDomain;
    int rID;
    if (pst->nLeaves > 1) {
	irf.nLower = pst->nLowerStore;
	irf.nUpper = pst->nUpperStore;
	irf.bnd = in->bnd;
	irf.dFraction = 1.0 * pst->nLower / pst->nLeaves;
	irf.dReserveFraction = 1.0 - 1.0 / log10(8.0+1.0*pst->nLeaves);

	pstOrbRootFind(pst,&irf,sizeof(irf),&out,NULL);
	pst->iSplitDim = out.iDim;

	iDomain = pst->idUpper;
        rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_SPLIT,NULL,0);
        pstOrbSplit(pst->pstLower,&iDomain,sizeof(iDomain),NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);


	in2.bnd = in->bnd;
	in2.bnd.fMax[out.iDim] = 0.5*(in->bnd.fCenter[out.iDim]+in->bnd.fMax[out.iDim] - out.dSplit);
	in2.bnd.fCenter[out.iDim] = out.dSplit + in2.bnd.fMax[out.iDim];

	in1.bnd = in->bnd;
	in1.bnd.fMax[out.iDim] = 0.5*(out.dSplit - in->bnd.fCenter[out.iDim]+in->bnd.fMax[out.iDim]);
	in1.bnd.fCenter[out.iDim] = out.dSplit - in1.bnd.fMax[out.iDim];

        rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORB_DECOMP,&in2,sizeof(in2));
        pstOrbDecomp(pst->pstLower,&in1,sizeof(in1),vout,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        while(pkdOrbRootFind(plcl->pkd,0,0,0,0,NULL,NULL,NULL)) {
	    pstOrbSplit(pst,NULL,0.0,NULL,NULL);
	    }
        }
    if (pnOut) *pnOut = 0;
    }

void pstPeanoHilbertDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inPeanoHilbertDecomp *in = (struct inPeanoHilbertDecomp *)vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPeanoHilbertDecomp));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PEANOHILBERTDECOMP,vin,nIn);
	pstPeanoHilbertDecomp(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdPeanoHilbertDecomp(plcl->pkd, in->nRungs, in->iMethod);
	}
    if (pnOut) *pnOut = 0;
    }

void pstRungOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inRungOrder *in = (struct inRungOrder *)vin;
    struct outRungOrder *out = vout;
    struct outRungOrder outRung;
    mdlassert(pst->mdl,nIn == sizeof(struct inRungOrder));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RUNGORDER,vin,nIn);
	pstRungOrder(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outRung,pnOut);
	BND_COMBINE(&pst->bnd,&out->bnd,&outRung.bnd);
	out->nMoved += outRung.nMoved;
	out->bnd = pst->bnd;
	}
    else {
	pkdRungOrder(plcl->pkd,in->iRung,&outRung.nMoved);
        out->bnd = pst->bnd = plcl->pkd->bnd;
	out->nMoved = outRung.nMoved;
	}
    if (pnOut) *pnOut = sizeof(struct outRungOrder);
    }
#endif

#define NEWSPLITDIMCUT 0.707
#define NMINFORROOTFIND 16

void pstDomainDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    int d=0,j,nBndWrapd;
    double dimsize;
    struct inDomainDecomp *in = vin;
    double l,u;
    mdlTimer t;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inDomainDecomp));

    mdlZeroTimer(pst->mdl,&t);
    mdlprintf(pst->mdl,"Starting pstDomainDecomp\n");

    pst->bnd = in->bnd;
    if (pst->nLeaves > 1) {
#ifdef USE_ITT
	__itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
	__itt_string_handle* shMyTask = __itt_string_handle_create("Domain Decomposition");
	__itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
	if (pst->iSplitDim != -1 && in->nActive < NMINFORROOTFIND) {
	    mdlprintf(pst->mdl,"Aborting RootFind -- too few actives.\n");
	    in->bDoSplitDimFind = 0;
	    in->bDoRootFind = 0;
	    }

	/*
	** Next determine the longest axis based on the bounds.
	** Don't switch dimensions unless the bounds have changed significantly.
	**
	** NB: Standard bnds don't work for Wrapped dimensions
	*/
	if (in->bDoSplitDimFind) {
	    d = pst->iSplitDim;
	    if (d==-1) {
		dimsize = -1;
		nBndWrapd = 0;
		}
	    else {
		dimsize = pst->bnd.fMax[d]*NEWSPLITDIMCUT;
		nBndWrapd = in->nBndWrap[d];
		}

	    for (j=0;j<3;++j) {
		if (in->nBndWrap[j] < nBndWrapd || pst->bnd.fMax[j] > dimsize) {
		    d=j;
		    dimsize = pst->bnd.fMax[d];
		    nBndWrapd = in->nBndWrap[d];
		    }
		}
	    }

	mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
	_pstRootSplit(pst,d,in->bDoRootFind,in->bDoSplitDimFind,in->bSplitVA);
	mdlPrintTimer(pst->mdl,"TIME RootSplit done in pstDomainDecomp",&t);

	mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
	/*
	** Now go on to DD of next levels, but pass correct wrapping bounds.
	*/
	d = pst->iSplitDim;
	nBndWrapd = in->nBndWrap[d];

	l = pst->bnd.fCenter[d] - pst->bnd.fMax[d];
	u = pst->bnd.fCenter[d] + pst->bnd.fMax[d];
	in->nBndWrap[d] = nBndWrapd;
	if (pst->fSplitInactive <= l || pst->fSplitInactive >= u) {
	    l = pst->fSplit;
	    }
	else if (pst->fSplitInactive > pst->fSplit) {
	    l = pst->fSplit;
	    u = pst->fSplitInactive;
	    }
	else
	    in->nBndWrap[d]++;

	in->bnd.fMax[d] = 0.5*(u - l);
	in->bnd.fCenter[d] = 0.5*(u + l);
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINDECOMP,
		      vin,sizeof(*in));
	    
	l = pst->bnd.fCenter[d] - pst->bnd.fMax[d];
	u = pst->bnd.fCenter[d] + pst->bnd.fMax[d];
	in->nBndWrap[d] = nBndWrapd;
	if (pst->fSplitInactive <= l || pst->fSplitInactive >= u) {
	    u = pst->fSplit;
	    }
	else if (pst->fSplitInactive < pst->fSplit) {
	    u = pst->fSplit;
	    l = pst->fSplitInactive;
	    }
	else
	    in->nBndWrap[d]++;

	in->bnd.fMax[d] = 0.5*(u - l);
	in->bnd.fCenter[d] = 0.5*(u + l);
	pstDomainDecomp(pst->pstLower,vin,sizeof(*in),NULL,NULL);

	mdlGetReply(pst->mdl,rID,NULL,NULL);
#ifdef USE_ITT
	__itt_task_end(domain);
#endif
	}
    else {
	float offs;
	/*
	** We always set plcl->pkd->bnd from pst->bnd.
	*/
	plcl->pkd->bnd = pst->bnd;   /* This resets the local bounding box, but doesn't squeeze! */
        offs= 0.5f / (plcl->pkd->nLocal*1.0f - 1.0f);
	for (j=0; j < 3; j++) {
	    pst->bnd.fMax[j] += offs;
	    }
	}
    if (pnOut) *pnOut = 0;
    }


void pstCalcBound(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    BND *out = vout;
    BND outBnd;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCBOUND,NULL,0);
	pstCalcBound(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,rID,&outBnd,NULL);
	BND_COMBINE(out,out,&outBnd);
	}
    else {
	pkdCalcBound(plcl->pkd,out);
	}
    if (pnOut) *pnOut = sizeof(BND);
    }

void pstCalcVBound(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    BND *out = vout;
    BND outBnd;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCVBOUND,NULL,0);
	pstCalcVBound(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,rID,&outBnd,NULL);
	BND_COMBINE(out,out,&outBnd);
	}
    else {
	pkdCalcVBound(plcl->pkd,out);
	}
    if (pnOut) *pnOut = sizeof(BND);
    }

void pstCombineBound(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    BND *out = vout;
    BND outBnd;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COMBINEBOUND,NULL,0);
	pstCombineBound(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,rID,&outBnd,NULL);
	BND_COMBINE(out,out,&outBnd);
	}
    else {
        *out = plcl->pkd->bnd;
	}
    if (pnOut) *pnOut = sizeof(BND);
    }


/*
** Make sure that the local particles are split into active and inactive
** when passing pFlag != 0.
** pFlag == 0 => weight all particles.
** pFlag > 0 => weight active particles.
** pFlag < 0 => weight inactive particles.
*/
void pstWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inWeight *in = vin;
    struct outWeight *out = vout;
    struct outWeight outWt;
    FLOAT fSplit,fLow,fHigh;
    int iSplitSide;
    int nLow,nHigh;

    mdlassert(pst->mdl,nIn == sizeof(struct inWeight));
    /*
      pkdStartTimer(plcl->pkd,7);
    */
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,in,nIn);
	pstWeight(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&outWt,NULL);
	out->nLow += outWt.nLow;
	out->nHigh += outWt.nHigh;
	out->fLow += outWt.fLow;
	out->fHigh += outWt.fHigh;
	}
    else {
	fSplit = in->fSplit;
	iSplitSide = in->iSplitSide;
	if (in->ittr == 0) {
	    /*
	    ** Initialize.
	    */
	    plcl->fSplit = fSplit;
	    if (in->pFlag == 0) {
		plcl->iWtFrom = 0;
		plcl->iWtTo = pkdLocal(plcl->pkd)-1;
		}
	    else if (in->pFlag > 0) {
		/*
		** Particles must be in the active-inactive order here!
		*/
		plcl->iWtFrom = 0;
		plcl->iWtTo = pkdActive(plcl->pkd)-1;
		}
	    else {
		/*
		** Particles must be in the active-inactive order here!
		*/
		plcl->iWtFrom = pkdActive(plcl->pkd);
		plcl->iWtTo = pkdLocal(plcl->pkd)-1;
		}
	    plcl->fWtLow = 0.0;
	    plcl->fWtHigh = 0.0;
	    }
	else {
	    /*
	    ** Update the Weight Sums and use smaller weight region.
	    */
	    if (fSplit < plcl->fSplit) {
		plcl->fWtHigh += plcl->fHigh;
		if (iSplitSide) plcl->iWtFrom = plcl->iPart;
		else plcl->iWtTo = plcl->iPart-1;
		}
	    else {
		plcl->fWtLow += plcl->fLow;
		if (iSplitSide) plcl->iWtTo = plcl->iPart-1;
		else plcl->iWtFrom = plcl->iPart;
		}
	    plcl->fSplit = fSplit;
	    }
	plcl->iPart = pkdWeight(plcl->pkd,in->iSplitDim,fSplit,iSplitSide,
				plcl->iWtFrom,plcl->iWtTo,
				&nLow,&nHigh,&fLow,&fHigh);
	out->nLow = nLow;
	out->nHigh = nHigh;
	out->fLow = fLow + plcl->fWtLow;
	out->fHigh = fHigh + plcl->fWtHigh;
	plcl->fLow = fLow;
	plcl->fHigh = fHigh;
	if (in->pFlag > 0) {
	    if (iSplitSide) out->nLow -= pkdInactive(plcl->pkd);
	    else out->nHigh -= pkdInactive(plcl->pkd);
	    }
	if (in->pFlag < 0) {
	    if (iSplitSide) out->nHigh -= pkdActive(plcl->pkd);
	    else out->nLow -= pkdActive(plcl->pkd);
	    }
	}
    if (pnOut) *pnOut = sizeof(struct outWeight);
    /*
      pkdStopTimer(plcl->pkd,7);
    */
    }



void pstCountVA(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCountVA *in = vin;
    struct outCountVA *out = vout;
    struct outCountVA outCt;

    mdlassert(pst->mdl,nIn == sizeof(struct inCountVA));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COUNTVA,in,nIn);
	pstCountVA(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&outCt,NULL);
	out->nLow += outCt.nLow;
	out->nHigh += outCt.nHigh;
	}
    else {
	pkdCountVA(plcl->pkd,in->iSplitDim,in->fSplit,&out->nLow,&out->nHigh);
	}
    if (pnOut) *pnOut = sizeof(struct outCountVA);
    }



/*
** Make sure that the local particles are split into active and inactive
** when passing pFlag != 0.
** pFlag == 0 => weight all particles.
** pFlag > 0 => weight active particles.
** pFlag < 0 => weight inactive particles.
*/
void pstWeightWrap(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inWeightWrap *in = vin;
    struct outWeightWrap *out = vout;
    struct outWeightWrap outWt;
    int nLow,nHigh;

    mdlassert(pst->mdl,nIn == sizeof(struct inWeightWrap));
    /*
      pkdStartTimer(plcl->pkd,7);
    */
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,in,nIn);
	pstWeightWrap(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&outWt,NULL);
	out->nLow += outWt.nLow;
	out->nHigh += outWt.nHigh;
	}
    else {
	if (in->ittr == 0) {
	    /*
	    ** Initialize.
	    */
	    plcl->fSplit = in->fSplit;
	    plcl->iWtFrom = 0;
	    plcl->iWtTo = pkdLocal(plcl->pkd)-1;
	    }
	else {
	    /*
	    ** Update the Weight Sums and use smaller weight region.
	    */
	    if ((in->fSplit < plcl->fSplit && in->fSplit<in->fSplit2 && plcl->fSplit<in->fSplit2) ||
		    (in->fSplit < plcl->fSplit && in->fSplit>in->fSplit2 && plcl->fSplit>in->fSplit2) ||
		    (in->fSplit > plcl->fSplit && in->fSplit>in->fSplit2 && plcl->fSplit<in->fSplit2)) {
		if (!in->iSplitSide) plcl->iWtFrom = plcl->iPart;
		else plcl->iWtTo = plcl->iPart-1;
		}
	    else {
		if (!in->iSplitSide) plcl->iWtTo = plcl->iPart-1;
		else plcl->iWtFrom = plcl->iPart;
		}
	    plcl->fSplit = in->fSplit;
	    }
	plcl->iPart = pkdWeightWrap(plcl->pkd,in->iSplitDim,in->fSplit,in->fSplit2,in->iSplitSide,
				    in->iVASplitSide,plcl->iWtFrom,plcl->iWtTo,&nLow,&nHigh);
	out->nLow = nLow;
	out->nHigh = nHigh;
	/* For collect rejects */
	plcl->nSplit = plcl->iPart;
	}
    if (pnOut) *pnOut = sizeof(struct outWeightWrap);
    /*
      pkdStopTimer(plcl->pkd,7);
    */
    }


/*
** Weight request for splitting into iOrder order.
*/
void pstOrdWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inOrdWeight *in = vin;
    struct outOrdWeight *out = vout;
    struct outOrdWeight outWt;
    int nLow,nHigh;

    mdlassert(pst->mdl,nIn == sizeof(struct inOrdWeight));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORDWEIGHT,in,nIn);
	pstOrdWeight(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&outWt,NULL);
	out->nLow += outWt.nLow;
	out->nHigh += outWt.nHigh;
	}
    else {
	if (in->ittr == 0) {
	    /*
	    ** Initialize.
	    */
	    plcl->iOrdSplit = in->iOrdSplit;
	    plcl->iWtFrom = 0;
	    plcl->iWtTo = pkdLocal(plcl->pkd)-1;
	    }
	else {
	    /*
	    ** Update the Weight Sums and use smaller weight region.
	    */
	    if (in->iOrdSplit < plcl->iOrdSplit) {
		if (in->iSplitSide) plcl->iWtFrom = plcl->iPart;
		else plcl->iWtTo = plcl->iPart-1;
		}
	    else {
		if (in->iSplitSide) plcl->iWtTo = plcl->iPart-1;
		else plcl->iWtFrom = plcl->iPart;
		}
	    plcl->iOrdSplit = in->iOrdSplit;
	    }
	plcl->iPart = pkdOrdWeight(plcl->pkd,in->iOrdSplit,in->iSplitSide,
				   plcl->iWtFrom,plcl->iWtTo,
				   &nLow,&nHigh);
	out->nLow = nLow;
	out->nHigh = nHigh;
	}
    if (pnOut) *pnOut = sizeof(struct outOrdWeight);
    }


void pstFreeStore(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outFreeStore *out = vout;
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
	pstFreeStore(pst->pstLower,NULL,0,out,NULL);
	pst->nLowerStore = out->nFreeStore;
	mdlGetReply(pst->mdl,rID,out,NULL);
	pst->nUpperStore = out->nFreeStore;
	out->nFreeStore = pst->nLowerStore + pst->nUpperStore;
	}
    else {
	out->nFreeStore = pkdFreeStore(plcl->pkd);
	}
    if (pnOut) *pnOut = sizeof(struct outFreeStore);
    }


void pstColRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper;

    mdlassert(pst->mdl,nIn == 0);
    /*
      pkdStartTimer(plcl->pkd,9);
    */
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,vin,nIn);
	pstColRejects(pst->pstLower,vin,nIn,&pOutRej[0],&nLower);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,rID,&pOutRej[iUpper],&nUpper);
	if (pnOut) *pnOut = nLower + nUpper;
	}
    else {
	pOutRej->nRejects = pkdColRejects(plcl->pkd, plcl->nSplit);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	pOutRej->nLocal = pkdLocal(plcl->pkd);
	if (pnOut) *pnOut = sizeof(OREJ);
	}
    /*
      pkdStopTimer(plcl->pkd,9);
    */
    }


void pstColOrdRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inColOrdRejects *in = vin;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inColOrdRejects));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,in,nIn);
	pstColOrdRejects(pst->pstLower,in,nIn,&pOutRej[0],&nLower);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,rID,&pOutRej[iUpper],&nUpper);
	if (pnOut) *pnOut = nLower + nUpper;
	}
    else {
	pOutRej->nRejects = pkdColOrdRejects(plcl->pkd,in->iOrdSplit,
					     in->iSplitSide);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	if (pnOut) *pnOut = sizeof(OREJ);
	}
    }


void pstSwapRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    int *pidSwap = vin;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper,idSwap;

    /*	pkdStartTimer(plcl->pkd,8);*/
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,vin,nIn);
	pstSwapRejects(pst->pstLower,vin,nIn,&pOutRej[0],&nLower);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,rID,&pOutRej[iUpper],&nUpper);
	if (pnOut) *pnOut = nLower + nUpper;
	}
    else {
	idSwap = pidSwap[pst->idSelf];
	pOutRej->nRejects = pkdSwapRejects(plcl->pkd,idSwap);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	if (pnOut) *pnOut = sizeof(OREJ);
	}
    /*	pkdStopTimer(plcl->pkd,8); */
    }

/*
 * Routine to swap all particles.  Note that this does not walk the pst
 * but simply works with one other processor.
 */
void pstSwapAll(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl;
    int *pidSwap = vin;
    PST lpst;

    mdlassert(pst->mdl,nIn == sizeof(*pidSwap));
    lpst = pst;
    while (lpst->nLeaves > 1) lpst = lpst->pstLower;
    plcl = lpst->plcl;
    pkdSwapAll(plcl->pkd, *pidSwap);
    }


uint64_t _pstOrdSplit(PST pst,uint64_t iMinOrder,uint64_t iMaxOrder) {
    struct outFreeStore outFree;
    struct inOrdWeight inWt;
    struct outOrdWeight outWtLow,outWtHigh;
    uint64_t im,imm,il,iu;
    uint64_t nLowerStore,nUpperStore,nLow,nHigh;
    struct inColOrdRejects inCol;
    OREJ *pLowerRej,*pUpperRej;
    int *pidSwap,iRet,nOut,ittr;
    int rID;

    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    rID = mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
    pstFreeStore(pst->pstLower,NULL,0,&outFree,NULL);
    nLowerStore = outFree.nFreeStore;
    mdlGetReply(pst->mdl,rID,&outFree,NULL);
    nUpperStore = outFree.nFreeStore;
    /*
    ** Find the correct iOrdSplit, such that all processors will
    ** have close to the same number of particles in the end.
    ** Start the ROOT finder based on balancing number of particles.
    */ 
    il = iMinOrder;
    iu = iMaxOrder;
    im = 0;        /* just initialization */
    nLow = 0;      /* just initialization */
    nHigh = iu+1;  /* just initialization */
    imm = (il + iu + 1)/2;
    ittr = 0;
    while (il < imm && imm < iu && ittr < MAX_ITTR) {
	im = imm;
	inWt.iOrdSplit = im;
	inWt.ittr = ittr;
	inWt.iSplitSide = 1;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORDWEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	pstOrdWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	mdlGetReply(pst->mdl,rID,&outWtHigh,NULL);
	/*
	** Add lower and Upper subsets weights and numbers
	*/
	nLow = outWtLow.nLow + outWtHigh.nLow;
	nHigh = outWtLow.nHigh + outWtHigh.nHigh;
	/*
	  printf("ittr:%d l:%d u:%d : %llu < %llu < %llu\n",
	      ittr,nLow,nHigh, il, imm, iu);
	*/
	if (nLow == 1 && nHigh == 1) /* break on trivial case */
	    break;
	else {		/* split on number */
	    if (nLow/(double)pst->nLower >
		    nHigh/(double)pst->nUpper) iu = im;
	    else if (nLow/(double)pst->nLower <
		     nHigh/(double)pst->nUpper) il = im;
	    else break;
	    }
	imm = (il + iu)/2;
	++ittr;
	}
    mdlassert(pst->mdl,nLow <= nLowerStore);
    mdlassert(pst->mdl,nHigh <= nUpperStore);
    pst->iOrdSplit = im;
    /*
    ** Collect rejects.
    */
    pLowerRej = malloc(pst->nLower*sizeof(OREJ));
    mdlassert(pst->mdl,pLowerRej != NULL);
    pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
    mdlassert(pst->mdl,pUpperRej != NULL);
    pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
    mdlassert(pst->mdl,pidSwap != NULL);
    inCol.iOrdSplit = pst->iOrdSplit;
    inCol.iSplitSide = 1;
    rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,&inCol,
		  sizeof(inCol));
    inCol.iSplitSide = 0;
    pstColOrdRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
    mdlGetReply(pst->mdl,rID,pUpperRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
    while (1) {
	iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,pst->nUpper,
			    pUpperRej,pidSwap);
	if (!iRet) break;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
		      mdlThreads(pst->mdl)*sizeof(int));
	pstSwapRejects(pst->pstLower,pidSwap,mdlThreads(pst->mdl)*sizeof(int),
		       pLowerRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,rID,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
	}
    free(pLowerRej);
    free(pUpperRej);
    free(pidSwap);
    return pst->iOrdSplit;
    }


void pstDomainOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inDomainOrder *in = vin;
    int rID=0;

    mdlassert(pst->mdl,nIn == sizeof(struct inDomainOrder));
    if (pst->nLeaves > 1) {
	uint64_t iMinOrder = in->iMinOrder;
	uint64_t iMaxOrder = in->iMaxOrder;
	uint64_t iMidOrder;
	iMidOrder = _pstOrdSplit(pst,iMinOrder,iMaxOrder);
	/*
	** Now go on to Domain Order of next levels.
	*/
	in->iMinOrder = iMidOrder;
	if (pst->nUpper > 1) rID = mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINORDER,in,nIn);
	in->iMinOrder = iMinOrder;
	in->iMaxOrder = iMidOrder-1;
	if (pst->nLower > 1) pstDomainOrder(pst->pstLower,in,nIn,NULL,NULL);
	if (pst->nUpper > 1) mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    if (pnOut) *pnOut = 0;
    }


void pstLocalOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inDomainOrder *in = vin;
    LCL *plcl = pst->plcl;
    uint64_t iMinOrder = in->iMinOrder;
    uint64_t iMaxOrder = in->iMaxOrder;

    mdlassert(pst->mdl,nIn == sizeof(struct inDomainOrder));
    if (pst->nLeaves > 1) {
	uint64_t iMinOrder = in->iMinOrder;
	uint64_t iMaxOrder = in->iMaxOrder;
	uint64_t iMidOrder = pst->iOrdSplit;
	assert(iMidOrder >= iMinOrder && iMidOrder <= iMaxOrder);
	in->iMinOrder = iMidOrder;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LOCALORDER,in,nIn);
	in->iMinOrder = iMinOrder;
	in->iMaxOrder = iMidOrder-1;
	pstLocalOrder(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdLocalOrder(plcl->pkd,iMinOrder,iMaxOrder);
	}
    if (pnOut) *pnOut = 0;
    }


void pstActiveOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    uint64_t *pnActive = vout;
    uint64_t nActiveLeaf;

    mdlassert(pst->mdl,nIn == 0);
    /*
      pkdStartTimer(pst->plcl->pkd,5);
    */
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEORDER,NULL,0);
	pstActiveOrder(pst->pstLower,NULL,0,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveOrder(plcl->pkd);
	}
    if (pnOut) *pnOut = sizeof(uint64_t);
    /*
      pkdStopTimer(pst->plcl->pkd,5);
    */
    }

void pstAddWriteStart(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inAddWriteStart *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inAddWriteStart));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ADDWRITESTART,in,nIn);
	pstAddWriteStart(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	plcl->nWriteStart += in->nWriteStart;
	}
    if (pnOut) *pnOut = 0;
    }

void pstCompressASCII(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCompressASCII *in = vin;
    struct outCompressASCII *out = vout;
    struct outCompressASCII outUp;
    struct inAddWriteStart inAdd;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inCompressASCII));
    if (pst->nLeaves > 1) {
	/* Instruct all processors to compress their data, and count it */
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_COMPRESSASCII,in,nIn);
	pstCompressASCII(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUp,NULL);
	/* Add the new offset from the left children to the right children */
	inAdd.nWriteStart = out->nBytes;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_ADDWRITESTART,&inAdd,sizeof(inAdd));
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	out->nBytes += outUp.nBytes;
	}
    else {
	plcl->pkdout = pkdStartOutASCII(plcl->pkd,in->iFile,in->iType);
	if (pst->idSelf==0 && in->iDim==0) {
	    pkdOutHdr(plcl->pkd,plcl->pkdout,in->nTotal);
	    }
	pkdOutASCII(plcl->pkd,plcl->pkdout,in->iType,in->iDim);
	pkdFinishOutASCII(plcl->pkd,plcl->pkdout);
	out->nBytes = pkdCountOutASCII(plcl->pkd,plcl->pkdout);
	plcl->nWriteStart = 0;
	}
    if (pnOut) *pnOut = sizeof(struct outCompressASCII);
    }

void pstWriteASCII(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inWriteASCII *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inWriteASCII));
    if (pst->nLeaves > 1) {
	/* Instruct all processors to compress their data, and count it */
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WRITEASCII,in,nIn);
	pstWriteASCII(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	FILE *fp = fopen (in->achOutFile,"r+");

	mdlassert(plcl->pkd->mdl,fp != NULL);
	setvbuf(fp,NULL,_IOFBF,PKDOUT_BUFFER_SIZE);
#ifdef HAVE_FSEEKO
	fseeko(fp,in->nFileOffset+plcl->nWriteStart,SEEK_SET);
#else
	assert(0);
#endif
	pkdDumpOutASCII(plcl->pkd,plcl->pkdout,fp);
	pkdFreeOutASCII(plcl->pkd,plcl->pkdout);
	fclose(fp);

	}
    if (pnOut) *pnOut = 0;
    }

static void makeName( char *achOutName, const char *inName, int iIndex ) {
    char *p;

    strcpy( achOutName, inName );
    p = strstr( achOutName, "&I" );
    if ( p ) {
	int n = p - achOutName;
	sprintf( p, "%d", iIndex );
	strcat( p, inName + n + 2 );
	}
    else {
	p = achOutName + strlen(achOutName);
	sprintf(p,".%d", iIndex);
	}
    }

void pstWrite(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    char achOutFile[PST_FILENAME_SIZE];
    struct inWrite *in = vin;
    FIO fio;
    uint32_t nCount;
    int i;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inWrite));

    if (pst->nLeaves > 1) {
	int nProcessors = in->nProcessors;
	/* If we are the writer (nProcessors==1) or a producer (==1) keep requesting */
	if (nProcessors<=1) {
	    in->nProcessors = 0;
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WRITE,in,nIn);
	    in->nProcessors = nProcessors;
	    pstWrite(pst->pstLower,in,nIn,vout,pnOut);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	/* Split writing tasks between child nodes */
	else {
	    int nLeft = nProcessors / 2;
	    int iUpper = in->iUpper;
	    assert(in->iLower == mdlSelf(pst->mdl));
	    in->iLower = pst->idUpper;
	    in->nProcessors -= nLeft;
	    in->iIndex += nLeft;
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WRITE,in,nIn);
	    in->iLower = mdlSelf(pst->mdl);
	    in->iUpper = pst->idUpper;
	    in->iIndex -= nLeft;
	    in->nProcessors = nLeft;
	    pstWrite(pst->pstLower,in,nIn,vout,pnOut);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	}
    else {
	LCL *plcl = pst->plcl;
	if (in->nProcessors==0) {
	    pkdWriteViaNode(plcl->pkd, in->iLower);
	    }
	else {
//	    if (in->bStandard==2) {
//		makeName(achOutFile,in->achOutFile,in->iIndex);
//		fio = fioGadgetCreate(achOutFile,in->mFlags,in->dTime,in->dBoxSize,
//		    in->Omega0,in->OmegaLambda,in->HubbleParam,
//		    int nTypes, const uint64_t *nPart,
//		    int nFiles, const uint64_t *nAll,
//		    const double *dMass );
//		}
	    if (in->bHDF5) {
#ifdef USE_HDF5
		makeName(achOutFile,in->achOutFile,in->iIndex);
		fio = fioHDF5Create(achOutFile,in->mFlags);
#else
		fio = NULL; /* Should never happen */
#endif
		}
	    else {
		if (strstr(in->achOutFile, "&I" )) {
		    makeName(achOutFile,in->achOutFile,in->iIndex);
		    fio = fioTipsyCreatePart(achOutFile,0,in->mFlags&FIO_FLAG_CHECKPOINT,
			in->bStandard, in->dTime, 
			in->nSph, in->nDark, in->nStar, plcl->nWriteStart);
		    }
		else {
		    fio = fioTipsyAppend(in->achOutFile,in->mFlags&FIO_FLAG_CHECKPOINT,in->bStandard);
		    if (fio) {
			fioSeek(fio,plcl->nWriteStart,FIO_SPECIES_ALL);
			}
		    }
		}
	    if (fio==NULL) {
		fprintf(stderr,"ERROR: unable to create file for output\n");
		perror(in->achOutFile);
		mdlassert(pst->mdl,fio!=NULL);
		}
	    fioSetAttr(fio, "dTime",    FIO_TYPE_DOUBLE, &in->dTime);
	    /* Restart information */
	    fioSetAttr(fio, "dEcosmo",  FIO_TYPE_DOUBLE, &in->dEcosmo );
	    fioSetAttr(fio, "dTimeOld", FIO_TYPE_DOUBLE, &in->dTimeOld );
	    fioSetAttr(fio, "dUOld",    FIO_TYPE_DOUBLE, &in->dUOld );

	    pkdWriteFIO(plcl->pkd,fio,in->dvFac,&in->bnd);
	    for(i=in->iLower+1; i<in->iUpper; ++i ) {
		pkdWriteFromNode(plcl->pkd,i,fio,in->dvFac,&in->bnd);
		}
	    fioClose(fio);
	    }
	}

    if (pnOut) *pnOut = 0;
    }


void pstSetSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSetSoft *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSetSoft));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETSOFT,in,nIn);
	pstSetSoft(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSetSoft(plcl->pkd,in->dSoft);
	}
    if (pnOut) *pnOut = 0;
    }

void pstDumpTrees(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    mdlassert(pst->mdl,nIn == 0);

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DUMPTREES,vin,nIn);
	pstDumpTrees(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdDumpTrees(pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstBuildTree(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inBuildTree *in = vin;
    TREESPEC *pSpec = (TREESPEC *)(in+1);
    KDN *pkdn = vout;
    KDN *pCell1, *pCell2, *pCell;
    FLOAT minside;
    uint32_t *puCells;
    int iCell,iLower;
    int i;

    /* We need to save our cells so we can update them later */
    puCells = stack_alloc( in->nTrees * sizeof(uint32_t));
    for(i=0; i<in->nTrees; ++i) puCells[i] = pSpec[i].uCell;

    if (pst->nLeaves > 1) {
	/* Get our cell ready */
	for(i=0; i<in->nTrees; ++i) {
	    iCell = pSpec[i].uCell;
	    pCell = pkdTreeNode(pkd,iCell);
	    pCell->bRemote = 1;                          /* Lower sibling is remote */
	    pCell->iLower = pkdTreeAllocNode(pkd);       /* Lower cell */
	    pCell->pLower = pSpec[i].uRoot;              /* Remote root node */
	    pCell->pUpper = pst->idUpper;                /* Remote processor */
	    }
	pCell1 = stack_alloc(in->nTrees*pkdNodeSize(plcl->pkd));
	pCell2 = stack_alloc(in->nTrees*pkdNodeSize(plcl->pkd));

	/* The right of the trees starts in a new domain */
	for(i=0; i<in->nTrees; ++i) pSpec[i].uCell = pSpec[i].uRoot;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,in,nIn);

	/* We descend down the left -- we will get back the combined cell */
	for(i=0; i<in->nTrees; ++i) pSpec[i].uCell = pkdTreeNode(pkd,puCells[i])->iLower;
	pstBuildTree(pst->pstLower,in,nIn,pCell1,pnOut);
	mdlGetReply(pst->mdl,rID,pCell2,pnOut);

	/*
	** Combine Cell1 and Cell2 into pCell
	** to find cell CoM, bounds and multipoles.
	** This also computes the opening radius for gravity.
	*/
	for(i=0; i<in->nTrees; ++i) {
	    pCell = pkdTreeNode(pkd,puCells[i]);
	    pCell->bMax = HUGE_VAL;  /* initialize bMax for CombineCells */
	    MINSIDE(pst->bnd.fMax,minside);
	    pkdCombineCells1(pkd,pCell,pkdNode(pkd,pCell1,i),pkdNode(pkd,pCell2,i));
	    CALCOPEN(pCell,minside);
	    pkdCombineCells2(pkd,pCell,pkdNode(pkd,pCell1,i),pkdNode(pkd,pCell2,i));
	    }
	stack_free(pCell1);
	stack_free(pCell2);
	}
    else {
	pkdTreeAlignNode(pkd);
	pkdTreeBuild(plcl->pkd,in->nBucket,in->nTrees,pSpec);
	/* This cell will now have bRemote=0. pLower and pUpper are now particle indexes. */
	}
    /*
    ** Calculated all cell properties, now pass up this cell info.
    */
    if (pkdn != NULL) {
	for(i=0; i<in->nTrees; ++i) {
	    pkdCopyNode(pkd,pkdNode(pkd,pkdn,i),pkdTreeNode(pkd,puCells[i]));
	    }
	if (pnOut) *pnOut = in->nTrees * pkdNodeSize(plcl->pkd);
	}
    else if (pnOut) *pnOut = 0;
    stack_free(puCells);
    }

void pstCalcRoot(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct ioCalcRoot *out = vout;
    struct ioCalcRoot temp;

    mdlassert(pst->mdl,nIn == 3*sizeof(double) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCROOT,vin,nIn);
	pstCalcRoot(pst->pstLower,vin,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&temp,NULL);
	momAddMomc(&out->momc,&temp.momc);
	}
    else {
	double *com = vin;
	pkdCalcRoot(plcl->pkd,com,&out->momc);
	}
    if (pnOut) *pnOut = sizeof(struct ioCalcRoot);
    }

void pstDistribRoot(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct ioDistribRoot *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct ioDistribRoot));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBROOT,vin,nIn);
	pstDistribRoot(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdDistribRoot(plcl->pkd,in->r,&in->momc);
	}
    if (pnOut) *pnOut = 0;
    }


void pstEnforcePeriodic(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    BND *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(BND));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ENFORCEPERIODIC,vin,nIn);
	pstEnforcePeriodic(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdEnforcePeriodic(plcl->pkd,in);
	}
    if (pnOut) *pnOut = 0;
    }

void pstPhysicalSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inPhysicalSoft *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPhysicalSoft));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PHYSICALSOFT,in,nIn);
	pstPhysicalSoft(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdPhysicalSoft(plcl->pkd,in->dSoftMax,in->dFac,in->bSoftMaxMul);
	}
    if (pnOut) *pnOut = 0;
    }

void pstHopLink(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inHopLink *in = vin;
    uint64_t *nOutGroups = (uint64_t *)vout;
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopLink));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_LINK,in,nIn);
	pstHopLink(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&nOutUpper,NULL);
	*nOutGroups += nOutUpper;
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;
	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
	*nOutGroups = smHopLink(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = sizeof(uint64_t);
    }

void pstHopJoin(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inHopLink *in = vin;
    struct outHopJoin *out = vout;
    struct outHopJoin outUpper, outLower;
    int nOut, nLocal;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopLink));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_JOIN,in,nIn);
        pstHopJoin(pst->pstLower,in,nIn,&outLower,pnOut);
        mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
        assert(nOut==sizeof(struct outHopJoin));
	out->bDone = outLower.bDone && outUpper.bDone;
	out->nGroups = outLower.nGroups + outUpper.nGroups;
        }
    else {
	LCL *plcl = pst->plcl;
	SMX smx;
	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
        out->bDone = smHopJoin(smx,&in->smf,in->dHopTau,&nLocal);
	smFinish(smx,&in->smf);
	out->nGroups = nLocal;
        }
    if (pnOut) *pnOut = sizeof(struct outHopJoin);
    }

void pstHopFinishUp(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inHopFinishUp *in = (struct inHopFinishUp *)vin;
    uint64_t *nOutGroups = (uint64_t *)vout;
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopFinishUp));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_FINISH_UP,vin,nIn);
        pstHopFinishUp(pst->pstLower,vin,nIn,vout,pnOut);
        mdlGetReply(pst->mdl,rID,&nOutUpper,pnOut);
	*nOutGroups += nOutUpper;
        }
    else {
	LCL *plcl = pst->plcl;
        *nOutGroups = pkdHopFinishUp(plcl->pkd,in->nMinGroupSize,in->bPeriodic,in->fPeriod);
        }
    if (pnOut) *pnOut = sizeof(uint64_t);
    }

void pstHopTreeBuild(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inHopTreeBuild *in = (struct inHopTreeBuild *)vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inHopTreeBuild));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_TREE_BUILD,vin,nIn);
        pstHopTreeBuild(pst->pstLower,vin,nIn,NULL,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,pnOut);
        }
    else {
	LCL *plcl = pst->plcl;
        pkdHopTreeBuild(plcl->pkd,in->nBucket);
        }
    if (pnOut) *pnOut = 0;
    }

void pstHopGravity(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inHopGravity *in = (struct inHopGravity *)vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inHopGravity));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_GRAVITY,vin,nIn);
        pstHopGravity(pst->pstLower,vin,nIn,NULL,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,pnOut);
        }
    else {
	LCL *plcl = pst->plcl;
	double dFlop=0, dPartSum=0, dCellSum=0;
	pkdGravWalkHop(plcl->pkd,in->dTime,in->nGroup,in->dThetaMin,&dFlop,&dPartSum,&dCellSum);
        }
    if (pnOut) *pnOut = 0;
    }

void pstHopUnbind(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inHopUnbind *in = (struct inHopUnbind *)vin;
    struct outHopUnbind *out = (struct outHopUnbind *)vout;
    struct outHopUnbind outUpper, outLower;
    mdlassert(pst->mdl,nIn == sizeof(struct inHopUnbind));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_UNBIND,vin,nIn);
        pstHopUnbind(pst->pstLower,vin,nIn,&outLower,pnOut);
        mdlGetReply(pst->mdl,rID,&outUpper,pnOut);
	out->nEvaporated = outLower.nEvaporated + outUpper.nEvaporated;
	out->nGroups = outLower.nGroups + outUpper.nGroups;
        }
    else {
	LCL *plcl = pst->plcl;
        out->nEvaporated = pkdHopUnbind(plcl->pkd,in->dTime,in->nMinGroupSize,in->bPeriodic,in->fPeriod);
	out->nGroups = plcl->pkd->nLocalGroups;
        }
    if (pnOut) *pnOut = sizeof(struct outHopUnbind);
    }

void pstHopAssignGID(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_ASSIGN_GID,vin,nIn);
        pstHopAssignGID(pst->pstLower,vin,nIn,NULL,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,pnOut);
        }
    else {
	LCL *plcl = pst->plcl;
        pkdHopAssignGID(plcl->pkd);
        }
    if (pnOut) *pnOut = 0;
    }

void pstHopSendStats(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    pkdHopSendStats(pst->plcl->pkd);
    if (pnOut) *pnOut = 0;
    }

void pstSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SMOOTH,in,nIn);
	pstSmooth(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;

	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
	smSmooth(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }


void pstFastGasPhase1(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FASTGASPHASE1,in,nIn);
	pstFastGasPhase1(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;

	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
	smFastGasPhase1(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }


void pstFastGasPhase2(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FASTGASPHASE2,in,nIn);
	pstFastGasPhase2(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;

	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
	smFastGasPhase2(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }


void pstFastGasCleanup(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FASTGASCLEANUP,NULL,0);
	pstFastGasCleanup(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdFastGasCleanup(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }


void pstReSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTH,in,nIn);
	pstReSmooth(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;

	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
	smReSmooth(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }

void pstGravity(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inGravity *in = vin;
    struct outGravity *out = vout;
    struct outGravity *outUp = out + pst->idUpper-pst->idSelf;

    mdlassert(pst->mdl,nIn == sizeof(struct inGravity));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GRAVITY,in,nIn);
	pstGravity(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,outUp,NULL);
	}
    else {
	pkdGravAll(plcl->pkd,in->uRungLo,in->uRungHi,in->bKickClose,in->bKickOpen,
	    in->dtClose,in->dtOpen,in->dAccFac,in->dTime,in->nReps,in->bPeriodic,
	    4,in->bEwald,in->nGroup,in->iRoot1,in->iRoot2,in->dEwCut,in->dEwhCut, in->dThetaMin,&out->nActive,
	    &out->dPartSum,&out->dCellSum,&out->cs,&out->dFlop,&out->uRungMax);
	out->nLocal = plcl->pkd->nLocal;
	out->dWalkTime = pkdGetWallClockTimer(plcl->pkd,1);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	out->dComputing     = mdlTimeComputing(pst->mdl);
	out->dWaiting       = mdlTimeWaiting(pst->mdl);
	out->dSynchronizing = mdlTimeSynchronizing(pst->mdl);
#endif
	}
    if (pnOut) *pnOut = pst->nLeaves*sizeof(struct outGravity);
    }


void pstCalcEandL(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outCalcEandL *out = vout;
    struct outCalcEandL outLcl;
    int k;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCEANDL,NULL,0);
	pstCalcEandL(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,rID,&outLcl,NULL);
	out->T += outLcl.T;
	out->U += outLcl.U;
	out->Eth += outLcl.Eth;
	for (k=0;k<3;k++) out->L[k] = outLcl.L[k];
	for (k=0;k<3;k++) out->F[k] = outLcl.F[k];
	out->W += outLcl.W;
	}
    else {
	pkdCalcEandL(plcl->pkd,&out->T,&out->U,&out->Eth,out->L,out->F,&out->W);
	}
    if (pnOut) *pnOut = sizeof(struct outCalcEandL);
    }


void pstDrift(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inDrift *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DRIFT,in,nIn);
	pstDrift(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdDrift(plcl->pkd,in->dDelta,in->dDeltaVPred,in->dDeltaUPred,in->uRungLo,in->uRungHi);
	}
    if (pnOut) *pnOut = 0;
    }

void pstScaleVel(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inScaleVel *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inScaleVel));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SCALEVEL,in,nIn);
	pstScaleVel(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdScaleVel(plcl->pkd,in->dvFac);
	}
    if (pnOut) *pnOut = 0;
    }

void pstCacheBarrier(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	pstCacheBarrier(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	mdlCacheBarrier(pst->mdl,CID_CELL);
	}
    if (pnOut) *pnOut = 0;
    }

void pstStepVeryActiveKDK(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inStepVeryActive *in = vin;
    struct outStepVeryActive *out = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inStepVeryActive));
    if (pst->nLeaves > 1) {
	if (pst->iVASplitSide > 0) {
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	    pstStepVeryActiveKDK(pst->pstLower,in,nIn,out,NULL);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	else if (pst->iVASplitSide < 0) {
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_STEPVERYACTIVE,in,nIn);
	    pstCacheBarrier(pst->pstLower,NULL,0,NULL,NULL);
	    mdlGetReply(pst->mdl,rID,out,NULL);
	    }
	else {
	    mdlassert(pst->mdl,pst->iVASplitSide != 0);
	    }
	}
    else {
//	assert(plcl->pkd->nVeryActive > 0);

	out->nMaxRung = in->nMaxRung;
	pkdStepVeryActiveKDK(plcl->pkd,in->uRungLo,in->uRungHi,in->dStep,in->dTime,in->dDelta,
			     in->iRung, in->iRung, in->iRung, 0, in->dThetaMin,
			     &out->nMaxRung, in->aSunInact, in->adSunInact,
			     in->dSunMass);
	mdlCacheBarrier(pst->mdl,CID_CELL);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

void pstROParticleCache(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ROPARTICLECACHE,vin,nIn);
	pstROParticleCache(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	/*
	** Start particle caching space.
	*/
	PKD pkd = plcl->pkd;
	mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),
		   pkdLocal(pkd));

	}
    if (pnOut) *pnOut = 0;
    }

void pstParticleCacheFinish(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PARTICLECACHEFINISH,vin,nIn);
	pstParticleCacheFinish(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	/*
	** Start particle caching space.
	*/
	PKD pkd = plcl->pkd;
	mdlFinishCache(pkd->mdl,CID_PARTICLE);
	}
    if (pnOut) *pnOut = 0;
    }

void pstKick(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inKick *in = vin;
    struct outKick *out = vout;
    struct outKick outUp;

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_KICK,in,nIn);
	pstKick(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&outUp,NULL);

	out->SumTime += outUp.SumTime;
	out->nSum += outUp.nSum;
	if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
	}
    else {
	pkdKick(plcl->pkd,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaU,in->dDeltaUPred,in->uRungLo,in->uRungHi);
	out->Time = pkdGetTimer(plcl->pkd,1);
	out->MaxTime = out->Time;
	out->SumTime = out->Time;
	out->nSum = 1;
	}
    if (pnOut) *pnOut = sizeof(struct outKick);
    }

void pstKickTree(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inKickTree *in = vin;
    struct outKickTree *out = vout;
    struct outKickTree outUp;

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_KICKTREE,in,nIn);
	pstKickTree(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&outUp,NULL);

	out->SumTime += outUp.SumTime;
	out->nSum += outUp.nSum;
	if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
	}
    else {
	pkdKickTree(plcl->pkd,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaU,in->dDeltaUPred,in->iRoot);
	out->Time = pkdGetTimer(plcl->pkd,1);
	out->MaxTime = out->Time;
	out->SumTime = out->Time;
	out->nSum = 1;
	}
    if (pnOut) *pnOut = sizeof(struct outKickTree);
    }

void pstSetTotal(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outSetTotal *out = vout;
    struct outSetTotal oute;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETTOTAL,NULL,0);
	pstSetTotal(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,rID,&oute,NULL);
	out->nTotal += oute.nTotal;
	pst->nTotal = out->nTotal;
	}
    else {
	/*pst->nTotal = pkdLocal(plcl->pkd);*/
	pst->nTotal = pkdNumSrcActive(plcl->pkd,0,MAX_RUNG);
	out->nTotal = pst->nTotal;
	}
    if (pnOut) *pnOut = sizeof(struct outSetTotal);
    }


void pstSetWriteStart(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSetWriteStart *in = vin;
    uint64_t nWriteStart;

    mdlassert(pst->mdl,nIn == sizeof(struct inSetWriteStart));
    nWriteStart = in->nWriteStart;
    if (pst->nLeaves > 1) {
	in->nWriteStart = nWriteStart + pst->pstLower->nTotal;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETWRITESTART,in,nIn);
	in->nWriteStart = nWriteStart;
	pstSetWriteStart(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	plcl->nWriteStart = nWriteStart;
	}
    if (pnOut) *pnOut = 0;
    }


void
pstSetRung(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSetRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETRUNG,vin,nIn);
	pstSetRung(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSetRung(plcl->pkd,in->uRungLo,in->uRungHi,in->uRung);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstZeroNewRung(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inZeroNewRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ZERONEWRUNG,vin,nIn);
	pstZeroNewRung(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdZeroNewRung(plcl->pkd,in->uRungLo,in->uRungHi,in->uRung);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstInitStep(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inInitStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITSTEP,vin,nIn);
	pstInitStep(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdInitStep(plcl->pkd,&in->param,&in->csm);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstActiveRung(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inActiveRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVERUNG,vin,nIn);
	pstActiveRung(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdActiveRung(plcl->pkd, in->iRung, in->bGreater);
	}
    }

void
pstCurrRung(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCurrRung *in = vin;
    struct outCurrRung *out = vout;
    int iCurrent;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CURRRUNG,vin,nIn);
	pstCurrRung(pst->pstLower,vin,nIn,vout,pnOut);
	iCurrent = out->iCurrent;
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	if (iCurrent)
	    out->iCurrent = iCurrent;
	}
    else {
	out->iCurrent = pkdCurrRung(plcl->pkd, in->iRung);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

void
pstAccelStep(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inAccelStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inAccelStep));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ACCELSTEP,vin,nIn);
	pstAccelStep(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdAccelStep(plcl->pkd,in->uRungLo,in->uRungHi,in->dEta,in->dVelFac,in->dAccFac,
		     in->bDoGravity,in->bEpsAcc,in->dhMinOverSoft);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstSphStep(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSphStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSphStep));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SPHSTEP,vin,nIn);
	pstSphStep(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSphStep(plcl->pkd,in->uRungLo,in->uRungHi,in->dAccFac);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstStarForm(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inStarForm *in = vin;
    struct outStarForm *out = vout;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inStarForm));
    if (pst->nLeaves > 1) {
	struct outStarForm fsStats;
	
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_STARFORM,in,nIn);
	pstStarForm(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&fsStats,NULL);
	out->nFormed += fsStats.nFormed;
	out->nDeleted += fsStats.nDeleted;
	out->dMassFormed += fsStats.dMassFormed;
		}
    else {
	pkdStarForm(pst->plcl->pkd, in->dRateCoeff, in->dTMax, in->dDenMin, in->dDelta, 
		    in->dTime,
		     in->dInitStarMass, in->dESNPerStarMass, in->dtCoolingShutoff,
		    in->dtFeedbackDelay,    in->dMassLossPerStarMass,    
		    in->dZMassPerStarMass,    in->dMinGasMass,
		    in->bdivv,
		     &out->nFormed, &out->dMassFormed, &out->nDeleted);
	}
    
    if (pnOut) *pnOut = sizeof(struct outStarForm);
    }

void
pstDensityStep(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inDensityStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DENSITYSTEP,vin,nIn);
	pstDensityStep(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdDensityStep(plcl->pkd,in->uRungLo,in->uRungHi,in->dEta,in->dRhoFac);
	}
    if (pnOut) *pnOut = 0;
    }

#ifdef COOLING
void pstCoolSetup(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inCoolSetup *in = vin;
    
    if (pst->nLeaves > 1) {
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_COOLSETUP,in,nIn);
	pstCoolSetup(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
#if defined(COOLDEBUG)
	(plcl->pkd->Cool)->mdl = plcl->pkd->mdl;
#endif
	CoolSetup((plcl->pkd->Cool),in->dGmPerCcUnit, in->dComovingGmPerCcUnit, in->dErgPerGmUnit, in->dSecUnit, in->dKpcUnit, in->dOmega0, in->dHubble0, in->dLambda, in->dOmegab, in->dOmegaRad, in->a, in->z, in->dTime, in->CoolParam);
    }
    if (pnOut) *pnOut = 0;
    }

void pstCooling(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inCooling *in = vin;
    struct outCooling *out = vout;
    struct outCooling outUp;
    
    mdlassert(pst->mdl,nIn == sizeof(struct inCooling));
    
    if (pst->nLeaves > 1) {
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_COOLING,in,nIn);
	pstCooling(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&outUp,NULL);
	
	out->SumTime += outUp.SumTime;
	out->nSum += outUp.nSum;
	if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
	}
    else {
	pkdCooling(plcl->pkd,in->dTime,in->z,in->bUpdateState,in->bUpdateTable,in->bIterateDt,in->bIsothermal);
	out->Time = pkdGetTimer(plcl->pkd,1);
	out->MaxTime = out->Time;
	out->SumTime = out->Time;
	out->nSum = 1;
	}

    if (pnOut) *pnOut = sizeof(struct outCooling);
    }
#endif

void pstCorrectEnergy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inCorrectEnergy *in = vin;
    
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CORRECTENERGY,vin,nIn);
	pstCorrectEnergy(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdCorrectEnergy(plcl->pkd,in->dTuFac,in->z,in->dTime,in->iDirection);
	}
    if (pnOut) *pnOut = 0;
    }

void pstUpdateRung(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outUpdateRung outTemp;
    struct inUpdateRung *in = vin;
    struct outUpdateRung *out = vout;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_UPDATERUNG,vin,nIn);
	pstUpdateRung(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outTemp,pnOut);
	for (i=0;i<in->uMaxRung;++i) {
	    out->nRungCount[i] += outTemp.nRungCount[i];
	    }
	}
    else {
	pkdUpdateRung(plcl->pkd,in->uRungLo,in->uRungHi,
		      in->uMinRung,in->uMaxRung,out->nRungCount);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

void pstUpdateRungByTree(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inUpdateRungByTree *in = vin;
    struct outUpdateRung *out = vout;
    struct outUpdateRung outTemp;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_UPDATERUNGBYTREE,vin,nIn);
	pstUpdateRungByTree(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outTemp,pnOut);
	for (i=0;i<in->uMaxRung;++i) {
	    out->nRungCount[i] += outTemp.nRungCount[i];
	    }
	}
    else {
	pkdUpdateRungByTree(plcl->pkd,in->iRoot,in->uMinRung,in->uMaxRung,out->nRungCount);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

void pstSetRungVeryActive(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSetRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETRUNGVERYACTIVE,vin,nIn);
	pstSetRungVeryActive(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSetRungVeryActive(plcl->pkd,in->uRung);
	}
    }


void
pstColNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outColNParts *out = vout;
    struct outColNParts *outUp = out + pst->idUpper-pst->idSelf;
    int i;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLNPARTS,vin,nIn);
	pstColNParts(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,outUp,pnOut);
	}
    else {
	pkdColNParts(plcl->pkd, &out->nNew,
		     &out->nDeltaGas,
		     &out->nDeltaDark,
		     &out->nDeltaStar);
	}
    if (pnOut) *pnOut = pst->nLeaves*sizeof(*out);
    }

void
pstNewOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    uint64_t *in = vin;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_NEWORDER,vin,nIn);
	pstNewOrder(pst->pstLower, vin, nIn, vout, pnOut);
	mdlGetReply(pst->mdl, rID, vout, pnOut);
	}
    else {
	pkdNewOrder(plcl->pkd, in[pst->idSelf]);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstGetNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct outGetNParts *out = vout;
    
    if(pst->nLeaves > 1) {
	struct outGetNParts outtmp;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GETNPARTS,vin,nIn);
	pstGetNParts(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,(void *) &outtmp,pnOut);
	
	out->n += outtmp.n;
	out->nGas += outtmp.nGas;
	out->nDark += outtmp.nDark;
	out->nStar += outtmp.nStar;
	if (outtmp.nMaxOrder > out->nMaxOrder) out->nMaxOrder = outtmp.nMaxOrder;
	}
    else {
	pkdGetNParts(pst->plcl->pkd, out);
	}
    if(pnOut) *pnOut = sizeof(struct outGetNParts);
}

void
pstSetNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inSetNParts *in = vin;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETNPARTS,vin,nIn);
	pstSetNParts(pst->pstLower, vin, nIn, vout, pnOut);
	mdlGetReply(pst->mdl, rID, vout, pnOut);
	}
    else {
	pkdSetNParts(pst->plcl->pkd, in->nGas, in->nDark, in->nStar);
	}
    if (pnOut) *pnOut = 0;
    }


void
pstClearTimer(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inClearTimer *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inClearTimer));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CLEARTIMER,in,nIn);
	pstClearTimer(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdClearTimer(pst->plcl->pkd,in->iTimer);
	}
    if (pnOut) *pnOut = 0;
    }


void pstFof(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inFof *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inFof));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FOF,in,nIn);
	pstFof(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;
	(&in->smf)->pkd = pst->plcl->pkd;
	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
	smFof(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }

void pstGroupMerge(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inGroupMerge *in = vin;
    int *nGroups = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inGroupMerge));
    if (pst->nLeaves > 1) {
	int nGroupsLeaf;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GROUPMERGE,in,nIn);
	pstGroupMerge(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&nGroupsLeaf,pnOut);
	*nGroups += nGroupsLeaf;
	}
    else {      
	(&in->smf)->pkd = pst->plcl->pkd;
	*nGroups = smGroupMerge(&in->smf,in->bPeriodic);
	}
    if (pnOut) *pnOut = sizeof(int);
    }

void pstGroupProfiles(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inGroupProfiles *in = vin;
    int *nBins = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inGroupProfiles));
    if (pst->nLeaves > 1) {
	int nBinsLeaf;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GROUPPROFILES,in,nIn);
	pstGroupProfiles(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&nBinsLeaf,pnOut);
	*nBins = nBinsLeaf;
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;
	(&in->smf)->pkd = pst->plcl->pkd;
	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType);
	*nBins = smGroupProfiles(smx,&in->smf,in->nTotalGroups);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = sizeof(int);
    }

void pstInitRelaxation(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITRELAXATION,vin,nIn);
	pstInitRelaxation(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdInitRelaxation(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

#ifdef MDL_FFTW
void pstInitializePStore(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inInitializePStore *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inInitializePStore));
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITIALIZEPSTORE,in,nIn);
	pstInitializePStore(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdInitialize(
	    &plcl->pkd,pst->mdl,in->nStore,in->nBucket,in->nGroup,
	    in->nTreeBitsLo,in->nTreeBitsHi,
	    in->iCacheSize,in->iWorkQueueSize,in->iCUDAQueueSize,in->fPeriod,
	    in->nDark,in->nGas,in->nStar,in->mMemoryModel, in->nDomainRungs);
	}
    }

void pstGetFFTMaxSizes(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inGetFFTMaxSizes *in = vin;
    struct outGetFFTMaxSizes *out = vout;
    struct outGetFFTMaxSizes outUp;
    uint64_t nTotal, nLocal, nStore;

    mdlassert(pst->mdl,nIn == sizeof(struct inGetFFTMaxSizes));
    mdlassert(pst->mdl,vout != NULL);
    assert(mdlCore(pst->mdl)==0);

    if (pstOffNode(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GETFFTMAXSIZES,in,nIn);
	pstGetFFTMaxSizes(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUp,pnOut);
	if (outUp.nMaxLocal > out->nMaxLocal) out->nMaxLocal = outUp.nMaxLocal;
	if (outUp.nMaxZ > out->nMaxZ) out->nMaxZ = outUp.nMaxZ;
	if (outUp.nMaxY > out->nMaxY) out->nMaxY = outUp.nMaxZ;
	}
    else {
	assert(pstAmNode(pst));
	out->nMaxLocal = mdlFFTlocalCount(pst->mdl,in->nx,in->ny,in->nz,
	    &out->nMaxZ,0,&out->nMaxY,0);
	}
    if (pnOut) *pnOut = 0;
    }

void pltConstructIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inConstructIC *in = vin;
    struct outConstructIC *out = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inConstructIC));
    assert(pstOnNode(pst)); /* We pass around pointers! */
    if (pstNotCore(pst)) {
	int iBegYr = in->iBegYr;
	int iBegZk = in->iBegZk;

	in->iBegYr = (in->iBegYr + in->iEndYr) / 2;
	in->iBegZk = (in->iBegZk + in->iEndZk) / 2;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CONSTRUCTIC,in,nIn);
	in->iEndYr = in->iBegYr;
	in->iEndZk = in->iBegZk;
	in->iBegYr = iBegYr;
	in->iBegZk = iBegZk;
	pltConstructIC(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdGenerateIC(plcl->pkd,in->iSeed,in->dBoxSize,in->omegam,in->omegav,in->a,
	    in->fft,in->iBegYr,in->iEndYr,in->iBegZk,in->iEndZk,in->dic,in->pos);
//	out->nLocal = mdlFFTlocalCount(pst->mdl,in->nGrid,in->nGrid,in->nGrid,0,0,0,0);
	}
    if (pnOut) *pnOut = 0;
    }

void pstGenerateIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inGenerateIC *in = vin;
    struct outGenerateIC *out = vout;
    uint64_t nTotal, nLocal;
    struct inConstructIC cic;

    mdlassert(pst->mdl,nIn == sizeof(struct inGenerateIC));
    mdlassert(pst->mdl,vout != NULL);

    if (pstOffNode(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GENERATEIC,in,nIn);
	pstGenerateIC(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	assert(pstAmNode(pst)); /* We are "on node" here */
	nTotal = in->nGrid + 2; /* Careful: 32 bit integer cubed => 64 bit integer */
	nTotal *= in->nGrid;
	nTotal *= in->nGrid;
	nLocal = nTotal / mdlThreads(pst->mdl);
	in->ps.nStore = nLocal + (int)ceil(nLocal*in->fExtraStore);

	/* Adjust upward if necessary */
	if (in->ps.nStore*mdlCores(pst->mdl) < in->nPerNode) {
	    in->ps.nStore = (in->nPerNode + mdlCores(pst->mdl) - 1) / mdlCores(pst->mdl);
	    }

	/*
	** Allocate memory for particle/grid store. We require 9 double arrays
	** and a PARTICLE is normally at least 10 when calculating gravity.
	*/
	pstInitializePStore(pst,&in->ps,sizeof(in->ps),NULL,NULL);
	assert( 9*sizeof(double) <= pkdParticleSize(plcl->pkd) );
	assert(3*sizeof(double) == sizeof(gridpos));

	cic.pos = (gridpos *)pkdParticleBase(plcl->pkd);
	cic.dic[0].r = (double *)(cic.pos + in->nPerNode);
	cic.dic[1].r = cic.dic[0].r + in->nPerNode;
	cic.dic[2].r = cic.dic[1].r + in->nPerNode;
	cic.dic[3].r = cic.dic[2].r + in->nPerNode;
	cic.dic[4].r = cic.dic[3].r + in->nPerNode;
	cic.dic[5].r = cic.dic[4].r + in->nPerNode;
	cic.vel = (gridpos *)cic.dic[3].r; /* We overlap here */

	assert(cic.dic[8].r+in->nPerNode <= (double *)pkdParticle(plcl->pkd,in->nPerNode));

	mdlFFTInitialize(pst->mdl,&cic.fft,in->nGrid,in->nGrid,in->nGrid,0,cic.dic[0].r);

	cic.iBegYr = 0;
	cic.iEndYr = cic.fft->rgrid->n2;
	cic.iBegZk = 0;
	cic.iEndZk = cic.fft->rgrid->n3;
	cic.iSeed = in->iSeed;
	cic.dBoxSize = in->dBoxSize;
	cic.omegav = in->omegav;
	cic.omegam = in->omegab + in->omegac;
	cic.a = in->dExpansion;

	pltConstructIC(pst,&cic,sizeof(cic),NULL,NULL);

	// generate noise into dic1, dic2, dic3

	// copy to dic4, dic5, dic6

	// fft to real dic4, dic5, dic6 -> this is the first order






	out->dExpansion = in->dExpansion;
	}
#if 0
	/* Okay, here we set it to 1/50 of the interparticle separation */
	fSoft = 1.0 / (50.0 * in->nGrid);
	/* Mass is easy */
	fMass = in->omegac / (1.0 * in->nGrid * in->nGrid * in->nGrid );

	dx = in->dBoxSize / in->nGrid;
//	graficInitialize( &gctx, dx, 1, in->iSeed, in->h*100,
//			  in->omegac, in->omegab, in->omegav,
//			  in->nGrid, in->nGrid, in->nGrid );
//	out->dExpansion = graficGetExpansionFactor(gctx);

	for ( d=1; d<=3; d++ ) {
//	    pkdGenerateIC( plcl->pkd, gctx, d,
//			   fSoft, fMass, in->bComove );
	    }
//	graficFinish(gctx);
	}
#endif
    if (pnOut) *pnOut = sizeof(struct outGenerateIC);
    }
#endif

void pstHostname(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct outHostname *out = vout;
    struct outHostname *outUp = out + pst->idUpper-pst->idSelf;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOSTNAME,vin,nIn);
	pstHostname(pst->pstLower,vin,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,outUp,NULL);
	}
    else {
	char *p;
	out->iMpiID = mdlSelf(pst->mdl);
	strncpy(out->szHostname,mdlName(pst->mdl),sizeof(out->szHostname));
	out->szHostname[sizeof(out->szHostname)-1] = 0;
	p = strchr(out->szHostname,'.');
	if (p) *p = 0;
	}
    if (pnOut) *pnOut = pst->nLeaves*sizeof(struct outHostname);
    }

void pstMemStatus(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outMemStatus *out = vout;
    struct outMemStatus *outUp = out + pst->idUpper-pst->idSelf;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEMSTATUS,vin,nIn);
	pstMemStatus(pst->pstLower,vin,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,outUp,NULL);
	}
    else {
	FILE *fp;
	char buffer[512], *save, *f;
	int i;

#ifdef __linux__
	fp = fopen("/proc/self/stat","r");
	if ( fp != NULL ) {
	    fgets(buffer,sizeof(buffer),fp);
	    fclose(fp);
	    f = strtok_r(buffer," ",&save);
	    for ( i=0; i<= 36 && f; i++ ) {
		switch (i) {
		case  9: out->minflt = atol(f); break;
		case 10: out->minflt+= atol(f); break;
		case 11: out->majflt = atol(f); break;
		case 12: out->majflt+= atol(f); break;
		case 22:
		    out->vsize  = atol(f)/1024/1024;
		    break;
		case 23:
#ifdef HAVE_GETPAGESIZE
		    out->rss    = atol(f)*getpagesize()/1024/1024;
#else
		    out->rss    = atol(f)/1024;
#endif
		    break;
		default: break;
		    }
		f = strtok_r(NULL," ",&save);
		}
	    }
#endif
	out->nBytesTree = pkdTreeMemory(plcl->pkd);
	out->nBytesCl   = pkdClMemory(plcl->pkd);
	out->nBytesIlp  = pkdIlpMemory(plcl->pkd);
	out->nBytesIlc  = pkdIlcMemory(plcl->pkd);

	}
    if (pnOut) *pnOut = pst->nLeaves*sizeof(struct outMemStatus);
    }


void pstGetClasses(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PARTCLASS *out = vout;
    PARTCLASS *outUp;
    int nUp;
    int n, i, j;

    mdlassert(pst->mdl,nIn==0);
    if (pst->nLeaves > 1) {
	outUp = malloc(PKD_MAX_CLASSES*sizeof(PARTCLASS));
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GETCLASSES,vin,nIn);
	pstGetClasses(pst->pstLower,vin,nIn,out,pnOut);
	mdlGetReply(pst->mdl,rID,outUp,&nUp);
	n = nUp / sizeof(PARTCLASS);
	mdlassert(pst->mdl,n*sizeof(PARTCLASS)==nUp);
	nUp = n;
	n = *pnOut / sizeof(PARTCLASS);
	mdlassert(pst->mdl,n*sizeof(PARTCLASS)== *pnOut);
	for ( i=0; i<nUp; i++ ) {
	    for ( j=0; j<n; j++ ) {
		if ( outUp[i].fMass==out[j].fMass && outUp[i].fSoft==out[j].fSoft && outUp[i].eSpecies==out[j].eSpecies)
		    break;
		}
	    if ( j == n ) {
		out[n++] = outUp[i];
		}
	    }
	free(outUp);
	*pnOut = n * sizeof(PARTCLASS);
	}
    else {
	n = pkdGetClasses(plcl->pkd,PKD_MAX_CLASSES,vout);
	*pnOut = n*sizeof(PARTCLASS);
	}
    }

void pstSetClasses(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PARTCLASS *in = vin;
    int n;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETCLASSES,vin,nIn);
	pstSetClasses(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	n = nIn / sizeof(PARTCLASS);
	mdlassert(pst->mdl,n*sizeof(PARTCLASS)==nIn);
	pkdSetClasses(plcl->pkd,n,in,1);
	}
    if (pnOut) *pnOut = 0;
    }

/*
 * Routine to swap the class table .  Note that this does not walk
 * the pst but simply set's the given table and returns the old table.
 */
void pstSwapClasses(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl;
    PARTCLASS *in = vin;
    PARTCLASS *out = vout;
    int n;
    PST lpst;

    lpst = pst;
    while (lpst->nLeaves > 1)
	lpst = lpst->pstLower;
    plcl = lpst->plcl;

    n = pkdGetClasses( plcl->pkd, PKD_MAX_CLASSES, out );
    *pnOut = n * sizeof(PARTCLASS);

    n = nIn / sizeof(PARTCLASS);
    mdlassert(pst->mdl,n*sizeof(PARTCLASS) == nIn);
    pkdSetClasses( plcl->pkd, n, in, 0 );
    }


void pstSelSrcAll(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCALL,vin,nIn);
	pstSelSrcAll(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelSrcAll(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }
void pstSelDstAll(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTALL,vin,nIn);
	pstSelDstAll(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelDstAll(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelSrcGas(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCGAS,vin,nIn);
	pstSelSrcGas(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelSrcGas(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelDstGas(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTGAS,vin,nIn);
	pstSelDstGas(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelDstGas(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelSrcStar(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCSTAR,vin,nIn);
	pstSelSrcStar(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelSrcStar(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelDstStar(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelDstStar *in = vin;
    assert(nIn==sizeof(struct inSelDstStar));

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTSTAR,vin,nIn);
	pstSelDstStar(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelDstStar(plcl->pkd,in->bFB,in->dTimeFB);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelSrcDeleted(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCDELETED,vin,nIn);
	pstSelSrcDeleted(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelSrcDeleted(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelDstDeleted(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTDELETED,vin,nIn);
	pstSelDstDeleted(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelDstDeleted(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelSrcById(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelById *in = vin;
    struct outSelById *out = vout;
    struct outSelById outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelById) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCBYID,vin,nIn);
	pstSelSrcById(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelById));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSrcById(plcl->pkd,in->idStart,in->idEnd,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelById);
    }
void pstSelDstById(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelById *in = vin;
    struct outSelById *out = vout;
    struct outSelById outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelById) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTBYID,vin,nIn);
	pstSelDstById(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelById));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelDstById(plcl->pkd,in->idStart,in->idEnd,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelById);
    }
void pstSelSrcMass(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelMass *in = vin;
    struct outSelMass *out = vout;
    struct outSelMass outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelMass) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCMASS,vin,nIn);
	pstSelSrcMass(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelMass));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSrcMass(plcl->pkd,in->dMinMass,in->dMaxMass,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelMass);
    }
void pstSelDstMass(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelMass *in = vin;
    struct outSelMass *out = vout;
    struct outSelMass outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelMass) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTMASS,vin,nIn);
	pstSelDstMass(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelMass));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelDstMass(plcl->pkd,in->dMinMass,in->dMaxMass,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelMass);
    }

void pstSelSrcPhaseDensity(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelPhaseDensity *in = vin;
    struct outSelPhaseDensity *out = vout;
    struct outSelPhaseDensity outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelPhaseDensity) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCPHASEDENSITY,vin,nIn);
	pstSelSrcPhaseDensity(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelPhaseDensity));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSrcPhaseDensity(plcl->pkd,in->dMinDensity,in->dMaxDensity,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelPhaseDensity);
    }

void pstSelDstPhaseDensity(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelPhaseDensity *in = vin;
    struct outSelPhaseDensity *out = vout;
    struct outSelPhaseDensity outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelPhaseDensity) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTPHASEDENSITY,vin,nIn);
	pstSelDstPhaseDensity(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelPhaseDensity));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelDstPhaseDensity(plcl->pkd,in->dMinDensity,in->dMaxDensity,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelPhaseDensity);
    }

void pstSelSrcBox(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelBox *in = vin;
    struct outSelBox *out = vout;
    struct outSelBox outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelBox) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCBOX,vin,nIn);
	pstSelSrcBox(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelBox));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSrcBox(
	    plcl->pkd,in->dCenter,in->dSize,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelBox);
    }
void pstSelDstBox(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelBox *in = vin;
    struct outSelBox *out = vout;
    struct outSelBox outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelBox) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTBOX,vin,nIn);
	pstSelDstBox(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelBox));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelDstBox(
	    plcl->pkd,in->dCenter,in->dSize,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelBox);
    }



void pstSelSrcSphere(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelSphere *in = vin;
    struct outSelSphere *out = vout;
    struct outSelSphere outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelSphere) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCSPHERE,vin,nIn);
	pstSelSrcSphere(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelSphere));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSrcSphere(
	    plcl->pkd,in->r,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelSphere);
    }

void pstSelDstSphere(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelSphere *in = vin;
    struct outSelSphere *out = vout;
    struct outSelSphere outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelSphere) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTSPHERE,vin,nIn);
	pstSelDstSphere(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelSphere));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelDstSphere(
	    plcl->pkd,in->r,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelSphere);
    }

void pstSelSrcCylinder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelCylinder *in = vin;
    struct outSelCylinder *out = vout;
    struct outSelCylinder outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelCylinder) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCCYLINDER,vin,nIn);
	pstSelSrcCylinder(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelCylinder));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSrcCylinder(
	    plcl->pkd,in->dP1,in->dP2,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelCylinder);
    }

void pstSelDstCylinder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelCylinder *in = vin;
    struct outSelCylinder *out = vout;
    struct outSelCylinder outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelCylinder) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTCYLINDER,vin,nIn);
	pstSelDstCylinder(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelCylinder));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelDstCylinder(
	    plcl->pkd,in->dP1,in->dP2,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelCylinder);
    }

void pstSelSrcGroup(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    int *in = vin;
    LCL *plcl = pst->plcl;
    assert( nIn==sizeof(int) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSRCGROUP,vin,nIn);
	pstSelSrcGroup(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelSrcGroup(plcl->pkd, *in);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelDstGroup(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    int *in = vin;
    LCL *plcl = pst->plcl;
    assert( nIn==sizeof(int) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDSTGROUP,vin,nIn);
	pstSelDstGroup(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelDstGroup(plcl->pkd, *in);
	}
    if (pnOut) *pnOut = 0;
    }

void pstDeepestPot(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inDeepestPot *in = vin;
    struct outDeepestPot *out = vout;
    struct outDeepestPot outUpper;
    int nOut;

    assert( nIn==sizeof(struct inDeepestPot) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DEEPESTPOT,vin,nIn);
	pstDeepestPot(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outDeepestPot));
	if ( out->nChecked==0 || (outUpper.nChecked && outUpper.fPot < out->fPot) ) {
	    out->r[0] = outUpper.r[0];
	    out->r[1] = outUpper.r[1];
	    out->r[2] = outUpper.r[2];
	    out->fPot = outUpper.fPot;
	    out->nChecked += outUpper.nChecked;
	    }
	}
    else {
	out->nChecked = pkdDeepestPot(plcl->pkd,in->uRungLo,in->uRungHi,
	    out->r,&out->fPot);
	}
    if (pnOut) *pnOut = sizeof(struct outDeepestPot);
    }

void pstProfile(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inProfile *in = vin;
    /*assert( nIn==sizeof(struct inProfile) );*/
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PROFILE,vin,nIn);
	pstProfile(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdProfile(plcl->pkd,in->uRungLo,in->uRungHi,
		   in->dCenter, in->dRadii, in->nBins,
		   in->com, in->vcm, in->L);
	}
    if (pnOut) *pnOut = 0;
    }

void pstCalcDistance(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCalcDistance *in = vin;

    assert( nIn==sizeof(struct inCalcDistance) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCDISTANCE,vin,nIn);
	pstCalcDistance(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdCalcDistance(plcl->pkd,in->dCenter);
	}
    if (pnOut) *pnOut = 0;
    }

void pstCalcCOM(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCalcCOM *in = vin;
    struct outCalcCOM *out = vout;
    struct outCalcCOM outUpper;
    int i;

    assert( nIn==sizeof(struct inCalcCOM) );
    assert( vout != NULL );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCCOM,vin,nIn);
	pstCalcCOM(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,pnOut);
	out->N += outUpper.N;
	out->M += outUpper.M;
	for(i=0; i<3; i++ ) {
	    out->com[i] += outUpper.com[i];
	    out->vcm[i] += outUpper.vcm[i];
	    out->L[i] += outUpper.L[i];
	    }
	}
    else {
	pkdCalcCOM(plcl->pkd,in->dCenter,in->dRadius,
		   out->com, out->vcm, out->L, &out->M, &out->N);
	}
    if (pnOut) *pnOut = sizeof(struct outCalcCOM);
    }

void pstCountDistance(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCountDistance *in = vin;
    struct outCountDistance *out = vout;
    struct outCountDistance outUpper;
    int nOut;

    assert( nIn==sizeof(struct inCountDistance) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COUNTDISTANCE,vin,nIn);
	pstCountDistance(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut==sizeof(struct outCountDistance));
	out->nCount += outUpper.nCount;
	}
    else {
	out->nCount = pkdCountDistance(plcl->pkd,in->dRadius2Inner,in->dRadius2Outer);
	}
    if (pnOut) *pnOut = sizeof(struct outCountDistance);
    }

void pstInitGrid(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inInitGrid *in = vin;
    int rID;
    int s = in->s;
    int n = in->n;
    assert( sizeof(struct inInitGrid) == nIn );
    if (pst->nLeaves > 1) {

	in->s = s + pst->nLower*n/pst->nLeaves;
	in->n = s + n - in->s;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITGRID,vin,nIn);
	in->n = n - in->n;
	in->s = in->n ? s : 0;
	pstInitGrid(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdGridInitialize(plcl->pkd,in->n1,in->n2,in->n3,in->a1,in->s,in->n);
	}
    if (pnOut) *pnOut = 0;
    }

void pstGridProject(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    /*struct inGridProject *in = vin;*/
    assert( sizeof(struct inGridProject) == nIn );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GRIDPROJECT,vin,nIn);
	pstGridProject(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdGridProject(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

#ifdef MDL_FFTW
void pstMeasurePk(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inMeasurePk *in = vin;
    struct outMeasurePk *out = vout;
    struct outMeasurePk outUpper;
    int nOut;
    int i;

    assert( nIn==sizeof(struct inMeasurePk) );
    if (pstOffNode(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEASUREPK,vin,nIn);
	pstMeasurePk(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut==sizeof(struct outMeasurePk));

	for(i=0;i<=in->nGrid/2; i++) {
	    out->fPower[i] += outUpper.fPower[i];
	    out->nPower[i] += outUpper.nPower[i];
	    }
	}
    else {
	pkdMeasurePk(plcl->pkd, in->dCenter, in->dRadius,
		     in->nGrid, out->fPower, out->nPower);
	}

    if (pstAmNode(pst)) {
	}


    if (pnOut) *pnOut = sizeof(struct outMeasurePk);
    }
#endif

void pstTotalMass(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outTotalMass *out = vout;
    struct outTotalMass outUpper;
    int nOut;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_TOTALMASS,vin,nIn);
	pstTotalMass(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut==sizeof(struct outTotalMass));
	out->dMass += outUpper.dMass;
	}
    else {
	out->dMass = pkdTotalMass(plcl->pkd);
	}
    if (pnOut) *pnOut = sizeof(struct outTotalMass);
    }

#ifdef USE_PSD
void pstBuildPsdTree(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inPSD *in = vin;
    KDN *pkdn = vout;
    KDN *ptmp, *pCell;
    int i,iCell,iLower,iNext;
    BND *bnd, *p1bnd, *p2bnd;

    mdlassert(pst->mdl,nIn == sizeof(struct inPSD));
    iCell = in->iCell;
    pCell = pkdNode(pkd,pkdn,iCell);
    if (pst->nLeaves > 1) {
        in->iCell = UPPER(iCell);
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_BUILDPSDTREE,in,nIn);
        in->iCell = LOWER(iCell);
        pstBuildPsdTree(pst->pstLower,in,nIn,vout,pnOut);
        in->iCell = iCell;
        ptmp = malloc(in->nCell*pkdNodeSize(pkd));
        mdlassert(pst->mdl,ptmp != NULL);
        mdlGetReply(pst->mdl,rID,ptmp,pnOut);
        for (i=1;i<in->nCell;++i) {
	    KDN *pSrc = pkdNode(pkd,ptmp,i);
	    if (pSrc->pUpper) {
		pkdCopyNode(pkd,pkdNode(pkd,pkdn,i),pSrc);
		}
	    }
        free(ptmp);
        /*
        ** Combine to find cell bounds.
        */
        iLower = LOWER(iCell);
        iNext = UPPER(iCell);

        bnd = pkdNodeBnd(pkd, pCell);
        p1bnd = pkdNodeBnd(pkd, pkdNode(pkd,pkdn,iLower));
        p2bnd = pkdNodeBnd(pkd, pkdNode(pkd,pkdn,iNext));
        BND_COMBINE(bnd,p1bnd,p2bnd);

        bnd = pkdNodeVBnd(pkd, pCell);
        p1bnd = pkdNodeVBnd(pkd, pkdNode(pkd,pkdn,iLower));
        p2bnd = pkdNodeVBnd(pkd, pkdNode(pkd,pkdn,iNext));
        BND_COMBINE(bnd,p1bnd,p2bnd);


        /*
        ** Set all the pointers and flags.
        */
        pCell->iLower = iLower;
        pCell->iParent = 0;
        pCell->pLower = -1;
        pCell->pUpper = 1;
        pkdNode(pkd,pkdn,iLower)->iParent = iCell;
        pkdNode(pkd,pkdn,iNext)->iParent = iCell;
        }
    else {
        for (i=1;i<in->nCell;++i) pkdNode(pkd,pkdn,i)->pUpper = 0; /* used flag = unused */

        psdBuildTree(pkd,pkd->psx,in,pCell);

        pCell->iLower = 0;
        pCell->pLower = pst->idSelf;
        pCell->pUpper = 1;
        }
    /*
    ** Calculated all cell properties, now pass up this cell info.
    */
    if (pnOut) *pnOut = in->nCell*pkdNodeSize(plcl->pkd);
    }

void pstPSD(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inPSD *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPSD));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD,in,nIn);
        pstPSD(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        PKD pkd = pst->plcl->pkd;
        psdSmooth(pkd,pkd->psx);
        }
    if (pnOut) *pnOut = 0;
    }

void pstPSDLink(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inPSD *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPSD));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSDLINK,in,nIn);
        pstPSDLink(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        PKD pkd = pst->plcl->pkd;
        psdSmoothLink(pkd,pkd->psx);
        }
    if (pnOut) *pnOut = 0;
    }

void pstPSDInit(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inPSD *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPSD));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_INIT,in,nIn);
        pstPSDInit(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        PKD pkd = pst->plcl->pkd;
        pkd->psx = calloc(1, sizeof(*pkd->psx));
        assert(pkd->psx != NULL);
        psdInitialize(pkd, pkd->psx,in->nSmooth, in->bPeriodic);
        }
    if (pnOut) *pnOut = 0;
    }

void pstPSDJoinBridges(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    /*struct inPSD *in = vin;*/
    int *done = vout;
    int doneUpper;
    int nOut;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_JOINBRIDGES,NULL,0);
        pstPSDJoinBridges(pst->pstLower,NULL,0,vout,pnOut);
        mdlGetReply(pst->mdl,rID,&doneUpper,&nOut);
        assert(nOut==sizeof(doneUpper));
        *done = *done && doneUpper;
        }
    else {
        PKD pkd = pst->plcl->pkd;
        *done = psdJoinBridges(pkd, pkd->psx);
        }
    if (pnOut) *pnOut = sizeof(*done);
    }

void pstPSDCountLocalGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    int *in = vin;
    struct inAssignGlobalIds *out = vout;
    struct inAssignGlobalIds *tmp;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
#ifdef HAVE_ALLOCA
        tmp = alloca(mdlThreads(pst->mdl)*sizeof(*tmp));
#else
        tmp = malloc(mdlThreads(pst->mdl)*sizeof(*tmp));
#endif
        mdlassert(pst->mdl,tmp != NULL);
        pstPSDCountLocalGroups(pst->pstLower,in,nIn,vout,pnOut);
        *in += pst->nLower;
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_CLG,in,nIn);
        mdlGetReply(pst->mdl,rID,tmp,pnOut);
        for (i=0;i<pst->nUpper;++i) {
        out[*in+i] = tmp[*in+i];
        }
        *in -= pst->nLower;
#ifndef HAVE_ALLOCA
        free(tmp);
#endif
        }
    else {
        PKD pkd = pst->plcl->pkd;
        assert(*in == pst->idSelf);
        out[*in].offs = 0;
        out[*in].count = psdCountLocalGroups(pkd);
        }
    if (pnOut) *pnOut = mdlThreads(pst->mdl)*sizeof(struct inAssignGlobalIds);
    }

void pstPSDAssignGlobalIds(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inAssignGlobalIds *in = vin;

    mdlassert(pst->mdl,nIn == mdlThreads(pst->mdl)*sizeof(struct inAssignGlobalIds));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_ASSIGN_GLOBAL_IDS,in,nIn);
        pstPSDAssignGlobalIds(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        PKD pkd = pst->plcl->pkd;
        struct inAssignGlobalIds *d = in+pkd->idSelf;
        psdAssignGlobalIds(pkd, d->offs, d->count);
        }
    if (pnOut) *pnOut = 0;
    }

void pstPSDMergeNoisyGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut) {

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_MERGENOISYGROUPS,NULL,0);
        pstPSDMergeNoisyGroups(pst->pstLower,NULL,0,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        PKD pkd = pst->plcl->pkd;
        psdMergeNoisyGroups(pkd, pkd->psx);
        }
    if (pnOut) *pnOut = 0;
    }

void pstPSDJoinGroupBridges(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    /*struct inPSD *in = vin;*/
    int *done = vout;
    int doneUpper;
    int nOut;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_JOINGROUPBRIDGES,NULL,0);
        pstPSDJoinGroupBridges(pst->pstLower,NULL,0,vout,pnOut);
        mdlGetReply(pst->mdl,rID,&doneUpper,&nOut);
        assert(nOut==sizeof(doneUpper));
        *done = *done && doneUpper;
        }
    else {
        PKD pkd = pst->plcl->pkd;
        *done = psdJoinGroupBridges(pkd,pkd->psx);
        }
    if (pnOut) *pnOut = sizeof(*done);
    }

void pstPSDSetGlobalId(PST pst,void *vin,int nIn,void *vout,int *pnOut) {

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_SETGLOBALID,NULL,0);
        pstPSDSetGlobalId(pst->pstLower,NULL,0,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        LCL *plcl = pst->plcl;
        psdSetGlobalId(plcl->pkd);
        }
    if (pnOut) *pnOut = 0;
    }

void pstPSDFinish(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inPSD *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPSD));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PSD_FINISH,in,nIn);
        pstPSDFinish(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        PKD pkd = pst->plcl->pkd;
        psdFinish(pkd,pkd->psx);
        }
    if (pnOut) *pnOut = 0;
    }

void pstWritePsGroups(PST pst,void *vin,int nIn,void *vout,int *pnOut) 
{
    struct inWritePsGroups *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    LCL *plcl = pst->plcl;
    pkdOutPsGroup(plcl->pkd, in->achOutFile, in->iType);
    if (pnOut) *pnOut = 0;
}

void pstUnbind(PST pst,void *vin,int nIn,void *vout,int *pnOut) 
{
    int *in = vin;
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_UNBIND,in,nIn);
        pstUnbind(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
        LCL *plcl = pst->plcl;
        ubUnbind(plcl->pkd, 32);
        }
    if (pnOut) *pnOut = 0;
}
#endif
