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

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
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
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#ifdef USE_ITT
#include "ittnotify.h"
#endif
#include "mdl.h"
#include "pst.h"
#include "pkd.h"
#include "smooth.h"
#include "hop.h"
#include "fof.h"
#include "group.h"
#include "groupstats.h"

void pstAddServices(PST pst,MDL mdl) {
    int nThreads;

    nThreads = mdlThreads(mdl);

    mdlAddService(mdl,PST_SETADD,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetAdd,
		  sizeof(struct inSetAdd),0);
    mdlAddService(mdl,PST_READFILE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstReadFile,
	          sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE),0);
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
    mdlAddService(mdl,PST_SENDPARTICLES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSendParticles,
		  sizeof(int),0);
    mdlAddService(mdl,PST_CHECKPOINT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCheckpoint,
		  sizeof(struct inWrite),0);
    mdlAddService(mdl,PST_OUTPUT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOutput,
		  sizeof(struct inOutput),0);
    mdlAddService(mdl,PST_OUTPUT_SEND,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOutputSend,
		  sizeof(struct inOutputSend),0);
    mdlAddService(mdl,PST_RESTORE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstRestore,
		  sizeof(struct inRestore),0);
    /*
    ** Calculate the number of levels in the top tree and use it to
    ** define the size of the messages.
    */
    mdlAddService(mdl,PST_BUILDTREE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstBuildTree,
	sizeof(struct inBuildTree),
	(nThreads==1?1:2*nThreads-1)*pkdMaxNodeSize());
    mdlAddService(mdl,PST_DISTRIBTOPTREE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDistribTopTree,
	sizeof(struct inDistribTopTree) + (nThreads==1?1:2*nThreads-1)*pkdMaxNodeSize(),0);
    mdlAddService(mdl,PST_DUMPTREES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDumpTrees,
	          sizeof(struct inDumpTrees),0);
    mdlAddService(mdl,PST_TREEINITMARKED,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstTreeInitMarked,
	          0,0);
    mdlAddService(mdl,PST_CALCROOT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcRoot,
	          sizeof(struct inCalcRoot),sizeof(struct outCalcRoot));
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
    mdlAddService(mdl,PST_GROUP_RELOCATE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGroupRelocate,
		  0,0);
    mdlAddService(mdl,PST_GROUP_STATS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGroupStats,
	          sizeof(struct inGroupStats),0);
    mdlAddService(mdl,PST_HOP_SEND_STATS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstHopSendStats,
		  0,0);
    mdlAddService(mdl,PST_SMOOTH,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSmooth,
		  sizeof(struct inSmooth),0);
#ifdef FAST_GAS
    mdlAddService(mdl,PST_FASTGASPHASE1,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFastGasPhase1,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_FASTGASPHASE2,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFastGasPhase2,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_FASTGASCLEANUP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFastGasCleanup,
		  0,0);
#endif
    mdlAddService(mdl,PST_GRAVITY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGravity,
	          sizeof(struct inGravity),
	          nThreads*sizeof(struct outGravityPerProc) + sizeof(struct outGravityReduct));
    mdlAddService(mdl,PST_LIGHTCONE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstLightCone,
		  sizeof(struct inLightCone),0);
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
    mdlAddService(mdl,PST_COUNTRUNGS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCountRungs,
		  0,sizeof(struct outCountRungs));
    mdlAddService(mdl,PST_DENSITYSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDensityStep,
		  sizeof(struct inDensityStep),0);
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
    mdlAddService(mdl,PST_NEW_FOF,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstNewFof,
		  sizeof(struct inNewFof),0);
    mdlAddService(mdl,PST_FOF_PHASES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFofPhases,
	          0,sizeof(struct outFofPhases));
    mdlAddService(mdl,PST_FOF_FINISH_UP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFofFinishUp,
	          sizeof(struct inFofFinishUp),sizeof(uint64_t));
    mdlAddService(mdl,PST_INITRELAXATION,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitRelaxation,0,0);
    mdlAddService(mdl,PST_INITIALIZEPSTORE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitializePStore,
		  sizeof(struct inInitializePStore),0);
#ifdef MDL_FFTW
    mdlAddService(mdl,PST_GETFFTMAXSIZES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetFFTMaxSizes,
		  sizeof(struct inGetFFTMaxSizes),sizeof(struct outGetFFTMaxSizes));
    mdlAddService(mdl,PST_GENERATEIC,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGenerateIC,
		  sizeof(struct inGenerateIC),sizeof(struct outGenerateIC));
    mdlAddService(mdl,PLT_GENERATEIC,pst,
		  (void (*)(void *,void *,int,void *,int *)) pltGenerateIC,
		  sizeof(struct inGenerateICthread),sizeof(struct outGenerateIC));
    mdlAddService(mdl,PLT_MOVEIC,pst,
		  (void (*)(void *,void *,int,void *,int *)) pltMoveIC,
		  sizeof(struct inMoveIC),0);
    mdlAddService(mdl,PST_MOVEIC,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstMoveIC,
		  sizeof(struct inGenerateIC),0);
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
    mdlAddService(mdl,PST_SELALL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelAll,
		  0, 0 );
    mdlAddService(mdl,PST_SELGAS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelGas,
		  0, 0 );
    mdlAddService(mdl,PST_SELSTAR,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelStar,
		  0, 0 );
    mdlAddService(mdl,PST_SELDELETED,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelDeleted,
		  0, 0 );
    mdlAddService(mdl,PST_SELBYID,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelById,
		  sizeof(struct inSelById), sizeof(struct outSelById));
    mdlAddService(mdl,PST_SELMASS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelMass,
		  sizeof(struct inSelMass), sizeof(struct outSelMass));
    mdlAddService(mdl,PST_SELPHASEDENSITY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelPhaseDensity,
		  sizeof(struct inSelPhaseDensity), sizeof(struct outSelPhaseDensity));
    mdlAddService(mdl,PST_SELBOX,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelBox,
		  sizeof(struct inSelBox), sizeof(struct outSelBox));
    mdlAddService(mdl,PST_SELSPHERE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelSphere,
		  sizeof(struct inSelSphere), sizeof(struct outSelSphere));
    mdlAddService(mdl,PST_SELCYLINDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelCylinder,
		  sizeof(struct inSelCylinder), sizeof(struct outSelCylinder));
    mdlAddService(mdl,PST_SELGROUP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSelGroup,
		  sizeof(int), 0);
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
    mdlAddService(mdl,PST_GRID_CREATE_FFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGridCreateFFT,
		  sizeof(struct inGridCreateFFT), 0);
    mdlAddService(mdl,PST_GRID_DELETE_FFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGridDeleteFFT,
		  0, 0);
    mdlAddService(mdl,PST_MEASUREPK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstMeasurePk,
		  sizeof(struct inMeasurePk), sizeof(struct outMeasurePk));
    mdlAddService(mdl,PST_ASSIGN_MASS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstAssignMass,
		  sizeof(struct inAssignMass), 0);
    mdlAddService(mdl,PST_SETLINGRID, pst,
           (void (*)(void*, void*, int, void*, int*)) pstSetLinGrid,
           sizeof(struct inSetLinGrid), 0);
    mdlAddService(mdl,PST_MEASURELINPK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstMeasureLinPk,
		  sizeof(struct inMeasureLinPk), sizeof(struct outMeasureLinPk));
#endif
    mdlAddService(mdl,PST_TOTALMASS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstTotalMass,
		  0, sizeof(struct outTotalMass));
    mdlAddService(mdl,PST_LIGHTCONE_OPEN,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstLightConeOpen,
		  sizeof(struct inLightConeOpen), 0);
    mdlAddService(mdl,PST_LIGHTCONE_CLOSE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstLightConeClose,
		  sizeof(struct inLightConeClose), 0);
    mdlAddService(mdl,PST_LIGHTCONEVEL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstLightConeVel,
		  0,0);
    mdlAddService(mdl,PST_INFLATE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInflate,
	          sizeof(struct inInflate), 0);
    mdlAddService(mdl,PST_GET_PARTICLES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetParticles,
	          sizeof(uint64_t)*GET_PARTICLES_MAX,
	          sizeof(struct outGetParticles)*GET_PARTICLES_MAX );
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

static void initializePStore(PKD *ppkd,MDL mdl,struct inInitializePStore *in) {
    pkdInitialize(
	ppkd,mdl,in->nStore,in->nMinTotalStore,in->nMinEphemeral,in->nEphemeralBytes,
	in->nTreeBitsLo,in->nTreeBitsHi,
	in->iCacheSize,in->iWorkQueueSize,in->iCUDAQueueSize,in->fPeriod,
	in->nSpecies[FIO_SPECIES_DARK],in->nSpecies[FIO_SPECIES_SPH],in->nSpecies[FIO_SPECIES_STAR],
	in->mMemoryModel,in->bLightCone,in->bLightConeParticles);
    }

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
	if (plcl->pkd) pkdFinish(plcl->pkd);
	initializePStore(&plcl->pkd,pst->mdl,in);
	}
    }

void pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inReadFile *in = vin;
    int *pout = vout;
    uint64_t nFileStart,nFileEnd,nFileTotal,nFileSplit;
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
    uint64_t nNodeStart,nNodeEnd,nNodeTotal,nNodeSplit;
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
    double fLow,fHigh;
    double fl,fu,fm=-1,fmm;
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
    double fSplit,fLow,fHigh;
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

static void makeName( char *achOutName, const char *inName, int iIndex,const char *prefix ) {
    char *p;

    strcpy( achOutName, inName );
    p = strstr( achOutName, "&I" );
    if ( p ) {
	int n = p - achOutName;
	sprintf( p, "%s%d", prefix, iIndex );
	strcat( p, inName + n + 2 );
	}
    else {
	p = achOutName + strlen(achOutName);
	sprintf(p,".%s%d", prefix, iIndex);
	}
    }

void pstRestore(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inRestore *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inRestore));
    if (pstNotCore(pst)) {
	uint64_t nProcessors = in->nProcessors;
	if (nProcessors>1) { /* Keep going parallel */
	    int nLower, nUpper;
	    nLower = nProcessors * pst->nLower / pst->nLeaves;
	    if (nLower==0) nLower=1;
	    nUpper = nProcessors - nLower;
	    in->nProcessors = nUpper;
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RESTORE,in,nIn);
	    in->nProcessors = nLower;
	    pstRestore(pst->pstLower,in,nIn,vout,pnOut);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	else { /* Serialize these processors now */
	    pstRestore(pst->pstLower,in,nIn,vout,pnOut);
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RESTORE,in,nIn);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	}
    else {
	PKD pkd = pst->plcl->pkd;
	char achInFile[PST_FILENAME_SIZE];
	makeName(achInFile,in->achInFile,mdlSelf(pkd->mdl),"chk.");
	pkdRestore(pkd,achInFile);
	}
    }

void pstCheckpoint(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inWrite *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inWrite));
    if (pstNotCore(pst)) {
	uint64_t nProcessors = in->nProcessors;
	if (nProcessors>1) { /* Keep going parallel */
	    int nLower, nUpper;
	    nLower = nProcessors * pst->nLower / pst->nLeaves;
	    if (nLower==0) nLower=1;
	    nUpper = nProcessors - nLower;
	    in->nProcessors = nUpper;
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CHECKPOINT,in,nIn);
	    in->nProcessors = nLower;
	    pstCheckpoint(pst->pstLower,in,nIn,vout,pnOut);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	else { /* Serialize these processors now */
	    pstCheckpoint(pst->pstLower,in,nIn,vout,pnOut);
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CHECKPOINT,in,nIn);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	}
    else {
	PKD pkd = pst->plcl->pkd;
	char achOutFile[PST_FILENAME_SIZE];
	makeName(achOutFile,in->achOutFile,mdlSelf(pkd->mdl),"chk.");
	pkdCheckpoint(pkd,achOutFile);
	}
    }

void pstSendParticles(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    int iTo = *(int *)vin;
    pkdWriteViaNode(pst->plcl->pkd, iTo);
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
//	    in->nProcessors = 0;
//	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WRITE,in,nIn);
//	    in->nProcessors = nProcessors;
	    pstWrite(pst->pstLower,in,nIn,vout,pnOut);
//	    mdlGetReply(pst->mdl,rID,NULL,NULL);
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
//	if (in->nProcessors==0) {
//	    pkdWriteViaNode(plcl->pkd, in->iLower);
//	    }
//	else {
	if (in->nProcessors!=0) {
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
		makeName(achOutFile,in->achOutFile,in->iIndex,"");
		fio = fioHDF5Create(achOutFile,in->mFlags);
#else
		fio = NULL; /* Should never happen */
#endif
		}
	    else {
		if (strstr(in->achOutFile, "&I" )) {
		    makeName(achOutFile,in->achOutFile,in->iIndex,"");
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
		int rID = mdlReqService(pst->mdl,i,PST_SENDPARTICLES,&pst->idSelf,sizeof(pst->idSelf));
		pkdWriteFromNode(plcl->pkd,i,fio,in->dvFac,&in->bnd);
		mdlGetReply(pst->mdl,rID,NULL,NULL);
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
    struct inDumpTrees *in = (struct inDumpTrees *)vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inDumpTrees));

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DUMPTREES,vin,nIn);
	pstDumpTrees(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdDumpTrees(pkd,in->bOnlyVA,in->uRungDD);
	}
    if (pnOut) *pnOut = 0;
    }

void pstTreeInitMarked(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    mdlassert(pst->mdl,nIn == 0);

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_TREEINITMARKED,vin,nIn);
	pstTreeInitMarked(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdTreeInitMarked(pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstDistribTopTree(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    KDN *pTop = vin;
    struct inDistribTopTree *in = vin;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBTOPTREE,vin,nIn);
	pstDistribTopTree(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdDistribTopTree(pkd,in->uRoot,in->nTop,(KDN *)(in+1));
	}
    if (pnOut) *pnOut = 0;
    }

void pstBuildTree(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inBuildTree *in = vin;
    uint32_t uRoot = in->uRoot;
    KDN *pTop = vout;
    KDN *pCell1, *pCell2;
    double minside;
    int nOutUpper;
    int iLower;
    int i;

    /* We need to save our cells so we can update them later */
    if (pst->nLeaves > 1) {
	pCell1 = pkdNode(pkd,pTop,1);
	pCell2 = pkdNode(pkd,pTop,pst->nLower*2);

	/* We will accumulate the top tree here */
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,vin,nIn);
	pstBuildTree(pst->pstLower,vin,nIn,pCell1,pnOut);
	assert(*pnOut == (pst->nLower*2-1) * pkdNodeSize(pkd));
	mdlGetReply(pst->mdl,rID,pCell2,&nOutUpper);
	assert(nOutUpper == (pst->nUpper*2-1) * pkdNodeSize(pkd));
	*pnOut += nOutUpper + pkdNodeSize(pkd);

	/*
	** Combine Cell1 and Cell2 into pCell
	** to find cell CoM, bounds and multipoles.
	** This also computes the opening radius for gravity.
	*/
	pTop->bMax = HUGE_VAL;  /* initialize bMax for CombineCells */
	MINSIDE(pst->bnd.fMax,minside);
	pkdCombineCells1(pkd,pTop,pCell1,pCell2);
	CALCOPEN(pTop,minside);
	pkdCombineCells2(pkd,pTop,pCell1,pCell2);

	/* Get our cell ready */
	pTop->bTopTree = 1;                         /* The replicated top tree */
	pTop->bGroup = 0;                           /* top tree can never be a group */
	pTop->bRemote = 0;                          /* top tree is not remote */
	pTop->iLower = 1;                           /* Relative index to lower cell */
	pTop->pUpper = pst->nLower*2;               /* Relative index to upper cell */
	pTop->pLower = 0;                           /* Not used */
	}
    else {
	KDN *pRoot = pkdTreeNode(pkd,uRoot);
	pkdTreeAlignNode(pkd);
	pkdTreeBuild(plcl->pkd,in->nBucket,in->nGroup,in->uRoot,in->utRoot);
	pkdCopyNode(pkd,pTop,pRoot);
	/* Get our cell ready */
	pTop->bTopTree = 1;
	pTop->bGroup = 0;
	pTop->bRemote = 1;
	pTop->pUpper = pTop->pLower = pst->idSelf;
	/* iLower is valid = ROOT */
	*pnOut = pkdNodeSize(pkd);
	}
    }

void pstCalcRoot(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCalcRoot *in = vin;
    struct outCalcRoot *out = vout;
    struct outCalcRoot temp;

    mdlassert(pst->mdl,nIn == sizeof(struct inCalcRoot) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCROOT,vin,nIn);
	pstCalcRoot(pst->pstLower,vin,nIn,out,NULL);
	mdlGetReply(pst->mdl,rID,&temp,NULL);
	momAddMomc(&out->momc,&temp.momc);
	}
    else {
	pkdCalcRoot(plcl->pkd,in->uRoot,in->com,&out->momc);
	}
    if (pnOut) *pnOut = sizeof(struct outCalcRoot);
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
        pkdHopTreeBuild(plcl->pkd,in->nBucket,in->nGroup);
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

void pstGroupRelocate(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GROUP_RELOCATE,NULL,0);
        pstGroupRelocate(pst->pstLower,vin,nIn,NULL,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,pnOut);
        }
    else {
	LCL *plcl = pst->plcl;
	PKD pkd = plcl->pkd;
        pkdGroupRelocate(pkd,pkd->nGroups,pkd->ga);
        }
    if (pnOut) *pnOut = 0;
    }

void pstGroupStats(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inGroupStats *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inGroupStats));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GROUP_STATS,vin,nIn);
        pstGroupStats(pst->pstLower,vin,nIn,NULL,pnOut);
        mdlGetReply(pst->mdl,rID,NULL,pnOut);
        }
    else {
	LCL *plcl = pst->plcl;
        pkdCalculateGroupStats(plcl->pkd,in->bPeriodic,in->dPeriod,in->rEnvironment);
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


#ifdef FAST_GAS
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
#endif

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


void pstInitStat(STAT *ps,int id) {
    ps->idMax = id;
    ps->dMax = ps->dSum;
    ps->dSum2 = ps->dSum*ps->dSum;
    ps->n = 1;  /* sometimes we will need to zero this */
    }
    
void pstCombStat(STAT *ps,STAT *pa) {
    if (pa->dMax > ps->dMax) {
	ps->idMax = pa->idMax;
	ps->dMax = pa->dMax;
	}
    ps->dSum += pa->dSum;
    ps->dSum2 += pa->dSum2;
    ps->n += pa->n;
    }
    
void pstGravity(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inGravity *in = vin;
    struct outGravityPerProc *out = vout;
    struct outGravityPerProc *outup = out + pst->idUpper - pst->idSelf;
    struct outGravityPerProc *outend = out + mdlThreads(pst->mdl) - pst->idSelf;
    struct outGravityReduct *outr = (struct outGravityReduct *)outend;
    struct outGravityReduct tmp;
    int i;


    mdlassert(pst->mdl,nIn == sizeof(struct inGravity));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GRAVITY,in,nIn);
	pstGravity(pst->pstLower,in,nIn,out,NULL);
	/*
	** Make a temporary copy of the reduct part of the out buffer before setting it as the 
	** reply buffer. The reduct part follows at the end of all the outGravityPerProc entries.
	*/
	tmp = *outr;  /* copy the whole GravityReduct structure */
	mdlGetReply(pst->mdl,rID,outup,NULL);
	/*
	** Now combine in the tempory copy of the lower branch reduct part.
	*/
#ifdef __linux__
	pstCombStat(&outr->sRSS, &tmp.sRSS);
	pstCombStat(&outr->sFreeMemory,&tmp.sFreeMemory);
#endif
	pstCombStat(&outr->sLocal,&tmp.sLocal);
	pstCombStat(&outr->sActive,&tmp.sActive);
	pstCombStat(&outr->sPart,&tmp.sPart);
	pstCombStat(&outr->sPartNumAccess,&tmp.sPartNumAccess);
	pstCombStat(&outr->sPartMissRatio,&tmp.sPartMissRatio);
	pstCombStat(&outr->sCell,&tmp.sCell);
	pstCombStat(&outr->sCellNumAccess,&tmp.sCellNumAccess);
	pstCombStat(&outr->sCellMissRatio,&tmp.sCellMissRatio);
	pstCombStat(&outr->sFlop,&tmp.sFlop);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	pstCombStat(&outr->sComputing,&tmp.sComputing);
	pstCombStat(&outr->sWaiting,&tmp.sWaiting);
	pstCombStat(&outr->sSynchronizing,&tmp.sSynchronizing);
#endif
	/*
	** Combine the rung counts for the next timestep...
	*/
	for (i=in->uRungLo;i<=IRUNGMAX;++i) outr->nRung[i] += tmp.nRung[i];
	/*
	** and the number of actives processed in this gravity call.
	*/
	outr->nActive += tmp.nActive;
	outr->dFlopSingleCPU += tmp.dFlopSingleCPU;
	outr->dFlopDoubleCPU += tmp.dFlopDoubleCPU;
	outr->dFlopSingleGPU += tmp.dFlopSingleGPU;
	outr->dFlopDoubleGPU += tmp.dFlopDoubleGPU;
	}
    else {
#ifdef __linux__
	FILE *fp;
	char buffer[512], *save, *f, *v;
#endif
	PKD pkd = plcl->pkd;
	pkdGravAll(pkd,in->uRungLo,in->uRungHi,in->bKickClose,in->bKickOpen,
	    in->dtClose,in->dtOpen,in->dtLCDrift,in->dtLCKick,in->dLookbackFac,in->dLookbackFacLCP,
	    in->dAccFac,in->dTime,in->nReps,in->bPeriodic,
	    in->bEwald,in->nGroup,in->iRoot1,in->iRoot2,in->dEwCut,in->dEwhCut,in->dThetaMin,
	    in->bLinearSpecies,
	    &outr->nActive,
	    &outr->sPart.dSum,&outr->sPartNumAccess.dSum,&outr->sPartMissRatio.dSum,
	    &outr->sCell.dSum,&outr->sCellNumAccess.dSum,&outr->sCellMissRatio.dSum,
	    &outr->sFlop.dSum,outr->nRung);
	outr->dFlopSingleCPU = 1e-9*pkd->dFlopSingleCPU;
	outr->dFlopDoubleCPU = 1e-9*pkd->dFlopDoubleCPU;
	outr->dFlopSingleGPU = 1e-9*pkd->dFlopSingleGPU;
	outr->dFlopDoubleGPU = 1e-9*pkd->dFlopDoubleGPU;
	outr->sLocal.dSum = plcl->pkd->nLocal;
	outr->sActive.dSum = (double)outr->nActive;
#ifdef __linux__
	fp = fopen("/proc/self/stat","r");
	if ( fp != NULL ) {
	    if (fgets(buffer,sizeof(buffer),fp)==NULL) buffer[0] = '\0';
	    fclose(fp);
	    f = strtok_r(buffer," ",&save);
	    for ( i=0; i<= 36 && f; i++ ) {
		switch (i) {
		case 23:
#ifdef HAVE_GETPAGESIZE
		    outr->sRSS.dSum = atol(f)*1.0*getpagesize()/1024/1024/1024;
#else
		    outr->sRSS.dSum = atol(f)*1.0/1024/1024;
#endif
		    break;
		default: break;
		    }
		f = strtok_r(NULL," ",&save);
		}
	    }
	else outr->sRSS.dSum = 0;

	fp = fopen("/proc/meminfo","r");
	outr->sFreeMemory.dSum = 0;
	if ( fp != NULL ) {
	    while(fgets(buffer,sizeof(buffer),fp)) {
		f = strtok_r(buffer,":",&save);
		v = strtok_r(NULL," ",&save);
		if (strcmp(f,"MemFree")==0 || strcmp(f,"Buffers")==0 || strcmp(f,"Cached")==0) {
		    outr->sFreeMemory.dSum -= atol(v) * 1.0 / 1024 / 1024;
		    }
		}
	    fclose(fp);
	    }
#endif
	/*
	** Initialize statistics tracking variables assuming dSum is set.
	*/
#ifdef __linux__
	pstInitStat(&outr->sRSS, pst->idSelf);
	pstInitStat(&outr->sFreeMemory,pst->idSelf);
#endif
	pstInitStat(&outr->sLocal,pst->idSelf);
	pstInitStat(&outr->sActive,pst->idSelf);
	pstInitStat(&outr->sPart,pst->idSelf);
	pstInitStat(&outr->sPartNumAccess,pst->idSelf);
	pstInitStat(&outr->sPartMissRatio,pst->idSelf);
	pstInitStat(&outr->sCell,pst->idSelf);
	pstInitStat(&outr->sCellNumAccess,pst->idSelf);
	pstInitStat(&outr->sCellMissRatio,pst->idSelf);
	pstInitStat(&outr->sFlop,pst->idSelf);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	outr->sComputing.dSum     = mdlTimeComputing(pst->mdl);
	outr->sWaiting.dSum       = mdlTimeWaiting(pst->mdl);
	outr->sSynchronizing.dSum = mdlTimeSynchronizing(pst->mdl);
	pstInitStat(&outr->sComputing,pst->idSelf);
	pstInitStat(&outr->sWaiting,pst->idSelf);
	pstInitStat(&outr->sSynchronizing,pst->idSelf);
#endif
	/*
	** If there were no actives then we need to zero statistics counts
	** for those quantities which are based on per active particles.
	*/
	if (!outr->nActive) {
	    outr->sPart.n = 0;
	    outr->sPartNumAccess.n = 0;
	    outr->sPartMissRatio.n = 0;
	    outr->sCell.n = 0;
	    outr->sCellNumAccess.n = 0;
	    outr->sCellMissRatio.n = 0;	    
	    }
	out->dWalkTime = pkdGetWallClockTimer(plcl->pkd,1);
	}
    if (pnOut) *pnOut = (mdlThreads(pst->mdl) - pst->idSelf)*sizeof(struct outGravityPerProc) + sizeof(struct outGravityReduct);
    }


void pstLightCone(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inLightCone *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inLightCone));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONE,in,nIn);
	pstLightCone(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdLightCone(plcl->pkd,in->uRungLo,in->uRungHi,in->dLookbackFac,in->dLookbackFacLCP,in->dtLCDrift,in->dtLCKick);
	}
    if (pnOut) *pnOut = 0;
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
	pkdDrift(plcl->pkd,in->iRoot,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaUPred);
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
	pst->nTotal = pkdLocal(plcl->pkd);
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
	pkdInitStep(plcl->pkd,&in->param,&in->cosmo);
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
pstCountRungs(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outCountRungs *out = vout, outUpper;
    int i;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COUNTRUNGS,vin,nIn);
	pstCountRungs(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,pnOut);
	for(i=0; i<=MAX_RUNG; ++i) out->nRungs[i] += outUpper.nRungs[i];
	}
    else {
	pkdCountRungs(plcl->pkd, out->nRungs);
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


void pstNewFof(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inNewFof *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inNewFof));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_NEW_FOF,in,nIn);
	pstNewFof(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	pkdNewFof(plcl->pkd,in->dTau2,in->nMinMembers);
	}
    if (pnOut) *pnOut = 0;
    }

void pstFofPhases(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct outFofPhases *out = vout;
    int bMadeProgress;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FOF_PHASES,vin,nIn);
	pstFofPhases(pst->pstLower,vin,nIn,out,NULL);
	bMadeProgress = out->bMadeProgress;
	mdlGetReply(pst->mdl,rID,out,NULL);
	if (!out->bMadeProgress) out->bMadeProgress = bMadeProgress;
	}
    else {
	LCL *plcl = pst->plcl;
	out->bMadeProgress = pkdFofPhases(plcl->pkd);
	}
    if (pnOut) *pnOut = sizeof(struct outFofPhases);
    }


/*
** This is an almost identical copy of HopFinishUp.
*/
void pstFofFinishUp(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inFofFinishUp *in = (struct inFofFinishUp *)vin;
    uint64_t *nOutGroups = (uint64_t *)vout;
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inFofFinishUp));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FOF_FINISH_UP,vin,nIn);
        pstFofFinishUp(pst->pstLower,vin,nIn,vout,pnOut);
        mdlGetReply(pst->mdl,rID,&nOutUpper,pnOut);
	*nOutGroups += nOutUpper;
        }
    else {
	LCL *plcl = pst->plcl;
        *nOutGroups = pkdFofFinishUp(plcl->pkd,in->nMinGroupSize);
        }
    if (pnOut) *pnOut = sizeof(uint64_t);
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
    if (pnOut) *pnOut = sizeof(struct outGetFFTMaxSizes);
    }

void pltMoveIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inMoveIC *in = vin;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(struct inMoveIC));
    assert(pstOnNode(pst)); /* We pass around pointers! */
    if (pstNotCore(pst)) {
	struct inMoveIC icUp;
	
	icUp.pBase = in->pBase;
	icUp.nMove = pst->nUpper * in->nMove / pst->nLeaves;
	in->nMove -= icUp.nMove;
	icUp.iStart = in->iStart + in->nMove;
	icUp.fMass = in->fMass;
	icUp.fSoft = in->fSoft;
	icUp.nGrid = in->nGrid;
	icUp.nInflateFactor = in->nInflateFactor;

	int rID = mdlReqService(pst->mdl,pst->idUpper,PLT_MOVEIC,&icUp,nIn);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	pltMoveIC(pst->pstLower,in,nIn,vout,pnOut);
	}
    else {
	PKD pkd = plcl->pkd;
	assert(in->nInflateFactor>0);
	if (in->nMove <= (pkd->nStore/in->nInflateFactor)) {}
	else {
	    printf("nMove=%"PRIu64" nStore=%d nInflateFactor=%d\n", in->nMove, pkd->nStore, in->nInflateFactor);
	}
	assert(in->nMove <= (pkd->nStore/in->nInflateFactor));
	double inGrid = 1.0 / in->nGrid;
	for(i=in->nMove-1; i>=0; --i) {
	    expandParticle *b = ((expandParticle *)in->pBase) + in->iStart + i;
	    PARTICLE *p = pkdParticle(pkd,i);
	    vel_t *pVel = pkdVel(pkd,p);
	    expandParticle temp;
	    memcpy(&temp,b,sizeof(temp));
	    pVel[2] = temp.v[2];
	    pVel[1] = temp.v[1];
	    pVel[0] = temp.v[0];
	    pkdSetPos(pkd,p,2,temp.dr[2] + (temp.iz+0.5) * inGrid - 0.5);
	    pkdSetPos(pkd,p,1,temp.dr[1] + (temp.iy+0.5) * inGrid - 0.5);
	    pkdSetPos(pkd,p,0,temp.dr[0] + (temp.ix+0.5) * inGrid - 0.5);
	    pkdSetClass(pkd,in->fMass,in->fSoft,FIO_SPECIES_DARK,p);
	    p->bMarked = 1;
	    p->uRung = 0;
	    if (pkd->bNoParticleOrder) {
		((UPARTICLE *)p)->iGroup = 0;
		}
	    else {
		p->iOrder = temp.ix + in->nGrid*(temp.iy + 1ul*in->nGrid*temp.iz);
		p->uNewRung = 0;
		}
	    }
	pkd->nLocal = pkd->nActive = in->nMove;
	}
    if (pnOut) *pnOut = 0;
    }

void pstMoveIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inGenerateIC *in = vin;

    if (pstOffNode(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MOVEIC,vin,nIn);
	pstMoveIC(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	MDLFFT fft = pkd->fft;
	int myProc = mdlProc(pst->mdl);
	int iProc;
	uint64_t iUnderBeg=0, iOverBeg=0;
	uint64_t iUnderEnd, iOverEnd;
	uint64_t iBeg, iEnd;

	int *scount = malloc(sizeof(int)*mdlProcs(pst->mdl)); assert(scount!=NULL);
	int *sdisps = malloc(sizeof(int)*mdlProcs(pst->mdl)); assert(sdisps!=NULL);
	int *rcount = malloc(sizeof(int)*mdlProcs(pst->mdl)); assert(rcount!=NULL);
	int *rdisps = malloc(sizeof(int)*mdlProcs(pst->mdl)); assert(rdisps!=NULL);

	assert(pstAmNode(pst));
	assert(fft != NULL);

	assert(in->nInflateFactor>0);
	uint64_t nPerNode = (uint64_t)mdlCores(pst->mdl) * (pkd->nStore / in->nInflateFactor);
	uint64_t nLocal = (int64_t)fft->rgrid->rn[myProc] * in->nGrid*in->nGrid;

	/* Calculate how many slots are free (under) and how many need to be sent (over) before my rank */
	iUnderBeg = iOverBeg = 0;
	for(iProc=0; iProc<myProc; ++iProc) {
	    uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
	    if (nOnNode>nPerNode) iOverBeg += nOnNode - nPerNode;
	    else iUnderBeg += nPerNode - nOnNode;
	    }
	expandParticle *pBase = (expandParticle *)pkdParticleBase(pkd);
	expandParticle *pRecv = pBase + nLocal;
	expandParticle *eBase;
	if (nLocal > nPerNode) {      /* Too much here: send extra particles to other nodes */
	    eBase = pBase + nPerNode;
	    iUnderBeg = 0;
	    iOverEnd = iOverBeg + nLocal - nPerNode;
	    for(iProc=0; iProc<mdlProcs(pst->mdl); ++iProc) {
		rcount[iProc] = rdisps[iProc] = 0; // We cannot receive anything
		uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
		if (nOnNode<nPerNode) {
		    iUnderEnd = iUnderBeg + nPerNode - nOnNode;
		    /* The transfer condition */
		    if (iUnderEnd>iOverBeg && iUnderBeg<iOverEnd) {
			iBeg = iOverBeg>iUnderBeg ? iOverBeg : iUnderBeg;
			iEnd = iOverEnd<iUnderEnd ? iOverEnd : iUnderEnd;
			scount[iProc] = (iEnd-iBeg);
			sdisps[iProc] = (iBeg-iOverBeg);
			nLocal -= iEnd-iBeg;
			}
		    else scount[iProc] = sdisps[iProc] = 0;
		    iUnderBeg = iUnderEnd;
		    }
		else scount[iProc] = sdisps[iProc] = 0;
		}
	    assert(nLocal == nPerNode);
	    }
	else if (nLocal < nPerNode) { /* We have room: *maybe* receive particles from other nodes */
	    eBase = pBase + nLocal;
	    iOverBeg = 0;
	    iUnderEnd = iUnderBeg + nPerNode - nLocal;
	    for(iProc=0; iProc<mdlProcs(pst->mdl); ++iProc) {
		scount[iProc] = sdisps[iProc] = 0; // We have nothing to send
		uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
		if (nOnNode>nPerNode) {
		    iOverEnd = iOverBeg + nOnNode - nPerNode;
		    if (iOverEnd>iUnderBeg && iOverBeg<iUnderEnd) {
			iBeg = iOverBeg>iUnderBeg ? iOverBeg : iUnderBeg;
			iEnd = iOverEnd<iUnderEnd ? iOverEnd : iUnderEnd;
			rcount[iProc] = (iEnd-iBeg);
			rdisps[iProc] = (iBeg-iUnderBeg);
			nLocal += iEnd-iBeg;
			}
		    else rcount[iProc] = rdisps[iProc] = 0;
		    iOverBeg = iOverEnd;
		    }
		else rcount[iProc] = rdisps[iProc] = 0;
		}
	    assert(nLocal <= nPerNode);
	    }
	else {
	    for(iProc=0; iProc<mdlProcs(pst->mdl); ++iProc) {
		rcount[iProc] = rdisps[iProc] = 0; // We cannot receive anything
		scount[iProc] = sdisps[iProc] = 0; // We have nothing to send
		}
	    }
	mdlAlltoallv(pst->mdl, sizeof(expandParticle),
	    pBase + nPerNode, scount, sdisps,
	    pRecv,            rcount, rdisps);
	free(scount);
	free(sdisps);
	free(rcount);
	free(rdisps);
	mdlFFTNodeFinish(pst->mdl,fft);
	pkd->fft = NULL;

	/* We need to relocate the particles */
	struct inMoveIC move;
	uint64_t nTotal;
	nTotal = in->nGrid; /* Careful: 32 bit integer cubed => 64 bit integer */
	nTotal *= in->nGrid;
	nTotal *= in->nGrid;
	move.pBase = (overlayedParticle *)pkdParticleBase(pkd);
	move.iStart = 0;
	move.nMove = nLocal;
	move.fMass = in->dBoxMass;
	move.fSoft = 1.0 / (50.0*in->nGrid);
	move.nGrid = in->nGrid;
	move.nInflateFactor = in->nInflateFactor;
	pltMoveIC(pst,&move,sizeof(move),NULL,0);
	}
    }

/* NOTE: only called when on-node -- pointers are passed around. */
void pltGenerateIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inGenerateICthread *tin = vin;
    struct inGenerateIC *in = tin->ic;
    struct outGenerateIC *out = vout, outUp;
    mdlassert(pst->mdl,nIn == sizeof(struct inGenerateICthread));
    mdlassert(pst->mdl,vout != NULL);
    assert(pstOnNode(pst)); /* We pass around pointers! */

    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PLT_GENERATEIC,vin,nIn);
	pltGenerateIC(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUp,pnOut);
	out->N += outUp.N;
	out->noiseMean += outUp.noiseMean;
	out->noiseCSQ += outUp.noiseCSQ;
	}
    else {
	if (in->bClass)
	    out->N = pkdGenerateClassICm(plcl->pkd,tin->fft,in->iSeed, in->bFixed,in->fPhase,
	        in->nGrid, in->dBoxSize,&in->cosmo,in->dExpansion,&out->noiseMean,&out->noiseCSQ);
	else
	    out->N = pkdGenerateIC(plcl->pkd,tin->fft,in->iSeed,in->bFixed,in->fPhase,
	        in->nGrid,in->b2LPT,in->dBoxSize, &in->cosmo,in->dExpansion,in->nTf,
	        in->k, in->tf,&out->noiseMean,&out->noiseCSQ);
	out->dExpansion = in->dExpansion;
	}

    if (pnOut) *pnOut = sizeof(struct outGenerateIC);
    }

void pstGenerateIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inGenerateIC *in = vin;
    struct outGenerateIC *out = vout, outUp;
    int64_t i;

    mdlassert(pst->mdl,nIn == sizeof(struct inGenerateIC));
    mdlassert(pst->mdl,vout != NULL);

    if (pstAmNode(pst)) {
	struct inGenerateICthread tin;
	MDLFFT fft = mdlFFTNodeInitialize(pst->mdl,in->nGrid,in->nGrid,in->nGrid,0,0);
	tin.ic = vin;
	tin.fft = fft;
	pltGenerateIC(pst,&tin,sizeof(tin),vout,pnOut);

	int myProc = mdlProc(pst->mdl);
	uint64_t nLocal = (int64_t)fft->rgrid->rn[myProc] * in->nGrid*in->nGrid;

	/* Expand the particles by adding an iOrder */
	assert(sizeof(expandParticle) >= sizeof(basicParticle));
	overlayedParticle  * pbBase = (overlayedParticle *)pkdParticleBase(pkd);
	int iz = fft->rgrid->rs[myProc] + fft->rgrid->rn[myProc];
	int iy=0, ix=0;
	float inGrid = 1.0 / in->nGrid;
	for(i=nLocal-1; i>=0; --i) {
	    basicParticle  *b = &pbBase->b + i;
	    expandParticle *p = &pbBase->e + i;
	    basicParticle temp;
	    if (ix>0) --ix;
	    else {
		ix = in->nGrid-1;
		if (iy>0) --iy;
		else {
		    iy = in->nGrid-1;
		    --iz;
		    assert(iz>=0);
		    }
		}
	    memcpy(&temp,b,sizeof(temp));
	    p->v[2] = temp.v[2];
	    p->v[1] = temp.v[1];
	    p->v[0] = temp.v[0];
	    p->dr[2] = temp.dr[2];
	    p->dr[1] = temp.dr[1];
	    p->dr[0] = temp.dr[0];
	    p->ix = ix;
	    p->iy = iy;
	    p->iz = iz;
	    }
	assert(ix==0 && iy==0 && iz==fft->rgrid->rs[myProc]);
	/* Now we need to move excess particles between nodes so nStore is obeyed. */
	pkd->fft = fft; /* This is freed in pstMoveIC() */
	}
    else if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GENERATEIC,in,nIn);
	pstGenerateIC(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUp,pnOut);
	out->N += outUp.N;
	out->noiseMean += outUp.noiseMean;
	out->noiseCSQ += outUp.noiseCSQ;
	}
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
	char buffer[512], *save, *f, *v;
	int i;

#ifdef __linux__
	fp = fopen("/proc/self/stat","r");
	if ( fp != NULL ) {
	    if (fgets(buffer,sizeof(buffer),fp)==NULL) buffer[0] = '\0';
	    fclose(fp);
	    f = strtok_r(buffer," ",&save);
	    for ( i=0; i<= 36 && f; i++ ) {
		switch (i) {
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
	out->freeMemory = 0;
	fp = fopen("/proc/meminfo","r");
	if ( fp != NULL ) {
	    while(fgets(buffer,sizeof(buffer),fp)) {
		f = strtok_r(buffer,":",&save);
		v = strtok_r(NULL," ",&save);
		if (strcmp(f,"MemFree")==0 || strcmp(f,"Buffers")==0 || strcmp(f,"Cached")==0) {
		    out->freeMemory += atol(v) / 1024;
		    }
		}
	    fclose(fp);
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


void pstSelAll(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELALL,vin,nIn);
	pstSelAll(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelAll(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelGas(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELGAS,vin,nIn);
	pstSelGas(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelGas(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelStar(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSTAR,vin,nIn);
	pstSelStar(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelStar(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelDeleted(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDELETED,vin,nIn);
	pstSelDeleted(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelDeleted(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSelById(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelById *in = vin;
    struct outSelById *out = vout;
    struct outSelById outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelById) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELBYID,vin,nIn);
	pstSelById(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelById));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelById(plcl->pkd,in->idStart,in->idEnd,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelById);
    }
void pstSelMass(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelMass *in = vin;
    struct outSelMass *out = vout;
    struct outSelMass outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelMass) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELMASS,vin,nIn);
	pstSelMass(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelMass));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelMass(plcl->pkd,in->dMinMass,in->dMaxMass,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelMass);
    }

void pstSelPhaseDensity(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelPhaseDensity *in = vin;
    struct outSelPhaseDensity *out = vout;
    struct outSelPhaseDensity outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelPhaseDensity) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELPHASEDENSITY,vin,nIn);
	pstSelPhaseDensity(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelPhaseDensity));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelPhaseDensity(plcl->pkd,in->dMinDensity,in->dMaxDensity,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelPhaseDensity);
    }

void pstSelBox(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelBox *in = vin;
    struct outSelBox *out = vout;
    struct outSelBox outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelBox) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELBOX,vin,nIn);
	pstSelBox(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelBox));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelBox(
	    plcl->pkd,in->dCenter,in->dSize,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelBox);
    }

void pstSelSphere(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelSphere *in = vin;
    struct outSelSphere *out = vout;
    struct outSelSphere outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelSphere) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSPHERE,vin,nIn);
	pstSelSphere(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelSphere));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSphere(
	    plcl->pkd,in->r,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelSphere);
    }

void pstSelCylinder(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inSelCylinder *in = vin;
    struct outSelCylinder *out = vout;
    struct outSelCylinder outUpper;
    int nOut;

    assert( nIn==sizeof(struct inSelCylinder) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELCYLINDER,vin,nIn);
	pstSelCylinder(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelCylinder));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelCylinder(
	    plcl->pkd,in->dP1,in->dP2,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    if (pnOut) *pnOut = sizeof(struct outSelCylinder);
    }

void pstSelGroup(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    int *in = vin;
    LCL *plcl = pst->plcl;
    assert( nIn==sizeof(int) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELGROUP,vin,nIn);
	pstSelGroup(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdSelGroup(plcl->pkd, *in);
	}
    if (pnOut) *pnOut = 0;
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
    struct outMeasurePk *outUpper;
    int nOut;
    int i;

    assert( nIn==sizeof(struct inMeasurePk) );
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEASUREPK,vin,nIn);
	pstMeasurePk(pst->pstLower,vin,nIn,vout,pnOut);
	outUpper = malloc(sizeof(struct outMeasurePk));
	assert(outUpper != NULL);
	mdlGetReply(pst->mdl,rID,outUpper,&nOut);
	assert(nOut==sizeof(struct outMeasurePk));

	for(i=0;i<in->nBins; i++) {
	    out->fK[i] += outUpper->fK[i];
	    out->fPower[i] += outUpper->fPower[i];
	    out->nPower[i] += outUpper->nPower[i];
	    }
	free(outUpper);
	}
    else {
	pkdMeasurePk(plcl->pkd, in->dTotalMass, in->iAssignment, in->bInterlace,
	    in->nGrid, in->nBins, out->fK, out->fPower, out->nPower);
	}
    if (pnOut) *pnOut = sizeof(struct outMeasurePk);
    }

void pstMeasureLinPk(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inMeasureLinPk *in = vin;
    struct outMeasureLinPk *out = vout;
    struct outMeasureLinPk *outUpper;
    int nOut;
    int i;

    assert( nIn==sizeof(struct inMeasureLinPk) );
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEASURELINPK,vin,nIn);
	pstMeasureLinPk(pst->pstLower,vin,nIn,vout,pnOut);
	outUpper = malloc(sizeof(struct outMeasureLinPk));
	assert(outUpper != NULL);
	mdlGetReply(pst->mdl,rID,outUpper,&nOut);
	assert(nOut==sizeof(struct outMeasureLinPk));

	for(i=0;i<in->nBins; i++) {
	    out->fK[i] += outUpper->fK[i];
	    out->fPower[i] += outUpper->fPower[i];
	    out->nPower[i] += outUpper->nPower[i];
	    }
	free(outUpper);
	}
    else {
	pkdMeasureLinPk(plcl->pkd, in->nGrid, in->dA, in->dBoxSize,
                        in->nBins, in->iSeed, in->bFixed, in->fPhase, 
                        out->fK, out->fPower, out->nPower);
	}
    if (pnOut) *pnOut = sizeof(struct outMeasureLinPk);
    }

void pstSetLinGrid(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
        LCL *plcl = pst->plcl;
        struct inSetLinGrid *in = vin;
        assert (nIn==sizeof(struct inSetLinGrid) );
        if (pstNotCore(pst)) {
            int rID = mdlReqService(pst->mdl, pst->idUpper, PST_SETLINGRID, vin, nIn);
            pstSetLinGrid(pst->pstLower, vin, nIn, vout, pnOut);
            mdlGetReply(pst->mdl,rID, vout,pnOut);
        }
        else {
            plcl->pkd->Linfft = mdlFFTInitialize(pst->mdl, in->nGrid, in->nGrid, in->nGrid, 0,0);
            pkdSetLinGrid(plcl->pkd, in->dTime,
                in->dBSize, in->nGrid, 
                in ->iSeed, in->bFixed, in->fPhase);
        }
    }

void pstGridCreateFFT(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
        LCL *plcl = pst->plcl;
        struct inGridCreateFFT *in = vin;
        assert (nIn==sizeof(struct inGridCreateFFT) );
        if (pstNotCore(pst)) {
            int rID = mdlReqService(pst->mdl, pst->idUpper, PST_GRID_CREATE_FFT, vin, nIn);
            pstGridCreateFFT(pst->pstLower, vin, nIn, vout, pnOut);
            mdlGetReply(pst->mdl,rID, vout,pnOut);
        }
        else {
            PKD pkd = plcl->pkd;
            assert(pkd->fft == NULL);
            pkd->fft = mdlFFTInitialize(pst->mdl, in->nGrid, in->nGrid, in->nGrid, 0, 0);
        }
    }

void pstGridDeleteFFT(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
        LCL *plcl = pst->plcl;
        if (pstNotCore(pst)) {
            int rID = mdlReqService(pst->mdl, pst->idUpper, PST_GRID_DELETE_FFT, vin, nIn);
            pstGridDeleteFFT(pst->pstLower, vin, nIn, vout, pnOut);
            mdlGetReply(pst->mdl,rID, vout,pnOut);
        }
        else {
            PKD pkd = plcl->pkd;
            assert(pkd->fft != NULL);
            mdlFFTFinish(pst->mdl,plcl->pkd->fft);
            pkd->fft = NULL;
        }
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

void pstLightConeOpen(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inLightConeOpen *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inLightConeOpen));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONE_OPEN,in,nIn);
        pstLightConeOpen(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	PKD pkd = pst->plcl->pkd;
	char achOutFile[PST_FILENAME_SIZE];
	if (in->achOutFile[0]) makeName(achOutFile,in->achOutFile,mdlSelf(pkd->mdl),"lcp.");
	else achOutFile[0] = 0;
        pkdLightConeOpen(pkd, achOutFile, in->nSideHealpix);
        }
    if (pnOut) *pnOut = 0;
}

void pstLightConeClose(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    struct inLightConeClose *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inLightConeClose));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONE_CLOSE,in,nIn);
        pstLightConeClose(pst->pstLower,in,nIn,NULL,NULL);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	PKD pkd = pst->plcl->pkd;
	char achOutFile[PST_FILENAME_SIZE];
	if (in->achOutFile[0]) makeName(achOutFile,in->achOutFile,mdlSelf(pkd->mdl),"hpb.");
	else achOutFile[0] = 0;
        pkdLightConeClose(pst->plcl->pkd,achOutFile);
        }
    if (pnOut) *pnOut = 0;
}

void pstLightConeVel(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inScaleVel *in = vin;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONEVEL,in,nIn);
	pstLightConeVel(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdLightConeVel(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstInflate(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INFLATE,vin,nIn);
	pstInflate(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	struct inInflate *in = vin;
	pkdInflate(plcl->pkd,in->nInflateReps);
	}
    if (pnOut) *pnOut = 0;
    }

void pstGetParticles(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct outGetParticles *out = vout;
    uint64_t *ID = vin;
    int nOutUpper;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GET_PARTICLES,vin,nIn);
	pstGetParticles(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,out + (*pnOut / sizeof(struct outGetParticles)),&nOutUpper);
	*pnOut += nOutUpper;
	}
    else {
	int nParticles = nIn / sizeof(uint64_t);
	int n = pkdGetParticles(plcl->pkd,nParticles, ID, out );
	*pnOut = n * sizeof(struct outGetParticles);
	}
    }
