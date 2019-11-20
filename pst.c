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

    mdlAddService(mdl,PST_SETADD,pst,(fcnService_t*)pstSetAdd,
		  sizeof(struct inSetAdd),0);
    mdlAddService(mdl,PST_READFILE,pst,(fcnService_t*)pstReadFile,
	          sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE),0);
    mdlAddService(mdl,PST_DOMAINDECOMP,pst,(fcnService_t*)pstDomainDecomp,
		  sizeof(struct inDomainDecomp),0);
    mdlAddService(mdl,PST_CALCBOUND,pst,(fcnService_t*)pstCalcBound,
		  0,sizeof(BND));
    mdlAddService(mdl,PST_CALCVBOUND,pst,(fcnService_t*)pstCalcVBound,
		  0,sizeof(BND));
    mdlAddService(mdl,PST_COMBINEBOUND,pst,(fcnService_t*)pstCombineBound,
		  0,sizeof(BND));
    mdlAddService(mdl,PST_WEIGHT,pst,(fcnService_t*)pstWeight,
		  sizeof(struct inWeight),sizeof(struct outWeight));
    mdlAddService(mdl,PST_COUNTVA,pst,(fcnService_t*)pstCountVA,
		  sizeof(struct inCountVA),sizeof(struct outCountVA));
    mdlAddService(mdl,PST_WEIGHTWRAP,pst,(fcnService_t*)pstWeightWrap,
		  sizeof(struct inWeightWrap),sizeof(struct outWeightWrap));
    mdlAddService(mdl,PST_ORDWEIGHT,pst,(fcnService_t*)pstOrdWeight,
		  sizeof(struct inOrdWeight),sizeof(struct outOrdWeight));
    mdlAddService(mdl,PST_FREESTORE,pst,(fcnService_t*)pstFreeStore,
		  0,sizeof(struct outFreeStore));
    mdlAddService(mdl,PST_COLREJECTS,pst,(fcnService_t*)pstColRejects,
		  0,nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_SWAPREJECTS,pst,(fcnService_t*)pstSwapRejects,
		  nThreads*sizeof(int),nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_COLORDREJECTS,pst,(fcnService_t*)pstColOrdRejects,
		  sizeof(struct inColOrdRejects),nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_DOMAINORDER,pst,(fcnService_t*)pstDomainOrder,
		  sizeof(struct inDomainOrder),0);
    mdlAddService(mdl,PST_LOCALORDER,pst,(fcnService_t*)pstLocalOrder,
		  sizeof(struct inDomainOrder),0);
    mdlAddService(mdl,PST_COMPRESSASCII,pst,(fcnService_t*)pstCompressASCII,
		  sizeof(struct inCompressASCII),sizeof(struct outCompressASCII));
    mdlAddService(mdl,PST_WRITEASCII,pst,(fcnService_t*)pstWriteASCII,
		  sizeof(struct inWriteASCII),0);
    mdlAddService(mdl,PST_WRITE,pst,(fcnService_t*)pstWrite,
		  sizeof(struct inWrite),0);
    mdlAddService(mdl,PST_SENDPARTICLES,pst,(fcnService_t*)pstSendParticles,
		  sizeof(int),0);
    mdlAddService(mdl,PST_CHECKPOINT,pst,(fcnService_t*)pstCheckpoint,
		  sizeof(struct inWrite),0);
    mdlAddService(mdl,PST_OUTPUT,pst,(fcnService_t*)pstOutput,
		  sizeof(struct inOutput),0);
    mdlAddService(mdl,PST_OUTPUT_SEND,pst,(fcnService_t*)pstOutputSend,
		  sizeof(struct inOutputSend),0);
    mdlAddService(mdl,PST_RESTORE,pst,(fcnService_t*)pstRestore,
		  sizeof(struct inRestore),0);
    /*
    ** Calculate the number of levels in the top tree and use it to
    ** define the size of the messages.
    */
    mdlAddService(mdl,PST_BUILDTREE,pst,(fcnService_t*)pstBuildTree,
	sizeof(struct inBuildTree),
	(nThreads==1?1:2*nThreads-1)*pkdMaxNodeSize());
    mdlAddService(mdl,PST_DISTRIBTOPTREE,pst,(fcnService_t*)pstDistribTopTree,
	sizeof(struct inDistribTopTree) + (nThreads==1?1:2*nThreads-1)*pkdMaxNodeSize(),0);
    mdlAddService(mdl,PST_DUMPTREES,pst,(fcnService_t*)pstDumpTrees,
	          sizeof(struct inDumpTrees),0);
    mdlAddService(mdl,PST_TREEINITMARKED,pst,(fcnService_t*)pstTreeInitMarked,
	          0,0);
    mdlAddService(mdl,PST_CALCROOT,pst,(fcnService_t*)pstCalcRoot,
	          sizeof(struct inCalcRoot),sizeof(struct outCalcRoot));
    mdlAddService(mdl,PST_DISTRIBROOT,pst,(fcnService_t*)pstDistribRoot,
		  sizeof(struct ioDistribRoot),0);
    mdlAddService(mdl,PST_ENFORCEPERIODIC,pst,(fcnService_t*)pstEnforcePeriodic,
		  sizeof(BND),0);
    mdlAddService(mdl,PST_HOP_LINK,pst,(fcnService_t*)pstHopLink,
	          sizeof(struct inHopLink),sizeof(uint64_t));
    mdlAddService(mdl,PST_HOP_JOIN,pst,(fcnService_t*)pstHopJoin,
		  sizeof(struct inHopLink),sizeof(struct outHopJoin));
    mdlAddService(mdl,PST_HOP_FINISH_UP,pst,(fcnService_t*)pstHopFinishUp,
	          sizeof(struct inHopFinishUp),sizeof(uint64_t));
    mdlAddService(mdl,PST_HOP_TREE_BUILD,pst,(fcnService_t*)pstHopTreeBuild,
	          sizeof(struct inHopTreeBuild),0);
    mdlAddService(mdl,PST_HOP_GRAVITY,pst,(fcnService_t*)pstHopGravity,
	          sizeof(struct inHopGravity),0);
    mdlAddService(mdl,PST_HOP_UNBIND,pst,(fcnService_t*)pstHopUnbind,
	          sizeof(struct inHopUnbind),sizeof(struct outHopUnbind));
    mdlAddService(mdl,PST_GROUP_RELOCATE,pst,(fcnService_t*)pstGroupRelocate,
		  0,0);
    mdlAddService(mdl,PST_GROUP_STATS,pst,(fcnService_t*)pstGroupStats,
	          sizeof(struct inGroupStats),0);
    mdlAddService(mdl,PST_SMOOTH,pst,(fcnService_t*)pstSmooth,
		  sizeof(struct inSmooth),0);
#ifdef FAST_GAS
    mdlAddService(mdl,PST_FASTGASPHASE1,pst,(fcnService_t*)pstFastGasPhase1,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_FASTGASPHASE2,pst,(fcnService_t*)pstFastGasPhase2,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_FASTGASCLEANUP,pst,(fcnService_t*)pstFastGasCleanup,
		  0,0);
#endif
    mdlAddService(mdl,PST_GRAVITY,pst,(fcnService_t*)pstGravity,
	          sizeof(struct inGravity),
	          nThreads*sizeof(struct outGravityPerProc) + sizeof(struct outGravityReduct));
    mdlAddService(mdl,PST_CALCEANDL,pst,(fcnService_t*)pstCalcEandL,
		  0,sizeof(struct outCalcEandL));
    mdlAddService(mdl,PST_DRIFT,pst,(fcnService_t*)pstDrift,
		  sizeof(struct inDrift),0);
    mdlAddService(mdl,PST_SCALEVEL,pst,(fcnService_t*)pstScaleVel,
		  sizeof(struct inScaleVel),0);
    mdlAddService(mdl,PST_CACHEBARRIER,pst,(fcnService_t*)pstCacheBarrier,
		  0,0);
    mdlAddService(mdl,PST_STEPVERYACTIVE,pst,(fcnService_t*)pstStepVeryActiveKDK,
		  sizeof(struct inStepVeryActive),
		  sizeof(struct outStepVeryActive));
    mdlAddService(mdl,PST_ROPARTICLECACHE,pst,(fcnService_t*)pstROParticleCache,
		  0,0);
    mdlAddService(mdl,PST_PARTICLECACHEFINISH,pst,(fcnService_t*)pstParticleCacheFinish,
		  0,0);
    mdlAddService(mdl,PST_KICK,pst,(fcnService_t*)pstKick,
		  sizeof(struct inKick),sizeof(struct outKick));
    mdlAddService(mdl,PST_KICKTREE,pst,(fcnService_t*)pstKickTree,
		  sizeof(struct inKickTree),sizeof(struct outKickTree));
    mdlAddService(mdl,PST_SETSOFT,pst,(fcnService_t*)pstSetSoft,
		  sizeof(struct inSetSoft),0);
    mdlAddService(mdl,PST_PHYSICALSOFT,pst,(fcnService_t*)pstPhysicalSoft,
		  sizeof(struct inPhysicalSoft),0);
    mdlAddService(mdl,PST_SETTOTAL,pst,(fcnService_t*)pstSetTotal,
		  0,sizeof(struct outSetTotal));
    mdlAddService(mdl,PST_SETWRITESTART,pst,(fcnService_t*)pstSetWriteStart,
		  sizeof(struct inSetWriteStart),0);
    mdlAddService(mdl,PST_ADDWRITESTART,pst,(fcnService_t*)pstAddWriteStart,
		  sizeof(struct inAddWriteStart),0);
    mdlAddService(mdl,PST_ONENODEREADINIT,pst,(fcnService_t*)pstOneNodeReadInit,
	          sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE),
	          nThreads*sizeof(int));
    mdlAddService(mdl,PST_SWAPALL,pst,(fcnService_t*)pstSwapAll,
		  sizeof(int),0);
    mdlAddService(mdl,PST_ACTIVEORDER,pst,(fcnService_t*)pstActiveOrder,
		  0,sizeof(uint64_t));
    mdlAddService(mdl,PST_INITCOSMOLOGY,pst,(fcnService_t*)pstInitCosmology,
		  sizeof(struct csmVariables),0);
    mdlAddService(mdl,PST_SETPARAMETERS,pst,(fcnService_t*)pstSetParameters,
		  sizeof(struct parameters),0);
    mdlAddService(mdl,PST_SETRUNG,pst,(fcnService_t*)pstSetRung,
		  sizeof(struct inSetRung),0);
    mdlAddService(mdl,PST_ZERONEWRUNG,pst,(fcnService_t*)pstZeroNewRung,
		  sizeof(struct inZeroNewRung),0);
    mdlAddService(mdl,PST_ACTIVERUNG,pst,(fcnService_t*)pstActiveRung,
		  sizeof(struct inActiveRung),0);
    mdlAddService(mdl,PST_COUNTRUNGS,pst,(fcnService_t*)pstCountRungs,
		  0,sizeof(struct outCountRungs));
    mdlAddService(mdl,PST_DENSITYSTEP,pst,(fcnService_t*)pstDensityStep,
		  sizeof(struct inDensityStep),0);
    mdlAddService(mdl,PST_CORRECTENERGY,pst,(fcnService_t*)pstCorrectEnergy,
		  sizeof(struct inCorrectEnergy),0);
    mdlAddService(mdl,PST_ACCELSTEP,pst,(fcnService_t*)pstAccelStep,
		  sizeof(struct inAccelStep), 0);
    mdlAddService(mdl,PST_SPHSTEP,pst,(fcnService_t*)pstSphStep,
		  sizeof(struct inSphStep), 0);
    mdlAddService(mdl,PST_STARFORM,pst,(fcnService_t*)pstStarForm,
		  sizeof(struct inStarForm),sizeof(struct outStarForm));
    mdlAddService(mdl,PST_SETRUNGVERYACTIVE,pst,(fcnService_t*)pstSetRungVeryActive,
		  sizeof(struct inSetRung),0);
    mdlAddService(mdl,PST_RESMOOTH,pst,(fcnService_t*)pstReSmooth,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_UPDATERUNG,pst,(fcnService_t*)pstUpdateRung,
		  sizeof(struct inUpdateRung),sizeof(struct outUpdateRung));
    mdlAddService(mdl,PST_UPDATERUNGBYTREE,pst,(fcnService_t*)pstUpdateRungByTree,
		  sizeof(struct inUpdateRungByTree),sizeof(struct outUpdateRung));
    mdlAddService(mdl,PST_COLNPARTS,pst,(fcnService_t*)pstColNParts,
		  0,nThreads*sizeof(struct outColNParts));
    mdlAddService(mdl,PST_NEWORDER,pst,(fcnService_t*)pstNewOrder,
		  nThreads*sizeof(uint64_t),0);
    mdlAddService(mdl,PST_GETNPARTS,pst,(fcnService_t*)pstGetNParts,
		  0,sizeof(struct outGetNParts));
    mdlAddService(mdl,PST_SETNPARTS,pst,(fcnService_t*)pstSetNParts,
		  sizeof(struct inSetNParts),0);
    mdlAddService(mdl,PST_CLEARTIMER,pst,(fcnService_t*)pstClearTimer,
		  sizeof(struct inClearTimer),0);
    mdlAddService(mdl,PST_NEW_FOF,pst,(fcnService_t*)pstNewFof,
		  sizeof(struct inNewFof),0);
    mdlAddService(mdl,PST_FOF_PHASES,pst,(fcnService_t*)pstFofPhases,
	          0,sizeof(struct outFofPhases));
    mdlAddService(mdl,PST_FOF_FINISH_UP,pst,(fcnService_t*)pstFofFinishUp,
	          sizeof(struct inFofFinishUp),sizeof(uint64_t));
    mdlAddService(mdl,PST_INITRELAXATION,pst,(fcnService_t*)pstInitRelaxation,0,0);
    mdlAddService(mdl,PST_INITIALIZEPSTORE,pst,(fcnService_t*)pstInitializePStore,
		  sizeof(struct inInitializePStore),0);
#ifdef MDL_FFTW
    mdlAddService(mdl,PST_GETFFTMAXSIZES,pst,(fcnService_t*)pstGetFFTMaxSizes,
		  sizeof(struct inGetFFTMaxSizes),sizeof(struct outGetFFTMaxSizes));
    mdlAddService(mdl,PST_GENERATEIC,pst,(fcnService_t*)pstGenerateIC,
		  sizeof(struct inGenerateIC),sizeof(struct outGenerateIC));
    mdlAddService(mdl,PLT_GENERATEIC,pst,(fcnService_t*)pltGenerateIC,
		  sizeof(struct inGenerateICthread),sizeof(struct outGenerateIC));
    mdlAddService(mdl,PLT_MOVEIC,pst,(fcnService_t*)pltMoveIC,
		  sizeof(struct inMoveIC),0);
    mdlAddService(mdl,PST_MOVEIC,pst,(fcnService_t*)pstMoveIC,
		  sizeof(struct inGenerateIC),0);
#endif
    mdlAddService(mdl,PST_HOSTNAME,pst,(fcnService_t*)pstHostname,
		  0,nThreads*sizeof(struct outHostname));
    mdlAddService(mdl,PST_MEMSTATUS,pst,(fcnService_t*)pstMemStatus,
		  0,nThreads*sizeof(struct outMemStatus));
    mdlAddService(mdl,PST_GETCLASSES,pst,(fcnService_t*)pstGetClasses,
		  0, PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlAddService(mdl,PST_SETCLASSES,pst,(fcnService_t*)pstSetClasses,
		  PKD_MAX_CLASSES*sizeof(PARTCLASS), 0);
    mdlAddService(mdl,PST_SWAPCLASSES,pst,(fcnService_t*)pstSwapClasses,
		  PKD_MAX_CLASSES*sizeof(PARTCLASS),
		  PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlAddService(mdl,PST_SELALL,pst,(fcnService_t*)pstSelAll,
		  0, 0 );
    mdlAddService(mdl,PST_SELGAS,pst,(fcnService_t*)pstSelGas,
		  0, 0 );
    mdlAddService(mdl,PST_SELSTAR,pst,(fcnService_t*)pstSelStar,
		  0, 0 );
    mdlAddService(mdl,PST_SELDELETED,pst,(fcnService_t*)pstSelDeleted,
		  0, 0 );
    mdlAddService(mdl,PST_SELBYID,pst,(fcnService_t*)pstSelById,
		  sizeof(struct inSelById), sizeof(struct outSelById));
    mdlAddService(mdl,PST_SELMASS,pst,(fcnService_t*)pstSelMass,
		  sizeof(struct inSelMass), sizeof(struct outSelMass));
    mdlAddService(mdl,PST_SELPHASEDENSITY,pst,(fcnService_t*)pstSelPhaseDensity,
		  sizeof(struct inSelPhaseDensity), sizeof(struct outSelPhaseDensity));
    mdlAddService(mdl,PST_SELBOX,pst,(fcnService_t*)pstSelBox,
		  sizeof(struct inSelBox), sizeof(struct outSelBox));
    mdlAddService(mdl,PST_SELSPHERE,pst,(fcnService_t*)pstSelSphere,
		  sizeof(struct inSelSphere), sizeof(struct outSelSphere));
    mdlAddService(mdl,PST_SELCYLINDER,pst,(fcnService_t*)pstSelCylinder,
		  sizeof(struct inSelCylinder), sizeof(struct outSelCylinder));
    mdlAddService(mdl,PST_SELGROUP,pst,(fcnService_t*)pstSelGroup,
		  sizeof(int), 0);
    mdlAddService(mdl,PST_PROFILE,pst,(fcnService_t*)pstProfile,
		  sizeof(struct inProfile), 0); 
    mdlAddService(mdl,PST_CALCDISTANCE,pst,(fcnService_t*)pstCalcDistance,
		  sizeof(struct inCalcDistance), 0);
    mdlAddService(mdl,PST_CALCCOM,pst,(fcnService_t*)pstCalcCOM,
		  sizeof(struct inCalcCOM), sizeof(struct outCalcCOM));
    mdlAddService(mdl,PST_COUNTDISTANCE,pst,(fcnService_t*)pstCountDistance,
		  sizeof(struct inCountDistance), sizeof(struct outCountDistance));
    mdlAddService(mdl,PST_INITGRID,pst,(fcnService_t*)pstInitGrid,
		  sizeof(struct inInitGrid), 0);
    mdlAddService(mdl,PST_GRIDPROJECT,pst,(fcnService_t*)pstGridProject,
		  sizeof(struct inGridProject), 0);
#ifdef MDL_FFTW
    mdlAddService(mdl,PST_GRID_CREATE_FFT,pst,(fcnService_t*)pstGridCreateFFT,
		  sizeof(struct inGridCreateFFT), 0);
    mdlAddService(mdl,PST_GRID_DELETE_FFT,pst,(fcnService_t*)pstGridDeleteFFT,
		  0, 0);
    mdlAddService(mdl,PST_MEASUREPK,pst,(fcnService_t*)pstMeasurePk,
		  sizeof(struct inMeasurePk), sizeof(struct outMeasurePk));
    mdlAddService(mdl,PST_ASSIGN_MASS,pst,(fcnService_t*)pstAssignMass,
		  sizeof(struct inAssignMass), 0);
    mdlAddService(mdl,PST_LINEARKICK, pst,(fcnService_t*)pstLinearKick,
           sizeof(struct inLinearKick), 0);
    mdlAddService(mdl,PST_SETLINGRID, pst,(fcnService_t*)pstSetLinGrid,
           sizeof(struct inSetLinGrid), 0);
    mdlAddService(mdl,PST_MEASURELINPK,pst,(fcnService_t*)pstMeasureLinPk,
		  sizeof(struct inMeasureLinPk), sizeof(struct outMeasureLinPk));
#endif
    mdlAddService(mdl,PST_TOTALMASS,pst,(fcnService_t*)pstTotalMass,
		  0, sizeof(struct outTotalMass));
    mdlAddService(mdl,PST_LIGHTCONE_OPEN,pst,(fcnService_t*)pstLightConeOpen,
		  sizeof(struct inLightConeOpen), 0);
    mdlAddService(mdl,PST_LIGHTCONE_CLOSE,pst,(fcnService_t*)pstLightConeClose,
		  sizeof(struct inLightConeClose), 0);
    mdlAddService(mdl,PST_LIGHTCONEVEL,pst,(fcnService_t*)pstLightConeVel,
		  0,0);
    mdlAddService(mdl,PST_INFLATE,pst,(fcnService_t*)pstInflate,
	          sizeof(struct inInflate), 0);
    mdlAddService(mdl,PST_GET_PARTICLES,pst,(fcnService_t*)pstGetParticles,
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

int pstSetAdd(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstSetAdd(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    return 0;
    }

static void initializePStore(PKD *ppkd,MDL mdl,struct inInitializePStore *in) {
    pkdInitialize(
	ppkd,mdl,in->nStore,in->nMinTotalStore,in->nMinEphemeral,in->nEphemeralBytes,
	in->nTreeBitsLo,in->nTreeBitsHi,
	in->iCacheSize,in->iWorkQueueSize,in->iCUDAQueueSize,in->fPeriod,
	in->nSpecies[FIO_SPECIES_DARK],in->nSpecies[FIO_SPECIES_SPH],in->nSpecies[FIO_SPECIES_STAR],
	in->mMemoryModel,in->bLightCone,in->bLightConeParticles);
    }

int pstInitializePStore(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inInitializePStore *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inInitializePStore));
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITIALIZEPSTORE,in,nIn);
	pstInitializePStore(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	if (plcl->pkd) pkdFinish(plcl->pkd);
	initializePStore(&plcl->pkd,pst->mdl,in);
	}
    return 0;
    }

int pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstOneNodeReadInit(pst->pstLower,in,nIn,pout,nOut);
	in->nNodeEnd = nFileEnd;
	mdlGetReply(pst->mdl,rID,pout+pst->nLower,NULL);
	}
    else {
	/*
	** Determine the size of the local particle store.
	*/
	*pout = nFileTotal; /* Truncated: okay */
	}
    return sizeof(*pout) * pst->nLeaves;
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

int pstReadFile(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstReadFile(pst->pstLower,in,nIn,NULL,0);

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
	int i;
	uint64_t nStart;
	PKD pkd;
	MDL mdl;
	PST pst0;

	assert(nParts!=NULL);
	pstOneNodeReadInit(pst,in,sizeof(*in),nParts,pst->nLeaves * sizeof(nParts));
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
    return 0;
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
    pstFreeStore(pst->pstLower,NULL,0,&outFree,sizeof(outFree));
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
	pstActiveOrder(pst, NULL, 0, &nActiveOrder, sizeof(nActiveOrder)); /* SOON NO MORE ACTIVE ORDER */

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
	pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
	mdlGetReply(pst->mdl,rID,&outWtHigh,NULL);
	nTotalActive = outWtLow.nLow + outWtHigh.nLow
		       + outWtLow.nHigh + outWtHigh.nHigh;
	mdlassert(pst->mdl,nActiveOrder == nTotalActive);
	pFlag = 1;
	if (nTotalActive <=1) {
	    pFlag = 0;			/* Divide them all */
	    rID = mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
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
	    pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
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
	pstCountVA(pst,&inCtVA,sizeof(inCtVA),&outCtVA,sizeof(outCtVA));
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
    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
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
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
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
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
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
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
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
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
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
    nOut = pstColRejects(pst->pstLower,NULL,0,pLowerRej,mdlThreads(pst->mdl)*sizeof(OREJ));
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
	nOut = pstSwapRejects(pst->pstLower,pidSwap,
		       mdlThreads(pst->mdl)*sizeof(int),pLowerRej,mdlThreads(pst->mdl)*sizeof(OREJ));
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

int pstDomainDecomp(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstDomainDecomp(pst->pstLower,vin,sizeof(*in),NULL,0);

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
    return 0;
    }


int pstCalcBound(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    BND *out = vout;
    BND outBnd;

    mdlassert(pst->mdl,nIn == 0);
    mdlassert(pst->mdl,nOut >= sizeof(BND));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCBOUND,NULL,0);
	pstCalcBound(pst->pstLower,NULL,0,out,nOut);
	mdlGetReply(pst->mdl,rID,&outBnd,NULL);
	BND_COMBINE(out,out,&outBnd);
	}
    else {
	pkdCalcBound(plcl->pkd,out);
	}
    return sizeof(BND);
    }

int pstCalcVBound(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    BND *out = vout;
    BND outBnd;

    mdlassert(pst->mdl,nIn == 0);
    mdlassert(pst->mdl,nOut >= sizeof(BND));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCVBOUND,NULL,0);
	pstCalcVBound(pst->pstLower,NULL,0,out,nOut);
	mdlGetReply(pst->mdl,rID,&outBnd,NULL);
	BND_COMBINE(out,out,&outBnd);
	}
    else {
	pkdCalcVBound(plcl->pkd,out);
	}
    return sizeof(BND);
    }

int pstCombineBound(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    BND *out = vout;
    BND outBnd;

    mdlassert(pst->mdl,nIn == 0);
    mdlassert(pst->mdl,nOut >= sizeof(BND));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COMBINEBOUND,NULL,0);
	pstCombineBound(pst->pstLower,NULL,0,out,nOut);
	mdlGetReply(pst->mdl,rID,&outBnd,NULL);
	BND_COMBINE(out,out,&outBnd);
	}
    else {
        *out = plcl->pkd->bnd;
	}
    return sizeof(BND);
    }


/*
** Make sure that the local particles are split into active and inactive
** when passing pFlag != 0.
** pFlag == 0 => weight all particles.
** pFlag > 0 => weight active particles.
** pFlag < 0 => weight inactive particles.
*/
int pstWeight(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstWeight(pst->pstLower,in,nIn,out,nOut);
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
    return sizeof(struct outWeight);
    /*
      pkdStopTimer(plcl->pkd,7);
    */
    }



int pstCountVA(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inCountVA *in = vin;
    struct outCountVA *out = vout;
    struct outCountVA outCt;

    mdlassert(pst->mdl,nIn == sizeof(struct inCountVA));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COUNTVA,in,nIn);
	pstCountVA(pst->pstLower,in,nIn,out,nOut);
	mdlGetReply(pst->mdl,rID,&outCt,NULL);
	out->nLow += outCt.nLow;
	out->nHigh += outCt.nHigh;
	}
    else {
	pkdCountVA(plcl->pkd,in->iSplitDim,in->fSplit,&out->nLow,&out->nHigh);
	}
    return sizeof(struct outCountVA);
    }



/*
** Make sure that the local particles are split into active and inactive
** when passing pFlag != 0.
** pFlag == 0 => weight all particles.
** pFlag > 0 => weight active particles.
** pFlag < 0 => weight inactive particles.
*/
int pstWeightWrap(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstWeightWrap(pst->pstLower,in,nIn,out,nOut);
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
    return sizeof(struct outWeightWrap);
    /*
      pkdStopTimer(plcl->pkd,7);
    */
    }


/*
** Weight request for splitting into iOrder order.
*/
int pstOrdWeight(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inOrdWeight *in = vin;
    struct outOrdWeight *out = vout;
    struct outOrdWeight outWt;
    int nLow,nHigh;

    mdlassert(pst->mdl,nIn == sizeof(struct inOrdWeight));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ORDWEIGHT,in,nIn);
	pstOrdWeight(pst->pstLower,in,nIn,out,nOut);
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
    return sizeof(struct outOrdWeight);
    }


int pstFreeStore(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outFreeStore *out = vout;
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
	pstFreeStore(pst->pstLower,NULL,0,out,nOut);
	pst->nLowerStore = out->nFreeStore;
	mdlGetReply(pst->mdl,rID,out,NULL);
	pst->nUpperStore = out->nFreeStore;
	out->nFreeStore = pst->nLowerStore + pst->nUpperStore;
	}
    else {
	out->nFreeStore = pkdFreeStore(plcl->pkd);
	}
    return sizeof(struct outFreeStore);
    }


int pstColRejects(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper;

    mdlassert(pst->mdl,nIn == 0);
    /*
      pkdStartTimer(plcl->pkd,9);
    */
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,vin,nIn);
	nLower = pstColRejects(pst->pstLower,vin,nIn,&pOutRej[0],nOut);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,rID,&pOutRej[iUpper],&nUpper);
	return nLower + nUpper;
	}
    else {
	pOutRej->nRejects = pkdColRejects(plcl->pkd, plcl->nSplit);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	pOutRej->nLocal = pkdLocal(plcl->pkd);
	return sizeof(OREJ);
	}
    /*
      pkdStopTimer(plcl->pkd,9);
    */
    }


int pstColOrdRejects(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inColOrdRejects *in = vin;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inColOrdRejects));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,in,nIn);
	nLower = pstColOrdRejects(pst->pstLower,in,nIn,&pOutRej[0],nOut);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,rID,&pOutRej[iUpper],&nUpper);
	return nLower + nUpper;
	}
    else {
	pOutRej->nRejects = pkdColOrdRejects(plcl->pkd,in->iOrdSplit,
					     in->iSplitSide);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	return sizeof(OREJ);
	}
    }


int pstSwapRejects(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    int *pidSwap = vin;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper,idSwap;

    /*	pkdStartTimer(plcl->pkd,8);*/
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,vin,nIn);
	nLower = pstSwapRejects(pst->pstLower,vin,nIn,&pOutRej[0],nOut);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,rID,&pOutRej[iUpper],&nUpper);
	return nLower + nUpper;
	}
    else {
	idSwap = pidSwap[pst->idSelf];
	pOutRej->nRejects = pkdSwapRejects(plcl->pkd,idSwap);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	return sizeof(OREJ);
	}
    /*	pkdStopTimer(plcl->pkd,8); */
    }

/*
 * Routine to swap all particles.  Note that this does not walk the pst
 * but simply works with one other processor.
 */
int pstSwapAll(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl;
    int *pidSwap = vin;
    PST lpst;

    mdlassert(pst->mdl,nIn == sizeof(*pidSwap));
    lpst = pst;
    while (lpst->nLeaves > 1) lpst = lpst->pstLower;
    plcl = lpst->plcl;
    pkdSwapAll(plcl->pkd, *pidSwap);
    return 0;
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
    pstFreeStore(pst->pstLower,NULL,0,&outFree,sizeof(outFree));
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
	pstOrdWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
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
    nOut = pstColOrdRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,mdlThreads(pst->mdl)*sizeof(OREJ));
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
    mdlGetReply(pst->mdl,rID,pUpperRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
    while (1) {
	iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,pst->nUpper,
			    pUpperRej,pidSwap);
	if (!iRet) break;
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
		      mdlThreads(pst->mdl)*sizeof(int));
	nOut = pstSwapRejects(pst->pstLower,pidSwap,mdlThreads(pst->mdl)*sizeof(int),
		       pLowerRej,mdlThreads(pst->mdl)*sizeof(OREJ));
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,rID,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
	}
    free(pLowerRej);
    free(pUpperRej);
    free(pidSwap);
    return pst->iOrdSplit;
    }


int pstDomainOrder(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	if (pst->nLower > 1) pstDomainOrder(pst->pstLower,in,nIn,NULL,0);
	if (pst->nUpper > 1) mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    return 0;
    }


int pstLocalOrder(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstLocalOrder(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdLocalOrder(plcl->pkd,iMinOrder,iMaxOrder);
	}
    return 0;
    }


int pstActiveOrder(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    uint64_t *pnActive = vout;
    uint64_t nActiveLeaf;

    mdlassert(pst->mdl,nIn == 0);
    /*
      pkdStartTimer(pst->plcl->pkd,5);
    */
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEORDER,NULL,0);
	pstActiveOrder(pst->pstLower,NULL,0,vout,nOut);
	mdlGetReply(pst->mdl,rID,&nActiveLeaf,NULL);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveOrder(plcl->pkd);
	}
    return sizeof(uint64_t);
    /*
      pkdStopTimer(pst->plcl->pkd,5);
    */
    }

int pstAddWriteStart(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inAddWriteStart *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inAddWriteStart));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ADDWRITESTART,in,nIn);
	pstAddWriteStart(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	plcl->nWriteStart += in->nWriteStart;
	}
    return 0;
    }

int pstCompressASCII(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstCompressASCII(pst->pstLower,in,nIn,vout,nOut);
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
    return sizeof(struct outCompressASCII);
    }

int pstWriteASCII(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inWriteASCII *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inWriteASCII));
    if (pst->nLeaves > 1) {
	/* Instruct all processors to compress their data, and count it */
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WRITEASCII,in,nIn);
	pstWriteASCII(pst->pstLower,in,nIn,NULL,0);
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
    return 0;
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

int pstRestore(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	    pstRestore(pst->pstLower,in,nIn,NULL,0);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	else { /* Serialize these processors now */
	    pstRestore(pst->pstLower,in,nIn,NULL,0);
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
    return 0;
    }

int pstCheckpoint(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	    pstCheckpoint(pst->pstLower,in,nIn,NULL,0);
	    mdlGetReply(pst->mdl,rID,NULL,NULL);
	    }
	else { /* Serialize these processors now */
	    pstCheckpoint(pst->pstLower,in,nIn,NULL,0);
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
    return 0;
    }

int pstSendParticles(PST pst,void *vin,int nIn,void *vout,int nOut) {
    int iTo = *(int *)vin;
    pkdWriteViaNode(pst->plcl->pkd, iTo);
    return 0;
    }

int pstWrite(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	    pstWrite(pst->pstLower,in,nIn,NULL,0);
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
	    pstWrite(pst->pstLower,in,nIn,NULL,0);
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

    return 0;
    }


int pstSetSoft(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSetSoft *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSetSoft));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETSOFT,in,nIn);
	pstSetSoft(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSetSoft(plcl->pkd,in->dSoft);
	}
    return 0;
    }

int pstDumpTrees(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inDumpTrees *in = (struct inDumpTrees *)vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inDumpTrees));

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DUMPTREES,vin,nIn);
	pstDumpTrees(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdDumpTrees(pkd,in->bOnlyVA,in->uRungDD);
	}
    return 0;
    }

int pstTreeInitMarked(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    mdlassert(pst->mdl,nIn == 0);

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_TREEINITMARKED,vin,nIn);
	pstTreeInitMarked(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdTreeInitMarked(pkd);
	}
    return 0;
    }

int pstDistribTopTree(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    KDN *pTop = vin;
    struct inDistribTopTree *in = vin;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBTOPTREE,vin,nIn);
	pstDistribTopTree(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdDistribTopTree(pkd,in->uRoot,in->nTop,(KDN *)(in+1));
	}
    return 0;
    }

int pstBuildTree(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	nOut = pstBuildTree(pst->pstLower,vin,nIn,pCell1,(pst->nLower*2-1) * pkdNodeSize(pkd));
	assert(nOut == (pst->nLower*2-1) * pkdNodeSize(pkd));
	mdlGetReply(pst->mdl,rID,pCell2,&nOutUpper);
	assert(nOutUpper == (pst->nUpper*2-1) * pkdNodeSize(pkd));
	nOut += nOutUpper + pkdNodeSize(pkd);

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
	nOut = pkdNodeSize(pkd);
	}
    return nOut;
    }

int pstCalcRoot(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inCalcRoot *in = vin;
    struct outCalcRoot *out = vout;
    struct outCalcRoot temp;

    mdlassert(pst->mdl,nIn == sizeof(struct inCalcRoot) );
    mdlassert(pst->mdl,nOut == sizeof(struct outCalcRoot) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCROOT,vin,nIn);
	pstCalcRoot(pst->pstLower,vin,nIn,out,nOut);
	mdlGetReply(pst->mdl,rID,&temp,NULL);
	momAddMomc(&out->momc,&temp.momc);
	}
    else {
	pkdCalcRoot(plcl->pkd,in->uRoot,in->com,&out->momc);
	}
    return sizeof(struct outCalcRoot);
    }

int pstDistribRoot(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct ioDistribRoot *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct ioDistribRoot));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBROOT,vin,nIn);
	pstDistribRoot(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdDistribRoot(plcl->pkd,in->r,&in->momc);
	}
    return 0;
    }


int pstEnforcePeriodic(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    BND *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(BND));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ENFORCEPERIODIC,vin,nIn);
	pstEnforcePeriodic(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdEnforcePeriodic(plcl->pkd,in);
	}
    return 0;
    }

int pstPhysicalSoft(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inPhysicalSoft *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPhysicalSoft));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PHYSICALSOFT,in,nIn);
	pstPhysicalSoft(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdPhysicalSoft(plcl->pkd,in->dSoftMax,in->dFac,in->bSoftMaxMul);
	}
    return 0;
    }

int pstHopLink(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inHopLink *in = vin;
    uint64_t *nOutGroups = (uint64_t *)vout;
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopLink));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_LINK,in,nIn);
	pstHopLink(pst->pstLower,in,nIn,vout,nOut);
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
    return sizeof(uint64_t);
    }

int pstHopJoin(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inHopLink *in = vin;
    struct outHopJoin *out = vout;
    struct outHopJoin outUpper, outLower;
    int nLocal;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopLink));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_JOIN,in,nIn);
        pstHopJoin(pst->pstLower,in,nIn,&outLower,nOut);
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
    return sizeof(struct outHopJoin);
    }

int pstHopFinishUp(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inHopFinishUp *in = (struct inHopFinishUp *)vin;
    uint64_t *nOutGroups = (uint64_t *)vout;
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopFinishUp));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_FINISH_UP,vin,nIn);
        pstHopFinishUp(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&nOutUpper,NULL);
	*nOutGroups += nOutUpper;
        }
    else {
	LCL *plcl = pst->plcl;
        *nOutGroups = pkdHopFinishUp(plcl->pkd,in->nMinGroupSize,in->bPeriodic,in->fPeriod);
        }
    return sizeof(uint64_t);
    }

int pstHopTreeBuild(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inHopTreeBuild *in = (struct inHopTreeBuild *)vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inHopTreeBuild));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_TREE_BUILD,vin,nIn);
        pstHopTreeBuild(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	LCL *plcl = pst->plcl;
        pkdHopTreeBuild(plcl->pkd,in->nBucket,in->nGroup);
        }
    return 0;
    }

int pstHopGravity(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inHopGravity *in = (struct inHopGravity *)vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inHopGravity));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_GRAVITY,vin,nIn);
        pstHopGravity(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	LCL *plcl = pst->plcl;
	double dFlop=0, dPartSum=0, dCellSum=0;
	pkdGravWalkHop(plcl->pkd,in->dTime,in->nGroup,in->dThetaMin,&dFlop,&dPartSum,&dCellSum);
        }
    return 0;
    }

int pstHopUnbind(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inHopUnbind *in = (struct inHopUnbind *)vin;
    struct outHopUnbind *out = (struct outHopUnbind *)vout;
    struct outHopUnbind outUpper, outLower;
    mdlassert(pst->mdl,nIn == sizeof(struct inHopUnbind));
    mdlassert(pst->mdl,nOut == sizeof(struct outHopUnbind));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOP_UNBIND,vin,nIn);
        pstHopUnbind(pst->pstLower,vin,nIn,&outLower,sizeof(outLower));
        mdlGetReply(pst->mdl,rID,&outUpper,NULL);
	out->nEvaporated = outLower.nEvaporated + outUpper.nEvaporated;
	out->nGroups = outLower.nGroups + outUpper.nGroups;
        }
    else {
	LCL *plcl = pst->plcl;
        out->nEvaporated = pkdHopUnbind(plcl->pkd,in->dTime,in->nMinGroupSize,in->bPeriodic,in->fPeriod);
	out->nGroups = plcl->pkd->nLocalGroups;
        }
    return sizeof(struct outHopUnbind);
    }

int pstGroupRelocate(PST pst,void *vin,int nIn,void *vout,int nOut) {
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GROUP_RELOCATE,NULL,0);
        pstGroupRelocate(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	LCL *plcl = pst->plcl;
	PKD pkd = plcl->pkd;
        pkdGroupRelocate(pkd,pkd->nGroups,pkd->ga);
        }
    return 0;
    }

int pstGroupStats(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inGroupStats *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inGroupStats));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GROUP_STATS,vin,nIn);
        pstGroupStats(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	LCL *plcl = pst->plcl;
        pkdCalculateGroupStats(plcl->pkd,in->bPeriodic,in->dPeriod,in->rEnvironment);
        }
    return 0;
    }

int pstSmooth(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SMOOTH,in,nIn);
	pstSmooth(pst->pstLower,in,nIn,NULL,0);
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
    return 0;
    }


#ifdef FAST_GAS
int pstFastGasPhase1(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FASTGASPHASE1,in,nIn);
	pstFastGasPhase1(pst->pstLower,in,nIn,NULL,0);
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
    return 0;
    }


int pstFastGasPhase2(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FASTGASPHASE2,in,nIn);
	pstFastGasPhase2(pst->pstLower,in,nIn,NULL,0);
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
    return 0;
    }


int pstFastGasCleanup(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FASTGASCLEANUP,NULL,0);
	pstFastGasCleanup(pst->pstLower,NULL,0,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdFastGasCleanup(plcl->pkd);
	}
    return 0;
    }
#endif

int pstReSmooth(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTH,in,nIn);
	pstReSmooth(pst->pstLower,in,nIn,NULL,0);
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
    return 0;
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
    
int pstGravity(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstGravity(pst->pstLower,in,nIn,out,nOut);
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
    return (mdlThreads(pst->mdl) - pst->idSelf)*sizeof(struct outGravityPerProc) + sizeof(struct outGravityReduct);
    }

int pstCalcEandL(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outCalcEandL *out = vout;
    struct outCalcEandL outLcl;
    int k;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCEANDL,NULL,0);
	pstCalcEandL(pst->pstLower,NULL,0,out,nOut);
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
    return sizeof(struct outCalcEandL);
    }


int pstDrift(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inDrift *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DRIFT,in,nIn);
	pstDrift(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdDrift(plcl->pkd,in->iRoot,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaUPred);
	}
    return 0;
    }

int pstScaleVel(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inScaleVel *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inScaleVel));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SCALEVEL,in,nIn);
	pstScaleVel(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdScaleVel(plcl->pkd,in->dvFac);
	}
    return 0;
    }

int pstCacheBarrier(PST pst,void *vin,int nIn,void *vout,int nOut) {
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	pstCacheBarrier(pst->pstLower,NULL,0,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	mdlCacheBarrier(pst->mdl,CID_CELL);
	}
    return 0;
    }

int pstStepVeryActiveKDK(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inStepVeryActive *in = vin;
    struct outStepVeryActive *out = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inStepVeryActive));
    if (pst->nLeaves > 1) {
	if (pst->iVASplitSide > 0) {
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	    pstStepVeryActiveKDK(pst->pstLower,in,nIn,out,nOut);
	    mdlGetReply(pst->mdl,rID,NULL,NULL); // FIXME: Shouldn't we look at nMaxRung?
	    }
	else if (pst->iVASplitSide < 0) {
	    int rID = mdlReqService(pst->mdl,pst->idUpper,PST_STEPVERYACTIVE,in,nIn);
	    pstCacheBarrier(pst->pstLower,NULL,0,NULL,0);
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
    return sizeof(*out);
    }

int pstROParticleCache(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ROPARTICLECACHE,vin,nIn);
	pstROParticleCache(pst->pstLower,NULL,0,NULL,0);
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
    return 0;
    }

int pstParticleCacheFinish(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PARTICLECACHEFINISH,vin,nIn);
	pstParticleCacheFinish(pst->pstLower,NULL,0,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	/*
	** Start particle caching space.
	*/
	PKD pkd = plcl->pkd;
	mdlFinishCache(pkd->mdl,CID_PARTICLE);
	}
    return 0;
    }

int pstKick(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inKick *in = vin;
    struct outKick *out = vout;
    struct outKick outUp;

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_KICK,in,nIn);
	pstKick(pst->pstLower,in,nIn,out,nOut);
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
    return sizeof(struct outKick);
    }

int pstKickTree(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inKickTree *in = vin;
    struct outKickTree *out = vout;
    struct outKickTree outUp;

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_KICKTREE,in,nIn);
	pstKickTree(pst->pstLower,in,nIn,out,nOut);
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
    return sizeof(struct outKickTree);
    }

int pstSetTotal(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outSetTotal *out = vout;
    struct outSetTotal oute;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETTOTAL,NULL,0);
	pstSetTotal(pst->pstLower,NULL,0,out,nOut);
	mdlGetReply(pst->mdl,rID,&oute,NULL);
	out->nTotal += oute.nTotal;
	pst->nTotal = out->nTotal;
	}
    else {
	/*pst->nTotal = pkdLocal(plcl->pkd);*/
	pst->nTotal = pkdLocal(plcl->pkd);
	out->nTotal = pst->nTotal;
	}
    return sizeof(struct outSetTotal);
    }


int pstSetWriteStart(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSetWriteStart *in = vin;
    uint64_t nWriteStart;

    mdlassert(pst->mdl,nIn == sizeof(struct inSetWriteStart));
    nWriteStart = in->nWriteStart;
    if (pst->nLeaves > 1) {
	in->nWriteStart = nWriteStart + pst->pstLower->nTotal;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETWRITESTART,in,nIn);
	in->nWriteStart = nWriteStart;
	pstSetWriteStart(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	plcl->nWriteStart = nWriteStart;
	}
    return 0;
    }

int pstSetRung(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSetRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETRUNG,vin,nIn);
	pstSetRung(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSetRung(plcl->pkd,in->uRungLo,in->uRungHi,in->uRung);
	}
    return 0;
    }

int pstZeroNewRung(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inZeroNewRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ZERONEWRUNG,vin,nIn);
	pstZeroNewRung(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdZeroNewRung(plcl->pkd,in->uRungLo,in->uRungHi,in->uRung);
	}
    return 0;
    }

int pstInitCosmology(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == sizeof(struct csmVariables));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITCOSMOLOGY,vin,nIn);
	pstInitCosmology(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	struct csmVariables *cosmo = vin;
	pkdInitCosmology(plcl->pkd,cosmo);
	}
    return 0;
    }

int pstSetParameters(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == sizeof(struct parameters));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETPARAMETERS,vin,nIn);
	pstSetParameters(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
    	struct parameters *param = vin;
	pkdSetParameters(plcl->pkd,param);
	}
    return 0;
    }

int pstActiveRung(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inActiveRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVERUNG,vin,nIn);
	pstActiveRung(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdActiveRung(plcl->pkd, in->iRung, in->bGreater);
	}
    return 0;
    }

int pstCountRungs(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outCountRungs *out = vout, outUpper;
    int i;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COUNTRUNGS,vin,nIn);
	pstCountRungs(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,NULL);
	for(i=0; i<=MAX_RUNG; ++i) out->nRungs[i] += outUpper.nRungs[i];
	}
    else {
	pkdCountRungs(plcl->pkd, out->nRungs);
	}
    return sizeof(*out);
    }

int pstAccelStep(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inAccelStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inAccelStep));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ACCELSTEP,vin,nIn);
	pstAccelStep(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdAccelStep(plcl->pkd,in->uRungLo,in->uRungHi,in->dEta,in->dVelFac,in->dAccFac,
		     in->bDoGravity,in->bEpsAcc,in->dhMinOverSoft);
	}
    return 0;
    }

int pstSphStep(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSphStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSphStep));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SPHSTEP,vin,nIn);
	pstSphStep(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSphStep(plcl->pkd,in->uRungLo,in->uRungHi,in->dAccFac);
	}
    return 0;
    }

int pstStarForm(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inStarForm *in = vin;
    struct outStarForm *out = vout;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inStarForm));
    if (pst->nLeaves > 1) {
	struct outStarForm fsStats;
	
	rID = mdlReqService(pst->mdl,pst->idUpper,PST_STARFORM,in,nIn);
	pstStarForm(pst->pstLower,in,nIn,vout,nOut);
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
    return sizeof(struct outStarForm);
    }

int pstDensityStep(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inDensityStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DENSITYSTEP,vin,nIn);
	pstDensityStep(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdDensityStep(plcl->pkd,in->uRungLo,in->uRungHi,in->dEta,in->dRhoFac);
	}
    return 0;
    }

int pstCorrectEnergy(PST pst,void *vin,int nIn,void *vout,int nOut)
    {
    LCL *plcl = pst->plcl;
    struct inCorrectEnergy *in = vin;
    
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CORRECTENERGY,vin,nIn);
	pstCorrectEnergy(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdCorrectEnergy(plcl->pkd,in->dTuFac,in->z,in->dTime,in->iDirection);
	}
    return 0;
    }

int pstUpdateRung(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outUpdateRung outTemp;
    struct inUpdateRung *in = vin;
    struct outUpdateRung *out = vout;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_UPDATERUNG,vin,nIn);
	pstUpdateRung(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outTemp,NULL);
	for (i=0;i<in->uMaxRung;++i) {
	    out->nRungCount[i] += outTemp.nRungCount[i];
	    }
	}
    else {
	pkdUpdateRung(plcl->pkd,in->uRungLo,in->uRungHi,
		      in->uMinRung,in->uMaxRung,out->nRungCount);
	}
    return sizeof(*out);
    }

int pstUpdateRungByTree(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inUpdateRungByTree *in = vin;
    struct outUpdateRung *out = vout;
    struct outUpdateRung outTemp;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_UPDATERUNGBYTREE,vin,nIn);
	pstUpdateRungByTree(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outTemp,NULL);
	for (i=0;i<in->uMaxRung;++i) {
	    out->nRungCount[i] += outTemp.nRungCount[i];
	    }
	}
    else {
	pkdUpdateRungByTree(plcl->pkd,in->iRoot,in->uMinRung,in->uMaxRung,out->nRungCount);
	}
    return sizeof(*out);
    }

int pstSetRungVeryActive(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSetRung *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETRUNGVERYACTIVE,vin,nIn);
	pstSetRungVeryActive(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSetRungVeryActive(plcl->pkd,in->uRung);
	}
    return 0;
    }


int pstColNParts(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outColNParts *out = vout;
    struct outColNParts *outUp = out + pst->idUpper-pst->idSelf;
    int i;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COLNPARTS,vin,nIn);
	pstColNParts(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,outUp,NULL);
	}
    else {
	pkdColNParts(plcl->pkd, &out->nNew,
		     &out->nDeltaGas,
		     &out->nDeltaDark,
		     &out->nDeltaStar);
	}
    return pst->nLeaves*sizeof(*out);
    }

int pstNewOrder(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    uint64_t *in = vin;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_NEWORDER,vin,nIn);
	pstNewOrder(pst->pstLower, vin, nIn, NULL, 0);
	mdlGetReply(pst->mdl, rID, NULL, NULL);
	}
    else {
	pkdNewOrder(plcl->pkd, in[pst->idSelf]);
	}
    return  0;
    }

int pstGetNParts(PST pst,void *vin,int nIn,void *vout,int nOut)
    {
    struct outGetNParts *out = vout;
    
    if(pst->nLeaves > 1) {
	struct outGetNParts outtmp;
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GETNPARTS,vin,nIn);
	pstGetNParts(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,(void *) &outtmp,NULL);
	
	out->n += outtmp.n;
	out->nGas += outtmp.nGas;
	out->nDark += outtmp.nDark;
	out->nStar += outtmp.nStar;
	if (outtmp.nMaxOrder > out->nMaxOrder) out->nMaxOrder = outtmp.nMaxOrder;
	}
    else {
	pkdGetNParts(pst->plcl->pkd, out);
	}
    return sizeof(struct outGetNParts);
}

int pstSetNParts(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inSetNParts *in = vin;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETNPARTS,vin,nIn);
	pstSetNParts(pst->pstLower, vin, nIn, NULL, 0);
	mdlGetReply(pst->mdl, rID, NULL, NULL);
	}
    else {
	pkdSetNParts(pst->plcl->pkd, in->nGas, in->nDark, in->nStar);
	}
    return 0;
    }

int pstClearTimer(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inClearTimer *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inClearTimer));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CLEARTIMER,in,nIn);
	pstClearTimer(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdClearTimer(pst->plcl->pkd,in->iTimer);
	}
    return 0;
    }


int pstNewFof(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inNewFof *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inNewFof));
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_NEW_FOF,in,nIn);
	pstNewFof(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	pkdNewFof(plcl->pkd,in->dTau2,in->nMinMembers);
	}
    return 0;
    }

int pstFofPhases(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct outFofPhases *out = vout;
    int bMadeProgress;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FOF_PHASES,vin,nIn);
	pstFofPhases(pst->pstLower,vin,nIn,out,nOut);
	bMadeProgress = out->bMadeProgress;
	mdlGetReply(pst->mdl,rID,out,NULL);
	if (!out->bMadeProgress) out->bMadeProgress = bMadeProgress;
	}
    else {
	LCL *plcl = pst->plcl;
	out->bMadeProgress = pkdFofPhases(plcl->pkd);
	}
    return sizeof(struct outFofPhases);
    }


/*
** This is an almost identical copy of HopFinishUp.
*/
int pstFofFinishUp(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inFofFinishUp *in = (struct inFofFinishUp *)vin;
    uint64_t *nOutGroups = (uint64_t *)vout;
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inFofFinishUp));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FOF_FINISH_UP,vin,nIn);
        pstFofFinishUp(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&nOutUpper,NULL);
	*nOutGroups += nOutUpper;
        }
    else {
	LCL *plcl = pst->plcl;
        *nOutGroups = pkdFofFinishUp(plcl->pkd,in->nMinGroupSize);
        }
    return sizeof(uint64_t);
    }


int pstInitRelaxation(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITRELAXATION,vin,nIn);
	pstInitRelaxation(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdInitRelaxation(plcl->pkd);
	}
    return 0;
    }

#ifdef MDL_FFTW
int pstGetFFTMaxSizes(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstGetFFTMaxSizes(pst->pstLower,in,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUp,NULL);
	if (outUp.nMaxLocal > out->nMaxLocal) out->nMaxLocal = outUp.nMaxLocal;
	if (outUp.nMaxZ > out->nMaxZ) out->nMaxZ = outUp.nMaxZ;
	if (outUp.nMaxY > out->nMaxY) out->nMaxY = outUp.nMaxZ;
	}
    else {
	assert(pstAmNode(pst));
	out->nMaxLocal = mdlFFTlocalCount(pst->mdl,in->nx,in->ny,in->nz,
	    &out->nMaxZ,0,&out->nMaxY,0);
	}
    return sizeof(struct outGetFFTMaxSizes);
    }
#endif

int pstHostname(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct outHostname *out = vout;
    struct outHostname *outUp = out + pst->idUpper-pst->idSelf;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_HOSTNAME,vin,nIn);
	pstHostname(pst->pstLower,vin,nIn,out,nOut);
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
    return pst->nLeaves*sizeof(struct outHostname);
    }

int pstMemStatus(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outMemStatus *out = vout;
    struct outMemStatus *outUp = out + pst->idUpper-pst->idSelf;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEMSTATUS,vin,nIn);
	pstMemStatus(pst->pstLower,vin,nIn,out,nOut);
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
    return pst->nLeaves*sizeof(struct outMemStatus);
    }


int pstGetClasses(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PARTCLASS *out = vout;
    PARTCLASS *outUp;
    int nUp;
    int n, i, j;

    mdlassert(pst->mdl,nIn==0);
    if (pst->nLeaves > 1) {
	outUp = malloc(PKD_MAX_CLASSES*sizeof(PARTCLASS));
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GETCLASSES,vin,nIn);
	nOut = pstGetClasses(pst->pstLower,vin,nIn,out,nOut);
	mdlGetReply(pst->mdl,rID,outUp,&nUp);
	n = nUp / sizeof(PARTCLASS);
	mdlassert(pst->mdl,n*sizeof(PARTCLASS)==nUp);
	nUp = n;
	n = nOut / sizeof(PARTCLASS);
	mdlassert(pst->mdl,n*sizeof(PARTCLASS)== nOut);
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
	}
    else {
	n = pkdGetClasses(plcl->pkd,PKD_MAX_CLASSES,vout);
	}
    return n * sizeof(PARTCLASS);
    }

int pstSetClasses(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PARTCLASS *in = vin;
    int n;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETCLASSES,vin,nIn);
	pstSetClasses(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	n = nIn / sizeof(PARTCLASS);
	mdlassert(pst->mdl,n*sizeof(PARTCLASS)==nIn);
	pkdSetClasses(plcl->pkd,n,in,1);
	}
    return 0;
    }

/*
 * Routine to swap the class table .  Note that this does not walk
 * the pst but simply set's the given table and returns the old table.
 */
int pstSwapClasses(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
    nOut = n * sizeof(PARTCLASS);

    n = nIn / sizeof(PARTCLASS);
    mdlassert(pst->mdl,n*sizeof(PARTCLASS) == nIn);
    pkdSetClasses( plcl->pkd, n, in, 0 );
    return nOut;
    }


int pstSelAll(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELALL,vin,nIn);
	pstSelAll(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSelAll(plcl->pkd);
	}
    return 0;
    }

int pstSelGas(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELGAS,vin,nIn);
	pstSelGas(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSelGas(plcl->pkd);
	}
    return 0;
    }

int pstSelStar(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSTAR,vin,nIn);
	pstSelStar(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSelStar(plcl->pkd);
	}
    return 0;
    }

int pstSelDeleted(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELDELETED,vin,nIn);
	pstSelDeleted(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSelDeleted(plcl->pkd);
	}
    return 0;
    }

int pstSelById(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSelById *in = vin;
    struct outSelById *out = vout;
    struct outSelById outUpper;

    assert( nIn==sizeof(struct inSelById) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELBYID,vin,nIn);
	pstSelById(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelById));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelById(plcl->pkd,in->idStart,in->idEnd,in->setIfTrue,in->clearIfFalse);
	}
    return sizeof(struct outSelById);
    }
int pstSelMass(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSelMass *in = vin;
    struct outSelMass *out = vout;
    struct outSelMass outUpper;

    assert( nIn==sizeof(struct inSelMass) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELMASS,vin,nIn);
	pstSelMass(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelMass));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelMass(plcl->pkd,in->dMinMass,in->dMaxMass,in->setIfTrue,in->clearIfFalse);
	}
    return sizeof(struct outSelMass);
    }

int pstSelPhaseDensity(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSelPhaseDensity *in = vin;
    struct outSelPhaseDensity *out = vout;
    struct outSelPhaseDensity outUpper;

    assert( nIn==sizeof(struct inSelPhaseDensity) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELPHASEDENSITY,vin,nIn);
	pstSelPhaseDensity(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelPhaseDensity));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelPhaseDensity(plcl->pkd,in->dMinDensity,in->dMaxDensity,in->setIfTrue,in->clearIfFalse);
	}
    return sizeof(struct outSelPhaseDensity);
    }

int pstSelBox(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSelBox *in = vin;
    struct outSelBox *out = vout;
    struct outSelBox outUpper;

    assert( nIn==sizeof(struct inSelBox) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELBOX,vin,nIn);
	pstSelBox(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelBox));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelBox(
	    plcl->pkd,in->dCenter,in->dSize,in->setIfTrue,in->clearIfFalse);
	}
    return sizeof(struct outSelBox);
    }

int pstSelSphere(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSelSphere *in = vin;
    struct outSelSphere *out = vout;
    struct outSelSphere outUpper;

    assert( nIn==sizeof(struct inSelSphere) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELSPHERE,vin,nIn);
	pstSelSphere(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelSphere));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelSphere(
	    plcl->pkd,in->r,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    return sizeof(struct outSelSphere);
    }

int pstSelCylinder(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSelCylinder *in = vin;
    struct outSelCylinder *out = vout;
    struct outSelCylinder outUpper;

    assert( nIn==sizeof(struct inSelCylinder) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELCYLINDER,vin,nIn);
	pstSelCylinder(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut == sizeof(struct outSelCylinder));
	out->nSelected += outUpper.nSelected;
	}
    else {
	out->nSelected = pkdSelCylinder(
	    plcl->pkd,in->dP1,in->dP2,in->dRadius,in->setIfTrue,in->clearIfFalse);
	}
    return sizeof(struct outSelCylinder);
    }

int pstSelGroup(PST pst,void *vin,int nIn,void *vout,int nOut) {
    int *in = vin;
    LCL *plcl = pst->plcl;
    assert( nIn==sizeof(int) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SELGROUP,vin,nIn);
	pstSelGroup(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdSelGroup(plcl->pkd, *in);
	}
    return 0;
    }

int pstProfile(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inProfile *in = vin;
    /*assert( nIn==sizeof(struct inProfile) );*/
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_PROFILE,vin,nIn);
	pstProfile(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdProfile(plcl->pkd,in->uRungLo,in->uRungHi,
		   in->dCenter, in->dRadii, in->nBins,
		   in->com, in->vcm, in->L);
	}
    return 0;
    }

int pstCalcDistance(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inCalcDistance *in = vin;

    assert( nIn==sizeof(struct inCalcDistance) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCDISTANCE,vin,nIn);
	pstCalcDistance(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdCalcDistance(plcl->pkd,in->dCenter);
	}
    return 0;
    }

int pstCalcCOM(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inCalcCOM *in = vin;
    struct outCalcCOM *out = vout;
    struct outCalcCOM outUpper;
    int i;

    assert( nIn==sizeof(struct inCalcCOM) );
    assert( vout != NULL );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCCOM,vin,nIn);
	pstCalcCOM(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,NULL);
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
    return sizeof(struct outCalcCOM);
    }

int pstCountDistance(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inCountDistance *in = vin;
    struct outCountDistance *out = vout;
    struct outCountDistance outUpper;

    assert( nIn==sizeof(struct inCountDistance) );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COUNTDISTANCE,vin,nIn);
	pstCountDistance(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut==sizeof(struct outCountDistance));
	out->nCount += outUpper.nCount;
	}
    else {
	out->nCount = pkdCountDistance(plcl->pkd,in->dRadius2Inner,in->dRadius2Outer);
	}
    return sizeof(struct outCountDistance);
    }

int pstInitGrid(PST pst,void *vin,int nIn,void *vout,int nOut) {
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
	pstInitGrid(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdGridInitialize(plcl->pkd,in->n1,in->n2,in->n3,in->a1,in->s,in->n);
	}
    return 0;
    }

int pstGridProject(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    /*struct inGridProject *in = vin;*/
    assert( sizeof(struct inGridProject) == nIn );
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GRIDPROJECT,vin,nIn);
	pstGridProject(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdGridProject(plcl->pkd);
	}
    return 0;
    }

#ifdef MDL_FFTW
int pstMeasurePk(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inMeasurePk *in = vin;
    struct outMeasurePk *out = vout;
    struct outMeasurePk *outUpper;
    int i;

    assert( nIn==sizeof(struct inMeasurePk) );
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEASUREPK,vin,nIn);
	pstMeasurePk(pst->pstLower,vin,nIn,vout,nOut);
	outUpper = malloc(sizeof(struct outMeasurePk));
	assert(outUpper != NULL);
	mdlGetReply(pst->mdl,rID,outUpper,&nOut);
	assert(nOut==sizeof(struct outMeasurePk));

	for(i=0;i<in->nBins; i++) {
	    out->fK[i] += outUpper->fK[i];
	    out->fPower[i] += outUpper->fPower[i];
	    out->fPowerAll[i] += outUpper->fPowerAll[i];
	    out->nPower[i] += outUpper->nPower[i];
	    }
	free(outUpper);
	}
    else {
	pkdMeasurePk(plcl->pkd, in->dTotalMass, in->iAssignment, in->bInterlace,
	    in->bLinear, in->iSeed, in->bFixed, in->fPhase, in->Lbox, in->a,
	    in->nGrid, in->nBins, out->fK, out->fPower, out->nPower, out->fPowerAll);
	}
    return sizeof(struct outMeasurePk);
    }

int pstGridCreateFFT(PST pst,void *vin,int nIn,void *vout,int nOut) {
        LCL *plcl = pst->plcl;
        struct inGridCreateFFT *in = vin;
        assert (nIn==sizeof(struct inGridCreateFFT) );
        if (pstNotCore(pst)) {
            int rID = mdlReqService(pst->mdl, pst->idUpper, PST_GRID_CREATE_FFT, vin, nIn);
            pstGridCreateFFT(pst->pstLower, vin, nIn, NULL, 0);
            mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
        else {
            PKD pkd = plcl->pkd;
            assert(pkd->fft == NULL);
            pkd->fft = mdlFFTInitialize(pst->mdl, in->nGrid, in->nGrid, in->nGrid, 0, 0);
        }
    return 0;
    }

int pstGridDeleteFFT(PST pst,void *vin,int nIn,void *vout,int nOut) {
        LCL *plcl = pst->plcl;
        if (pstNotCore(pst)) {
            int rID = mdlReqService(pst->mdl, pst->idUpper, PST_GRID_DELETE_FFT, vin, nIn);
            pstGridDeleteFFT(pst->pstLower, vin, nIn, NULL, 0);
            mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
        else {
            PKD pkd = plcl->pkd;
            assert(pkd->fft != NULL);
            mdlFFTFinish(pst->mdl,plcl->pkd->fft);
            pkd->fft = NULL;
        }
    return 0;
    }
#endif

int pstTotalMass(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outTotalMass *out = vout;
    struct outTotalMass outUpper;

    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_TOTALMASS,vin,nIn);
	pstTotalMass(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
	assert(nOut==sizeof(struct outTotalMass));
	out->dMass += outUpper.dMass;
	}
    else {
	out->dMass = pkdTotalMass(plcl->pkd);
	}
    return sizeof(struct outTotalMass);
    }

int pstLightConeOpen(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inLightConeOpen *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inLightConeOpen));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONE_OPEN,in,nIn);
        pstLightConeOpen(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	PKD pkd = pst->plcl->pkd;
	char achOutFile[PST_FILENAME_SIZE];
	if (in->achOutFile[0]) makeName(achOutFile,in->achOutFile,mdlSelf(pkd->mdl),"lcp.");
	else achOutFile[0] = 0;
        pkdLightConeOpen(pkd, achOutFile, in->nSideHealpix);
        }
    return 0;
}

int pstLightConeClose(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inLightConeClose *in = vin;
    mdlassert(pst->mdl,nIn == sizeof(struct inLightConeClose));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONE_CLOSE,in,nIn);
        pstLightConeClose(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
        }
    else {
	PKD pkd = pst->plcl->pkd;
	char achOutFile[PST_FILENAME_SIZE];
	if (in->achOutFile[0]) makeName(achOutFile,in->achOutFile,mdlSelf(pkd->mdl),"hpb.");
	else achOutFile[0] = 0;
        pkdLightConeClose(pst->plcl->pkd,achOutFile);
        }
    return 0;
}

int pstLightConeVel(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inScaleVel *in = vin;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONEVEL,in,nIn);
	pstLightConeVel(pst->pstLower,in,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	pkdLightConeVel(plcl->pkd);
	}
    return 0;
    }

int pstInflate(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INFLATE,vin,nIn);
	pstInflate(pst->pstLower,vin,nIn,NULL,0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
	struct inInflate *in = vin;
	pkdInflate(plcl->pkd,in->nInflateReps);
	}
    return 0;
    }

int pstGetParticles(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct outGetParticles *out = vout;
    uint64_t *ID = vin;
    int nOutUpper;
    if (pst->nLeaves > 1) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GET_PARTICLES,vin,nIn);
	nOut = pstGetParticles(pst->pstLower,vin,nIn,vout,nOut);
	mdlGetReply(pst->mdl,rID,out + (nOut / sizeof(struct outGetParticles)),&nOutUpper);
	nOut += nOutUpper;
	}
    else {
	int nParticles = nIn / sizeof(uint64_t);
	int n = pkdGetParticles(plcl->pkd,nParticles, ID, out );
	nOut = n * sizeof(struct outGetParticles);
	}
    return nOut;
    }
