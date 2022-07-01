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
#include "smooth/smooth.h"
#include "group/hop.h"
#include "group/fof.h"
#include "group/group.h"
#include "group/groupstats.h"

#ifdef COOLING
    #include "cooling/cooling.h"
#endif
#ifdef GRACKLE
    #include "cooling_grackle/cooling_grackle.h"
#endif
#ifdef BLACKHOLES
    #include "blackhole/merger.h"
    #include "blackhole/seed.h"
    #include "blackhole/init.h"
#endif
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#ifdef STAR_FORMATION
    #include "starformation/starformation.h"
#endif

void pstAddServices(PST pst,MDL mdl) {
    int nThreads;

    nThreads = mdlThreads(mdl);

    mdlAddService(mdl,PST_READFILE,pst,(fcnService_t *)pstReadFile,
                  sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE),0);
    mdlAddService(mdl,PST_COMPRESSASCII,pst,(fcnService_t *)pstCompressASCII,
                  sizeof(struct inCompressASCII),sizeof(struct outCompressASCII));
    mdlAddService(mdl,PST_WRITEASCII,pst,(fcnService_t *)pstWriteASCII,
                  sizeof(struct inWriteASCII),0);
    mdlAddService(mdl,PST_WRITE,pst,(fcnService_t *)pstWrite,
                  sizeof(struct inWrite),0);
    mdlAddService(mdl,PST_SENDPARTICLES,pst,(fcnService_t *)pstSendParticles,
                  sizeof(int),0);
    mdlAddService(mdl,PST_SENDARRAY,pst,(fcnService_t *)pstSendArray,
                  sizeof(struct inSendArray),0);
    mdlAddService(mdl,PST_CHECKPOINT,pst,(fcnService_t *)pstCheckpoint,
                  sizeof(struct inWrite),0);
    mdlAddService(mdl,PST_OUTPUT,pst,(fcnService_t *)pstOutput,
                  sizeof(struct inOutput),0);
    mdlAddService(mdl,PST_OUTPUT_SEND,pst,(fcnService_t *)pstOutputSend,
                  sizeof(struct inOutputSend),0);
    mdlAddService(mdl,PST_RESTORE,pst,(fcnService_t *)pstRestore,
                  sizeof(struct inRestore),0);
    /*
    ** Calculate the number of levels in the top tree and use it to
    ** define the size of the messages.
    */
    mdlAddService(mdl,PST_BUILDTREE,pst,(fcnService_t *)pstBuildTree,
                  sizeof(struct inBuildTree),
                  (nThreads==1?1:2*nThreads-1)*pkdContext::MaxNodeSize());
    mdlAddService(mdl,PST_TREEINITMARKED,pst,(fcnService_t *)pstTreeInitMarked,
                  0,0);
    mdlAddService(mdl,PST_HOP_LINK,pst,(fcnService_t *)pstHopLink,
                  sizeof(struct inHopLink),sizeof(uint64_t));
    mdlAddService(mdl,PST_HOP_JOIN,pst,(fcnService_t *)pstHopJoin,
                  sizeof(struct inHopLink),sizeof(struct outHopJoin));
    mdlAddService(mdl,PST_HOP_FINISH_UP,pst,(fcnService_t *)pstHopFinishUp,
                  sizeof(struct inHopFinishUp),sizeof(uint64_t));
    mdlAddService(mdl,PST_HOP_TREE_BUILD,pst,(fcnService_t *)pstHopTreeBuild,
                  sizeof(struct inHopTreeBuild),0);
    mdlAddService(mdl,PST_HOP_GRAVITY,pst,(fcnService_t *)pstHopGravity,
                  sizeof(struct inHopGravity),0);
    mdlAddService(mdl,PST_HOP_UNBIND,pst,(fcnService_t *)pstHopUnbind,
                  sizeof(struct inHopUnbind),sizeof(struct outHopUnbind));
    mdlAddService(mdl,PST_GROUP_RELOCATE,pst,(fcnService_t *)pstGroupRelocate,
                  0,0);
    mdlAddService(mdl,PST_GROUP_STATS,pst,(fcnService_t *)pstGroupStats,
                  sizeof(struct inGroupStats),0);
    mdlAddService(mdl,PST_SMOOTH,pst,(fcnService_t *)pstSmooth,
                  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_GRAVITY,pst,(fcnService_t *)pstGravity,
                  sizeof(struct inGravity),
                  sizeof(struct outGravityReduct));
    mdlAddService(mdl,PST_CALCEANDL,pst,(fcnService_t *)pstCalcEandL,
                  0,sizeof(struct outCalcEandL));
    mdlAddService(mdl,PST_DRIFT,pst,(fcnService_t *)pstDrift,
                  sizeof(struct inDrift),0);
    // IA: New PST functions
    mdlAddService(mdl,PST_RESETFLUXES,pst,
                  (fcnService_t *) pstResetFluxes,
                  sizeof(struct inDrift),0);
    mdlAddService(mdl,PST_COMPUTEPRIMVARS,pst,
                  (fcnService_t *) pstEndTimestepIntegration,
                  sizeof(struct inEndTimestep),0);
    mdlAddService(mdl,PST_WAKEPARTICLES,pst,
                  (fcnService_t *) pstWakeParticles,
                  sizeof(struct inDrift),0);
#ifdef DEBUG_CACHED_FLUXES
    mdlAddService(mdl,PST_FLUXSTATS,pst,
                  (fcnService_t *) pstFluxStats,
                  sizeof(struct inFluxStats), sizeof(struct outFluxStats));
#endif
    mdlAddService(mdl,PST_SETGLOBALDT,pst,
                  (fcnService_t *) pstSetGlobalDt,
                  sizeof(struct outGetMinDt),0);
#ifdef COOLING
    mdlAddService(mdl,PST_COOLINGUPDATE,pst,
                  (fcnService_t *) pstCoolingUpdate,
                  sizeof(struct inCoolUpdate),0);
    mdlAddService(mdl,PST_COOLINGUPDATEZ,pst,
                  (fcnService_t *) pstCoolingUpdateZ,
                  sizeof(float),0);
    mdlAddService(mdl,PST_COOLINGINIT,pst,
                  (fcnService_t *) pstCoolingInit,
                  sizeof(struct inCoolInit),0);
    mdlAddService(mdl,PST_COOLINGHYDREION,pst,
                  (fcnService_t *) pstCoolingHydReion,
                  0,0);
#endif
#ifdef GRACKLE
    mdlAddService(mdl,PST_GRACKLEINIT,pst,
                  (fcnService_t *) pstGrackleInit,
                  sizeof(struct inGrackleInit),0);
#endif
    mdlAddService(mdl,PST_CHEMCOMPINIT,pst,
                  (fcnService_t *) pstChemCompInit,
                  sizeof(struct inChemCompInit),0);
#ifdef BLACKHOLES
    mdlAddService(mdl,PST_BH_PLACESEED,pst,
                  (fcnService_t *) pstPlaceBHSeed,
                  sizeof(struct inPlaceBHSeed), sizeof(struct outPlaceBHSeed));
    mdlAddService(mdl,PST_BH_INIT,pst,
                  (fcnService_t *) pstBHInit,
                  sizeof(struct inPlaceBHSeed), 0);
    mdlAddService(mdl,PST_BH_REPOSITION,pst,
                  (fcnService_t *) pstRepositionBH,
                  0, 0);
#endif
    mdlAddService(mdl,PST_MOVEDELETED,pst,
                  (fcnService_t *)pstMoveDeletedParticles,
                  0, sizeof(struct outGetNParts) );
    mdlAddService(mdl,PST_GETMINDT,pst,
                  (fcnService_t *) pstGetMinDt,
                  0, sizeof(struct outGetMinDt));
    //
    mdlAddService(mdl,PST_CACHEBARRIER,pst,(fcnService_t *)pstCacheBarrier,
                  0,0);
    mdlAddService(mdl,PST_ROPARTICLECACHE,pst,(fcnService_t *)pstROParticleCache,
                  0,0);
    mdlAddService(mdl,PST_PARTICLECACHEFINISH,pst,(fcnService_t *)pstParticleCacheFinish,
                  0,0);
    mdlAddService(mdl,PST_KICK,pst,(fcnService_t *)pstKick,
                  sizeof(struct inKick),0);
    mdlAddService(mdl,PST_KICKTREE,pst,(fcnService_t *)pstKickTree,
                  sizeof(struct inKickTree),0);
    mdlAddService(mdl,PST_PHYSICALSOFT,pst,(fcnService_t *)pstPhysicalSoft,
                  sizeof(struct inPhysicalSoft),0);
    mdlAddService(mdl,PST_SETTOTAL,pst,(fcnService_t *)pstSetTotal,
                  0,sizeof(struct outSetTotal));
    mdlAddService(mdl,PST_SETWRITESTART,pst,(fcnService_t *)pstSetWriteStart,
                  sizeof(struct inSetWriteStart),0);
    mdlAddService(mdl,PST_ADDWRITESTART,pst,(fcnService_t *)pstAddWriteStart,
                  sizeof(struct inAddWriteStart),0);
    mdlAddService(mdl,PST_ONENODEREADINIT,pst,(fcnService_t *)pstOneNodeReadInit,
                  sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE),
                  nThreads*sizeof(int));
    mdlAddService(mdl,PST_ACTIVEORDER,pst,(fcnService_t *)pstActiveOrder,
                  0,sizeof(uint64_t));
    mdlAddService(mdl,PST_DENSITYSTEP,pst,(fcnService_t *)pstDensityStep,
                  sizeof(struct inDensityStep),0);
    mdlAddService(mdl,PST_CORRECTENERGY,pst,(fcnService_t *)pstCorrectEnergy,
                  sizeof(struct inCorrectEnergy),0);
    mdlAddService(mdl,PST_ACCELSTEP,pst,(fcnService_t *)pstAccelStep,
                  sizeof(struct inAccelStep), 0);
    mdlAddService(mdl,PST_SPHSTEP,pst,(fcnService_t *)pstSphStep,
                  sizeof(struct inSphStep), 0);
#ifdef STAR_FORMATION
    mdlAddService(mdl,PST_STARFORM,pst,
                  (fcnService_t *) pstStarForm,
                  sizeof(struct inStarForm),sizeof(struct outStarForm));
#endif
#if defined(STAR_FORMATION) || defined(FEEDBACK)
    mdlAddService(mdl,PST_STARFORMINIT,pst,
                  (fcnService_t *) pstStarFormInit,
                  sizeof(struct inStarFormInit),sizeof(struct outStarForm));
#endif
    mdlAddService(mdl,PST_RESMOOTH,pst,(fcnService_t *) pstReSmooth,
                  sizeof(struct inSmooth),sizeof(struct outSmooth));
#ifdef OPTIM_SMOOTH_NODE
    mdlAddService(mdl,PST_RESMOOTHNODE,pst,(fcnService_t *) pstReSmoothNode,
                  sizeof(struct inSmooth),sizeof(struct outSmooth));
#endif
#ifdef OPTIM_REORDER_IN_NODES
    mdlAddService(mdl,PST_REORDERINNODES,pst,(fcnService_t *) pstReorderWithinNodes,
                  0,0);
#endif
#ifdef STELLAR_EVOLUTION
    mdlAddService(mdl,PST_STELLAREVOLUTIONINIT,pst,
                  (fcnService_t *) pstStellarEvolutionInit,
                  sizeof(struct inStellarEvolutionInit),0);
#endif
    mdlAddService(mdl,PST_UPDATERUNG,pst,(fcnService_t *)pstUpdateRung,
                  sizeof(struct inUpdateRung),sizeof(struct outUpdateRung));
    mdlAddService(mdl,PST_COLNPARTS,pst,(fcnService_t *)pstColNParts,
                  0,nThreads*sizeof(struct outColNParts));
    mdlAddService(mdl,PST_NEWORDER,pst,(fcnService_t *)pstNewOrder,
                  nThreads*sizeof(uint64_t),0);
    mdlAddService(mdl,PST_GETNPARTS,pst,(fcnService_t *)pstGetNParts,
                  0,sizeof(struct outGetNParts));
    mdlAddService(mdl,PST_SETNPARTS,pst,(fcnService_t *)pstSetNParts,
                  sizeof(struct inSetNParts),0);
    mdlAddService(mdl,PST_NEW_FOF,pst,(fcnService_t *)pstNewFof,
                  sizeof(struct inNewFof),0);
    mdlAddService(mdl,PST_FOF_PHASES,pst,(fcnService_t *)pstFofPhases,
                  0,sizeof(struct outFofPhases));
    mdlAddService(mdl,PST_FOF_FINISH_UP,pst,(fcnService_t *)pstFofFinishUp,
                  sizeof(struct inFofFinishUp),sizeof(uint64_t));
    mdlAddService(mdl,PST_INITRELAXATION,pst,(fcnService_t *)pstInitRelaxation,0,0);
    mdlAddService(mdl,PST_INITIALIZEPSTORE,pst,(fcnService_t *)pstInitializePStore,
                  sizeof(struct inInitializePStore),0);
#ifdef MDL_FFTW
    mdlAddService(mdl,PST_GETFFTMAXSIZES,pst,(fcnService_t *)pstGetFFTMaxSizes,
                  sizeof(struct inGetFFTMaxSizes),sizeof(struct outGetFFTMaxSizes));
    mdlAddService(mdl,PST_GENERATEIC,pst,(fcnService_t *)pstGenerateIC,
                  sizeof(struct inGenerateIC),sizeof(struct outGenerateIC));
    mdlAddService(mdl,PLT_GENERATEIC,pst,(fcnService_t *)pltGenerateIC,
                  sizeof(struct inGenerateICthread),sizeof(struct outGenerateIC));
    mdlAddService(mdl,PLT_MOVEIC,pst,(fcnService_t *)pltMoveIC,
                  sizeof(struct inMoveIC),0);
    mdlAddService(mdl,PST_MOVEIC,pst,(fcnService_t *)pstMoveIC,
                  sizeof(struct inGenerateIC),0);
#endif
    mdlAddService(mdl,PST_MEMSTATUS,pst,(fcnService_t *)pstMemStatus,
                  0,nThreads*sizeof(struct outMemStatus));
    mdlAddService(mdl,PST_GETCLASSES,pst,(fcnService_t *)pstGetClasses,
                  0, PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlAddService(mdl,PST_SETCLASSES,pst,(fcnService_t *)pstSetClasses,
                  PKD_MAX_CLASSES*sizeof(PARTCLASS), 0);
    mdlAddService(mdl,PST_SWAPCLASSES,pst,(fcnService_t *)pstSwapClasses,
                  PKD_MAX_CLASSES*sizeof(PARTCLASS),
                  PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlAddService(mdl,PST_PROFILE,pst,(fcnService_t *)pstProfile,
                  sizeof(struct inProfile), 0);
    mdlAddService(mdl,PST_CALCDISTANCE,pst,(fcnService_t *)pstCalcDistance,
                  sizeof(struct inCalcDistance), 0);
    mdlAddService(mdl,PST_CALCCOM,pst,(fcnService_t *)pstCalcCOM,
                  sizeof(struct inCalcCOM), sizeof(struct outCalcCOM));
    mdlAddService(mdl,PST_CALCMTOT,pst,(fcnService_t *)pstCalcMtot,
                  sizeof(struct inCalcMtot), sizeof(struct outCalcMtot));
    mdlAddService(mdl,PST_SETSPHOPTIONS,pst,(fcnService_t *)pstSetSPHoptions,
                  sizeof(struct inSetSPHoptions), 0);
    mdlAddService(mdl,PST_TREEUPDATEFLAGBOUNDS,pst,(fcnService_t *)pstTreeUpdateFlagBounds,
                  sizeof(struct inTreeUpdateFlagBounds),
                  (nThreads==1?1:2*nThreads-1)*pkdContext::MaxNodeSize());
    mdlAddService(mdl,PST_COUNTDISTANCE,pst,(fcnService_t *)pstCountDistance,
                  sizeof(struct inCountDistance), sizeof(struct outCountDistance));
#ifdef MDL_FFTW
    mdlAddService(mdl,PST_GRID_CREATE_FFT,pst,(fcnService_t *)pstGridCreateFFT,
                  sizeof(struct inGridCreateFFT), 0);
    mdlAddService(mdl,PST_GRID_DELETE_FFT,pst,(fcnService_t *)pstGridDeleteFFT,
                  0, 0);
    mdlAddService(mdl,PST_ADD_LINEAR_SIGNAL,pst,(fcnService_t *)pstAddLinearSignal,
                  sizeof(struct inAddLinearSignal), 0);
    mdlAddService(mdl,PST_ASSIGN_MASS,pst,(fcnService_t *)pstAssignMass,
                  sizeof(struct inAssignMass), 0);
    mdlAddService(mdl,PST_DENSITY_CONTRAST,pst,(fcnService_t *)pstDensityContrast,
                  sizeof(struct inDensityContrast), 0);
    mdlAddService(mdl,PST_WINDOW_CORRECTION,pst,(fcnService_t *)pstWindowCorrection,
                  sizeof(struct inWindowCorrection), 0);
    mdlAddService(mdl,PST_INTERLACE,pst,(fcnService_t *)pstInterlace,
                  sizeof(struct inInterlace), 0);
    mdlAddService(mdl,PST_GRID_BIN_K,pst,(fcnService_t *)pstGridBinK,
                  sizeof(struct inGridBinK), sizeof(struct outGridBinK));
    mdlAddService(mdl,PST_BISPECTRUM_SELECT,pst,(fcnService_t *)pstBispectrumSelect,
                  sizeof(struct inBispectrumSelect), 0);
    mdlAddService(mdl,PST_BISPECTRUM_CALCULATE,pst,(fcnService_t *)pstBispectrumCalculate,
                  sizeof(struct inBispectrumCalculate), sizeof(double));
    mdlAddService(mdl,PST_LINEARKICK, pst,(fcnService_t *)pstLinearKick,
                  sizeof(struct inLinearKick), 0);
    mdlAddService(mdl,PST_SETLINGRID, pst,(fcnService_t *)pstSetLinGrid,
                  sizeof(struct inSetLinGrid), 0);
    mdlAddService(mdl,PST_MEASURELINPK,pst,(fcnService_t *)pstMeasureLinPk,
                  sizeof(struct inMeasureLinPk), sizeof(struct outMeasureLinPk));
#endif
    mdlAddService(mdl,PST_TOTALMASS,pst,(fcnService_t *)pstTotalMass,
                  0, sizeof(struct outTotalMass));
    mdlAddService(mdl,PST_LIGHTCONE_OPEN,pst,(fcnService_t *)pstLightConeOpen,
                  sizeof(struct inLightConeOpen), 0);
    mdlAddService(mdl,PST_LIGHTCONE_CLOSE,pst,(fcnService_t *)pstLightConeClose,
                  sizeof(struct inLightConeClose), 0);
    mdlAddService(mdl,PST_LIGHTCONEVEL,pst,(fcnService_t *)pstLightConeVel,
                  sizeof(struct inLightConeVel),0);
    mdlAddService(mdl,PST_GET_PARTICLES,pst,(fcnService_t *)pstGetParticles,
                  sizeof(uint64_t)*GET_PARTICLES_MAX,
                  sizeof(struct outGetParticles)*GET_PARTICLES_MAX );
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
    pst->idUpper = -1;  /* invalidate upper 'id' */
    pst->nLeaves = 1;
    pst->nLower = 0;
    pst->nUpper = 0;
    pst->iSplitDim = -1;
}


void pstFinish(PST pst) {
    PST pstKill;

    while (pst) {
        pstKill = pst;
        if (pst->nLeaves == 1 && pst->plcl->pkd)
            delete pst->plcl->pkd;
        pst = pst->pstLower;
        free(pstKill);
    }
}

static void initializePStore(PKD *ppkd,MDL mdl,struct inInitializePStore *in) {
    *ppkd = new pkdContext(
        static_cast<mdl::mdlClass *>(mdl),in->nStore,in->nMinTotalStore,in->nMinEphemeral,in->nEphemeralBytes,
        in->nTreeBitsLo,in->nTreeBitsHi,
        in->iCacheSize,in->iWorkQueueSize,in->iCUDAQueueSize,in->fPeriod,
        in->nSpecies[FIO_SPECIES_DARK],in->nSpecies[FIO_SPECIES_SPH],in->nSpecies[FIO_SPECIES_STAR], in->nSpecies[FIO_SPECIES_BH],
        in->mMemoryModel,in->bLightCone,in->bLightConeParticles);
}

int pstInitializePStore(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inInitializePStore *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inInitializePStore));
    if (pstNotCore(pst)) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_INITIALIZEPSTORE,in,nIn);
        pstInitializePStore(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        if (plcl->pkd) delete plcl->pkd;
        initializePStore(&plcl->pkd,pst->mdl,in);
    }
    return 0;
}

int pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inReadFile *>(vin);
    auto pout = static_cast<int *>(vout);
    uint64_t nFileStart,nFileEnd,nFileTotal,nFileSplit;

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

    pClass = new PARTCLASS[PKD_MAX_CLASSES];
    assert(pClass!=NULL);

    n = pkd->particles.getClasses( PKD_MAX_CLASSES, pClass );
    rID = mdlReqService(pkd->mdl,id,PST_SWAPCLASSES,pClass,n*sizeof(PARTCLASS));
    mdlGetReply(pkd->mdl,rID,pClass,&n);
    n = n / sizeof(PARTCLASS);
    pkd->particles.setClasses( n, pClass, 0 );
    delete [] pClass;
}

int pstReadFile(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto plcl = pst->plcl;
    auto in = static_cast<struct inReadFile *>(vin);
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
        auto nParts = new int[pst->nLeaves];
        int i;
        uint64_t nStart;
        PKD pkd;
        MDL mdl;

        assert(nParts!=NULL);
        pstOneNodeReadInit(pst,in,sizeof(*in),nParts,pst->nLeaves * sizeof(nParts));
        pkd = plcl->pkd;
        mdl = pkd->mdl;

        fio = fioLoad(in+1,in->dOmega0,in->dOmegab);
        assert(fio!=NULL);

        nStart = nNodeStart + nParts[0];
        for (i=1; i<pst->nLeaves; ++i) {
            int id = mdlSelf(mdl) + i;
            int inswap;
            /*
             * Read particles into the local storage.
             */
            assert(pkd->FreeStore() >= nParts[i]);
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
        delete [] nParts;
        fioClose(fio);
    }
    return 0;
}

int pstActiveOrder(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto plcl = pst->plcl;
    auto pnActive = static_cast<uint64_t *>(vout);
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
    auto plcl = pst->plcl;
    auto in = static_cast<struct inAddWriteStart *>(vin);

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
    auto in = static_cast<struct inCompressASCII *>(vin);
    auto out = static_cast<struct outCompressASCII *>(vout);
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
    auto in = static_cast<struct inWriteASCII *>(vin);

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
    auto in = static_cast<struct inRestore *>(vin);

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
        makeName(achInFile,in->achInFile,mdlSelf(pkd->mdl),"");
        pkdRestore(pkd,achInFile);
    }
    return 0;
}

int pstCheckpoint(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inWrite *>(vin);

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
        makeName(achOutFile,in->achOutFile,mdlSelf(pkd->mdl),"");
        pkdCheckpoint(pkd,achOutFile);
    }
    return 0;
}

int pstSendParticles(PST pst,void *vin,int nIn,void *vout,int nOut) {
    int iTo = *(int *)vin;
    pkdWriteViaNode(pst->plcl->pkd, iTo);
    return 0;
}

int pstSendArray(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inSendArray *>(vin);
    pkdSendArray(pst->plcl->pkd, in->iTo, in->field, in->iUnitSize, in->dvFac, in->bMarked);
    return 0;
}

int pstWrite(PST pst,void *vin,int nIn,void *vout,int nOut) {
    char achOutFile[PST_FILENAME_SIZE];
    auto in = static_cast<struct inWrite *>(vin);
    FIO fio;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(struct inWrite));

    if (pst->nLeaves > 1) {
        int nProcessors = in->nProcessors;
        /* If we are the writer (nProcessors==1) or a producer (==1) keep requesting */
        if (nProcessors<=1) {
            pstWrite(pst->pstLower,in,nIn,NULL,0);
        }
        /* Split writing tasks between child nodes */
        else {
            int nLeft = nProcessors / 2;
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
        if (in->nProcessors!=0) {
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
                                             in->bStandard, pst->plcl->pkd->particles.present(PKD_FIELD::oNewSph) ? in->dTime : in->dExp,
                                             in->nGas, in->nDark, in->nStar, plcl->nWriteStart);
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

            pkdWriteHeaderFIO(plcl->pkd, fio, 1./sqrt(in->dvFac), in->dTime,
                              in->nDark, in->nGas, in->nStar, in->nBH,
                              in->dBoxSize, in->HubbleParam, in->nProcessors, in->units);
            pkdWriteFIO(plcl->pkd,fio,in->dvFac,in->dTuFac,&in->bnd);
            for (i=in->iLower+1; i<in->iUpper; ++i ) {
                int rID = mdlReqService(pst->mdl,i,PST_SENDPARTICLES,&pst->idSelf,sizeof(pst->idSelf));
                pkdWriteFromNode(plcl->pkd,i,fio,in->dvFac,in->dTuFac,&in->bnd);
                mdlGetReply(pst->mdl,rID,NULL,NULL);
            }
            fioClose(fio);
        }
    }

    return 0;
}

int pstSetSoft(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inSetSoft *>(vin);

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

int pstBuildTree(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    auto in = static_cast<struct inBuildTree *>(vin);
    uint32_t uRoot = in->uRoot;
    auto pTop = static_cast<KDN *>(vout);
    KDN *pCell1, *pCell2;
    double minside;
    int nOutUpper;

    /* We need to save our cells so we can update them later */
    if (pst->nLeaves > 1) {
        pCell1 = pkd->Node(pTop,1);
        pCell2 = pkd->Node(pTop,pst->nLower*2);

        /* We will accumulate the top tree here */
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,vin,nIn);
        nOut = pstBuildTree(pst->pstLower,vin,nIn,pCell1,(pst->nLower*2-1) * pkd->NodeSize());
        assert(nOut == (pst->nLower*2-1) * pkd->NodeSize());
        mdlGetReply(pst->mdl,rID,pCell2,&nOutUpper);
        assert(nOutUpper == (pst->nUpper*2-1) * pkd->NodeSize());
        nOut += nOutUpper + pkd->NodeSize();

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
        KDN *pRoot = pkd->TreeNode(uRoot);
        pkd->TreeAlignNode();
        pkdTreeBuild(plcl->pkd,in->nBucket,in->nGroup,in->uRoot,in->utRoot,in->ddHonHLimit);
        pkdCopyNode(pkd,pTop,pRoot);
        /* Get our cell ready */
        pTop->bTopTree = 1;
        pTop->bGroup = 0;
        pTop->bRemote = 1;
        pTop->pUpper = pTop->pLower = pst->idSelf;
        /* iLower is valid = ROOT */
        nOut = pkd->NodeSize();
    }
    return nOut;
}

int pstPhysicalSoft(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inPhysicalSoft *>(vin);

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
    auto in = static_cast<struct inHopLink *>(vin);
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
    auto in = static_cast<struct inHopLink *>(vin);
    auto out = static_cast<struct outHopJoin *>(vout);
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
        pkdGravWalkHop(plcl->pkd,in->dTime,in->nGroup,in->dTheta,&dFlop,&dPartSum,&dCellSum);
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
    auto in = static_cast<struct inGroupStats *>(vin);
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

#ifdef BLACKHOLES
int pstPlaceBHSeed(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inPlaceBHSeed *>(vin);
    auto out = static_cast<struct outPlaceBHSeed *>(vout);
    struct outPlaceBHSeed outUpper;
    mdlassert(pst->mdl,nIn == sizeof(struct inPlaceBHSeed));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_BH_PLACESEED,in,nIn);
        pstPlaceBHSeed(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
        assert(nOut==sizeof(struct outPlaceBHSeed));
        out->nBHs += outUpper.nBHs;
    }
    else {
        LCL *plcl = pst->plcl;
        out->nBHs = pkdPlaceBHSeed(plcl->pkd,in->dTime, in->dScaleFactor, in->uRungMax, in->dDenMin,
                                   in->dBHMhaloMin, in->dTau, in->dInitialH, in->dBHSeedMass);
    }
    return sizeof(struct outPlaceBHSeed);

}
int pstBHInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inPlaceBHSeed *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inPlaceBHSeed));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_BH_INIT,in,nIn);
        pstBHInit(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdBHInit(plcl->pkd,in->uRungMax);
    }
    return 0;

}
int pstRepositionBH(PST pst,void *vin,int nIn,void *vout,int nOut) {
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_BH_REPOSITION,NULL,0);
        pstRepositionBH(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdRepositionBH(plcl->pkd);
    }
    return 0;

}
#endif

int pstMoveDeletedParticles(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outGetNParts *>(vout);
    struct outGetNParts outUpper;

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MOVEDELETED,vin,nIn);
        pstMoveDeletedParticles(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
        out->n += outUpper.n;
        out->nGas += outUpper.nGas;
        out->nDark += outUpper.nDark;
        out->nStar += outUpper.nStar;
        out->nBH += outUpper.nBH;
    }
    else {
        pkdMoveDeletedParticles(plcl->pkd, &out->n, &out->nGas, &out->nDark, &out->nStar, &out->nBH);
    }
    return sizeof(struct outGetNParts);
}

int pstSmooth(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inSmooth *>(vin);

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

int pstReSmooth(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inSmooth *>(vin);
    auto out = static_cast<struct outSmooth *>(vout);
    struct outSmooth outUpper;

//    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTH,in,nIn);
        pstReSmooth(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
        assert(nOut==sizeof(struct outSmooth));
        out->nSmoothed += outUpper.nSmoothed;
    }
    else {
        LCL *plcl = pst->plcl;
        SMX smx;

        smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
                     in->bPeriodic,in->bSymmetric,in->iSmoothType);
        out->nSmoothed = smReSmooth(smx,&in->smf, in->iSmoothType);
        smFinish(smx,&in->smf);
    }
    return sizeof(struct outSmooth);
}

#ifdef OPTIM_SMOOTH_NODE
int pstReSmoothNode(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inSmooth *>(vin);
    auto out = static_cast<struct outSmooth *>(vout);
    struct outSmooth outUpper;

//    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTHNODE,in,nIn);
        pstReSmoothNode(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
        assert(nOut==sizeof(struct outSmooth));
        out->nSmoothed += outUpper.nSmoothed;
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER) && defined(DEBUG_FLUX_INFO)
        struct outSmooth tmp = *out;
        pstCombStat(&out->sComputing,&tmp.sComputing);
        pstCombStat(&out->sWaiting,&tmp.sWaiting);
        pstCombStat(&out->sSynchronizing,&tmp.sSynchronizing);
        pstCombStat(&out->sPartNumAccess,&tmp.sPartNumAccess);
        pstCombStat(&out->sPartMissRatio,&tmp.sPartMissRatio);
        pstCombStat(&out->sCellNumAccess,&tmp.sCellNumAccess);
        pstCombStat(&out->sCellMissRatio,&tmp.sCellMissRatio);
#endif
    }
    else {
        LCL *plcl = pst->plcl;
        SMX smx;

        smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,
                     in->bPeriodic,in->bSymmetric,in->iSmoothType);
#ifdef BLACKHOLES
        if (in->iSmoothType == SMX_BH_MERGER)
            out->nSmoothed = smReSmoothBHNode(smx,&in->smf, in->iSmoothType);
        else
#endif
            out->nSmoothed = smReSmoothNode(smx,&in->smf, in->iSmoothType);
        smFinish(smx,&in->smf);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER) && defined(DEBUG_FLUX_INFO)
        if (out->nSmoothed) {
            out->sCellNumAccess.dSum = mdlNumAccess(pst->mdl,CID_CELL)/out->nSmoothed;
            out->sPartNumAccess.dSum = mdlNumAccess(pst->mdl,CID_PARTICLE)/out->nSmoothed;
        }
        else {
            out->sCellNumAccess.dSum = 0;
            out->sPartNumAccess.dSum = 0;
        }
        out->sCellMissRatio.dSum = 100.0*mdlMissRatio(pst->mdl,CID_CELL);      /* as a percentage */
        out->sPartMissRatio.dSum = 100.0*mdlMissRatio(pst->mdl,CID_PARTICLE);  /* as a percentage */
        out->sComputing.dSum     = mdlTimeComputing(pst->mdl);
        out->sWaiting.dSum       = mdlTimeWaiting(pst->mdl);
        out->sSynchronizing.dSum = mdlTimeSynchronizing(pst->mdl);
        pstInitStat(&out->sComputing,pst->idSelf);
        pstInitStat(&out->sWaiting,pst->idSelf);
        pstInitStat(&out->sSynchronizing,pst->idSelf);
        pstInitStat(&out->sPartNumAccess,pst->idSelf);
        pstInitStat(&out->sPartMissRatio,pst->idSelf);
        pstInitStat(&out->sCellNumAccess,pst->idSelf);
        pstInitStat(&out->sCellMissRatio,pst->idSelf);
#endif
    }
    return sizeof(struct outSmooth);
}
#endif

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
    auto in = static_cast<struct inGravity *>(vin);
    auto outr = static_cast<struct outGravityReduct *>(vout);
    struct outGravityReduct tmp;
    int i;


    mdlassert(pst->mdl,nIn == sizeof(struct inGravity));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GRAVITY,in,nIn);
        pstGravity(pst->pstLower,in,nIn,vout,nOut);
        /*
        ** Make a temporary copy of the reduct part of the out buffer before setting it as the
        ** reply buffer. The reduct part follows at the end of all the outGravityPerProc entries.
        */
        mdlGetReply(pst->mdl,rID,&tmp,NULL);
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
        for (i=in->ts.uRungLo; i<=IRUNGMAX; ++i) outr->nRung[i] += tmp.nRung[i];
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
        pkdGravAll(pkd,&in->kick,&in->lc,&in->ts,
                   in->dTime,in->nReps,in->bPeriodic,
                   in->bEwald,in->nGroup,in->iRoot1,in->iRoot2,in->dEwCut,in->dEwhCut,in->dTheta,&in->SPHoptions,
                   &outr->nActive,
                   &outr->sPart.dSum,&outr->sPartNumAccess.dSum,&outr->sPartMissRatio.dSum,
                   &outr->sCell.dSum,&outr->sCellNumAccess.dSum,&outr->sCellMissRatio.dSum,
                   &outr->sFlop.dSum,outr->nRung);
        outr->dFlopSingleCPU = 1e-9*pkd->dFlopSingleCPU;
        outr->dFlopDoubleCPU = 1e-9*pkd->dFlopDoubleCPU;
        outr->dFlopSingleGPU = 1e-9*pkd->dFlopSingleGPU;
        outr->dFlopDoubleGPU = 1e-9*pkd->dFlopDoubleGPU;
        outr->sLocal.dSum = plcl->pkd->Local();
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
            while (fgets(buffer,sizeof(buffer),fp)) {
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
    }
    return sizeof(struct outGravityReduct);
}

int pstCalcEandL(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outCalcEandL *>(vout);
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
        for (k=0; k<3; k++) out->L[k] = outLcl.L[k];
        for (k=0; k<3; k++) out->F[k] = outLcl.F[k];
        out->W += outLcl.W;
    }
    else {
        pkdCalcEandL(plcl->pkd,&out->T,&out->U,&out->Eth,out->L,out->F,&out->W);
    }
    return sizeof(struct outCalcEandL);
}


int pstDrift(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inDrift *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DRIFT,in,nIn);
        pstDrift(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdDrift(plcl->pkd,in->iRoot,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaUPred,in->bDoGas);
    }
    return 0;
}


int pstSetGlobalDt(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct outGetMinDt *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct outGetMinDt));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETGLOBALDT,in,nIn);
        pstSetGlobalDt(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdSetGlobalDt(plcl->pkd, in->uMinDt);
    }
    return 0;
}


#ifdef OPTIM_REORDER_IN_NODES
int pstReorderWithinNodes(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_REORDERINNODES,NULL,0);
        pstReorderWithinNodes(pst->pstLower,NULL,0,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdReorderWithinNodes(plcl->pkd);
    }
    return 0;
}
#endif

int pstResetFluxes(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inDrift *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_RESETFLUXES,in,nIn);
        pstResetFluxes(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdResetFluxes(plcl->pkd,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaUPred);
    }
    return 0;
}

#ifdef DEBUG_CACHED_FLUXES
int pstFluxStats(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inFluxStats *>(vin);
    auto out = static_cast<struct outFluxStats *>(vout);
    struct outFluxStats outUpper;

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_FLUXSTATS,in,nIn);
        pstFluxStats(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,&nOut);
        assert(nOut==sizeof(struct outFluxStats));
        out->nAvoided += outUpper.nAvoided;
        out->nComputed += outUpper.nComputed;
    }
    else {
        LCL *plcl = pst->plcl;

        int avoided = 0;
        int computed = 0;
        pkdFluxStats(plcl->pkd, &computed, &avoided);
        out->nAvoided = avoided;
        out->nComputed = computed;
    }
    return sizeof(struct outFluxStats);
}
#endif

int pstEndTimestepIntegration(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inEndTimestep *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inEndTimestep));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COMPUTEPRIMVARS,in,nIn);
        pstEndTimestepIntegration(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdEndTimestepIntegration(plcl->pkd, *in);
    }
    return 0;
}

int pstWakeParticles(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inDrift *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_WAKEPARTICLES,in,nIn);
        pstWakeParticles(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdWakeParticles(plcl->pkd,in->iRoot, in->dTime, in->dDelta);
    }
    return 0;
}

#ifdef GRACKLE
int pstGrackleInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    auto in = static_cast<struct inGrackleInit *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inGrackleInit));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GRACKLEINIT,in,nIn);
        pstGrackleInit(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdGrackleInit(plcl->pkd, in->bComove, in->dScaleFactor, in->achCoolingTable,
                       in->units);
    }
    return 0;
}
#endif

#ifdef COOLING
int pstCoolingInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inCoolInit *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inCoolInit));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COOLINGINIT,in,nIn);
        pstCoolingInit(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkd_cooling_init_backend(plcl->pkd, in->in_cooling_data,
                                 in->Redshifts,
                                 in->nH,
                                 in->Temp,
                                 in->HeFrac,
                                 in->Therm,
                                 in->SolarAbundances,
                                 in->SolarAbundances_inv);
    }
    return 0;
}

int pstCoolingUpdate(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inCoolUpdate *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inCoolUpdate));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COOLINGUPDATE,in,nIn);
        pstCoolingUpdate(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkd_cooling_update(plcl->pkd, in);
    }
    return 0;
}
int pstCoolingUpdateZ(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<float *>(vin);

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COOLINGUPDATEZ,in,nIn);
        pstCoolingUpdateZ(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        plcl->pkd->cooling->dz = *in;
    }
    return 0;
}

int pstCoolingHydReion(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_COOLINGHYDREION,NULL,0);
        pstCoolingHydReion(pst->pstLower,NULL,0,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        cooling_Hydrogen_reionization(plcl->pkd);
    }
    return 0;
}
#endif

int pstChemCompInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inChemCompInit *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inChemCompInit));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CHEMCOMPINIT,in,nIn);
        pstChemCompInit(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdChemCompInit(plcl->pkd, *in);
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
        mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkd->ParticleBase(),pkd->ParticleSize(),
                   pkd->Local());

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
    auto in = static_cast<struct inKick *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_KICK,in,nIn);
        pstKick(pst->pstLower,in,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdKick(plcl->pkd,in->dTime,in->dDelta,in->bDoGas,in->dDeltaVPred,in->dDeltaU,in->dDeltaUPred,in->uRungLo,in->uRungHi);
    }
    return 0;
}

int pstKickTree(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inKickTree *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_KICKTREE,in,nIn);
        pstKickTree(pst->pstLower,in,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdKickTree(plcl->pkd,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaU,in->dDeltaUPred,in->iRoot);
    }
    return 0;
}

int pstSetTotal(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outSetTotal *>(vout);
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
        pst->nTotal = plcl->pkd->Local();
        out->nTotal = pst->nTotal;
    }
    return sizeof(struct outSetTotal);
}


int pstSetWriteStart(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inSetWriteStart *>(vin);
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

int pstAccelStep(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inAccelStep *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inAccelStep));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_ACCELSTEP,vin,nIn);
        pstAccelStep(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdAccelStep(plcl->pkd,in->uRungLo,in->uRungHi,
                     in->dDelta,in->iMaxRung,
                     in->dEta,in->dVelFac,in->dAccFac,
                     in->bDoGravity,in->bEpsAcc,in->dhMinOverSoft);
    }
    return 0;
}

int pstSphStep(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inSphStep *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inSphStep));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SPHSTEP,vin,nIn);
        pstSphStep(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdSphStep(plcl->pkd,in->uRungLo,in->uRungHi,in->dDelta,in->iMaxRung,in->dEta,in->dAccFac,in->dEtaUDot);
    }
    return 0;
}


int pstDensityStep(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inDensityStep *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_DENSITYSTEP,vin,nIn);
        pstDensityStep(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdDensityStep(plcl->pkd,in->uRungLo,in->uRungHi,in->iMaxRung,in->dDelta,in->dEta,in->dRhoFac);
    }
    return 0;
}

int pstCorrectEnergy(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inCorrectEnergy *>(vin);

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
    auto in = static_cast<struct inUpdateRung *>(vin);
    auto out = static_cast<struct outUpdateRung *>(vout);
    int i;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_UPDATERUNG,vin,nIn);
        pstUpdateRung(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outTemp,NULL);
        for (i=0; i<in->uMaxRung; ++i) {
            out->nRungCount[i] += outTemp.nRungCount[i];
        }
    }
    else {
        pkdUpdateRung(plcl->pkd,in->uRungLo,in->uRungHi,
                      in->uMinRung,in->uMaxRung,out->nRungCount);
    }
    return sizeof(*out);
}

int pstColNParts(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outColNParts *>(vout);
    struct outColNParts *outUp = out + pst->idUpper-pst->idSelf;

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
    auto in = static_cast<uint64_t *>(vin);

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

int pstGetNParts(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto out = static_cast<struct outGetNParts *>(vout);

    if (pst->nLeaves > 1) {
        struct outGetNParts outtmp;
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GETNPARTS,vin,nIn);
        pstGetNParts(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,(void *) &outtmp,NULL);

        out->n += outtmp.n;
        out->nGas += outtmp.nGas;
        out->nDark += outtmp.nDark;
        out->nStar += outtmp.nStar;
        out->nBH += outtmp.nBH;
        if (outtmp.nMaxOrder > out->nMaxOrder) out->nMaxOrder = outtmp.nMaxOrder;
    }
    else {
        pkdGetNParts(pst->plcl->pkd, out);
    }
    return sizeof(struct outGetNParts);
}

int pstSetNParts(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inSetNParts *>(vin);

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETNPARTS,vin,nIn);
        pstSetNParts(pst->pstLower, vin, nIn, NULL, 0);
        mdlGetReply(pst->mdl, rID, NULL, NULL);
    }
    else {
        pkdSetNParts(pst->plcl->pkd, in->nGas, in->nDark, in->nStar, in->nBH);
    }
    return 0;
}

int pstNewFof(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inNewFof *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inNewFof));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_NEW_FOF,in,nIn);
        pstNewFof(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdNewFof(plcl->pkd,in->dTau2,in->nMinMembers,in->bPeriodic,in->nReplicas,in->nBucket);
    }
    return 0;
}

int pstFofPhases(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto out = static_cast<struct outFofPhases *>(vout);
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
    auto in = static_cast<struct inGetFFTMaxSizes *>(vin);
    auto out = static_cast<struct outGetFFTMaxSizes *>(vout);
    struct outGetFFTMaxSizes outUp;

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

int pstMemStatus(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outMemStatus *>(vout);
    struct outMemStatus *outUp = out + pst->idUpper-pst->idSelf;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEMSTATUS,vin,nIn);
        pstMemStatus(pst->pstLower,vin,nIn,out,nOut);
        mdlGetReply(pst->mdl,rID,outUp,NULL);
    }
    else {
#ifdef __linux__
        FILE *fp;
        char buffer[512], *save, *f, *v;
        int i;

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
            while (fgets(buffer,sizeof(buffer),fp)) {
                f = strtok_r(buffer,":",&save);
                v = strtok_r(NULL," ",&save);
                if (strcmp(f,"MemFree")==0 || strcmp(f,"Buffers")==0 || strcmp(f,"Cached")==0) {
                    out->freeMemory += atol(v) / 1024;
                }
            }
            fclose(fp);
        }
#endif
        out->nBytesTree = plcl->pkd->TreeMemory();
        out->nBytesCl   = pkdClMemory(plcl->pkd);
        out->nBytesIlp  = pkdIlpMemory(plcl->pkd);
        out->nBytesIlc  = pkdIlcMemory(plcl->pkd);

    }
    return pst->nLeaves*sizeof(struct outMemStatus);
}


int pstGetClasses(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<PARTCLASS *>(vout);
    int nUp;
    int n, i, j;

    mdlassert(pst->mdl,nIn==0);
    if (pst->nLeaves > 1) {
        auto outUp = new PARTCLASS[PKD_MAX_CLASSES];
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
        delete [] outUp;
    }
    else {
        n = plcl->pkd->particles.getClasses(PKD_MAX_CLASSES,out);
    }
    return n * sizeof(PARTCLASS);
}

int pstSetClasses(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<PARTCLASS *>(vin);
    int n;

    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETCLASSES,vin,nIn);
        pstSetClasses(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        n = nIn / sizeof(PARTCLASS);
        mdlassert(pst->mdl,n*sizeof(PARTCLASS)==nIn);
        plcl->pkd->particles.setClasses(n,in,1);
    }
    return 0;
}

/*
 * Routine to swap the class table .  Note that this does not walk
 * the pst but simply set's the given table and returns the old table.
 */
int pstSwapClasses(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl;
    auto in = static_cast<PARTCLASS *>(vin);
    auto out = static_cast<PARTCLASS *>(vout);
    int n;
    PST lpst;

    lpst = pst;
    while (lpst->nLeaves > 1)
        lpst = lpst->pstLower;
    plcl = lpst->plcl;

    n = plcl->pkd->particles.getClasses( PKD_MAX_CLASSES, out );
    nOut = n * sizeof(PARTCLASS);

    n = nIn / sizeof(PARTCLASS);
    mdlassert(pst->mdl,n*sizeof(PARTCLASS) == nIn);
    plcl->pkd->particles.setClasses( n, in, 0 );
    return nOut;
}

int pstProfile(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inProfile *>(vin);
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
    auto in = static_cast<struct inCalcDistance *>(vin);

    assert( nIn==sizeof(struct inCalcDistance) );
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCDISTANCE,vin,nIn);
        pstCalcDistance(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdCalcDistance(plcl->pkd,in->dCenter,in->bPeriodic);
    }
    return 0;
}

int pstCalcCOM(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inCalcCOM *>(vin);
    auto out = static_cast<struct outCalcCOM *>(vout);
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
        for (i=0; i<3; i++ ) {
            out->com[i] += outUpper.com[i];
            out->vcm[i] += outUpper.vcm[i];
            out->L[i] += outUpper.L[i];
        }
    }
    else {
        pkdCalcCOM(plcl->pkd,in->dCenter,in->dRadius,in->bPeriodic,
                   out->com, out->vcm, out->L, &out->M, &out->N);
    }
    return sizeof(struct outCalcCOM);
}

int pstCalcMtot(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    //auto in = static_cast<struct inCalcMtot *>(vin);
    auto out = static_cast<struct outCalcMtot *>(vout);
    struct outCalcMtot outUpper;

    assert( nIn==sizeof(struct inCalcMtot) );
    assert( vout != NULL );
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_CALCMTOT,vin,nIn);
        pstCalcMtot(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,NULL);
        out->N += outUpper.N;
        out->M += outUpper.M;
    }
    else {
        pkdCalcMtot(plcl->pkd,&out->M, &out->N);
    }
    return sizeof(struct outCalcMtot);
}


int pstSetSPHoptions(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inSetSPHoptions *>(vin);

    assert( nIn==sizeof(struct inSetSPHoptions) );
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_SETSPHOPTIONS,vin,nIn);
        pstSetSPHoptions(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        copySPHOptions(&in->SPHoptions, &plcl->pkd->SPHoptions);
    }
    return 0;
}

int pstTreeUpdateFlagBounds(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    auto in = static_cast<struct inTreeUpdateFlagBounds *>(vin);
    uint32_t uRoot = in->uRoot;
    auto pTop = static_cast<KDN *>(vout);
    KDN *pCell1, *pCell2;
    double minside;
    int nOutUpper;

    /* We need to save our cells so we can update them later */
    if (pst->nLeaves > 1) {
        pCell1 = pkd->Node(pTop,1);
        pCell2 = pkd->Node(pTop,pst->nLower*2);

        /* We will accumulate the top tree here */
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_TREEUPDATEFLAGBOUNDS,vin,nIn);
        nOut = pstTreeUpdateFlagBounds(pst->pstLower,vin,nIn,pCell1,(pst->nLower*2-1) * pkd->NodeSize());
        assert(nOut == (pst->nLower*2-1) * pkd->NodeSize());
        mdlGetReply(pst->mdl,rID,pCell2,&nOutUpper);
        assert(nOutUpper == (pst->nUpper*2-1) * pkd->NodeSize());
        nOut += nOutUpper + pkd->NodeSize();

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
        KDN *pRoot = pkd->TreeNode(uRoot);
        pkd->TreeAlignNode();
        pkdTreeUpdateFlagBounds(plcl->pkd,uRoot,&in->SPHoptions);
        pkdCopyNode(pkd,pTop,pRoot);
        /* Get our cell ready */
        pTop->bTopTree = 1;
        pTop->bGroup = 0;
        pTop->bRemote = 1;
        pTop->pUpper = pTop->pLower = pst->idSelf;
        /* iLower is valid = ROOT */
        nOut = pkd->NodeSize();
    }
    return nOut;
}

int pstCountDistance(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inCountDistance *>(vin);
    auto out = static_cast<struct outCountDistance *>(vout);
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

#ifdef MDL_FFTW
int pstGridCreateFFT(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inGridCreateFFT *>(vin);
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
    auto out = static_cast<struct outTotalMass *>(vout);
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


int pstGetMinDt(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outGetMinDt *>(vout);
    struct outGetMinDt outUpper;



    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GETMINDT,vin,nIn);
        pstGetMinDt(pst->pstLower,vin,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,NULL);
        assert(nOut==sizeof(struct outGetMinDt));
        if (out->uMinDt < outUpper.uMinDt) out->uMinDt = outUpper.uMinDt;
    }
    else {
        out->uMinDt = pkdGetMinDt(plcl->pkd);
    }
    return sizeof(struct outGetMinDt);
}


int pstLightConeOpen(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inLightConeOpen *>(vin);
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
    auto in = static_cast<struct inLightConeClose *>(vin);
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
    auto in = static_cast<struct inLightConeVel *>(vin);

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LIGHTCONEVEL,in,nIn);
        pstLightConeVel(pst->pstLower,in,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdLightConeVel(plcl->pkd,in->dBoxSize);
    }
    return 0;
}

int pstGetParticles(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outGetParticles *>(vout);
    auto ID = static_cast<uint64_t *>(vin);
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
