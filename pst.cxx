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
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
    #include <malloc.h>
#endif
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <cinttypes>
#ifdef __linux__
    #include <unistd.h>
#endif
#ifdef HAVE_SYS_STAT_H
    #include <sys/stat.h>
#endif
#ifdef USE_ITT
    #include "ittnotify.h"
#endif
#include "fmt/format.h"
#include "mdl.h"
#include "pst.h"
#include "pkd.h"
#include "smooth/smooth.h"
#include "group/hop.h"
#include "group/fof.h"
#include "group/group.h"
#include "group/groupstats.h"
#include "io/output.h"
#include "core/swapall.h"

#include "cooling/cooling.h"
#ifdef GRACKLE
    #include "cooling_grackle/cooling_grackle.h"
#endif
#ifdef BLACKHOLES
    #include "blackhole/merger.h"
    #include "blackhole/evolve.h"
    #include "blackhole/seed.h"
    #include "blackhole/init.h"
    #include "blackhole/drift.h"
    #include "blackhole/accretion.h"
#endif
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif

#define GET_PARTICLES_MAX 20 /* We have a nested loop, so don't increase this */

void pstAddServices(PST pst,mdl::mdlClass *mdl) {
    int nThreads;

    nThreads = mdl->Threads();

    mdl->AddService(PST_READFILE,pst,(fcnService_t *)pstReadFile,
                    sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PATH_MAX),0,
                    Legacy::ReadFile);
    mdl->AddService(PST_COMPRESSASCII,pst,(fcnService_t *)pstCompressASCII,
                    sizeof(struct inCompressASCII),sizeof(struct outCompressASCII),
                    Legacy::CompressASCII);
    mdl->AddService(PST_WRITEASCII,pst,(fcnService_t *)pstWriteASCII,
                    sizeof(struct inWriteASCII),0,
                    Legacy::WriteASCII);
    mdl->AddService(PST_WRITE,pst,(fcnService_t *)pstWrite,
                    sizeof(struct inWrite),0,
                    Legacy::Write);
    mdl->AddService(PST_SENDPARTICLES,pst,(fcnService_t *)pstSendParticles,
                    sizeof(int),0,
                    Legacy::SendParticles);
    mdl->AddService(PST_SENDARRAY,pst,(fcnService_t *)pstSendArray,
                    sizeof(struct inSendArray),0,
                    Legacy::SendArray);
    mdl->AddService(PST_CHECKPOINT,pst,(fcnService_t *)pstCheckpoint,
                    sizeof(struct inWrite),0,
                    Legacy::Checkpoint);
    mdl->AddService(PST_OUTPUT,pst,(fcnService_t *)pstOutput,
                    sizeof(struct inOutput),0,
                    Legacy::Output);
    mdl->AddService(PST_OUTPUT_SEND,pst,(fcnService_t *)pstOutputSend,
                    sizeof(struct inOutputSend),0,
                    Legacy::OutputSend);
    /*
    ** Calculate the number of levels in the top tree and use it to
    ** define the size of the messages.
    */
    mdl->AddService(PST_BUILDTREE,pst,(fcnService_t *)pstBuildTree,
                    sizeof(struct inBuildTree),
                    (nThreads==1?1:2*nThreads-1)*pkdContext::MaxNodeSize(),
                    Legacy::BuildTree);
    mdl->AddService(PST_TREEINITMARKED,pst,(fcnService_t *)pstTreeInitMarked,
                    0,0,
                    Legacy::TreeInitMarked);
    mdl->AddService(PST_HOP_LINK,pst,(fcnService_t *)pstHopLink,
                    sizeof(struct inHopLink),sizeof(uint64_t),
                    Legacy::HopLink);
    mdl->AddService(PST_HOP_JOIN,pst,(fcnService_t *)pstHopJoin,
                    sizeof(struct inHopLink),sizeof(struct outHopJoin),
                    Legacy::HopJoin);
    mdl->AddService(PST_HOP_FINISH_UP,pst,(fcnService_t *)pstHopFinishUp,
                    sizeof(struct inHopFinishUp),sizeof(uint64_t),
                    Legacy::HopFinishUp);
    mdl->AddService(PST_HOP_TREE_BUILD,pst,(fcnService_t *)pstHopTreeBuild,
                    sizeof(struct inHopTreeBuild),0,
                    Legacy::HopTreeBuild);
    mdl->AddService(PST_HOP_GRAVITY,pst,(fcnService_t *)pstHopGravity,
                    sizeof(struct inHopGravity),0,
                    Legacy::HopGravity);
    mdl->AddService(PST_HOP_UNBIND,pst,(fcnService_t *)pstHopUnbind,
                    sizeof(struct inHopUnbind),sizeof(struct outHopUnbind),
                    Legacy::HopUnbind);
    mdl->AddService(PST_GROUP_RELOCATE,pst,(fcnService_t *)pstGroupRelocate,
                    0,0,
                    Legacy::GroupRelocate);
    mdl->AddService(PST_GROUP_STATS,pst,(fcnService_t *)pstGroupStats,
                    sizeof(struct inGroupStats),0,
                    Legacy::GroupStats);
    mdl->AddService(PST_SMOOTH,pst,(fcnService_t *)pstSmooth,
                    sizeof(struct inSmooth),0,
                    Legacy::Smooth);
    mdl->AddService(PST_GRAVITY,pst,(fcnService_t *)pstGravity,
                    sizeof(struct inGravity),
                    sizeof(struct outGravityReduct),
                    Legacy::Gravity);
    mdl->AddService(PST_CALCEANDL,pst,(fcnService_t *)pstCalcEandL,
                    0,sizeof(struct outCalcEandL),
                    Legacy::CalcEandL);
    mdl->AddService(PST_DRIFT,pst,(fcnService_t *)pstDrift,
                    sizeof(struct inDrift),0,
                    Legacy::Drift);
    mdl->AddService(PST_SETGLOBALDT,pst,
                    (fcnService_t *) pstSetGlobalDt,
                    sizeof(struct outGetMinDt),0,
                    Legacy::SetGlobalDt);
#ifdef GRACKLE
    mdl->AddService(PST_GRACKLEINIT,pst,
                    (fcnService_t *) pstGrackleInit,
                    sizeof(struct inGrackleInit),0,
                    Legacy::GrackleInit);
#endif
    mdl->AddService(PST_CHEMCOMPINIT,pst,
                    (fcnService_t *) pstChemCompInit,
                    sizeof(struct inChemCompInit),0,
                    Legacy::ChemCompInit);
#ifdef BLACKHOLES
    mdl->AddService(PST_BH_PLACESEED,pst,
                    (fcnService_t *) pstPlaceBHSeed,
                    sizeof(struct inPlaceBHSeed), sizeof(struct outPlaceBHSeed),
                    Legacy::PlaceBHSeed);
    mdl->AddService(PST_BH_INIT,pst,
                    (fcnService_t *) pstBHInit,
                    sizeof(struct inPlaceBHSeed), 0,
                    Legacy::BHInit);
    mdl->AddService(PST_BH_REPOSITION,pst,
                    (fcnService_t *) pstBHReposition,
                    0, 0,
                    Legacy::BHReposition);
    mdl->AddService(PST_BH_ACCRETION,pst,
                    (fcnService_t *) pstBHAccretion,
                    sizeof(struct inBHAccretion), 0,
                    Legacy::BHAccretion);
#endif
    mdl->AddService(PST_GETMINDT,pst,
                    (fcnService_t *) pstGetMinDt,
                    0, sizeof(struct outGetMinDt),
                    Legacy::GetMinDt);
    //
    mdl->AddService(PST_ROPARTICLECACHE,pst,(fcnService_t *)pstROParticleCache,
                    0,0,
                    Legacy::ROParticleCache);
    mdl->AddService(PST_PARTICLECACHEFINISH,pst,(fcnService_t *)pstParticleCacheFinish,
                    0,0,
                    Legacy::ParticleCacheFinish);
    mdl->AddService(PST_KICK,pst,(fcnService_t *)pstKick,
                    sizeof(struct inKick),0,
                    Legacy::Kick);
    mdl->AddService(PST_KICKTREE,pst,(fcnService_t *)pstKickTree,
                    sizeof(struct inKickTree),0,
                    Legacy::KickTree);
    mdl->AddService(PST_PHYSICALSOFT,pst,(fcnService_t *)pstPhysicalSoft,
                    sizeof(struct inPhysicalSoft),0,
                    Legacy::PhysicalSoft);
    mdl->AddService(PST_SETTOTAL,pst,(fcnService_t *)pstSetTotal,
                    0,sizeof(struct outSetTotal),
                    Legacy::SetTotal);
    mdl->AddService(PST_SETWRITESTART,pst,(fcnService_t *)pstSetWriteStart,
                    sizeof(struct inSetWriteStart),0,
                    Legacy::SetWriteStart);
    mdl->AddService(PST_ADDWRITESTART,pst,(fcnService_t *)pstAddWriteStart,
                    sizeof(struct inAddWriteStart),0,
                    Legacy::AddWriteStart);
    mdl->AddService(PST_ONENODEREADINIT,pst,(fcnService_t *)pstOneNodeReadInit,
                    sizeof(struct inReadFile) + PST_MAX_FILES*(sizeof(fioSpeciesList)+PATH_MAX),
                    nThreads*sizeof(fioSpeciesList),
                    Legacy::OneNodeReadInit);
    mdl->AddService(PST_ACTIVEORDER,pst,(fcnService_t *)pstActiveOrder,
                    0,sizeof(uint64_t),
                    Legacy::ActiveOrder);
    mdl->AddService(PST_DENSITYSTEP,pst,(fcnService_t *)pstDensityStep,
                    sizeof(struct inDensityStep),0,
                    Legacy::DensityStep);
    mdl->AddService(PST_CORRECTENERGY,pst,(fcnService_t *)pstCorrectEnergy,
                    sizeof(struct inCorrectEnergy),0,
                    Legacy::CorrectEnergy);
    mdl->AddService(PST_ACCELSTEP,pst,(fcnService_t *)pstAccelStep,
                    sizeof(struct inAccelStep), 0,
                    Legacy::AccelStep);
#ifdef STAR_FORMATION
    mdl->AddService(PST_STARFORM,pst,
                    (fcnService_t *) pstStarForm,
                    sizeof(struct inStarForm),sizeof(struct outStarForm),
                    Legacy::StarForm);
#endif
#if defined(STAR_FORMATION) || defined(FEEDBACK)
    mdl->AddService(PST_STARFORMINIT,pst,
                    (fcnService_t *) pstStarFormInit,
                    sizeof(struct inStarFormInit),sizeof(struct outStarForm),
                    Legacy::StarFormInit);
#endif
    mdl->AddService(PST_RESMOOTH,pst,(fcnService_t *) pstReSmooth,
                    sizeof(struct inSmooth),sizeof(struct outSmooth),
                    Legacy::ReSmooth);
#ifdef OPTIM_SMOOTH_NODE
    mdl->AddService(PST_RESMOOTHNODE,pst,(fcnService_t *) pstReSmoothNode,
                    sizeof(struct inSmooth),sizeof(struct outSmooth),
                    Legacy::ReSmoothNode);
#endif
#ifdef OPTIM_REORDER_IN_NODES
    mdl->AddService(PST_REORDERINNODES,pst,(fcnService_t *) pstReorderWithinNodes,
                    0,0,
                    Legacy::ReorderInNodes);
#endif
#ifdef STELLAR_EVOLUTION
    mdl->AddService(PST_STELLAREVOLUTIONINIT,pst,
                    (fcnService_t *) pstStellarEvolutionInit,
                    sizeof(struct inStellarEvolutionInit),0,
                    Legacy::StellarEvolutionInit);
#endif
    mdl->AddService(PST_COLNPARTS,pst,(fcnService_t *)pstColNParts,
                    0,nThreads*sizeof(struct outColNParts),
                    Legacy::ColNParts);
    mdl->AddService(PST_NEWORDER,pst,(fcnService_t *)pstNewOrder,
                    nThreads*sizeof(uint64_t),0,
                    Legacy::NewOrder);
    mdl->AddService(PST_SETNPARTS,pst,(fcnService_t *)pstSetNParts,
                    sizeof(struct inSetNParts),0,
                    Legacy::SetNParts);
    mdl->AddService(PST_NEW_FOF,pst,(fcnService_t *)pstNewFof,
                    sizeof(struct inNewFof),0,
                    Legacy::NewFof);
    mdl->AddService(PST_FOF_PHASES,pst,(fcnService_t *)pstFofPhases,
                    0,sizeof(struct outFofPhases),
                    Legacy::FofPhases);
    mdl->AddService(PST_FOF_FINISH_UP,pst,(fcnService_t *)pstFofFinishUp,
                    sizeof(struct inFofFinishUp),sizeof(uint64_t),
                    Legacy::FofFinishUp);
    mdl->AddService(PST_INITIALIZEPSTORE,pst,(fcnService_t *)pstInitializePStore,
                    sizeof(struct inInitializePStore),sizeof(struct outInitializePStore),
                    Legacy::InitializePStore);
#ifdef MDL_FFTW
    mdl->AddService(PST_GENERATEIC,pst,(fcnService_t *)pstGenerateIC,
                    sizeof(struct inGenerateIC),sizeof(struct outGenerateIC),
                    Legacy::GenerateIC);
    mdl->AddService(PLT_GENERATEIC,pst,(fcnService_t *)pltGenerateIC,
                    sizeof(struct inGenerateICthread),sizeof(struct outGenerateIC),
                    Legacy::LocalGenerateIC);
    mdl->AddService(PLT_MOVEIC,pst,(fcnService_t *)pltMoveIC,
                    sizeof(struct inMoveIC),0,
                    Legacy::MoveIC);
    mdl->AddService(PST_MOVEIC,pst,(fcnService_t *)pstMoveIC,
                    sizeof(struct inGenerateIC),0,
                    Legacy::LocalMoveIC);
#endif
    mdl->AddService(PST_MEMSTATUS,pst,(fcnService_t *)pstMemStatus,
                    0,nThreads*sizeof(struct outMemStatus),
                    Legacy::MemStatus);
    mdl->AddService(PST_GETCLASSES,pst,(fcnService_t *)pstGetClasses,
                    0, PKD_MAX_CLASSES*sizeof(PARTCLASS),
                    Legacy::GetClasses);
    mdl->AddService(PST_SETCLASSES,pst,(fcnService_t *)pstSetClasses,
                    PKD_MAX_CLASSES*sizeof(PARTCLASS), 0,
                    Legacy::SetClasses);
    mdl->AddService(PST_SWAPCLASSES,pst,(fcnService_t *)pstSwapClasses,
                    PKD_MAX_CLASSES*sizeof(PARTCLASS),
                    PKD_MAX_CLASSES*sizeof(PARTCLASS),
                    Legacy::SwapClasses);
    mdl->AddService(PST_PROFILE,pst,(fcnService_t *)pstProfile,
                    sizeof(struct inProfile), 0,
                    Legacy::Profile);
    mdl->AddService(PST_CALCDISTANCE,pst,(fcnService_t *)pstCalcDistance,
                    sizeof(struct inCalcDistance), 0,
                    Legacy::CalcDistance);
    mdl->AddService(PST_CALCCOM,pst,(fcnService_t *)pstCalcCOM,
                    sizeof(struct inCalcCOM), sizeof(struct outCalcCOM),
                    Legacy::CalcCOM);
    mdl->AddService(PST_CALCMTOT,pst,(fcnService_t *)pstCalcMtot,
                    sizeof(struct inCalcMtot), sizeof(struct outCalcMtot),
                    Legacy::CalcMtot);
    mdl->AddService(PST_SETSPHOPTIONS,pst,(fcnService_t *)pstSetSPHoptions,
                    sizeof(struct inSetSPHoptions), 0,
                    Legacy::SetSPHoptions);
    mdl->AddService(PST_RESETCOM,pst,(fcnService_t *)pstResetCOM,
                    sizeof(struct inResetCOM), 0,
                    Legacy::ResetCOM);
    mdl->AddService(PST_INITIALIZEEOS,pst,(fcnService_t *)pstInitializeEOS,
                    0, 0,
                    Legacy::InitializeEOS);
    mdl->AddService(PST_UPDATEGASVALUES,pst,(fcnService_t *)pstUpdateGasValues,
                    sizeof(struct inUpdateGasValues), 0,
                    Legacy::UpdateGasValues);
    mdl->AddService(PST_TREEUPDATEFLAGBOUNDS,pst,(fcnService_t *)pstTreeUpdateFlagBounds,
                    sizeof(struct inTreeUpdateFlagBounds),
                    (nThreads==1?1:2*nThreads-1)*pkdContext::MaxNodeSize(),
                    Legacy::TreeUpdateFlagBounds);
#ifdef MDL_FFTW
    mdl->AddService(PST_GRID_CREATE_FFT,pst,(fcnService_t *)pstGridCreateFFT,
                    sizeof(struct inGridCreateFFT), 0,
                    Legacy::GridCreateFFT);
    mdl->AddService(PST_GRID_DELETE_FFT,pst,(fcnService_t *)pstGridDeleteFFT,
                    0, 0,
                    Legacy::GridDeleteFFT);
    mdl->AddService(PST_ADD_LINEAR_SIGNAL,pst,(fcnService_t *)pstAddLinearSignal,
                    sizeof(struct inAddLinearSignal), 0,
                    Legacy::AddLinearSignal);
    mdl->AddService(PST_ASSIGN_MASS,pst,(fcnService_t *)pstAssignMass,
                    sizeof(struct inAssignMass), 0,
                    Legacy::AssignMass);
    mdl->AddService(PST_DENSITY_CONTRAST,pst,(fcnService_t *)pstDensityContrast,
                    sizeof(struct inDensityContrast), 0,
                    Legacy::DensityContrast);
    mdl->AddService(PST_WINDOW_CORRECTION,pst,(fcnService_t *)pstWindowCorrection,
                    sizeof(struct inWindowCorrection), 0,
                    Legacy::WindowCorrection);
    mdl->AddService(PST_INTERLACE,pst,(fcnService_t *)pstInterlace,
                    sizeof(struct inInterlace), 0,
                    Legacy::Interlace);
    mdl->AddService(PST_GRID_BIN_K,pst,(fcnService_t *)pstGridBinK,
                    sizeof(struct inGridBinK), sizeof(struct outGridBinK),
                    Legacy::GridBinK);
    mdl->AddService(PST_BISPECTRUM_SELECT,pst,(fcnService_t *)pstBispectrumSelect,
                    sizeof(struct inBispectrumSelect), 0,
                    Legacy::BispectrumSelect);
    mdl->AddService(PST_BISPECTRUM_CALCULATE,pst,(fcnService_t *)pstBispectrumCalculate,
                    sizeof(struct inBispectrumCalculate), sizeof(double),
                    Legacy::BispectrumCalculate);
    mdl->AddService(PST_LINEARKICK, pst,(fcnService_t *)pstLinearKick,
                    sizeof(struct inLinearKick), 0,
                    Legacy::LinearKick);
    mdl->AddService(PST_SETLINGRID, pst,(fcnService_t *)pstSetLinGrid,
                    sizeof(struct inSetLinGrid), 0,
                    Legacy::SetLinGrid);
    mdl->AddService(PST_MEASURELINPK,pst,(fcnService_t *)pstMeasureLinPk,
                    sizeof(struct inMeasureLinPk), sizeof(struct outMeasureLinPk),
                    Legacy::MeasureLinPk);
#endif
    mdl->AddService(PST_TOTALMASS,pst,(fcnService_t *)pstTotalMass,
                    0, sizeof(struct outTotalMass),
                    Legacy::TotalMass);
    mdl->AddService(PST_LIGHTCONE_OPEN,pst,(fcnService_t *)pstLightConeOpen,
                    sizeof(struct inLightConeOpen), 0,
                    Legacy::LightConeOpen);
    mdl->AddService(PST_LIGHTCONE_CLOSE,pst,(fcnService_t *)pstLightConeClose,
                    sizeof(struct inLightConeClose), 0,
                    Legacy::LightConeClose);
    mdl->AddService(PST_LIGHTCONEVEL,pst,(fcnService_t *)pstLightConeVel,
                    sizeof(struct inLightConeVel),0,
                    Legacy::LightConeVel);
    mdl->AddService(PST_GET_PARTICLES,pst,(fcnService_t *)pstGetParticles,
                    sizeof(uint64_t)*GET_PARTICLES_MAX,
                    sizeof(struct outGetParticles)*GET_PARTICLES_MAX,
                    Legacy::GetParticles);
}

void pstInitialize(PST *ppst,mdl::mdlClass *mdl,LCL *plcl) {
    PST pst;

    pst = (PST)malloc(sizeof(struct pstContext));
    mdlassert(mdl,pst != NULL);
    *ppst = pst;
    pst->plcl = plcl;
    pst->mdl = mdl;
    pst->idSelf = mdl->Self();
    pst->pstLower = NULL;
    pst->idUpper = -1;  /* invalidate upper 'id' */
    pst->nLeaves = 1;
    pst->nLower = 0;
    pst->nUpper = 0;
    pst->iSplitDim = -1;
    pst->iLvl = 0;
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
        static_cast<mdl::mdlClass *>(mdl),in->nStore,in->nDecoratorStore,in->nBucket,in->nMinTotalStore,in->nMinEphemeral,in->nEphemeralBytes,
        in->nTreeBitsLo,in->nTreeBitsHi,
        in->iCacheSize,in->iCacheMaxInflight,in->iWorkQueueSize,
        in->mpi_recv_mode,in->mpi_recv_count,in->mpi_RDMA,in->fPeriod,
        in->nSpecies[FIO_SPECIES_DARK],in->nSpecies[FIO_SPECIES_SPH],in->nSpecies[FIO_SPECIES_STAR], in->nSpecies[FIO_SPECIES_BH],
        in->mMemoryModel,in->nIntegerFactor);
}

int pstInitializePStore(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inInitializePStore *>(vin);
    auto out = static_cast<struct outInitializePStore *>(vout);
    mdlassert(pst->mdl,nIn == sizeof(struct inInitializePStore));
    if (pstNotCore(pst)) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_INITIALIZEPSTORE,in,nIn);
        pstInitializePStore(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID,nOut,vout);
    }
    else {
        if (plcl->pkd) delete plcl->pkd;
        initializePStore(&plcl->pkd,pst->mdl,in);
        out->nSizeParticle = plcl->pkd->particles.ElementSize();
        out->nSizeDecorator = plcl->pkd->decorators.ElementSize();
        out->nSizeNode = plcl->pkd->tree.ElementSize();
    }
    return sizeof(*out);
}

int pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inReadFile *>(vin);
    auto pout = static_cast<fioSpeciesList *>(vout);
    fioSpeciesList nFileStart,nFileEnd,nFileTotal,nFileSplit;

    mdlassert(pst->mdl,nIn == sizeof(struct inReadFile));
    for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
        nFileStart[s] = in->nSpeciesStart[s];
        nFileEnd[s] = in->nSpeciesEnd[s];
        nFileTotal[s] = nFileEnd[s] - nFileStart[s] + 1;
    }

    if (pst->nLeaves > 1) {
        int rID;
        for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
            nFileSplit[s] = nFileStart[s] + pst->nLower*(nFileTotal[s]/pst->nLeaves);
            in->nSpeciesStart[s] = nFileSplit[s];
        }
        rID = pst->mdl->ReqService(pst->idUpper,PST_ONENODEREADINIT,in,nIn);
        for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
            in->nSpeciesStart[s] = nFileStart[s];
            in->nSpeciesEnd[s] = nFileSplit[s] - 1;
        }
        pstOneNodeReadInit(pst->pstLower,in,nIn,pout,nOut);
        for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
            in->nSpeciesEnd[s] = nFileEnd[s];
        }
        auto nBytesLower = pst->nLower * sizeof(*pout);
        pst->mdl->GetReply(rID,nOut-nBytesLower,pout+pst->nLower);
    }
    else {
        /*
        ** Determine the size of the local particle store.
        */
        for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
            pout[0][s] = nFileTotal[s]; /* Truncated: okay */
        }
    }
    return sizeof(*pout) * pst->nLeaves;
}

static void swap_classes(PKD pkd, int id) {
    auto pClass = std::make_unique<PARTCLASS[]>(PKD_MAX_CLASSES);
    int n = pkd->particles.getClasses( PKD_MAX_CLASSES, pClass.get() );
    int rID = pkd->mdl->ReqService(id,PST_SWAPCLASSES,pClass.get(),n*sizeof(PARTCLASS));
    n = pkd->mdl->GetReply(rID,PKD_MAX_CLASSES*sizeof(PARTCLASS),pClass.get());
    n = n / sizeof(PARTCLASS);
    pkd->particles.setClasses( n, pClass.get(), 0 );
}

int pstReadFile(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto plcl = pst->plcl;
    auto in = static_cast<struct inReadFile *>(vin);
    FIO fio;
    fioSpeciesList nNodeSpeciesStart,nNodeSpeciesEnd,nNodeSpeciesTotal,nNodeSpeciesSplit;
    int rID;

    mdlassert(pst->mdl,nIn >= sizeof(struct inReadFile));

    for ( auto s=FIO_SPECIES(0); s<=FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
        if (in->nSpeciesEnd[s] == 0) continue;
        // First element to read
        nNodeSpeciesStart[s] = in->nSpeciesStart[s];
        // Last element to read
        nNodeSpeciesEnd[s] = in->nSpeciesEnd[s];
        // Total count of read elements
        nNodeSpeciesTotal[s] = nNodeSpeciesEnd[s] - nNodeSpeciesStart[s] + 1;
    }
    if (pst->nLeaves > 1 && in->nProcessors > 1) {
        int nProcessors = in->nProcessors;
        int nProcUpper = pst->nUpper * nProcessors / pst->nLeaves;
        int nProcLower = nProcessors - nProcUpper;

        // Index of the first element in the second half of the split
        for ( auto s=FIO_SPECIES(0); s<=FIO_SPECIES_LAST; s=FIO_SPECIES(s+1))
            nNodeSpeciesSplit[s] = nNodeSpeciesStart[s] + pst->nLower*(nNodeSpeciesTotal[s]/pst->nLeaves);

        // Upper split
        if ( nProcessors > 1 ) {
            in->nProcessors = nProcUpper;
            for ( auto s=FIO_SPECIES(0); s<=FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
                if (in->nSpeciesEnd[s] == 0) continue;
                in->nSpeciesStart[s] = nNodeSpeciesSplit[s];
                in->nSpeciesEnd[s] = nNodeSpeciesEnd[s];
            }
            rID = pst->mdl->ReqService(pst->idUpper,PST_READFILE,in,nIn);
        }

        // Lower split
        in->nProcessors = nProcLower;
        for ( auto s=FIO_SPECIES(0); s<=FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
            if (in->nSpeciesEnd[s] == 0) continue;
            in->nSpeciesStart[s] = nNodeSpeciesStart[s];
            in->nSpeciesEnd[s] = nNodeSpeciesSplit[s] - 1;
        }
        pstReadFile(pst->pstLower,in,nIn,NULL,0);

        if ( nProcessors <= 1 ) {
            in->nProcessors = nProcUpper;
            for ( auto s=FIO_SPECIES(0); s<=FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
                if (in->nSpeciesEnd[s] == 0) continue;
                in->nSpeciesStart[s] = nNodeSpeciesSplit[s];
                in->nSpeciesEnd[s] = nNodeSpeciesEnd[s];
            }
            rID = pst->mdl->ReqService(pst->idUpper,PST_READFILE,in,nIn);
        }

        pst->mdl->GetReply(rID);
    }
    else {
        // Number of particles to be read per core per species
        auto nSpeciesParts = new fioSpeciesList[pst->nLeaves];
        int i;
        PKD pkd = plcl->pkd;
        auto mdl = pkd->mdl;

        pstOneNodeReadInit(pst,in,sizeof(*in),nSpeciesParts,pst->nLeaves * sizeof(fioSpeciesList));

        fio = fioLoad(in+1,in->dOmega0,in->dOmegab);
        assert(fio!=NULL);

        fioSpeciesList nRead, nStart, nBasePart;
        for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
            // We take care of the first nSpeciesParts later
            nStart[s] = in->nSpeciesStart[s] + nSpeciesParts[0][s];
        }
        for (i=1; i<pst->nLeaves; ++i) {
            int id = mdl->Self() + i;
            int inswap;
            for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
                if (in->nSpeciesEnd[s] == 0) continue;
                if (nSpeciesParts[i][s] == 0) continue;

                // Read particles into the local storage.
                assert(pkd->FreeStore() >= nSpeciesParts[i][s]);
                pkdReadFIO(pkd, fio, nStart[s], nSpeciesParts[i][s], in->dvFac,in->dTuFac, FIO_SPECIES(s));
                nStart[s] += nSpeciesParts[i][s];
            }

            // Now shove them over to the remote processor.
            swap_classes(pkd,id);
            inswap = mdl->Self();
            rID = pst->mdl->ReqService<ServiceSwapAll>(id,&inswap,sizeof(inswap));
            pkdSwapAll(pkd, id);
            pst->mdl->GetReply(rID);
        }
        // Finish with our own particles
        for ( auto s=FIO_SPECIES(0); s<FIO_SPECIES_LAST; s=FIO_SPECIES(s+1)) {
            if (in->nSpeciesEnd[s] == 0) continue;
            if (nSpeciesParts[0][s] == 0) continue;
            pkdReadFIO(pkd, fio, in->nSpeciesStart[s], nSpeciesParts[0][s], in->dvFac,in->dTuFac, FIO_SPECIES(s));
        }
        delete [] nSpeciesParts;
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_ACTIVEORDER,NULL,0);
        pstActiveOrder(pst->pstLower,NULL,0,vout,nOut);
        pst->mdl->GetReply(rID,nActiveLeaf);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_ADDWRITESTART,in,nIn);
        pstAddWriteStart(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        rID = pst->mdl->ReqService(pst->idUpper,PST_COMPRESSASCII,in,nIn);
        pstCompressASCII(pst->pstLower,in,nIn,vout,nOut);
        pst->mdl->GetReply(rID,outUp);
        /* Add the new offset from the left children to the right children */
        inAdd.nWriteStart = out->nBytes;
        rID = pst->mdl->ReqService(pst->idUpper,PST_ADDWRITESTART,&inAdd,sizeof(inAdd));
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_WRITEASCII,in,nIn);
        pstWriteASCII(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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

static std::string makeName(std::string_view inName, int iIndex) {
    return fmt::format(inName, fmt::arg("i", iIndex));
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
            int rID = pst->mdl->ReqService(pst->idUpper,PST_CHECKPOINT,in,nIn);
            in->nProcessors = nLower;
            pstCheckpoint(pst->pstLower,in,nIn,NULL,0);
            pst->mdl->GetReply(rID);
        }
        else { /* Serialize these processors now */
            pstCheckpoint(pst->pstLower,in,nIn,NULL,0);
            int rID = pst->mdl->ReqService(pst->idUpper,PST_CHECKPOINT,in,nIn);
            pst->mdl->GetReply(rID);
        }
    }
    else {
        PKD pkd = pst->plcl->pkd;
        auto achOutFile = makeName(in->achOutFile,pkd->mdl->Self());
        auto achOutDecoratorFile = makeName(in->achOutFile,pkd->mdl->Self()) + ".dec";
        // Need to add prefix as with previous "makeName" implementation
        pkdCheckpoint(pkd,achOutFile.c_str(), achOutDecoratorFile.c_str());
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
            assert(in->iLower == pst->mdl->Self());
            in->iLower = pst->idUpper;
            in->nProcessors -= nLeft;
            in->iIndex += nLeft;
            int rID = pst->mdl->ReqService(pst->idUpper,PST_WRITE,in,nIn);
            in->iLower = pst->mdl->Self();
            in->iUpper = pst->idUpper;
            in->iIndex -= nLeft;
            in->nProcessors = nLeft;
            pstWrite(pst->pstLower,in,nIn,NULL,0);
            pst->mdl->GetReply(rID);
        }
    }
    else {
        LCL *plcl = pst->plcl;
        if (in->nProcessors!=0) {
            auto achOutFile = makeName(in->achOutFile,in->iIndex);
            if (in->bHDF5) {
#ifdef USE_HDF5
                fio = fioHDF5Create(achOutFile.c_str(),in->mFlags);
#else
                fio = NULL; /* Should never happen */
#endif
            }
            else {
                if (achOutFile != makeName(in->achOutFile,in->iIndex+1)) {
                    fio = fioTipsyCreatePart(achOutFile.c_str(),0,in->mFlags&FIO_FLAG_CHECKPOINT,
                                             in->bStandard, pst->plcl->pkd->particles.present(PKD_FIELD::oNewSph) ? in->dTime : in->dExp,
                                             in->nGas, in->nDark, in->nStar, plcl->nWriteStart);
                }
                else {
                    fio = fioTipsyAppend(achOutFile.c_str(),in->mFlags&FIO_FLAG_CHECKPOINT,in->bStandard);
                    if (fio) {
                        fioSeek(fio,plcl->nWriteStart,FIO_SPECIES_ALL);
                    }
                }
            }
            if (fio==NULL) {
                fprintf(stderr,"ERROR: unable to create file for output\n");
                perror(achOutFile.c_str());
                mdlassert(pst->mdl,fio!=NULL);
            }

            fioSpeciesList nTotal, nThisFile{0};
            nTotal[FIO_SPECIES_DARK] = in->nDark;
            nTotal[FIO_SPECIES_SPH]  = in->nGas;
            nTotal[FIO_SPECIES_STAR] = in->nStar;
            nTotal[FIO_SPECIES_BH]   = in->nBH;

            pkdWriteFIO(plcl->pkd,fio,in->dvFac,in->dTuFac,in->bnd,nThisFile);
            for (i=in->iLower+1; i<in->iUpper; ++i ) {
                int rID = pst->mdl->ReqService(i,PST_SENDPARTICLES,&pst->idSelf,sizeof(pst->idSelf));
                pkdWriteFromNode(plcl->pkd,i,fio,in->dvFac,in->dTuFac,in->bnd,nThisFile);
                pst->mdl->GetReply(rID);
            }

            pkdWriteHeaderFIO(plcl->pkd, fio, 1./sqrt(in->dvFac), in->dTime, nTotal, nThisFile,
                              in->dBoxSize, in->HubbleParam, in->nProcessorsTotal, in->units);

            fioClose(fio);
        }
    }

    return 0;
}

int pstTreeInitMarked(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    mdlassert(pst->mdl,nIn == 0);

    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_TREEINITMARKED,vin,nIn);
        pstTreeInitMarked(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
    auto pTop = pkd->tree[static_cast<KDN *>(vout)];
    // double minside;

    /* We need to save our cells so we can update them later */
    if (pst->nLeaves > 1) {
        auto pCell1 = pTop + 1;
        auto pCell2 = pTop + pst->nLower*2;

        /* We will accumulate the top tree here */
        int rID = pst->mdl->ReqService(pst->idUpper,PST_BUILDTREE,vin,nIn);
        nOut = pstBuildTree(pst->pstLower,vin,nIn,pCell1,(pst->nLower*2-1) * pkd->tree.ElementSize());
        assert(nOut == (pst->nLower*2-1) * pkd->tree.ElementSize());
        auto nOutUpper = pst->mdl->GetReply(rID,(pst->nUpper*2-1) * pkd->tree.ElementSize(),pCell2);
        assert(nOutUpper == (pst->nUpper*2-1) * pkd->tree.ElementSize());
        nOut += nOutUpper + pkd->tree.ElementSize();

        /*
        ** Combine Cell1 and Cell2 into pCell
        ** to find cell CoM, bounds and multipoles.
        ** This also computes the opening radius for gravity.
        */
        pTop->bMax() = HUGE_VAL;  /* initialize bMax for CombineCells */
        auto minside = pst->bnd.minside();
        pkdCombineCells1(pkd,pTop,pCell1,pCell2);
        pTop->calc_open(minside);
        pkdCombineCells2(pkd,pTop,pCell1,pCell2);

        /* Get our cell ready */
        pTop->set_top_tree(pst->nLower*2);
    }
    else {
        auto pRoot = pkd->tree[uRoot];
        pkd->TreeAlignNode();
        pkdTreeBuild(plcl->pkd,in->nBucket,in->nGroup,in->uRoot,in->utRoot,in->ddHonHLimit);
        *pTop = *pRoot;
        /* Get our cell ready */
        pTop->set_remote(pst->idSelf);
        nOut = pkd->tree.ElementSize();
    }
    return nOut;
}

int pstPhysicalSoft(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inPhysicalSoft *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inPhysicalSoft));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_PHYSICALSOFT,in,nIn);
        pstPhysicalSoft(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdPhysicalSoft(plcl->pkd,in->dSoftMax,in->dFac,in->bSoftMaxMul);
    }
    return 0;
}

int pstHopLink(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inHopLink *>(vin);
    auto nOutGroups = static_cast<uint64_t *>(vout);
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopLink));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_HOP_LINK,in,nIn);
        pstHopLink(pst->pstLower,in,nIn,vout,nOut);
        pst->mdl->GetReply(rID,nOutUpper);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_HOP_JOIN,in,nIn);
        pstHopJoin(pst->pstLower,in,nIn,&outLower,nOut);
        nOut = pst->mdl->GetReply(rID,outUpper);
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
    auto in = static_cast<struct inHopFinishUp *>(vin);
    auto nOutGroups = static_cast<uint64_t *>(vout);
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inHopFinishUp));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_HOP_FINISH_UP,vin,nIn);
        pstHopFinishUp(pst->pstLower,vin,nIn,vout,nOut);
        pst->nLowerGroups = *nOutGroups;
        pst->mdl->GetReply(rID,nOutUpper);
        *nOutGroups += nOutUpper;
    }
    else {
        LCL *plcl = pst->plcl;
        *nOutGroups = pkdHopFinishUp(plcl->pkd,in->nMinGroupSize,in->bPeriodic,in->fPeriod);
        pst->nLowerGroups = 0;
    }
    return sizeof(uint64_t);
}

int pstHopTreeBuild(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inHopTreeBuild *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inHopTreeBuild));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_HOP_TREE_BUILD,vin,nIn);
        pstHopTreeBuild(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdHopTreeBuild(plcl->pkd,in->nBucket,in->nGroup);
    }
    return 0;
}

int pstHopGravity(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inHopGravity *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inHopGravity));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_HOP_GRAVITY,vin,nIn);
        pstHopGravity(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        LCL *plcl = pst->plcl;
        double dFlop=0, dPartSum=0, dCellSum=0;
        pkdGravWalkHop(plcl->pkd,in->dTime,in->dTheta,&dFlop,&dPartSum,&dCellSum);
    }
    return 0;
}

int pstHopUnbind(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inHopUnbind *>(vin);
    auto out = static_cast<struct outHopUnbind *>(vout);
    struct outHopUnbind outUpper, outLower;
    mdlassert(pst->mdl,nIn == sizeof(struct inHopUnbind));
    mdlassert(pst->mdl,nOut == sizeof(struct outHopUnbind));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_HOP_UNBIND,vin,nIn);
        pstHopUnbind(pst->pstLower,vin,nIn,&outLower,sizeof(outLower));
        pst->mdl->GetReply(rID,outUpper);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_GROUP_RELOCATE,NULL,0);
        pstGroupRelocate(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        in->iGlobalStart += pst->nLowerGroups; /* for the upper subset the global id starting point is more by the number of lower groups */
        int rID = pst->mdl->ReqService(pst->idUpper,PST_GROUP_STATS,vin,nIn);
        in->iGlobalStart -= pst->nLowerGroups; /* set the value back for the lower subset traversal */
        pstGroupStats(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdCalculateGroupStats(plcl->pkd,in->bPeriodic,in->dPeriod,in->rEnvironment,in->iGlobalStart);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_BH_PLACESEED,in,nIn);
        pstPlaceBHSeed(pst->pstLower,vin,nIn,vout,nOut);
        nOut = pst->mdl->GetReply(rID,outUpper);
        assert(nOut==sizeof(struct outPlaceBHSeed));
        out->nBHs += outUpper.nBHs;
    }
    else {
        LCL *plcl = pst->plcl;
        out->nBHs = pkdPlaceBHSeed(plcl->pkd,in->dTime, in->dScaleFactor, in->uRungMax, in->dDenMin,
                                   in->dBHMhaloMin, in->dTau, in->dBHSeedMass);
    }
    return sizeof(struct outPlaceBHSeed);

}
int pstBHInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inPlaceBHSeed *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inPlaceBHSeed));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_BH_INIT,in,nIn);
        pstBHInit(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdBHInit(plcl->pkd,in->uRungMax);
    }
    return 0;

}
int pstBHReposition(PST pst,void *vin,int nIn,void *vout,int nOut) {
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_BH_REPOSITION,NULL,0);
        pstBHReposition(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdBHReposition(plcl->pkd);
    }
    return 0;

}
int pstBHAccretion(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inBHAccretion *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inBHAccretion));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_BH_ACCRETION,in,nIn);
        pstBHAccretion(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
    }
    else {
        LCL *plcl = pst->plcl;
        pkdBHAccretion(plcl->pkd, in->dScaleFactor);
    }
    return 0;

}
#endif

int pstSmooth(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inSmooth *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_SMOOTH,in,nIn);
        pstSmooth(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_RESMOOTH,in,nIn);
        pstReSmooth(pst->pstLower,vin,nIn,vout,nOut);
        nOut = pst->mdl->GetReply(rID,outUpper);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_RESMOOTHNODE,in,nIn);
        pstReSmoothNode(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID,outUpper);
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
        out->sComputing.dSum     = pst->mdl->TimeComputing();
        out->sWaiting.dSum       = pst->mdl->TimeWaiting();
        out->sSynchronizing.dSum = pst->mdl->TimeSynchronizing();
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_GRAVITY,in,nIn);
        pstGravity(pst->pstLower,in,nIn,vout,nOut);
        /*
        ** Make a temporary copy of the reduct part of the out buffer before setting it as the
        ** reply buffer. The reduct part follows at the end of all the outGravityPerProc entries.
        */
        pst->mdl->GetReply(rID,tmp);
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
        outr->nTilesTotal    += tmp.nTilesTotal;
        outr->nTilesCPU      += tmp.nTilesCPU;
    }
    else {
#ifdef __linux__
        FILE *fp;
        char buffer[512], *save, *f, *v;
#endif
        PKD pkd = plcl->pkd;
        pkdGravAll(pkd,&in->kick,&in->lc,&in->ts,
                   in->dTime,in->nReps,in->bPeriodic,in->bGPU,
                   in->bEwald,in->iRoot1,in->iRoot2,in->dEwCut,in->dEwhCut,in->dTheta,&in->SPHoptions,
                   &outr->nActive,
                   &outr->sPart.dSum,&outr->sPartNumAccess.dSum,&outr->sPartMissRatio.dSum,
                   &outr->sCell.dSum,&outr->sCellNumAccess.dSum,&outr->sCellMissRatio.dSum,
                   &outr->sFlop.dSum,outr->nRung);
        outr->dFlopSingleCPU = 1e-9*pkd->dFlopSingleCPU;
        outr->dFlopDoubleCPU = 1e-9*pkd->dFlopDoubleCPU;
        outr->dFlopSingleGPU = 1e-9*pkd->dFlopSingleGPU;
        outr->dFlopDoubleGPU = 1e-9*pkd->dFlopDoubleGPU;
        outr->nTilesTotal    = pkd->nTilesTotal;
        outr->nTilesCPU      = pkd->nTilesCPU;
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
        outr->sComputing.dSum     = pst->mdl->TimeComputing();
        outr->sWaiting.dSum       = pst->mdl->TimeWaiting();
        outr->sSynchronizing.dSum = pst->mdl->TimeSynchronizing();
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_CALCEANDL,NULL,0);
        pstCalcEandL(pst->pstLower,NULL,0,out,nOut);
        pst->mdl->GetReply(rID,outLcl);
        out->T += outLcl.T;
        out->U += outLcl.U;
        out->Eth += outLcl.Eth;
        for (k=0; k<3; k++) out->L[k] += outLcl.L[k];
        for (k=0; k<3; k++) out->F[k] += outLcl.F[k];
        out->W += outLcl.W;
    }
    else {
        pkdCalcEandL(plcl->pkd,out->T,out->U,out->Eth,out->L,out->F,out->W);
    }
    return sizeof(struct outCalcEandL);
}

int pstDrift(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inDrift *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_DRIFT,in,nIn);
        pstDrift(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdDrift(plcl->pkd,in->iRoot,in->dTime,in->dDelta,in->dDeltaVPred,in->dDeltaUPred,in->bDoGas,in->bGasEvolveDensity);
    }
    return 0;
}

int pstSetGlobalDt(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct outGetMinDt *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct outGetMinDt));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_SETGLOBALDT,in,nIn);
        pstSetGlobalDt(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_REORDERINNODES,NULL,0);
        pstReorderWithinNodes(pst->pstLower,NULL,0,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdReorderWithinNodes(plcl->pkd);
    }
    return 0;
}
#endif


#ifdef GRACKLE
int pstGrackleInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    auto in = static_cast<struct inGrackleInit *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inGrackleInit));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_GRACKLEINIT,in,nIn);
        pstGrackleInit(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdGrackleInit(plcl->pkd, in->bComove, in->dScaleFactor, in->achCoolingTable,
                       in->units);
    }
    return 0;
}
#endif

int pstChemCompInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inChemCompInit *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inChemCompInit));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_CHEMCOMPINIT,in,nIn);
        pstChemCompInit(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdChemCompInit(plcl->pkd, *in);
    }
    return 0;
}

int pstROParticleCache(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_ROPARTICLECACHE,vin,nIn);
        pstROParticleCache(pst->pstLower,NULL,0,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        /*
        ** Start particle caching space.
        */
        PKD pkd = plcl->pkd;
        pkd->mdl->CacheInitialize(CID_PARTICLE,std::make_shared<particleCACHEhelper>(pkd));
    }
    return 0;
}

int pstParticleCacheFinish(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_PARTICLECACHEFINISH,vin,nIn);
        pstParticleCacheFinish(pst->pstLower,NULL,0,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        /*
        ** Start particle caching space.
        */
        PKD pkd = plcl->pkd;
        pkd->mdl->FinishCache(CID_PARTICLE);
    }
    return 0;
}

int pstKick(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inKick *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_KICK,in,nIn);
        pstKick(pst->pstLower,in,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_KICKTREE,in,nIn);
        pstKickTree(pst->pstLower,in,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_SETTOTAL,NULL,0);
        pstSetTotal(pst->pstLower,NULL,0,out,nOut);
        pst->mdl->GetReply(rID,oute);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_SETWRITESTART,in,nIn);
        in->nWriteStart = nWriteStart;
        pstSetWriteStart(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_ACCELSTEP,vin,nIn);
        pstAccelStep(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdAccelStep(plcl->pkd,in->uRungLo,in->uRungHi,
                     in->dDelta,in->iMaxRung,
                     in->dEta,in->dVelFac,in->dAccFac,
                     in->bDoGravity,in->bEpsAcc,in->dhMinOverSoft);
    }
    return 0;
}

int pstDensityStep(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inDensityStep *>(vin);

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_DENSITYSTEP,vin,nIn);
        pstDensityStep(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_CORRECTENERGY,vin,nIn);
        pstCorrectEnergy(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdCorrectEnergy(plcl->pkd,in->dTuFac,in->z,in->dTime,in->iDirection);
    }
    return 0;
}

int pstColNParts(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outColNParts *>(vout);
    auto iUpper = pst->idUpper - pst->idSelf;
    auto nBytesLower = iUpper*sizeof(*out);
    struct outColNParts *outUp = out + iUpper;

    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_COLNPARTS,vin,nIn);
        pstColNParts(pst->pstLower,vin,nIn,vout,nBytesLower);
        pst->mdl->GetReply(rID,nOut-nBytesLower,outUp);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_NEWORDER,vin,nIn);
        pstNewOrder(pst->pstLower, vin, nIn, NULL, 0);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdNewOrder(plcl->pkd, in[pst->idSelf]);
    }
    return  0;
}

int pstSetNParts(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inSetNParts *>(vin);

    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_SETNPARTS,vin,nIn);
        pstSetNParts(pst->pstLower, vin, nIn, NULL, 0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_NEW_FOF,in,nIn);
        pstNewFof(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_FOF_PHASES,vin,nIn);
        pstFofPhases(pst->pstLower,vin,nIn,out,nOut);
        bMadeProgress = out->bMadeProgress;
        pst->mdl->GetReply(rID,nOut,out);
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
** JST: added the count of lower subtree number of local groups.
*/
int pstFofFinishUp(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inFofFinishUp *>(vin);
    auto nOutGroups = static_cast<uint64_t *>(vout);
    uint64_t nOutUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inFofFinishUp));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_FOF_FINISH_UP,vin,nIn);
        pstFofFinishUp(pst->pstLower,vin,nIn,vout,nOut);
        pst->nLowerGroups = *nOutGroups;
        pst->mdl->GetReply(rID,nOutUpper);
        *nOutGroups += nOutUpper;
    }
    else {
        LCL *plcl = pst->plcl;
        *nOutGroups = pkdFofFinishUp(plcl->pkd,in->nMinGroupSize);
        pst->nLowerGroups = 0;
    }
    return sizeof(uint64_t);
}

int pstMemStatus(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto out = static_cast<struct outMemStatus *>(vout);
    auto iUpper = pst->idUpper - pst->idSelf;
    auto nLower = iUpper * sizeof(struct outMemStatus);
    struct outMemStatus *outUp = out + iUpper;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_MEMSTATUS,vin,nIn);
        pstMemStatus(pst->pstLower,vin,nIn,out,nLower);
        pst->mdl->GetReply(rID,nOut-nLower,outUp);
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
        out->nBytesTree = plcl->pkd->tree.memory();
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
        auto outUp = std::make_unique<PARTCLASS[]>(PKD_MAX_CLASSES);
        int rID = pst->mdl->ReqService(pst->idUpper,PST_GETCLASSES,vin,nIn);
        nOut = pstGetClasses(pst->pstLower,vin,nIn,out,nOut);
        nUp = pst->mdl->GetReply(rID,PKD_MAX_CLASSES*sizeof(outUp[0]),outUp.get());
        n = nUp / sizeof(PARTCLASS);
        mdlassert(pst->mdl,n*sizeof(PARTCLASS)==nUp);
        nUp = n;
        n = nOut / sizeof(PARTCLASS);
        mdlassert(pst->mdl,n*sizeof(PARTCLASS)== nOut);
        for ( i=0; i<nUp; i++ ) {
            if (std::find(out, out + n, outUp[i]) == out + n) {
                out[n++] = outUp[i];
            }
        }
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_SETCLASSES,vin,nIn);
        pstSetClasses(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_PROFILE,vin,nIn);
        pstProfile(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_CALCDISTANCE,vin,nIn);
        pstCalcDistance(pst->pstLower,vin,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_CALCCOM,vin,nIn);
        pstCalcCOM(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID,outUpper);
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
                   out->com, out->vcm, out->L, out->M, out->N);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_CALCMTOT,vin,nIn);
        pstCalcMtot(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID,outUpper);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_SETSPHOPTIONS,vin,nIn);
        pstSetSPHoptions(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
    }
    else {
        copySPHOptions(&in->SPHoptions, &plcl->pkd->SPHoptions);
    }
    return 0;
}

int pstResetCOM(PST pst,void *vin, int nIn, void *vout, int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inResetCOM *>(vin);

    assert( nIn==sizeof(struct inResetCOM) );
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_RESETCOM,vin,nIn);
        pstResetCOM(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdResetCOM(plcl->pkd,in->r_com,in->v_com);
    }
    return 0;
}

int pstInitializeEOS(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_INITIALIZEEOS,vin,nIn);
        pstInitializeEOS(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdInitializeEOS(plcl->pkd);
    }
    return 0;
}

int pstUpdateGasValues(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inUpdateGasValues *>(vin);

    assert( nIn==sizeof(struct inUpdateGasValues) );
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_UPDATEGASVALUES,vin,nIn);
        pstUpdateGasValues(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID);
    }
    else {
        pkdUpdateGasValues(plcl->pkd, &in->kick, &in->SPHoptions);
    }
    return 0;
}

int pstTreeUpdateFlagBounds(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    auto in = static_cast<struct inTreeUpdateFlagBounds *>(vin);
    uint32_t uRoot = in->uRoot;
    auto pTop = pkd->tree[static_cast<KDN *>(vout)];

    /* We need to save our cells so we can update them later */
    if (pst->nLeaves > 1) {
        auto pCell1 = pTop + 1;
        auto pCell2 = pTop + pst->nLower*2;
        auto nRemoveLower = pst->nLower * 2 * pkd->tree.ElementSize();

        /* We will accumulate the top tree here */
        int rID = pst->mdl->ReqService(pst->idUpper,PST_TREEUPDATEFLAGBOUNDS,vin,nIn);
        auto nOutLower = pstTreeUpdateFlagBounds(pst->pstLower,vin,nIn,pCell1,(pst->nLower*2-1) * pkd->tree.ElementSize());
        assert(nOutLower == (pst->nLower*2-1) * pkd->tree.ElementSize());
        auto nOutUpper = pst->mdl->GetReply(rID,nOut-nRemoveLower,pCell2);
        assert(nOutUpper == (pst->nUpper*2-1) * pkd->tree.ElementSize());
        nOut = nOutLower + nOutUpper + pkd->tree.ElementSize();

        /*
        ** Combine Cell1 and Cell2 into pCell
        ** to find cell CoM, bounds and multipoles.
        ** This also computes the opening radius for gravity.
        */
        pTop->bMax() = HUGE_VAL;  /* initialize bMax for CombineCells */
        auto minside = pst->bnd.minside();
        pkdCombineCells1(pkd,pTop,pCell1,pCell2);
        pTop->calc_open(minside);
        pkdCombineCells2(pkd,pTop,pCell1,pCell2);

        /* Get our cell ready */
        pTop->set_top_tree(pst->nLower*2);
    }
    else {
        auto pRoot = pkd->tree[uRoot];
        pkd->TreeAlignNode();
        pkdTreeUpdateFlagBounds(plcl->pkd,uRoot,&in->SPHoptions);
        *pTop = *pRoot;
        /* Get our cell ready */
        pTop->set_remote(pst->idSelf);
        /* iLower is valid = ROOT */
        nOut = pkd->tree.ElementSize();
    }
    return nOut;
}

#ifdef MDL_FFTW
int pstGridCreateFFT(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inGridCreateFFT *>(vin);
    assert (nIn==sizeof(struct inGridCreateFFT) );
    if (pstNotCore(pst)) {
        int rID = pst->mdl->ReqService(pst->idUpper, PST_GRID_CREATE_FFT, vin, nIn);
        pstGridCreateFFT(pst->pstLower, vin, nIn, NULL, 0);
        pst->mdl->GetReply(rID);
    }
    else {
        PKD pkd = plcl->pkd;
        assert(pkd->fft == NULL);
        pkd->fft = pst->mdl->FFTInitialize(in->nGrid, in->nGrid, in->nGrid, 0, 0);
    }
    return 0;
}

int pstGridDeleteFFT(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    if (pstNotCore(pst)) {
        int rID = pst->mdl->ReqService(pst->idUpper, PST_GRID_DELETE_FFT, vin, nIn);
        pstGridDeleteFFT(pst->pstLower, vin, nIn, NULL, 0);
        pst->mdl->GetReply(rID);
    }
    else {
        PKD pkd = plcl->pkd;
        assert(pkd->fft != NULL);
        pst->mdl->FFTFinish(plcl->pkd->fft);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_TOTALMASS,vin,nIn);
        pstTotalMass(pst->pstLower,vin,nIn,vout,nOut);
        nOut = pst->mdl->GetReply(rID,outUpper);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_GETMINDT,vin,nIn);
        pstGetMinDt(pst->pstLower,vin,nIn,vout,nOut);
        pst->mdl->GetReply(rID,outUpper);
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
        int rID = pst->mdl->ReqService(pst->idUpper,PST_LIGHTCONE_OPEN,in,nIn);
        pstLightConeOpen(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        PKD pkd = pst->plcl->pkd;
        pkdLightConeOpen(pkd, in->achOutFileLcp, in->nSideHealpix);
    }
    return 0;
}

int pstLightConeClose(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<struct inLightConeClose *>(vin);
    mdlassert(pst->mdl,nIn == sizeof(struct inLightConeClose));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_LIGHTCONE_CLOSE,in,nIn);
        pstLightConeClose(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
    }
    else {
        PKD pkd = pst->plcl->pkd;
        auto achOutFile = makeName(in->achOutFileHpb,pkd->mdl->Self());
        pkdLightConeClose(pst->plcl->pkd,achOutFile.c_str());
    }
    return 0;
}

int pstLightConeVel(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    auto in = static_cast<struct inLightConeVel *>(vin);

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_LIGHTCONEVEL,in,nIn);
        pstLightConeVel(pst->pstLower,in,nIn,NULL,0);
        pst->mdl->GetReply(rID);
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
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper,PST_GET_PARTICLES,vin,nIn);
        auto nOutLower = pstGetParticles(pst->pstLower,vin,nIn,vout,nOut);
        auto nOutUpper = pst->mdl->GetReply(rID,nOut-nOutLower,out + (nOutLower / sizeof(struct outGetParticles)));
        nOut = nOutLower + nOutUpper;
    }
    else {
        int nParticles = nIn / sizeof(uint64_t);
        int n = pkdGetParticles(plcl->pkd,nParticles, ID, out );
        nOut = n * sizeof(struct outGetParticles);
    }
    return nOut;
}
