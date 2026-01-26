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

#include <cstddef>
#include <type_traits>
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
#include "io/classes.h"
#include "core/swapall.h"

#include "cooling/cooling.h"
#ifdef GRACKLE
    #include "cooling_grackle/cooling_grackle.h"
#endif
#include "blackhole/merger.h"
#include "blackhole/evolve.h"
#include "blackhole/seed.h"
#include "blackhole/init.h"
#include "blackhole/drift.h"
#include "blackhole/accretion.h"
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#include "ic/ic.h"

#define GET_PARTICLES_MAX 20 /* We have a nested loop, so don't increase this */

void pstAddServices(PST pst, mdl::mdlClass* mdl) {
    int nThreads;

    nThreads = mdl->Threads();

    mdl->AddService(PST_READFILE, pst, (fcnService_t*)pstReadFile, sizeof(struct inReadFile) + PST_MAX_FILES * (sizeof(fioSpeciesList) + PATH_MAX), 0,
                    Legacy::ReadFile);
    mdl->AddService(PST_WRITE, pst, (fcnService_t*)pstWrite, sizeof(struct inWrite), 0, Legacy::Write);
    mdl->AddService(PST_SENDARRAY, pst, (fcnService_t*)pstSendArray, sizeof(struct inSendArray), 0, Legacy::SendArray);
    mdl->AddService(PST_CHECKPOINT, pst, (fcnService_t*)pstCheckpoint, sizeof(struct inWrite), 0, Legacy::Checkpoint);
    mdl->AddService(PST_OUTPUT, pst, (fcnService_t*)pstOutput, sizeof(struct inOutput), 0, Legacy::Output);
    mdl->AddService(PST_OUTPUT_SEND, pst, (fcnService_t*)pstOutputSend, sizeof(struct inOutputSend), 0, Legacy::OutputSend);
    /*
    ** Calculate the number of levels in the top tree and use it to
    ** define the size of the messages.
    */
    mdl->AddService(PST_SMOOTH, pst, (fcnService_t*)pstSmooth, sizeof(struct inSmooth), 0, Legacy::Smooth);
    mdl->AddService(PST_GRAVITY, pst, (fcnService_t*)pstGravity, sizeof(struct inGravity), sizeof(struct outGravityReduct), Legacy::Gravity);
    mdl->AddService(PST_CALCEANDL, pst, (fcnService_t*)pstCalcEandL, 0, sizeof(struct outCalcEandL), Legacy::CalcEandL);
#ifdef GRACKLE
    mdl->AddService(PST_GRACKLEINIT, pst, (fcnService_t*)pstGrackleInit, sizeof(struct inGrackleInit), 0, Legacy::GrackleInit);
#endif
    mdl->AddService(PST_CHEMCOMPINIT, pst, (fcnService_t*)pstChemCompInit, sizeof(struct inChemCompInit), 0, Legacy::ChemCompInit);
    mdl->AddService(PST_BH_PLACESEED, pst, (fcnService_t*)pstPlaceBHSeed, sizeof(struct inPlaceBHSeed), sizeof(struct outPlaceBHSeed), Legacy::PlaceBHSeed);
    mdl->AddService(PST_BH_INIT, pst, (fcnService_t*)pstBHInit, sizeof(struct inPlaceBHSeed), 0, Legacy::BHInit);
    mdl->AddService(PST_BH_REPOSITION, pst, (fcnService_t*)pstBHReposition, 0, 0, Legacy::BHReposition);
    mdl->AddService(PST_BH_ACCRETION, pst, (fcnService_t*)pstBHAccretion, sizeof(struct inBHAccretion), 0, Legacy::BHAccretion);
    //
    mdl->AddService(PST_PHYSICALSOFT, pst, (fcnService_t*)pstPhysicalSoft, sizeof(struct inPhysicalSoft), 0, Legacy::PhysicalSoft);
    mdl->AddService(PST_ONENODEREADINIT, pst, (fcnService_t*)pstOneNodeReadInit,
                    sizeof(struct inReadFile) + PST_MAX_FILES * (sizeof(fioSpeciesList) + PATH_MAX), nThreads * sizeof(fioSpeciesList),
                    Legacy::OneNodeReadInit);
    mdl->AddService(PST_RESMOOTH, pst, (fcnService_t*)pstReSmooth, sizeof(struct inSmooth), sizeof(struct outSmooth), Legacy::ReSmooth);
#ifdef OPTIM_SMOOTH_NODE
    mdl->AddService(PST_RESMOOTHNODE, pst, (fcnService_t*)pstReSmoothNode, sizeof(struct inSmooth), sizeof(struct outSmooth), Legacy::ReSmoothNode);
#endif
#ifdef OPTIM_REORDER_IN_NODES
    mdl->AddService(PST_REORDERINNODES, pst, (fcnService_t*)pstReorderWithinNodes, 0, 0, Legacy::ReorderInNodes);
#endif
#ifdef STELLAR_EVOLUTION
    mdl->AddService(PST_STELLAREVOLUTIONINIT, pst, (fcnService_t*)pstStellarEvolutionInit, sizeof(struct inStellarEvolutionInit), 0,
                    Legacy::StellarEvolutionInit);
#endif
    mdl->AddService(PST_INITIALIZEPSTORE, pst, (fcnService_t*)pstInitializePStore, sizeof(struct inInitializePStore), sizeof(struct outInitializePStore),
                    Legacy::InitializePStore);
#ifdef MDL_FFTW
    mdl->AddService(PLT_MOVEIC, pst, (fcnService_t*)pltMoveIC, sizeof(struct inMoveIC), 0, Legacy::MoveIC);
    mdl->AddService(PST_MOVEIC, pst, (fcnService_t*)pstMoveIC, sizeof(struct ServiceGenerateIC::input), 0, Legacy::LocalMoveIC);
#endif
    mdl->AddService(PST_MEMSTATUS, pst, (fcnService_t*)pstMemStatus, 0, nThreads * sizeof(struct outMemStatus), Legacy::MemStatus);
    mdl->AddService(PST_SETSPHOPTIONS, pst, (fcnService_t*)pstSetSPHoptions, sizeof(struct inSetSPHoptions), 0, Legacy::SetSPHoptions);
    mdl->AddService(PST_INITIALIZEEOS, pst, (fcnService_t*)pstInitializeEOS, 0, 0, Legacy::InitializeEOS);
}

static void initializePStore(PKD* ppkd, MDL mdl, struct inInitializePStore* in) {
    *ppkd = new pkdContext(static_cast<mdl::mdlClass*>(mdl), in->nStore, in->nDecoratorStore, in->nBucket, in->nMinTotalStore, in->nMinEphemeral,
                           in->nEphemeralBytes, in->nTreeBitsLo, in->nTreeBitsHi, in->iCacheSize, in->iCacheMaxInflight, in->iWorkQueueSize, in->mpi_recv_mode,
                           in->mpi_recv_count, in->mpi_RDMA, in->fPeriod, in->nSpecies[FIO_SPECIES_DARK], in->nSpecies[FIO_SPECIES_SPH],
                           in->nSpecies[FIO_SPECIES_STAR], in->nSpecies[FIO_SPECIES_BH], in->mMemoryModel, in->nIntegerFactor);
}

int pstInitializePStore(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto in = static_cast<struct inInitializePStore*>(vin);
    auto out = static_cast<struct outInitializePStore*>(vout);
    assert(nIn == sizeof(struct inInitializePStore));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_INITIALIZEPSTORE, in, nIn);
        pstInitializePStore(pst->lower(), vin, nIn, vout, nOut);
        pst->mdl()->GetReply(rID, nOut, vout);
    }
    else {
        if (plcl->pkd) delete plcl->pkd;
        initializePStore(&plcl->pkd, pst->mdl(), in);
        out->nSizeParticle = plcl->pkd->particles.ElementSize();
        out->nSizeDecorator = plcl->pkd->decorators.ElementSize();
        out->nSizeNode = plcl->pkd->tree.ElementSize();
    }
    return sizeof(*out);
}

int pstOneNodeReadInit(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inReadFile*>(vin);
    auto pout = static_cast<fioSpeciesList*>(vout);
    fioSpeciesList nFileStart, nFileEnd, nFileTotal, nFileSplit;

    assert(nIn == sizeof(struct inReadFile));
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        nFileStart[s] = in->nSpeciesStart[s];
        nFileEnd[s] = in->nSpeciesEnd[s];
        nFileTotal[s] = nFileEnd[s] - nFileStart[s] + 1;
    }

    if (pst->NotCore()) {
        int rID;
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            nFileSplit[s] = nFileStart[s] + pst->n_lower() * (nFileTotal[s] / pst->n_total());
            in->nSpeciesStart[s] = nFileSplit[s];
        }
        rID = pst->mdl()->ReqService(pst->upper(), PST_ONENODEREADINIT, in, nIn);
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            in->nSpeciesStart[s] = nFileStart[s];
            in->nSpeciesEnd[s] = nFileSplit[s] - 1;
        }
        pstOneNodeReadInit(pst->lower(), in, nIn, pout, nOut);
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            in->nSpeciesEnd[s] = nFileEnd[s];
        }
        auto nBytesLower = pst->n_lower() * sizeof(*pout);
        pst->mdl()->GetReply(rID, nOut - nBytesLower, pout + pst->n_lower());
    }
    else {
        /*
        ** Determine the size of the local particle store.
        */
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            pout[0][s] = nFileTotal[s]; /* Truncated: okay */
        }
    }
    return sizeof(*pout) * pst->n_total();
}

int pstReadFile(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto plcl = pst->lcl();
    auto in = static_cast<struct inReadFile*>(vin);
    FIO fio;
    fioSpeciesList nNodeSpeciesStart, nNodeSpeciesEnd, nNodeSpeciesTotal, nNodeSpeciesSplit;
    int rID;

    assert(nIn >= sizeof(struct inReadFile));

    for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        if (in->nSpeciesEnd[s] == 0) continue;
        // First element to read
        nNodeSpeciesStart[s] = in->nSpeciesStart[s];
        // Last element to read
        nNodeSpeciesEnd[s] = in->nSpeciesEnd[s];
        // Total count of read elements
        nNodeSpeciesTotal[s] = nNodeSpeciesEnd[s] - nNodeSpeciesStart[s] + 1;
    }
    if (pst->NotCore() && in->nProcessors > 1) {
        int nProcessors = in->nProcessors;
        int nProcUpper = pst->n_upper() * nProcessors / pst->n_total();
        int nProcLower = nProcessors - nProcUpper;

        // Index of the first element in the second half of the split
        for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1))
            nNodeSpeciesSplit[s] = nNodeSpeciesStart[s] + pst->n_lower() * (nNodeSpeciesTotal[s] / pst->n_total());

        // Upper split
        if (nProcessors > 1) {
            in->nProcessors = nProcUpper;
            for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
                if (in->nSpeciesEnd[s] == 0) continue;
                in->nSpeciesStart[s] = nNodeSpeciesSplit[s];
                in->nSpeciesEnd[s] = nNodeSpeciesEnd[s];
            }
            rID = pst->mdl()->ReqService(pst->upper(), PST_READFILE, in, nIn);
        }

        // Lower split
        in->nProcessors = nProcLower;
        for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            if (in->nSpeciesEnd[s] == 0) continue;
            in->nSpeciesStart[s] = nNodeSpeciesStart[s];
            in->nSpeciesEnd[s] = nNodeSpeciesSplit[s] - 1;
        }
        pstReadFile(pst->lower(), in, nIn, NULL, 0);

        if (nProcessors <= 1) {
            in->nProcessors = nProcUpper;
            for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
                if (in->nSpeciesEnd[s] == 0) continue;
                in->nSpeciesStart[s] = nNodeSpeciesSplit[s];
                in->nSpeciesEnd[s] = nNodeSpeciesEnd[s];
            }
            rID = pst->mdl()->ReqService(pst->upper(), PST_READFILE, in, nIn);
        }

        pst->mdl()->GetReply(rID);
    }
    else {
        // Number of particles to be read per core per species
        auto nSpeciesParts = new fioSpeciesList[pst->n_total()];
        int i;
        PKD pkd = plcl->pkd;
        auto mdl = pkd->mdl;

        pstOneNodeReadInit(pst, in, sizeof(*in), nSpeciesParts, pst->n_total() * sizeof(fioSpeciesList));

        fio = fioLoad(in + 1, in->dOmega0, in->dOmegab);
        assert(fio != NULL);

        fioSpeciesList nRead, nStart, nBasePart;
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            // We take care of the first nSpeciesParts later
            nStart[s] = in->nSpeciesStart[s] + nSpeciesParts[0][s];
        }
        for (i = 1; i < pst->n_total(); ++i) {
            int id = mdl->Self() + i;
            int inswap;
            for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
                if (in->nSpeciesEnd[s] == 0) continue;
                if (nSpeciesParts[i][s] == 0) continue;

                // Read particles into the local storage.
                assert(pkd->FreeStore() >= nSpeciesParts[i][s]);
                pkdReadFIO(pkd, fio, nStart[s], nSpeciesParts[i][s], in->dvFac, in->dTuFac, FIO_SPECIES(s));
                nStart[s] += nSpeciesParts[i][s];
            }

            // Now shove them over to the remote processor.
            pkd->swap_classes(id);
            inswap = mdl->Self();
            rID = pst->mdl()->ReqService<ServiceSwapAll>(id, inswap);
            pkdSwapAll(pkd, id);
            pst->mdl()->GetReply(rID);
        }
        // Finish with our own particles
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            if (in->nSpeciesEnd[s] == 0) continue;
            if (nSpeciesParts[0][s] == 0) continue;
            pkdReadFIO(pkd, fio, in->nSpeciesStart[s], nSpeciesParts[0][s], in->dvFac, in->dTuFac, FIO_SPECIES(s));
        }
        delete[] nSpeciesParts;
        fioClose(fio);
    }
    return 0;
}

static std::string makeName(std::string_view inName, int iIndex) {
    return fmt::format(inName, fmt::arg("i", iIndex));
}

int pstCheckpoint(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inWrite*>(vin);

    assert(nIn == sizeof(struct inWrite));
    if (pst->NotCore()) {
        uint64_t nProcessors = in->nProcessors;
        if (nProcessors > 1) { /* Keep going parallel */
            int nLower, nUpper;
            nLower = nProcessors * pst->n_lower() / pst->n_total();
            if (nLower == 0) nLower = 1;
            nUpper = nProcessors - nLower;
            in->nProcessors = nUpper;
            int rID = pst->mdl()->ReqService(pst->upper(), PST_CHECKPOINT, in, nIn);
            in->nProcessors = nLower;
            pstCheckpoint(pst->lower(), in, nIn, NULL, 0);
            pst->mdl()->GetReply(rID);
        }
        else { /* Serialize these processors now */ pstCheckpoint(pst->lower(), in, nIn, NULL, 0);
            int rID = pst->mdl()->ReqService(pst->upper(), PST_CHECKPOINT, in, nIn);
            pst->mdl()->GetReply(rID);
        }
    }
    else {
        PKD pkd = pst->lcl()->pkd;
        auto achOutFile = makeName(in->achOutFile, pkd->mdl->Self());
        auto achOutDecoratorFile = makeName(in->achOutFile, pkd->mdl->Self()) + ".dec";
        // Need to add prefix as with previous "makeName" implementation
        pkdCheckpoint(pkd, achOutFile.c_str(), achOutDecoratorFile.c_str());
    }
    return 0;
}

int pstSendArray(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inSendArray*>(vin);
    pkdSendArray(pst->lcl()->pkd, in->iTo, in->field, in->iUnitSize, in->dvFac, in->bMarked);
    return 0;
}

int pstWrite(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inWrite*>(vin);
    FIO fio;
    int i;

    assert(nIn == sizeof(struct inWrite));

    if (pst->NotCore()) {
        int nProcessors = in->nProcessors;
        /* If we are the writer (nProcessors==1) or a producer (==1) keep requesting */
        if (nProcessors <= 1) {
            pstWrite(pst->lower(), in, nIn, NULL, 0);
        }
        /* Split writing tasks between child nodes */
        else {
            int nLeft = nProcessors / 2;
            assert(in->iLower == pst->mdl()->Self());
            in->iLower = pst->upper();
            in->nProcessors -= nLeft;
            in->iIndex += nLeft;
            int rID = pst->mdl()->ReqService(pst->upper(), PST_WRITE, in, nIn);
            in->iLower = pst->mdl()->Self();
            in->iUpper = pst->upper();
            in->iIndex -= nLeft;
            in->nProcessors = nLeft;
            pstWrite(pst->lower(), in, nIn, NULL, 0);
            pst->mdl()->GetReply(rID);
        }
    }
    else {
        LCL* plcl = pst->lcl();
        if (in->nProcessors != 0) {
            auto achOutFile = makeName(in->achOutFile, in->iIndex);
            if (in->bHDF5) {
#ifdef USE_HDF5
                fio = fioHDF5Create(achOutFile.c_str(), in->mFlags);
#else
                fio = NULL; /* Should never happen */
#endif
            }
            else {
                if (achOutFile != makeName(in->achOutFile, in->iIndex + 1)) {
                    fio = fioTipsyCreatePart(achOutFile.c_str(), 0, in->mFlags & FIO_FLAG_CHECKPOINT, in->bStandard,
                                             pst->lcl()->pkd->particles.present(PKD_FIELD::oNewSph) ? in->dTime : in->dExp, in->nGas, in->nDark, in->nStar,
                                             plcl->nWriteStart);
                }
                else {
                    fio = fioTipsyAppend(achOutFile.c_str(), in->mFlags & FIO_FLAG_CHECKPOINT, in->bStandard);
                    if (fio) {
                        fioSeek(fio, plcl->nWriteStart, FIO_SPECIES_ALL);
                    }
                }
            }
            if (fio == NULL) {
                fprintf(stderr, "ERROR: unable to create file for output\n");
                perror(achOutFile.c_str());
                assert(fio != NULL);
            }

            fioSpeciesList nTotal, nThisFile{0};
            nTotal[FIO_SPECIES_DARK] = in->nDark;
            nTotal[FIO_SPECIES_SPH] = in->nGas;
            nTotal[FIO_SPECIES_STAR] = in->nStar;
            nTotal[FIO_SPECIES_BH] = in->nBH;

            pkdWriteFIO(plcl->pkd, fio, in->dvFac, in->dTuFac, in->bnd, nThisFile);
            for (i = in->iLower + 1; i < in->iUpper; ++i) {
                auto me = pst->Self();
                auto rID = pst->mdl()->ReqService<ServiceSendParticles>(i, me);
                pkdWriteFromNode(plcl->pkd, i, fio, in->dvFac, in->dTuFac, in->bnd, nThisFile);
                pst->mdl()->GetReply(rID);
            }

            pkdWriteHeaderFIO(plcl->pkd, fio, 1. / sqrt(in->dvFac), in->dTime, nTotal, nThisFile, in->dBoxSize, in->HubbleParam, in->nProcessorsTotal,
                              in->units);

            fioClose(fio);
        }
    }

    return 0;
}

int pstPhysicalSoft(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto in = static_cast<struct inPhysicalSoft*>(vin);

    assert(nIn == sizeof(struct inPhysicalSoft));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_PHYSICALSOFT, in, nIn);
        pstPhysicalSoft(pst->lower(), in, nIn, NULL, 0);
        pst->mdl()->GetReply(rID);
    }
    else {
        pkdPhysicalSoft(plcl->pkd, in->dSoftMax, in->dFac, in->bSoftMaxMul);
    }
    return 0;
}

int pstPlaceBHSeed(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inPlaceBHSeed*>(vin);
    auto out = static_cast<struct outPlaceBHSeed*>(vout);
    struct outPlaceBHSeed outUpper;
    assert(nIn == sizeof(struct inPlaceBHSeed));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_BH_PLACESEED, in, nIn);
        pstPlaceBHSeed(pst->lower(), vin, nIn, vout, nOut);
        nOut = pst->mdl()->GetReply(rID, outUpper);
        assert(nOut == sizeof(struct outPlaceBHSeed));
        out->nBHs += outUpper.nBHs;
    }
    else {
        LCL* plcl = pst->lcl();
        out->nBHs = pkdPlaceBHSeed(plcl->pkd, in->dTime, in->dScaleFactor, in->uRungMax, in->dDenMin, in->dBHMhaloMin, in->dTau, in->dBHSeedMass);
    }
    return sizeof(struct outPlaceBHSeed);
}
int pstBHInit(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inPlaceBHSeed*>(vin);
    assert(nIn == sizeof(struct inPlaceBHSeed));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_BH_INIT, in, nIn);
        pstBHInit(pst->lower(), vin, nIn, vout, nOut);
        pst->mdl()->GetReply(rID);
    }
    else {
        LCL* plcl = pst->lcl();
        pkdBHInit(plcl->pkd, in->uRungMax);
    }
    return 0;
}
int pstBHReposition(PST pst, void* vin, int nIn, void* vout, int nOut) {
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_BH_REPOSITION, NULL, 0);
        pstBHReposition(pst->lower(), vin, nIn, vout, nOut);
        pst->mdl()->GetReply(rID);
    }
    else {
        LCL* plcl = pst->lcl();
        pkdBHReposition(plcl->pkd);
    }
    return 0;
}
int pstBHAccretion(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inBHAccretion*>(vin);
    assert(nIn == sizeof(struct inBHAccretion));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_BH_ACCRETION, in, nIn);
        pstBHAccretion(pst->lower(), vin, nIn, vout, nOut);
        pst->mdl()->GetReply(rID);
    }
    else {
        LCL* plcl = pst->lcl();
        pkdBHAccretion(plcl->pkd, in->dScaleFactor);
    }
    return 0;
}

int pstSmooth(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inSmooth*>(vin);

    assert(nIn == sizeof(struct inSmooth));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_SMOOTH, in, nIn);
        pstSmooth(pst->lower(), in, nIn, NULL, 0);
        pst->mdl()->GetReply(rID);
    }
    else {
        LCL* plcl = pst->lcl();
        SMX smx;

        smInitialize(&smx, plcl->pkd, &in->smf, in->nSmooth, in->bPeriodic, in->bSymmetric, in->iSmoothType);
        smSmooth(smx, &in->smf);
        smFinish(smx, &in->smf);
    }
    return 0;
}

int pstReSmooth(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inSmooth*>(vin);
    auto out = static_cast<struct outSmooth*>(vout);
    struct outSmooth outUpper;

    //    assert(nIn == sizeof(struct inSmooth));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_RESMOOTH, in, nIn);
        pstReSmooth(pst->lower(), vin, nIn, vout, nOut);
        nOut = pst->mdl()->GetReply(rID, outUpper);
        assert(nOut == sizeof(struct outSmooth));
        out->nSmoothed += outUpper.nSmoothed;
    }
    else {
        LCL* plcl = pst->lcl();
        SMX smx;

        smInitialize(&smx, plcl->pkd, &in->smf, in->nSmooth, in->bPeriodic, in->bSymmetric, in->iSmoothType);
        out->nSmoothed = smReSmooth(smx, &in->smf, in->iSmoothType);
        smFinish(smx, &in->smf);
    }
    return sizeof(struct outSmooth);
}

#ifdef OPTIM_SMOOTH_NODE
int pstReSmoothNode(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = static_cast<struct inSmooth*>(vin);
    auto out = static_cast<struct outSmooth*>(vout);
    struct outSmooth outUpper;

    //    assert(nIn == sizeof(struct inSmooth));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_RESMOOTHNODE, in, nIn);
        pstReSmoothNode(pst->lower(), vin, nIn, vout, nOut);
        nOut = pst->mdl()->GetReply(rID, outUpper);
        assert(nOut == sizeof(struct outSmooth));
        out->nSmoothed += outUpper.nSmoothed;
    #if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER) && defined(DEBUG_FLUX_INFO)
        struct outSmooth tmp = *out;
        pstCombStat(&out->sComputing, &tmp.sComputing);
        pstCombStat(&out->sWaiting, &tmp.sWaiting);
        pstCombStat(&out->sSynchronizing, &tmp.sSynchronizing);
        pstCombStat(&out->sPartNumAccess, &tmp.sPartNumAccess);
        pstCombStat(&out->sPartMissRatio, &tmp.sPartMissRatio);
        pstCombStat(&out->sCellNumAccess, &tmp.sCellNumAccess);
        pstCombStat(&out->sCellMissRatio, &tmp.sCellMissRatio);
    #endif
    }
    else {
        LCL* plcl = pst->lcl();
        SMX smx;

        smInitialize(&smx, plcl->pkd, &in->smf, in->nSmooth, in->bPeriodic, in->bSymmetric, in->iSmoothType);
        if (in->iSmoothType == SMX_BH_MERGER) {
            out->nSmoothed = smReSmoothBHNode(smx, &in->smf, in->iSmoothType);
        }
        else {
            out->nSmoothed = smReSmoothNode(smx, &in->smf, in->iSmoothType);
        }
        smFinish(smx, &in->smf);
    #if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER) && defined(DEBUG_FLUX_INFO)
        if (out->nSmoothed) {
            out->sCellNumAccess.dSum = mdlNumAccess(pst->mdl(), CID_CELL) / out->nSmoothed;
            out->sPartNumAccess.dSum = mdlNumAccess(pst->mdl(), CID_PARTICLE) / out->nSmoothed;
        }
        else {
            out->sCellNumAccess.dSum = 0;
            out->sPartNumAccess.dSum = 0;
        }
        out->sCellMissRatio.dSum = 100.0 * mdlMissRatio(pst->mdl(), CID_CELL);     /* as a percentage */
        out->sPartMissRatio.dSum = 100.0 * mdlMissRatio(pst->mdl(), CID_PARTICLE); /* as a percentage */
        out->sComputing.dSum = pst->mdl()->TimeComputing();
        out->sWaiting.dSum = pst->mdl()->TimeWaiting();
        out->sSynchronizing.dSum = pst->mdl()->TimeSynchronizing();
        pstInitStat(&out->sComputing, pst->Self());
        pstInitStat(&out->sWaiting, pst->Self());
        pstInitStat(&out->sSynchronizing, pst->Self());
        pstInitStat(&out->sPartNumAccess, pst->Self());
        pstInitStat(&out->sPartMissRatio, pst->Self());
        pstInitStat(&out->sCellNumAccess, pst->Self());
        pstInitStat(&out->sCellMissRatio, pst->Self());
    #endif
    }
    return sizeof(struct outSmooth);
}
#endif

void pstInitStat(STAT* ps, int id) {
    ps->idMax = id;
    ps->dMax = ps->dSum;
    ps->dSum2 = ps->dSum * ps->dSum;
    ps->n = 1; /* sometimes we will need to zero this */
}

void pstCombStat(STAT* ps, STAT* pa) {
    if (pa->dMax > ps->dMax) {
        ps->idMax = pa->idMax;
        ps->dMax = pa->dMax;
    }
    ps->dSum += pa->dSum;
    ps->dSum2 += pa->dSum2;
    ps->n += pa->n;
}

int pstGravity(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto in = static_cast<struct inGravity*>(vin);
    auto outr = static_cast<struct outGravityReduct*>(vout);
    struct outGravityReduct tmp;
    int i;

    assert(nIn == sizeof(struct inGravity));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_GRAVITY, in, nIn);
        pstGravity(pst->lower(), in, nIn, vout, nOut);
        /*
        ** Make a temporary copy of the reduct part of the out buffer before setting it as the
        ** reply buffer. The reduct part follows at the end of all the outGravityPerProc entries.
        */
        pst->mdl()->GetReply(rID, tmp);
        /*
        ** Now combine in the tempory copy of the lower branch reduct part.
        */
#ifdef __linux__
        pstCombStat(&outr->sRSS, &tmp.sRSS);
        pstCombStat(&outr->sFreeMemory, &tmp.sFreeMemory);
#endif
        pstCombStat(&outr->sLocal, &tmp.sLocal);
        pstCombStat(&outr->sActive, &tmp.sActive);
        pstCombStat(&outr->sPart, &tmp.sPart);
        pstCombStat(&outr->sPartNumAccess, &tmp.sPartNumAccess);
        pstCombStat(&outr->sPartMissRatio, &tmp.sPartMissRatio);
        pstCombStat(&outr->sCell, &tmp.sCell);
        pstCombStat(&outr->sCellNumAccess, &tmp.sCellNumAccess);
        pstCombStat(&outr->sCellMissRatio, &tmp.sCellMissRatio);
        pstCombStat(&outr->sFlop, &tmp.sFlop);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
        pstCombStat(&outr->sComputing, &tmp.sComputing);
        pstCombStat(&outr->sWaiting, &tmp.sWaiting);
        pstCombStat(&outr->sSynchronizing, &tmp.sSynchronizing);
#endif
        /*
        ** Combine the rung counts for the next timestep...
        */
        for (i = in->ts.uRungLo; i <= IRUNGMAX; ++i) outr->nRung[i] += tmp.nRung[i];
        /*
        ** and the number of actives processed in this gravity call.
        */
        outr->nActive += tmp.nActive;
        outr->dFlopSingleCPU += tmp.dFlopSingleCPU;
        outr->dFlopDoubleCPU += tmp.dFlopDoubleCPU;
        outr->dFlopSingleGPU += tmp.dFlopSingleGPU;
        outr->dFlopDoubleGPU += tmp.dFlopDoubleGPU;
        outr->nTilesTotal += tmp.nTilesTotal;
        outr->nTilesCPU += tmp.nTilesCPU;
    }
    else {
#ifdef __linux__
        FILE* fp;
        char buffer[512], *save, *f, *v;
#endif
        PKD pkd = plcl->pkd;
        pkdGravAll(pkd, &in->kick, &in->lc, &in->ts, in->dTime, in->nReps, in->bPeriodic, in->bGPU, in->bEwald, in->iRoot1, in->iRoot2, in->dEwCut, in->dEwhCut,
                   in->dTheta, &in->SPHoptions, &outr->nActive, &outr->sPart.dSum, &outr->sPartNumAccess.dSum, &outr->sPartMissRatio.dSum, &outr->sCell.dSum,
                   &outr->sCellNumAccess.dSum, &outr->sCellMissRatio.dSum, &outr->sFlop.dSum, outr->nRung);
        outr->dFlopSingleCPU = 1e-9 * pkd->dFlopSingleCPU;
        outr->dFlopDoubleCPU = 1e-9 * pkd->dFlopDoubleCPU;
        outr->dFlopSingleGPU = 1e-9 * pkd->dFlopSingleGPU;
        outr->dFlopDoubleGPU = 1e-9 * pkd->dFlopDoubleGPU;
        outr->nTilesTotal = pkd->nTilesTotal;
        outr->nTilesCPU = pkd->nTilesCPU;
        outr->sLocal.dSum = plcl->pkd->Local();
        outr->sActive.dSum = (double)outr->nActive;
#ifdef __linux__
        fp = fopen("/proc/self/stat", "r");
        if (fp != NULL) {
            if (fgets(buffer, sizeof(buffer), fp) == NULL) buffer[0] = '\0';
            fclose(fp);
            f = strtok_r(buffer, " ", &save);
            for (i = 0; i <= 36 && f; i++) {
                switch (i) {
                case 23:
    #ifdef HAVE_GETPAGESIZE
                    outr->sRSS.dSum = atol(f) * 1.0 * getpagesize() / 1024 / 1024 / 1024;
    #else
                    outr->sRSS.dSum = atol(f) * 1.0 / 1024 / 1024;
    #endif
                    break;
                default:
                    break;
                }
                f = strtok_r(NULL, " ", &save);
            }
        }
        else outr->sRSS.dSum = 0;

        fp = fopen("/proc/meminfo", "r");
        outr->sFreeMemory.dSum = 0;
        if (fp != NULL) {
            while (fgets(buffer, sizeof(buffer), fp)) {
                f = strtok_r(buffer, ":", &save);
                v = strtok_r(NULL, " ", &save);
                if (strcmp(f, "MemFree") == 0 || strcmp(f, "Buffers") == 0 || strcmp(f, "Cached") == 0) {
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
        pstInitStat(&outr->sRSS, pst->Self());
        pstInitStat(&outr->sFreeMemory, pst->Self());
#endif
        pstInitStat(&outr->sLocal, pst->Self());
        pstInitStat(&outr->sActive, pst->Self());
        pstInitStat(&outr->sPart, pst->Self());
        pstInitStat(&outr->sPartNumAccess, pst->Self());
        pstInitStat(&outr->sPartMissRatio, pst->Self());
        pstInitStat(&outr->sCell, pst->Self());
        pstInitStat(&outr->sCellNumAccess, pst->Self());
        pstInitStat(&outr->sCellMissRatio, pst->Self());
        pstInitStat(&outr->sFlop, pst->Self());
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
        outr->sComputing.dSum = pst->mdl()->TimeComputing();
        outr->sWaiting.dSum = pst->mdl()->TimeWaiting();
        outr->sSynchronizing.dSum = pst->mdl()->TimeSynchronizing();
        pstInitStat(&outr->sComputing, pst->Self());
        pstInitStat(&outr->sWaiting, pst->Self());
        pstInitStat(&outr->sSynchronizing, pst->Self());
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

int pstCalcEandL(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto out = static_cast<struct outCalcEandL*>(vout);
    struct outCalcEandL outLcl;
    int k;

    assert(nIn == 0);
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_CALCEANDL, NULL, 0);
        pstCalcEandL(pst->lower(), NULL, 0, out, nOut);
        pst->mdl()->GetReply(rID, outLcl);
        out->T += outLcl.T;
        out->U += outLcl.U;
        out->Eth += outLcl.Eth;
        for (k = 0; k < 3; k++) out->L[k] += outLcl.L[k];
        for (k = 0; k < 3; k++) out->F[k] += outLcl.F[k];
        out->W += outLcl.W;
    }
    else {
        pkdCalcEandL(plcl->pkd, out->T, out->U, out->Eth, out->L, out->F, out->W);
    }
    return sizeof(struct outCalcEandL);
}

#ifdef OPTIM_REORDER_IN_NODES
int pstReorderWithinNodes(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();

    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_REORDERINNODES, NULL, 0);
        pstReorderWithinNodes(pst->lower(), NULL, 0, NULL, 0);
        pst->mdl()->GetReply(rID);
    }
    else {
        plcl->pkd->ReorderWithinNodes();
    }
    return 0;
}
#endif

#ifdef GRACKLE
int pstGrackleInit(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();

    auto in = static_cast<struct inGrackleInit*>(vin);
    assert(nIn == sizeof(struct inGrackleInit));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_GRACKLEINIT, in, nIn);
        pstGrackleInit(pst->lower(), in, nIn, NULL, 0);
        pst->mdl()->GetReply(rID);
    }
    else {
        pkdGrackleInit(plcl->pkd, in->bComove, in->dScaleFactor, in->achCoolingTable, in->units);
    }
    return 0;
}
#endif

int pstChemCompInit(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto in = static_cast<struct inChemCompInit*>(vin);

    assert(nIn == sizeof(struct inChemCompInit));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_CHEMCOMPINIT, in, nIn);
        pstChemCompInit(pst->lower(), in, nIn, NULL, 0);
        pst->mdl()->GetReply(rID);
    }
    else {
        pkdChemCompInit(plcl->pkd, *in);
    }
    return 0;
}

int pstMemStatus(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto out = static_cast<struct outMemStatus*>(vout);
    auto iUpper = pst->upper() - pst->Self();
    auto nLower = iUpper * sizeof(struct outMemStatus);
    struct outMemStatus* outUp = out + iUpper;

    assert(nIn == 0);
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_MEMSTATUS, vin, nIn);
        pstMemStatus(pst->lower(), vin, nIn, out, nLower);
        pst->mdl()->GetReply(rID, nOut - nLower, outUp);
    }
    else {
#ifdef __linux__
        FILE* fp;
        char buffer[512], *save, *f, *v;
        int i;

        fp = fopen("/proc/self/stat", "r");
        if (fp != NULL) {
            if (fgets(buffer, sizeof(buffer), fp) == NULL) buffer[0] = '\0';
            fclose(fp);
            f = strtok_r(buffer, " ", &save);
            for (i = 0; i <= 36 && f; i++) {
                switch (i) {
                case 23:
    #ifdef HAVE_GETPAGESIZE
                    out->rss = atol(f) * getpagesize() / 1024 / 1024;
    #else
                    out->rss = atol(f) / 1024;
    #endif
                    break;
                default:
                    break;
                }
                f = strtok_r(NULL, " ", &save);
            }
        }
        out->freeMemory = 0;
        fp = fopen("/proc/meminfo", "r");
        if (fp != NULL) {
            while (fgets(buffer, sizeof(buffer), fp)) {
                f = strtok_r(buffer, ":", &save);
                v = strtok_r(NULL, " ", &save);
                if (strcmp(f, "MemFree") == 0 || strcmp(f, "Buffers") == 0 || strcmp(f, "Cached") == 0) {
                    out->freeMemory += atol(v) / 1024;
                }
            }
            fclose(fp);
        }
#endif
        out->nBytesTree = plcl->pkd->tree.memory();
        out->nBytesCl = pkdClMemory(plcl->pkd);
        out->nBytesIlp = pkdIlpMemory(plcl->pkd);
        out->nBytesIlc = pkdIlcMemory(plcl->pkd);
    }
    return pst->n_total() * sizeof(struct outMemStatus);
}

int pstSetSPHoptions(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto in = static_cast<struct inSetSPHoptions*>(vin);

    assert(nIn == sizeof(struct inSetSPHoptions));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_SETSPHOPTIONS, vin, nIn);
        pstSetSPHoptions(pst->lower(), vin, nIn, vout, nOut);
        pst->mdl()->GetReply(rID);
    }
    else {
        copySPHOptions(&in->SPHoptions, &plcl->pkd->SPHoptions);
    }
    return 0;
}

int pstInitializeEOS(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_INITIALIZEEOS, vin, nIn);
        pstInitializeEOS(pst->lower(), vin, nIn, vout, nOut);
        pst->mdl()->GetReply(rID);
    }
    else {
        pkdInitializeEOS(plcl->pkd);
    }
    return 0;
}
