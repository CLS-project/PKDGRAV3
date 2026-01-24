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

#include "hydro/hydro.h"
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
    mdl->AddService(PST_OUTPUT, pst, (fcnService_t*)pstOutput, sizeof(struct inOutput), 0, Legacy::Output);
    mdl->AddService(PST_OUTPUT_SEND, pst, (fcnService_t*)pstOutputSend, sizeof(struct inOutputSend), 0, Legacy::OutputSend);
    /*
    ** Calculate the number of levels in the top tree and use it to
    ** define the size of the messages.
    */
    mdl->AddService(PST_ONENODEREADINIT, pst, (fcnService_t*)pstOneNodeReadInit,
                    sizeof(struct inReadFile) + PST_MAX_FILES * (sizeof(fioSpeciesList) + PATH_MAX), nThreads * sizeof(fioSpeciesList),
                    Legacy::OneNodeReadInit);
#ifdef MDL_FFTW
    mdl->AddService(PLT_MOVEIC, pst, (fcnService_t*)pltMoveIC, sizeof(struct inMoveIC), 0, Legacy::MoveIC);
    mdl->AddService(PST_MOVEIC, pst, (fcnService_t*)pstMoveIC, sizeof(struct ServiceGenerateIC::input), 0, Legacy::LocalMoveIC);
#endif
    mdl->AddService(PST_SETSPHOPTIONS, pst, (fcnService_t*)pstSetSPHoptions, sizeof(struct inSetSPHoptions), 0, Legacy::SetSPHoptions);
    mdl->AddService(PST_INITIALIZEEOS, pst, (fcnService_t*)pstInitializeEOS, 0, 0, Legacy::InitializeEOS);
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
