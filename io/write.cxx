// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "write.h"
#include "io/service.h"
#include "io/output.h"
#include "fmt/format.h"

auto ServiceWrite::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<input*>(vin.data());
    /* If we are the writer (nProcessors==1) or a producer (==1) keep requesting */
    auto nProcessors = in->nProcessors;
    if (nProcessors <= 1) {
        Traverse(pst->lower(), vin, vout);
    }
    /* Split writing tasks between child nodes */
    else {
        int nLeft = nProcessors / 2;
        assert(in->iLower == pst->mdl()->Self());
        in->iLower = pst->upper();
        in->nProcessors -= nLeft;
        in->iIndex += nLeft;
        auto rID = ReqService(pst, vin);
        in->iLower = pst->mdl()->Self();
        in->iUpper = pst->upper();
        in->iIndex -= nLeft;
        in->nProcessors = nLeft;
        Traverse(pst->lower(), vin, vout);
        pst->mdl()->GetReply(rID, vout);
    }
    return 0;
}

static auto makeName(std::string_view inName, int iIndex) -> std::string {
    return fmt::format(inName, fmt::arg("i", iIndex));
}

auto ServiceWrite::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<input*>(vin.data());
    auto pkd = pst->pkd();

    LCL* plcl = pst->lcl();
    if (in->nProcessors != 0) {
        FIO fio;

        auto achOutFile = makeName(in->achOutFile, in->iIndex);
        if (in->bHDF5) {
#ifdef USE_HDF5
            fio = fioHDF5Create(achOutFile.c_str(), in->mFlags);
#else
            fio = NULL; /* Should never happen */
#endif
        }
        else {
            auto nWriteStart = plcl->data<ServiceSetWriteStart::InfoLCL>().nWriteStart;
            if (achOutFile != makeName(in->achOutFile, in->iIndex + 1)) {
                fio = fioTipsyCreatePart(achOutFile.c_str(), 0, in->mFlags & FIO_FLAG_CHECKPOINT, in->bStandard,
                                         (pkd->particles.present(PKD_FIELD::oNewSph) || pkd->particles.present(PKD_FIELD::oSSDEM)) ? in->dTime : in->dExp,
                                         in->nGas, in->nDark, in->nStar, nWriteStart);
            }
            else {
                fio = fioTipsyAppend(achOutFile.c_str(), in->mFlags & FIO_FLAG_CHECKPOINT, in->bStandard);
                if (fio) {
                    fioSeek(fio, nWriteStart, FIO_SPECIES_ALL);
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

        pkdWriteFIO(pkd, fio, in->dvFac, in->dTuFac, in->bnd, nThisFile);
        for (auto i = in->iLower + 1; i < in->iUpper; ++i) {
            auto me = pst->Self();
            auto rID = pst->mdl()->ReqService<ServiceSendParticles>(i, me);
            pkdWriteFromNode(pkd, i, fio, in->dvFac, in->dTuFac, in->bnd, nThisFile);
            pst->mdl()->GetReply(rID);
        }

        pkdWriteHeaderFIO(pkd, fio, 1. / sqrt(in->dvFac), in->dTime, nTotal, nThisFile, in->dBoxSize, in->HubbleParam, in->nProcessorsTotal, in->units);
        fioClose(fio);
    }

    return 0;
}
