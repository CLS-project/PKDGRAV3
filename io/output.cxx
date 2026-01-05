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

#include "pkd_config.h"
#include "output.h"
#include "io/iochunk.h"
#include "pst.h"
#include <complex>
#include <cstddef>
#define COMPLEX std::complex<float>

#include "fmt/format.h"
using namespace fmt::literals; // Gives us ""_a and ""_format literals

void ServiceSaveGroupStats::Write(PST pst, std::span<std::byte> vin, int iGroup, const std::string& filename, int iSegment, int nSegment) {
    pst->lcl()->pkd->GroupStatsSave(iGroup, filename, iSegment, nSegment);
}

void pkdContext::GroupStatsSave(int iGroup, const std::string& filename, int iSegment, int nSegment) {
    auto nBytes = nLocalGroups * sizeof(TinyGroupTable);
    io_chunk_write(filename.c_str(), tinyGroupTable + 1, nBytes, iSegment > 0);
}

/*
** We do not do the write, rather we send to another thread.
*/

template <typename F>
auto pack_grid(PKD pkd, MDLGRID grid, int iGrid, std::span<std::byte> buffer, size_t& iIndex) -> size_t {
    auto fftData = pkd->Ephemeral<F>() + iGrid * grid->nLocal; // Start of the grid data
    auto pOutput = reinterpret_cast<F*>(buffer.data());        // Output pointer
    size_t nLocal = 1ul * grid->a1 * grid->n2 * grid->nSlab;   // Total local size (including padding)
    size_t n = buffer.size() / sizeof(F);                      // Number of elements that fit in buffer
    size_t iOutput = 0;
    while (iIndex < nLocal && iOutput + grid->n1 <= n) {
        memcpy(pOutput + iOutput, fftData + iIndex, grid->n1 * sizeof(*fftData));
        iOutput += grid->n1;
        iIndex += grid->a1;
    }
    return iOutput * sizeof(F);
}

void pkdOutputSend(PKD pkd, outType eOutputType, int iPartner, int iGrid) {
    size_t iIndex = 0;
    switch (eOutputType) {
    case OUT_TINY_GROUP:
        pkd->mdl->Send(iPartner, [&](std::span<std::byte> buffer) -> size_t {
            auto nLeft = pkd->nLocalGroups - iIndex;
            auto n = buffer.size() / sizeof(TinyGroupTable);
            if (n > nLeft) n = nLeft;
            memcpy(buffer.data(), pkd->tinyGroupTable + 1 + iIndex, n * sizeof(TinyGroupTable));
            iIndex += n;
            return n * sizeof(TinyGroupTable);
        });
        break;
    case OUT_KGRID:
        pkd->mdl->Send(iPartner, [&](std::span<std::byte> buffer) -> size_t {
            if (pkd->mdl->Core() == 0) {
                return pack_grid<COMPLEX>(pkd, pkd->fft->kgrid, iGrid, buffer, iIndex);
            }
            else return 0;
        });
        break;
    case OUT_RGRID:
        pkd->mdl->Send(iPartner, [&](std::span<std::byte> buffer) -> size_t {
            if (pkd->mdl->Core() == 0) {
                return pack_grid<float>(pkd, pkd->fft->rgrid, iGrid, buffer, iIndex);
            }
            else return 0;
        });
        break;
    default:
        fprintf(stderr, "ERROR: invalid output type %d\n", eOutputType);
        abort();
    }
}

int pstOutputSend(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = reinterpret_cast<struct inOutputSend*>(vin);
    pkdOutputSend(pst->lcl()->pkd, in->eOutputType, in->iPartner, in->iGrid);
    return 0;
}

/*
** We are the writer. We may need to receive as well.
*/

template <typename F>
static void localWrite(PKD pkd, asyncFileInfo* info, MDLGRID grid, int iGrid) {
    const size_t SEND_BUFFER_SIZE = 64 * 1024;
    size_t iIndex = 0;
    auto vOut = new char[SEND_BUFFER_SIZE];
    while (auto n = pack_grid<F>(pkd, grid, iGrid, std::span<std::byte>(reinterpret_cast<std::byte*>(vOut), SEND_BUFFER_SIZE), iIndex)) {
        io_write(info, vOut, n);
    }
    delete[] vOut;
}

void pkdOutput(PKD pkd, outType eOutputType, int iProcessor, int nProcessor, int iPartner, int nPartner, const char* fname, int iGrid) {
    asyncFileInfo info;
    auto achOutFile = fmt::format(fname, "i"_a = iProcessor);
    io_init(&info, 4, 1024 * 1024, IO_AIO | IO_LIBAIO);
    if (io_create(&info, achOutFile.c_str()) < 0) {
        perror(fname);
        abort();
    }

    switch (eOutputType) {
    case OUT_TINY_GROUP:
        io_write(&info, pkd->tinyGroupTable + 1, sizeof(TinyGroupTable) * pkd->nLocalGroups);
        break;
    case OUT_KGRID:
        localWrite<COMPLEX>(pkd, &info, pkd->fft->kgrid, iGrid);
        break;
    case OUT_RGRID:
        localWrite<float>(pkd, &info, pkd->fft->rgrid, iGrid);
        break;
    default:
        fprintf(stderr, "ERROR: invalid output type %d\n", eOutputType);
        abort();
    }
    while (--nPartner) {
        struct inOutputSend send;
        send.iPartner = pkd->Self();
        send.eOutputType = eOutputType;
        send.iGrid = iGrid;
        ++iPartner;
        int rID = pkd->mdl->ReqService(iPartner, PST_OUTPUT_SEND, &send, sizeof(send));
        pkd->mdl->Recv(iPartner, [&](std::span<std::byte> buffer) -> size_t {
            if (buffer.empty()) return 0;
            io_write(&info, buffer.data(), buffer.size());
            return 1;
        });

        pkd->mdl->GetReply(rID);
    }
    io_close(&info);
    io_free(&info);
}

int pstOutput(PST pst, void* vin, int nIn, void* vout, int nOut) {
    auto in = reinterpret_cast<struct inOutput*>(vin);

    assert(nIn >= sizeof(struct inOutput));
    if (pst->NotCore()) {
        int nProcessor = in->nProcessor;
        int iProcessor = in->iProcessor;

        /* Still allowed to write more in parallel */
        if (nProcessor > 1) {
            int nLower, nUpper;
            nLower = nProcessor * pst->n_lower() / pst->n_total();
            if (nLower == 0) nLower = 1;
            nUpper = nProcessor - nLower;
            in->nProcessor = nUpper;
            in->iProcessor = iProcessor + nLower;
            int rID = pst->mdl()->ReqService(pst->upper(), PST_OUTPUT, in, nIn);
            in->nProcessor = nLower;
            in->iProcessor = iProcessor;
            pstOutput(pst->lower(), in, nIn, NULL, 0);
            pst->mdl()->GetReply(rID);
        }
        /* We are the node that will be the writer for all of the pst children */
        else if (nProcessor == 1) {
            in->iPartner = pst->Self();
            in->nPartner = pst->n_total();
            in->nProcessor = 0;
            pstOutput(pst->lower(), in, nIn, NULL, 0); /* Keep decending to write */
        }
        else {
            pstOutput(pst->lower(), in, nIn, NULL, 0); /* Keep decending to write */
        }
    }
    else {
        /* If it is fully parallel then there is just us writing. */
        if (in->nProcessor > 0) {
            in->iPartner = pst->Self();
            in->nPartner = 1;
        }
        PKD pkd = pst->lcl()->pkd;
        pkdOutput(pkd, in->eOutputType, in->iProcessor, in->nProcessor, in->iPartner, in->nPartner, in->achOutFile, in->iGrid);
    }
    return 0;
}

auto ServiceSendParticles::Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int {
    auto iTo = *reinterpret_cast<const uint32_t*>(in.data());
    pst->lcl()->pkd->WriteViaNode(iTo);
    return 0;
}

/*****************************************************************************\
* Send particles to be written
\*****************************************************************************/

void pkdContext::WriteViaNode(int iNode) {
#ifdef MPI_VERSION
    size_t iIndex = 0;
    size_t iDecoratorIndex = 0;

    mdl->Send(iNode, [&](std::span<std::byte> buffer) -> size_t {
        int nDecorator = 0;
        int nPart = 0;
        size_t i = iIndex;

        // Header: [nPart, nDecorator]
        size_t n = 2 * sizeof(uint64_t);

        const size_t pSize = particles.ParticleSize();
        const size_t dSize = decorators.DecoratorSize();

        while (n < buffer.size() - (pSize + dSize) && i < Local()) {
            n += pSize;
            if (particles[i].have_decorator()) {
                n += dSize;
                ++nDecorator;
            }
            ++nPart;
            ++i;
        }

        if (nPart == 0) return size_t{0};

        const size_t nPartSize = nPart * pSize;
        const size_t nDecoratorSize = nDecorator * dSize;

        auto* header = reinterpret_cast<uint64_t*>(buffer.data());
        header[0] = static_cast<uint64_t>(nPart);
        header[1] = static_cast<uint64_t>(nDecorator);
        std::memcpy(buffer.data() + 2 * sizeof(uint64_t), Particle(iIndex), nPartSize);
        std::memcpy(buffer.data() + 2 * sizeof(uint64_t) + nPartSize, Decorator(iDecoratorIndex), nDecoratorSize);

        iIndex += nPart;
        iDecoratorIndex += nDecorator;

        return n;
    });
#endif
}
