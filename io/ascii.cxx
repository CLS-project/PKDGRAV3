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
#include "ascii.h"

auto ServiceAddWriteStart::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<input*>(vin.data());
    pst->lcl()->nWriteStart += *in;
    return 0;
}

auto ServiceCompressASCII::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto plcl = pst->lcl();
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto& pkdout = plcl->data_emplace<PKDOUT>(plcl->pkd, in->iFile, in->iType);
    if (pst->Self() == 0 && in->iDim == 0) {
        pkdout.OutHdr(in->nTotal);
    }
    pkdout.OutASCII(in->iType, in->iDim);
    pkdout.FinishOutASCII();
    *out = pkdout.CountOutASCII();
    plcl->nWriteStart = 0;
    return sizeof(output);
}

auto ServiceCompressASCII::Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST pst) -> int {
    auto out1 = reinterpret_cast<output*>(vout1.data());
    auto out2 = reinterpret_cast<const output*>(vout2.data());
    mdl::mdlClass* mdl = pst->mdl();
    auto rID = mdl->ReqService<ServiceAddWriteStart>(pst->upper(), vout1);
    mdl->GetReply(rID);
    *out1 += *out2;
    return sizeof(output);
}

auto ServiceWriteASCII::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto plcl = pst->lcl();
    auto& pkdout = plcl->data<PKDOUT>();
    auto in = reinterpret_cast<const input*>(vin.data());
    FILE* fp = fopen(in->achOutFile, "r+");

    mdlassert(plcl->pkd->mdl, fp != NULL);
    setvbuf(fp, NULL, _IOFBF, PKDOUT_BUFFER_SIZE);
    fseeko(fp, in->nFileOffset + plcl->nWriteStart, SEEK_SET);
    pkdout.DumpOutASCII(fp);
    plcl->data_erase<PKDOUT>();
    fclose(fp);
    return 0;
}
