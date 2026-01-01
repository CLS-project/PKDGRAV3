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
#include "fftgrid.h"

// Ensure the communication structures are "standard" so that they can be moved around with "memcpy" (required by MDL)
static_assert(std::is_void_v<ServiceFftSizes::input> || std::is_standard_layout_v<ServiceFftSizes::input>);
static_assert(std::is_void_v<ServiceFftSizes::output> || std::is_standard_layout_v<ServiceFftSizes::output>);

int ServiceFftSizes::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == sizeof(output));
    out->nMaxLocal = pst->mdl()->FFTlocalCount(in->nx, in->ny, in->nz, &out->nMaxZ, 0, &out->nMaxY, 0);
    return sizeof(output);
}

auto ServiceGridCreateFFT::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == sizeof(output));
    auto pkd = pst->lcl()->pkd;
    auto [nMaxLocal, rs, rn, ks, kn] = pkd->GridCreateFFT(in->nGrid);
    *out = {.nMaxLocal = nMaxLocal, .nMaxZ = rn, .nMaxY = kn};
    return sizeof(output);
}

auto pkdContext::GridCreateFFT(int nGrid) -> std::tuple<size_t, int, int, int, int> {
    assert(fft == nullptr);
    fft = mdl->FFTInitialize(nGrid, nGrid, nGrid, 0, 0);
    return {fft->rgrid->nLocal, fft->rgrid->sSlab, fft->rgrid->nSlab, fft->kgrid->sSlab, fft->kgrid->nSlab};
}

auto ServiceGridDeleteFFT::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    assert(vin.size() == 0);
    assert(vout.size() == 0);
    pst->lcl()->pkd->GridDeleteFFT();
    return 0;
}

void pkdContext::GridDeleteFFT() {
    assert(fft != nullptr);
    mdl->FFTFinish(fft);
    fft = nullptr;
}
