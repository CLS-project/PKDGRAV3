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
#include "bispectrum.h"
#include "core/gridinfo.hpp"
using namespace gridinfo;

auto ServiceBispectrumSelect::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pkd = pst->pkd();
    if (in->iGridSource < 0) pkd->BispectrumNormalize(in->iGridTarget, in->kmin, in->kmax);
    else pkd->BispectrumSelect(in->iGridTarget, in->iGridSource, in->kmin, in->kmax);
    return 0;
}

void pkdContext::BispectrumNormalize(int iGridTarget, double kmin, double kmax) {
    auto ephemeral = Ephemeral<real_t>();
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;

    GridInfo G(mdl, fft);

    complex_array_t K1;
    auto data1 = mdl->SetArray(ephemeral, 0) + fft->rgrid->nLocal * iGridTarget;
    G.setupArray(data1, K1);
    double k2min = kmin * kmin;
    double k2max = kmax * kmax;
    for (auto index = K1.begin(); index != K1.end(); ++index) {
        auto pos = index.position();
        auto i = pos[0]; // i,j,k are all positive (absolute value)
        auto j = pos[1] > iNyquist ? nGrid - pos[1] : pos[1];
        auto k = pos[2] > iNyquist ? nGrid - pos[2] : pos[2];
        auto k2 = i * i + j * j + k * k;
        *index = (k2 >= k2min && k2 < k2max) ? 1.0 : 0.0;
    }
    mdl->IFFT(fft, (FFTW3(complex)*)K1.dataFirst());
}

void pkdContext::BispectrumSelect(int iGridTarget, int iGridSource, double kmin, double kmax) {
    auto ephemeral = Ephemeral<real_t>();
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;

    GridInfo G(mdl, fft);

    complex_array_t K1, K2;
    auto data1 = mdl->SetArray(ephemeral) + fft->rgrid->nLocal * iGridTarget;
    auto data2 = mdl->SetArray(ephemeral) + fft->rgrid->nLocal * iGridSource;
    G.setupArray(data1, K1);
    G.setupArray(data2, K2);

    double k2min = kmin * kmin;
    double k2max = kmax * kmax;
    for (auto index = K1.begin(); index != K1.end(); ++index) {
        auto pos = index.position();
        auto i = pos[0]; // i,j,k are all positive (absolute value)
        auto j = pos[1] > iNyquist ? nGrid - pos[1] : pos[1];
        auto k = pos[2] > iNyquist ? nGrid - pos[2] : pos[2];
        auto k2 = i * i + j * j + k * k;
        if (k2 < k2min || k2 >= k2max) *index = 0.0;
        else *index = K2(index.position());
    }
    mdl->IFFT(fft, (FFTW3(complex)*)K1.dataFirst());
}

auto ServiceBispectrumCalculate::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->pkd();
    *out = pkd->BispectrumCalculate(in->iGrid1, in->iGrid2, in->iGrid3);
    return sizeof(output);
}

auto pkdContext::BispectrumCalculate(int iGrid1, int iGrid2, int iGrid3) -> double {
    auto ephemeral = Ephemeral<real_t>();
    GridInfo G(mdl, fft);
    real_array_t R1, R2, R3;
    auto data1 = mdl->SetArray(ephemeral) + fft->rgrid->nLocal * iGrid1;
    auto data2 = mdl->SetArray(ephemeral) + fft->rgrid->nLocal * iGrid2;
    auto data3 = mdl->SetArray(ephemeral) + fft->rgrid->nLocal * iGrid3;
    G.setupArray(data1, R1);
    G.setupArray(data2, R2);
    G.setupArray(data3, R3);
    return sum(R1 * R2 * R3);
}
