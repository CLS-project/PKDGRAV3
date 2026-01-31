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
#include "measurepk.h"
#include "core/gridinfo.hpp"
#include "ic/whitenoise.hpp"
using namespace gridinfo;

class LinearSignal : public NoiseGenerator {
  private:
    CSM csm;
    gsl_interp_accel* acc;
    gsl_spline* spline;
    double iLbox;
    double *logk, *field;

  protected:
    virtual void update(gridinfo::complex_vector_t& pencil, gridinfo::complex_vector_t& noise, int j, int k);

  public:
    explicit LinearSignal(CSM csm, double a, double Lbox, unsigned long seed, bool bFixed = false, float fPhase = 0);
    virtual ~LinearSignal();
};

LinearSignal::LinearSignal(CSM csm, double a, double Lbox, unsigned long seed, bool bFixed, float fPhase) : NoiseGenerator(seed, bFixed, fPhase) {
    auto size = csm->val.classData.perturbations.size_k;
    this->csm = csm;
    iLbox = 2 * M_PI / Lbox;
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    logk = (double*)malloc(sizeof(double) * size);
    field = (double*)malloc(sizeof(double) * size);
    for (auto i = 0; i < size; i++) {
        auto k = csm->val.classData.perturbations.k[i];
        logk[i] = log(k);
        field[i] = csmDeltaRho_pk(csm, a, k);
        field[i] /= csmZeta(csm, k);
        assert(!std::isnan(field[i]));
    }
    gsl_spline_init(spline, logk, field, size);
}

LinearSignal::~LinearSignal() {
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
}

void LinearSignal::update(complex_vector_t& pencil, complex_vector_t& noise, int iy, int iz) {
    float k2jk = iy * iy + iz * iz;
    for (auto index = noise.begin(); index != noise.end(); ++index) {
        auto ix = index.position()[0];
        auto k = sqrt(k2jk + ix * ix) * iLbox;
        if (k > 0) {
            float signal = csmZeta(csm, k) * gsl_spline_eval(spline, log(k), acc);
            auto wnoise = *index;
            pencil(index.position()) += wnoise * signal;
        }
    }
}

auto ServiceAddLinearSignal::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == 0);
    pkd->AddLinearSignal(in->iGrid, in->iSeed, in->bFixed, in->fPhase, in->Lbox, in->a);
    return 0;
}

void pkdContext::AddLinearSignal(int iGrid, int iSeed, bool bFixed, float fPhase, double Lbox, double a) {
    assert(fft != nullptr);
    int nGrid = fft->rgrid->n1;
    GridInfo G(mdl, fft);

    complex_array_t K1;
    auto data1 = mdl->SetArray(Ephemeral<real_t>()) + fft->rgrid->nLocal * iGrid;
    G.setupArray(data1, K1);

    LinearSignal ng(csm, a, Lbox, iSeed, bFixed, fPhase);
    ng.FillNoise(K1, nGrid);
}

auto ServiceGridBinK::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == sizeof(output));
    pkd->GridBinK(in->nBins, in->iGrid, std::span{out->fK}, std::span{out->fPower}, std::span{out->nPower});
    return sizeof(output);
}

auto ServiceGridBinK::Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte> vin, PST) -> int {
    auto out1 = reinterpret_cast<output*>(vout1.data());
    auto out2 = reinterpret_cast<const output*>(vout2.data());
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    assert(vout1.size() == sizeof(output));
    assert(vout2.size() == sizeof(output));
    for (auto i = 0; i < in->nBins; i++) {
        out1->fK[i] += out2->fK[i];
        out1->fPower[i] += out2->fPower[i];
        out1->nPower[i] += out2->nPower[i];
    }
    return sizeof(output);
}

void pkdContext::GridBinK(int nBins, int iGrid, std::span<double> fK, std::span<double> fPower, std::span<uint64_t> nPower) {
    assert(fft != NULL);
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;
    GridInfo G(mdl, fft);

    complex_array_t K1;
    auto data1 = mdl->SetArray(Ephemeral<real_t>()) + fft->rgrid->nLocal * iGrid;
    G.setupArray(data1, K1);

    for (auto i = 0; i < nBins; i++) {
        fK[i] = 0.0;
        fPower[i] = 0.0;
        nPower[i] = 0;
    }
#ifdef LINEAR_PK
    double scale = nBins * 1.0 / iNyquist;
#else
    double scale = nBins * 1.0 / log(iNyquist + 1);
#endif
    for (auto index = K1.begin(); index != K1.end(); ++index) {
        auto pos = index.position();
        auto i = pos[0]; // i,j,k are all positive (absolute value)
        auto j = pos[1] > iNyquist ? nGrid - pos[1] : pos[1];
        auto k = pos[2] > iNyquist ? nGrid - pos[2] : pos[2];
        auto v1 = *index;

        auto ak = sqrt(i * i + j * j + k * k);
        auto ks = int(ak);
        if (ks >= 1 && ks <= iNyquist) {
#ifdef LINEAR_PK
            ks = floor((ks - 1.0) * scale);
#else
            ks = floor(log(ks) * scale);
#endif
            assert(ks >= 0 && ks < nBins);
            fK[ks] += log(ak);
            fPower[ks] += std::norm(v1);
            nPower[ks] += 1;
            if (i != 0 && i != iNyquist) { // Account for negative Kx values
                fK[ks] += log(ak);
                fPower[ks] += std::norm(v1);
                nPower[ks] += 1;
            }
        }
    }
}

auto ServiceDensityContrast::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    pkd->DensityContrast(in->dTotalMass, in->iGrid, in->k);
    return 0;
}

void pkdContext::DensityContrast(double dTotalMass, int iGrid, bool k) {
    GridInfo G(mdl, fft);
    int nGrid = fft->rgrid->n1;
    auto data = mdl->SetArray(Ephemeral<real_t>()) + fft->rgrid->nLocal * iGrid;
    real_array_t R;
    G.setupArray(data, R);
    if (k) { // Delta(k): same as Delta(r) below, but apply normalization for the FFT (divide by nGrid^3)
        real_t diTotalMass = 1.0 / dTotalMass;
        real_t fftNormalize = 1.0 / (1.0 * nGrid * nGrid * nGrid);
        R = R * diTotalMass - fftNormalize;
        mdl->FFT(fft, data);
    }
    else { // Delta(r) = rho / rho_bar - 1.0
        real_t diRhoBar = (1.0 * nGrid * nGrid * nGrid) / dTotalMass;
        R = R * diRhoBar - 1.0;
    }
}

auto ServiceInterlace::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    pkd->Interlace(in->iGridTarget, in->iGridSource);
    return 0;
}

void pkdContext::Interlace(int iGridTarget, int iGridSource) {
    auto ephemeral = Ephemeral<real_t>();
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;

    GridInfo G(mdl, fft);

    complex_array_t K1, K2;
    auto data1 = mdl->SetArray(ephemeral) + fft->rgrid->nLocal * iGridTarget;
    auto data2 = mdl->SetArray(ephemeral) + fft->rgrid->nLocal * iGridSource;
    G.setupArray(data1, K1);
    G.setupArray(data2, K2);

    for (auto index = K1.begin(); index != K1.end(); ++index) {
        const auto v1 = *index;
        auto pos = index.position();
        auto i = pos[0]; // i is positive; j and k can be negative
        auto j = pos[1] > iNyquist ? pos[1] - nGrid : pos[1];
        auto k = pos[2] > iNyquist ? pos[2] - nGrid : pos[2];
        float theta = M_PI / nGrid * (i + j + k);
        auto v2 = K2(index.position()) * complex_t(cosf(theta), sinf(theta));
        *index = complex_t(0.5) * (v1 + v2);
    }
}
