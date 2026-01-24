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
    auto pkd = pst->lcl()->pkd;
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

void pkdGridBinK(PKD pkd, int nBins, int iGrid, double* fK, double* fPower, uint64_t* nPower) {
    assert(pkd->fft != NULL);
    auto fft = pkd->fft;
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;
    GridInfo G(pkd->mdl, fft);

    complex_array_t K1;
    auto data1 = pkd->mdl->SetArray(pkd->Ephemeral<real_t>()) + fft->rgrid->nLocal * iGrid;
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

int pstGridBinK(PST pst, void* vin, int nIn, void* vout, int nOut) {
    LCL* plcl = pst->lcl();
    auto in = reinterpret_cast<struct inGridBinK*>(vin);
    auto out = reinterpret_cast<struct outGridBinK*>(vout);
    int i;

    assert(nIn == sizeof(struct inGridBinK));
    assert(nOut == sizeof(struct outGridBinK));
    if (pst->NotCore()) {
        auto outUpper = new struct outGridBinK;
        int rID = pst->mdl()->ReqService(pst->upper(), PST_GRID_BIN_K, vin, nIn);
        pstGridBinK(pst->lower(), vin, nIn, vout, nOut);
        nOut = pst->mdl()->GetReply(rID, sizeof(*outUpper), outUpper);
        assert(nOut == sizeof(struct outGridBinK));

        for (i = 0; i < in->nBins; i++) {
            out->fK[i] += outUpper->fK[i];
            out->fPower[i] += outUpper->fPower[i];
            out->nPower[i] += outUpper->nPower[i];
        }
        delete outUpper;
    }
    else {
        pkdGridBinK(plcl->pkd, in->nBins, in->iGrid, out->fK, out->fPower, out->nPower);
    }
    return sizeof(struct outGridBinK);
}
