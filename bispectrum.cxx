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
#include "master.h"
#include "core/gridinfo.hpp"
using namespace gridinfo;

void pkdBispectrumNormalize(PKD pkd, int iGridTarget, double kmin,double kmax) {
    auto fft = pkd->fft;
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;

    GridInfo G(pkd->mdl,fft);

    complex_array_t K1;
    auto data1 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGridTarget;
    G.setupArray(data1,K1);
    double k2min = kmin*kmin;
    double k2max = kmax*kmax;
    for ( auto index=K1.begin(); index!=K1.end(); ++index ) {
        auto pos = index.position();
        auto i = pos[0]; // i,j,k are all positive (absolute value)
        auto j = pos[1]>iNyquist ? nGrid - pos[1] : pos[1];
        auto k = pos[2]>iNyquist ? nGrid - pos[2] : pos[2];
        auto k2 = i*i + j*j + k*k;
        *index = (k2>=k2min && k2<k2max) ? 1.0 : 0.0;
    }
    mdlIFFT(pkd->mdl,fft,(FFTW3(complex) *)K1.dataFirst());
}

void pkdBispectrumSelect(PKD pkd, int iGridTarget, int iGridSource,double kmin,double kmax) {
    auto fft = pkd->fft;
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;

    GridInfo G(pkd->mdl,fft);

    complex_array_t K1,K2;
    auto data1 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGridTarget;
    auto data2 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGridSource;
    G.setupArray(data1,K1);
    G.setupArray(data2,K2);

    double k2min = kmin*kmin;
    double k2max = kmax*kmax;
    for ( auto index=K1.begin(); index!=K1.end(); ++index ) {
        auto pos = index.position();
        auto i = pos[0]; // i,j,k are all positive (absolute value)
        auto j = pos[1]>iNyquist ? nGrid - pos[1] : pos[1];
        auto k = pos[2]>iNyquist ? nGrid - pos[2] : pos[2];
        auto k2 = i*i + j*j + k*k;
        if (k2<k2min || k2 >=k2max) *index = 0.0;
        else *index = K2(index.position());
    }
    mdlIFFT(pkd->mdl,fft,(FFTW3(complex) *)K1.dataFirst());
}

int pstBispectrumSelect(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inBispectrumSelect *in = reinterpret_cast<struct inBispectrumSelect *>(vin);
    assert (nIn==sizeof(struct inBispectrumSelect) );
    if (pstNotCore(pst)) {
        int rID = mdlReqService(pst->mdl, pst->idUpper, PST_BISPECTRUM_SELECT, vin, nIn);
        pstBispectrumSelect(pst->pstLower, vin, nIn, NULL, 0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        if (in->iGridSource<0) pkdBispectrumNormalize(plcl->pkd,in->iGridTarget,in->kmin,in->kmax);
        else pkdBispectrumSelect(plcl->pkd,in->iGridTarget,in->iGridSource,in->kmin,in->kmax);
    }
    return 0;
}

void MSR::BispectrumSelect(int iGridTarget,int iGridSource,double kmin,double kmax) {
    struct inBispectrumSelect in;
    in.iGridTarget = iGridTarget;
    in.iGridSource = iGridSource;
    in.kmin = kmin;
    in.kmax = kmax;
    pstBispectrumSelect(pst, &in, sizeof(in), NULL, 0);
}

double pkdBispectrumCalculate(PKD pkd, int iGrid1, int iGrid2, int iGrid3 ) {
    auto fft = pkd->fft;
    GridInfo G(pkd->mdl,fft);
    real_array_t R1,R2,R3;
    auto data1 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGrid1;
    auto data2 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGrid2;
    auto data3 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGrid3;
    G.setupArray(data1,R1);
    G.setupArray(data2,R2);
    G.setupArray(data3,R3);
    return sum(R1 * R2 * R3);
}

int pstBispectrumCalculate(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inBispectrumCalculate *in = reinterpret_cast<struct inBispectrumCalculate *>(vin);
    auto out = reinterpret_cast<double *>(vout);
    assert (nIn==sizeof(struct inBispectrumCalculate) );
    assert (nOut==sizeof(double) );
    if (pstNotCore(pst)) {
        double outUpper;
        int nOutUpper;
        int rID = mdlReqService(pst->mdl, pst->idUpper, PST_BISPECTRUM_CALCULATE, vin, nIn);
        pstBispectrumCalculate(pst->pstLower, vin, nIn, vout, nOut);
        mdlGetReply(pst->mdl,rID,&outUpper,&nOutUpper);
        *out += outUpper;
    }
    else {
        *out = pkdBispectrumCalculate(plcl->pkd,in->iGrid1,in->iGrid2,in->iGrid3);
    }
    return nOut;
}

double MSR::BispectrumCalculate(int iGrid1, int iGrid2, int iGrid3) {
    struct inBispectrumCalculate in;
    double result;
    in.iGrid1 = iGrid1;
    in.iGrid2 = iGrid2;
    in.iGrid3 = iGrid3;
    pstBispectrumCalculate(pst, &in, sizeof(in), &result, sizeof(result));
    return result;
}
