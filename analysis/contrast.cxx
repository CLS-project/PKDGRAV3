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

static void pkdDensityContrast(PKD pkd,double dTotalMass,int iGrid,bool k=false) {
    auto fft = pkd->fft;
    GridInfo G(pkd->mdl,fft);
    int nGrid = fft->rgrid->n1;
    auto data = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGrid;
    real_array_t R;
    G.setupArray(data,R);
    if (k) { // Delta(k): same as Delta(r) below, but apply normalization for the FFT (divide by nGrid^3)
        real_t diTotalMass = 1.0 / dTotalMass;
        real_t fftNormalize = 1.0 / (1.0*nGrid*nGrid*nGrid);
        R = R * diTotalMass - fftNormalize;
        mdlFFT(pkd->mdl,fft,data);
    }
    else { // Delta(r) = rho / rho_bar - 1.0
        real_t diRhoBar = (1.0*nGrid*nGrid*nGrid) / dTotalMass;
        R = R * diRhoBar - 1.0;
    }
}

int pstDensityContrast(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inDensityContrast *in = reinterpret_cast<struct inDensityContrast *>(vin);
    assert (nIn==sizeof(struct inDensityContrast) );
    if (pstNotCore(pst)) {
        int rID = mdlReqService(pst->mdl, pst->idUpper, PST_DENSITY_CONTRAST, vin, nIn);
        pstDensityContrast(pst->pstLower, vin, nIn, NULL, 0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdDensityContrast(plcl->pkd,in->dTotalMass,in->iGrid,in->k);
    }
    return 0;
}

void MSR::DensityContrast(int iGrid,bool k) {
    struct inDensityContrast in;
    in.dTotalMass = TotalMass();
    in.iGrid = iGrid;
    in.k = k;
    pstDensityContrast(pst, &in, sizeof(in), NULL, 0);
}
