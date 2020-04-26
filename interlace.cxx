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
#include "gridinfo.hpp"
using namespace gridinfo;

void pkdInterlace(PKD pkd, int iGridTarget, int iGridSource) {
    auto fft = pkd->fft;
    int nGrid = fft->rgrid->n1;
    auto iNyquist = nGrid / 2;

    GridInfo G(pkd->mdl,fft);

    complex_array_t K1,K2;
    auto data1 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGridTarget;
    auto data2 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite)) + fft->rgrid->nLocal * iGridSource;
    G.setupArray(data1,K1);
    G.setupArray(data2,K2);

    for( auto index=K1.begin(); index!=K1.end(); ++index ) {
    	auto pos = index.position();
	auto i = pos[0]; // i,j,k are all positive (absolute value)
	auto j = pos[1]>iNyquist ? nGrid - pos[1] : pos[1];
	auto k = pos[2]>iNyquist ? nGrid - pos[2] : pos[2];
	auto v1 = *index;

	auto jj = pos[1]>iNyquist ? pos[1] - nGrid : pos[1];
	auto kk = pos[2]>iNyquist ? pos[2] - nGrid : pos[2];
	float theta = M_PI/nGrid * (i + jj + kk);
	auto v2 = K2(index.position()) * complex_t(cosf(theta),sinf(theta));
	v1 = complex_t(0.5) * (v1 + v2);
	*index = v1;
	}
    }

int pstInterlace(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inInterlace *in = reinterpret_cast<struct inInterlace *>(vin);
    assert (nIn==sizeof(struct inInterlace) );
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl, pst->idUpper, PST_INTERLACE, vin, nIn);
	pstInterlace(pst->pstLower, vin, nIn, NULL, 0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
    	pkdInterlace(plcl->pkd,in->iGridTarget,in->iGridSource);
	}
    return 0;
    }

void MSR::Interlace(int iGridTarget,int iGridSource) {
    struct inInterlace in;
    in.iGridTarget = iGridTarget;
    in.iGridSource = iGridSource;
    pstInterlace(pst, &in, sizeof(in), NULL, 0);
    }
