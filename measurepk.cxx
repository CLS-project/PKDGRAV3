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
#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif
#include "pst.h"
#include "master.h"
#include "gridinfo.hpp"
using namespace gridinfo;

#include <vector>

static void initPk(void *vpkd, void *g) {
    FFTW3(real) * r = (FFTW3(real) *)g;
    *r = 0.0;
    }
static void combPk(void *vpkd, void *g1, void *g2) {
    FFTW3(real) * r1 = (FFTW3(real) *)g1;
    FFTW3(real) * r2 = (FFTW3(real) *)g2;
    *r1 += *r2;
    }

static void assign_mass(PKD pkd, int iAssignment, double dTotalMass, double fDelta, MDLFFT fft, FFTW3(real) *fftData) {
    int nGrid = fft->rgrid->n1;
    double fftNormalize = 1.0 / (1.0*nGrid*nGrid*nGrid);
    mdlGridCoord first, last;
    int i, j;
    double r[3];

    mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&first,&last,1);
    for( i=first.i; i<last.i; ++i ) fftData[i] = 0.0;
    mdlCOcache(pkd->mdl,CID_PK,NULL,fftData,sizeof(FFTW3(real)),last.i,pkd,initPk,combPk);
    pkdAssignMass(pkd, ROOT, nGrid, fDelta, iAssignment);
    mdlFinishCache(pkd->mdl,CID_PK);

    for( i=first.i; i<last.i; ++i ) {
	assert(fftData[i] >= 0.0);
	}
    double dRhoMean = dTotalMass * fftNormalize;
    double diRhoMean = 1.0 / dRhoMean;

    /*printf( "Calculating density contrast\n" );*/
    for( i=first.i; i<last.i; ++i ) {
	fftData[i] = fftData[i]*diRhoMean - 1.0;
	}

    mdlFFT(pkd->mdl,fft,fftData);
    }

class Window : public std::vector<float> {
public:
    Window(int nGrid,int iAssignment);
    };

Window::Window(int nGrid,int iAssignment) {
    reserve(nGrid);
    for( auto i=0; i<nGrid; ++i) {
	float win = M_PI * i / nGrid;
	if(win>0.1) win = win / sinf(win);
	else win=1.0 / (1.0-win*win/6.0*(1.0-win*win/20.0*(1.0-win*win/76.0)));
	push_back(powf(win,iAssignment));
	}
    }

extern "C"
void pkdMeasurePk(PKD pkd, double dTotalMass, int iAssignment, int bInterlace,
    int nGrid, int nBins, double *fK, double *fPower, uint64_t *nPower) {
    mdlGridCoord first, last, index;
    assert(pkd->fft != NULL);
    auto iNyquist = nGrid / 2;
    auto fft = pkd->fft;
    GridInfo G(pkd->mdl,fft);
    Window W(nGrid,iAssignment);

    mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&first,&last,1);
    auto fftData1 = reinterpret_cast<FFTW3(real) *>(mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(real)),pkd->pLite));
    auto fftData2 = reinterpret_cast<FFTW3(real) *>(mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(real)),fftData1 + fft->rgrid->nLocal));

    assign_mass(pkd,iAssignment,dTotalMass,0.0f,fft,fftData1);
    if (bInterlace) {
	assign_mass(pkd,iAssignment,dTotalMass,0.5f,fft,fftData2);
	}

    complex_array_t K1,K2;
    auto data1 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite));
    auto data2 = data1 + fft->rgrid->nLocal;
    G.setupArray(data1,K1);
    G.setupArray(data2,K2);

    for( auto i=0; i<nBins; i++ ) {
	fK[i] = 0.0;
	fPower[i] = 0.0;
	nPower[i] = 0;
	}

    complex_t fftNormalize = 1.0 / (1.0*nGrid*nGrid*nGrid);
#ifdef LINEAR_PK
    double scale = nBins * 1.0 / iNyquist;
#else
    double scale = nBins * 1.0 / log(iNyquist+1);
#endif
    for( auto index=K1.begin(); index!=K1.end(); ++index ) {
    	auto pos = index.position();
	auto i = pos[0]; // i,j,k are all positive (absolute value)
	auto j = pos[1]>iNyquist ? nGrid - pos[1] : pos[1];
	auto k = pos[2]>iNyquist ? nGrid - pos[2] : pos[2];
	auto v1 = *index;
	if (bInterlace) { // jj,kk can be negative
	    auto jj = pos[1]>iNyquist ? pos[1] - nGrid : pos[1];
	    auto kk = pos[2]>iNyquist ? pos[2] - nGrid : pos[2];
	    float theta = M_PI/nGrid * (i + jj + kk);
	    auto v2 = K2(index.position()) * complex_t(cosf(theta),sinf(theta));
	    v1 = complex_t(0.5) * (v1 + v2);
	    }
	v1 *= fftNormalize;       // Normalize for FFT
	v1 *= W[i] * W[j] * W[k]; // Correction for mass assignment
	*index = v1;

	auto ak = sqrt(i*i + j*j + k*k);
	auto ks = int(ak);
	if ( ks >= 1 && ks <= iNyquist ) {
#ifdef LINEAR_PK
	    ks = floor((ks-1.0) * scale);
#else
	    ks = floor(log(ks) * scale);
#endif
	    assert(ks>=0 && ks<nBins);
	    fK[ks] += log(ak);
	    fPower[ks] += std::norm(v1);
	    nPower[ks] += 1;
	    }
	}
    }
