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
#include "whitenoise.hpp"
using namespace gridinfo;

#include <vector>

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

class LinearSignal : public NoiseGenerator {
private:
    CSM csm;
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double iLbox;
    double *logk, *field;
protected:
    virtual void update(gridinfo::complex_vector_t &pencil,gridinfo::complex_vector_t &noise,int j,int k);
public:
    explicit LinearSignal(CSM csm,double a,double Lbox,unsigned long seed,bool bFixed=false,float fPhase=0);
    virtual ~LinearSignal();
    };

LinearSignal::LinearSignal(CSM csm,double a,double Lbox,unsigned long seed,bool bFixed,float fPhase)
    : NoiseGenerator(seed,bFixed,fPhase) {
    auto size = csm->val.classData.perturbations.size_k;
    this->csm = csm;
    iLbox = 2*M_PI / Lbox;
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    logk = (double*)malloc(sizeof(double)*size);
    field = (double*)malloc(sizeof(double)*size);
    for (auto i = 0; i < size; i++){
        auto k = csm->val.classData.perturbations.k[i];
        logk[i] = log(k);
        field[i] = csmDeltaRho_pk(csm, a, k);
        field[i] /= csmZeta(csm, k);
    }
    gsl_spline_init(spline, logk, field, size);
    }

LinearSignal::~LinearSignal() {
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
    }

void LinearSignal::update(complex_vector_t &pencil,complex_vector_t &noise,int iy,int iz) {
    float k2jk = iy*iy + iz*iz;
    for( auto index=noise.begin(); index!=noise.end(); ++index ) {
	auto ix = index.position()[0];
	auto k = sqrt(k2jk + ix*ix) * iLbox;
	if (k>0) {
	    float signal = csmZeta(csm, k)*gsl_spline_eval(spline, log(k), acc);
	    auto wnoise = *index;
            pencil(index.position()) += wnoise*signal;
	    }
	}
    }

extern "C"
void pkdMeasurePk(PKD pkd, int iAssignment,
    int bLinear, int iSeed, int bFixed, float fPhase, double Lbox, double a,
    int nGrid, int nBins, double *fK, double *fPower, uint64_t *nPower, double *fPowerAll) {
    //mdlGridCoord first, last, index;
    assert(pkd->fft != NULL);
    auto iNyquist = nGrid / 2;
    auto fft = pkd->fft;
    GridInfo G(pkd->mdl,fft);
    Window W(nGrid,iAssignment);

    complex_array_t K1,K2;
    auto data1 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite));
    auto data2 = data1 + fft->rgrid->nLocal;
    G.setupArray(data1,K1);
    G.setupArray(data2,K2);

    for( auto i=0; i<nBins; i++ ) {
	fK[i] = 0.0;
	fPower[i] = 0.0;
	fPowerAll[i] = 0.0;
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
	    if (i!=0 && i!=iNyquist) { // Account for negative Kx values
		fK[ks] += log(ak);
		fPower[ks] += std::norm(v1);
		nPower[ks] += 1;
		}
	    }
	}
    if (bLinear) {
	LinearSignal ng(pkd->csm,a,Lbox,iSeed,bFixed,fPhase);
	ng.FillNoise(K1,nGrid);
	for( auto index=K1.begin(); index!=K1.end(); ++index ) {
    	 auto pos = index.position();
	    auto i = pos[0]; // i,j,k are all positive (absolute value)
	    auto j = pos[1]>iNyquist ? nGrid - pos[1] : pos[1];
	    auto k = pos[2]>iNyquist ? nGrid - pos[2] : pos[2];
	    auto ak = sqrt(i*i + j*j + k*k);
	    auto ks = int(ak);
	    if ( ks >= 1 && ks <= iNyquist ) {
    #ifdef LINEAR_PK
		ks = floor((ks-1.0) * scale);
    #else
		ks = floor(log(ks) * scale);
    #endif
		assert(ks>=0 && ks<nBins);
		fPowerAll[ks] += std::norm(*index);
		if (i!=0 && i!=iNyquist) { // Account for negative Kx values
		    fPowerAll[ks] += std::norm(*index);
		    }
		}
	    }
	}
    }
