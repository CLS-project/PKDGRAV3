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

#include "RngStream.h"
#include "whitenoise.hpp"
using namespace gridinfo;
using namespace blitz;

static const std::complex<float> I(0,1);

/* Gaussian noise in k-space. Note correction sqrt(2) because of FFT normalization. */
static complex_t pairc( RngStream g, int bFixed, float fPhase ) {
    float x1, x2, w;
    do {
	x1 = 2.0 * RngStream_RandU01(g) - 1.0;
	x2 = 2.0 * RngStream_RandU01(g) - 1.0;
	w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 || w == 0.0 ); /* Loop ~ 21.5% of the time */
    if (!bFixed) {
	/* Proper Gaussian Deviate */
	w = sqrt(-log(w)/w);
	return w * (x1 + I * x2);
	}
    else {
	float theta = atan2f(x2,x1) + fPhase;
	return cosf(theta) + I * sinf(theta);
	}
    }

/*
** Generate Gaussian white noise in k-space. The noise is in the proper form for
** an inverse FFT. The complex conjugates in the Nyquist planes are correct, and
** the normalization is such that that the inverse FFT needs to be normalized
** by sqrt(Ngrid^3) compared with Ngrid^3 with FFT followed by IFFT.
*/
void pkdGenerateNoise(PKD pkd,unsigned long seed,int bFixed, float fPhase,MDLFFT fft,complex_array_t &K,double *mean,double *csq) {
    MDL mdl = pkd->mdl;
    const int nGrid = fft->kgrid->n3;
    const int iNyquist = nGrid / 2;
    RngStream g;
    unsigned long fullKey[6];
    int i,j,k;
    complex_t v_ny,v_wn;

    fullKey[0] = seed;
    fullKey[1] = fullKey[0];
    fullKey[2] = fullKey[0];
    fullKey[3] = seed;
    fullKey[4] = fullKey[3];
    fullKey[5] = fullKey[3];

    /* Remember, we take elements from the stream and because we have increased the
    ** precision, we take pairs. Also, the method of determining Gaussian deviates
    ** also requires pairs (so four elements from the random stream), and some of
    ** these pairs are rejected.
    */
    g = RngStream_CreateStream ("IC");
#ifndef USE_SINGLE
    RngStream_IncreasedPrecis (g, 1);
#endif
    RngStream_SetSeed(g,fullKey);

    *mean = *csq = 0.0;

    j = k = nGrid; /* Start with invalid values so we advance the RNG correctly. */
    for( auto index=K.begin(); index!=K.end(); ++index ) {
    	auto pos = index.position();
	if (j!=pos[1] || k!=pos[2]) { // Start a new pencil
	    assert(pos[0]==0);
	    j = pos[1];
	    k = pos[2];
	    int jj = j<=iNyquist ? j*2 : (nGrid-j)*2 % nGrid + 1;
	    int kk = k<=iNyquist ? k*2 : (nGrid-k)*2 % nGrid + 1;

	    /* We need the sample for x==0 AND/OR x==iNyquist, usually both but at least one. */
	    RngStream_ResetStartStream (g);
	    if ( pos[2] <= iNyquist && (pos[2]%iNyquist!=0||pos[1]<=iNyquist) ) { /* Positive zone */
		RngStream_AdvanceState (g, 0, (1LL<<40)*jj + (1LL<<20)*kk );
		v_ny = pairc(g,bFixed,fPhase);
		v_wn = pairc(g,bFixed,fPhase);

		if ( (pos[1]==0 || pos[1]==iNyquist)  && (pos[2]==0 || pos[2]==iNyquist) ) {
		    /* These are real because they must be a complex conjugate of themselves. */
		    v_ny = std::real(v_ny);
		    v_wn = std::real(v_wn);
		    /* DC mode is zero */
		    if ( pos[2]==0 && pos[1]==0) v_wn = 0.0;
		    }
		}
	    /* We need to generate the correct complex conjugates */
	    else {
		int jjc = j<=iNyquist ? j*2 + 1 : (nGrid-j) % nGrid * 2;
		int kkc = k<=iNyquist ? k*2 + 1 : (nGrid-k) % nGrid * 2;
		if (k%iNyquist == 0) { kkc = kk; }
		if (j%iNyquist == 0) { jjc = jj; }
		RngStream_AdvanceState (g, 0, (1LL<<40)*jjc + (1LL<<20)*kkc );
		v_ny = conj(pairc(g,bFixed,fPhase));
		v_wn = conj(pairc(g,bFixed,fPhase));
		RngStream_ResetStartStream (g);
		RngStream_AdvanceState (g, 0, (1LL<<40)*jj + (1LL<<20)*kk );
		pairc(g,bFixed,fPhase); pairc(g,bFixed,fPhase); /* Burn the two samples we didn't use. */
		}
	    K(iNyquist,pos[1],pos[2]) = v_ny;
	    *index = v_wn;
	    }
	else if (pos[0]!=iNyquist) {
	    *index = pairc(g,bFixed,fPhase);
	    }
	*mean += std::real(*index) + std::imag(*index);
	*csq += std::norm(*index);
	}
    RngStream_DeleteStream(&g);
    }
