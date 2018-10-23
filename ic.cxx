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
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "pst.h"
#include "ic.h"
#include "RngStream.h"
#include "gridinfo.hpp"
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

static int wrap(int v,int h,int m) {
    return v - (v > h ? m : 0);
    }

// A blitz++ friendly wrap function. returns "ik" given array index
// Range: (-iNyquist,iNyquist] where iNyquist = m/2
static float fwrap(float v,float m) {
    return v - (v > m*0.5 ? m : 0);
    }
BZ_DEFINE_BINARY_FUNC(Fn_fwrap,fwrap)
BZ_DECLARE_ARRAY_ET_BINARY(fwrap,     Fn_fwrap)
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(fwrap,     Fn_fwrap, float)

/*
** Generate Gaussian white noise in k-space. The noise is in the proper form for
** an inverse FFT. The complex conjugates in the Nyquist planes are correct, and
** the normalization is such that that the inverse FFT needs to be normalized
** by sqrt(Ngrid^3) compared with Ngrid^3 with FFT followed by IFFT.
*/
static void pkdGenerateNoise(PKD pkd,unsigned long seed,int bFixed, float fPhase,MDLFFT fft,complex_array_t &K,double *mean,double *csq) {
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

class Power {
    typedef struct {
	class Power *power;
	double r;
	} varianceParameters;
    static double variance_integrand(double ak, void * params);
public:
    virtual double getAmplitude(double k) = 0;
    double variance(double dRadius,double k0,double k1);
    };

double Power::variance_integrand(double ak, void * params) {
    varianceParameters *vprm = reinterpret_cast<varianceParameters *>(params);
    double x, w;
    /* Window function for spherical tophat of given radius (e.g., 8 Mpc/h) */
    x = ak * vprm->r;
    w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
    return vprm->power->getAmplitude(ak)*ak*ak*w*w*4.0*M_PI;
    }

double Power::variance(double dRadius,double k0,double k1) {
    varianceParameters vprm;
    gsl_function F;
    double result, error;
    gsl_integration_workspace *W = gsl_integration_workspace_alloc (1000);
    vprm.power = this;
    vprm.r = dRadius; /* 8 Mpc/h for example */
    F.function = &variance_integrand;
    F.params = &vprm;
    gsl_integration_qag(&F, exp(k0), exp(k1),
	0.0, 1e-6, 1000, GSL_INTEG_GAUSS61, W, &result, &error);
    gsl_integration_workspace_free(W);
    return result;
    }

class PowerTransfer : public Power {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double *tk, *tf;
    double spectral;
    double normalization;
    int nTf;
    double rise0, rise1;
    double dSigma8;
public:
    virtual double getAmplitude(double k);
    PowerTransfer(CSM csm, double a,int nTf, double *tk, double *tf);
    double variance(double dRadius);
    };

double PowerTransfer::getAmplitude(double k) {
    double lk = log(k);
    double T;

    if (lk > tk[nTf-1]) /* Extrapolate beyond kmax */
	T =  tf[nTf-1] + (lk - tk[nTf-1]) * rise1;
    else if (lk < tk[0]) /* Extrapolate beyond kmin */
	T =  tf[0] + (lk - tk[0]) * rise0;
    else
	T = gsl_spline_eval(spline,lk,acc);
    T = exp(T);
    return pow(k,spectral) * normalization * T * T;
    }

double PowerTransfer::variance(double dRadius) {
    return Power::variance(dRadius,tk[0], tk[nTf-1]);
    }

PowerTransfer::PowerTransfer(CSM csm, double a,int nTf, double *tk, double *tf) {
    double D1_0, D2_0, D1_a, D2_a; 
    double f1_0, f2_0, f1_a, f2_a;
    csmComoveGrowth(csm, 1.0, &D1_0, &D2_0, &f1_0, &f2_0); 
    csmComoveGrowth(csm, a, &D1_a, &D2_a, &f1_a, &f2_a);
    normalization = 1.0;
    spectral = csm->val.dSpectral;
    this->nTf = nTf;
    this->tk = tk;
    this->tf = tf;
    rise0 = (tf[0] - tf[1]) / (tk[0] - tk[1]);
    rise1 = (tf[nTf-1] - tf[nTf-2]) / (tk[nTf-1] - tk[nTf-2]);
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc (gsl_interp_cspline, nTf);
    gsl_spline_init(spline, tk, tf, nTf);
    float dSigma8 = csm->val.dSigma8;
    if (dSigma8 > 0) {
	dSigma8 *= D1_a/D1_0;
	normalization *= dSigma8*dSigma8 / variance(8.0);
	}
    else if (csm->val.dNormalization > 0) {
	normalization = csm->val.dNormalization * D1_a/D1_0;
	dSigma8 = sqrt(variance(8.0));
	}
    }

#ifdef __cplusplus
extern "C"
#endif
int pkdGenerateIC(PKD pkd,MDLFFT fft,int iSeed,int bFixed,float fPhase,int nGrid,int b2LPT,double dBoxSize,
    struct csmVariables *cosmo,double a,int nTf, double *tk, double *tf,
    double *noiseMean, double *noiseCSQ) {
    MDL mdl = pkd->mdl;
    float twopi = 2.0 * 4.0 * atan(1.0);
    float itwopi = 1.0 / twopi;
    float inGrid = 1.0 / nGrid;
    float fftNormalize = inGrid*inGrid*inGrid;
    float ix, iy, iz;
    int iNyquist = nGrid / 2;
    float iLbox = twopi / dBoxSize;
    float iLbox3 = pow(iLbox,3.0);
    float ak, ak2, amp;
    float dOmega;
    double D1_0, D2_0, D1_a, D2_a; 
    double f1_0, f2_0, f1_a, f2_a;
    
    float velFactor; 
    basicParticle *p;
    int nLocal;
    float dSigma8 = cosmo->dSigma8;
    CSM csm;

    csmInitialize(&csm);
    csm->val = *cosmo;

    if (csm->val.classData.bClass){
        csmClassGslInitialize(csm);
    }

    PowerTransfer transfer(csm,a,nTf,tk,tf);
    csmComoveGrowth(csm, 1.0, &D1_0, &D2_0, &f1_0, &f2_0); 
    csmComoveGrowth(csm, a, &D1_a, &D2_a, &f1_a, &f2_a);
    dOmega = cosmo->dOmega0 / (a*a*a*pow(csmExp2Hub(csm, a)/cosmo->dHubble0,2.0));

    double f1 = pow(dOmega,5.0/9.0);
    double f2 = 2.0 * pow(dOmega,6.0/11.0);
    if (mdlSelf(mdl)==0) {
	printf("sigma8=%.15g\n",dSigma8); 
	printf("f1=%.12g (exact) or %.12g (approx)\n", f1_a, f1);
	printf("f2=%.12g (exact) or %.12g (approx)\n", f2_a, f2);//2.0 * f1_a, f2_a);
	fflush(stdout);
	}

    velFactor = csmExp2Hub(csm,a);
    velFactor *= a*a; /* Comoving */

    // Construct the arrays. "data" points to the same block for all threads!
    GridInfo G(pkd->mdl,fft);
    complex_array_t K[10];
    real_array_t R[10];
    auto data = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkdParticleBase(pkd)));
    for(auto i=0; i<10; ++i) {
	G.setupArray(data,K[i]);
	G.setupArray(data,R[i]);
	data += fft->rgrid->nLocal;
	}

    /* Particles will overlap K[0] through K[5] eventually */
    nLocal = blitz::product(R[7].shape());
    blitz::Array<basicParticle,3> output(
        reinterpret_cast<basicParticle *>(mdlSetArray(pkd->mdl,0,0,pkdParticleBase(pkd))),
        blitz::shape(G.n1r(),G.n2(),G.nz()), blitz::neverDeleteData,
        RegularArray(G.sz()));
    
    /* Generate white noise realization -> K[6] */
    if (mdlSelf(mdl)==0) {printf("Generating random noise\n"); fflush(stdout); }
    pkdGenerateNoise(pkd,iSeed,bFixed,fPhase,fft,K[6],noiseMean,noiseCSQ);

    if (mdlSelf(mdl)==0) {printf("Imprinting power\n"); fflush(stdout); }
    for( auto index=K[6].begin(); index!=K[6].end(); ++index ) {
	auto pos = index.position();
	iz = fwrap(pos[2],nGrid); // Range: (-iNyquist,iNyquist]
	iy = fwrap(pos[1],nGrid);
	ix = fwrap(pos[0],nGrid);
	ak2 = ix*ix + iy*iy + iz*iz;
	if (ak2>0) {
	    ak = sqrt(ak2) * iLbox;
	    amp = sqrt(transfer.getAmplitude(ak) * iLbox3) * itwopi / ak2;
	    }
	else amp = 0.0;
	K[9](pos) = *index * amp * iz * -I; // z
	K[8](pos) = *index * amp * iy * -I; // y
	K[7](pos) = *index * amp * ix * -I; // x
	if (b2LPT) {
	    K[0](pos) = K[7](pos) * twopi * ix * -I; // xx
	    K[1](pos) = K[8](pos) * twopi * iy * -I; // yy
	    K[2](pos) = K[9](pos) * twopi * iz * -I; // zz
	    K[3](pos) = K[7](pos) * twopi * iy * -I; // xy
	    K[4](pos) = K[8](pos) * twopi * iz * -I; // yz
	    K[5](pos) = K[9](pos) * twopi * ix * -I; // zx
	    }
	}

    if (mdlSelf(mdl)==0) {printf("Generating x displacements\n"); fflush(stdout); }
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[7].dataFirst() );
    if (mdlSelf(mdl)==0) {printf("Generating y displacements\n"); fflush(stdout); }
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[8].dataFirst() );
    if (mdlSelf(mdl)==0) {printf("Generating z displacements\n"); fflush(stdout); }
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[9].dataFirst() );
    if (b2LPT) {
	if (mdlSelf(mdl)==0) {printf("Generating xx term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[0].dataFirst() );
	if (mdlSelf(mdl)==0) {printf("Generating yy term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[1].dataFirst() );
	if (mdlSelf(mdl)==0) {printf("Generating zz term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[2].dataFirst() );
	if (mdlSelf(mdl)==0) {printf("Generating xy term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[3].dataFirst() );
	if (mdlSelf(mdl)==0) {printf("Generating yz term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[4].dataFirst() );
	if (mdlSelf(mdl)==0) {printf("Generating xz term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[5].dataFirst() );

	/* Calculate the source term */
	if (mdlSelf(mdl)==0) {printf("Generating source term\n"); fflush(stdout); }
	for( auto index=K[6].begin(); index!=K[6].end(); ++index ) {
	    auto pos = index.position();
	    *index = K[0](pos)*K[1](pos) + K[0](pos)*K[2](pos) + K[1](pos)*K[2](pos)
		- K[3](pos)*K[3](pos) - K[4](pos)*K[4](pos) - K[5](pos)*K[5](pos);
	    }
	mdlFFT(mdl, fft, R[6].dataFirst() );
	}

    /* Move the 1LPT positions/velocities to the particle area */
    if (mdlSelf(mdl)==0) {printf("Transfering 1LPT results to output area\n"); fflush(stdout); }
    for( auto index=R[7].begin(); index!=R[7].end(); ++index ) {
	auto pos = index.position();
	float x = R[7](pos);
	float y = R[8](pos);
	float z = R[9](pos);
	output(pos).dr[0] = x;
	output(pos).dr[1] = y;
	output(pos).dr[2] = z;
	output(pos).v[0] = f1_a * x * velFactor;
	output(pos).v[1] = f1_a * y * velFactor;
	output(pos).v[2] = f1_a * z * velFactor;
	}

    if (b2LPT) {
        float aeq = csmRadMatEquivalence(csm);
        float aeqDratio = aeq/D1_a;
	float D2 = D2_a/pow(D1_a,2);

	for( auto index=K[6].begin(); index!=K[6].end(); ++index ) {
	    auto pos = index.position();
	    iz = fwrap(pos[2],nGrid); // Range: (-iNyquist,iNyquist]
	    iy = fwrap(pos[1],nGrid);
	    ix = fwrap(pos[0],nGrid);
	    ak2 = ix*ix + iy*iy + iz*iz;
	    if (ak2>0.0) {
      		D2 = D2_a/pow(D1_a,2)/(ak2 * twopi);  // The source term contains phi^2 which in turn contains D1, we need to divide by D1^2 (cf Scoccimarro Transients paper, appendix D). Here we use normalized D1 values in contrast to there.
		K[7](pos) = D2 * *index * ix * -I;
		K[8](pos) = D2 * *index * iy * -I;
		K[9](pos) = D2 * *index * iz * -I;
		}
	    else K[7](pos) = K[8](pos) = K[9](pos) = 0.0;
	    }

	if (mdlSelf(mdl)==0) {printf("Generating x2 displacements\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[7].dataFirst() );
	if (mdlSelf(mdl)==0) {printf("Generating y2 displacements\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[8].dataFirst() );
	if (mdlSelf(mdl)==0) {printf("Generating z2 displacements\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)K[9].dataFirst() );

	/* Add the 2LPT positions/velocities corrections to the particle area */
	if (mdlSelf(mdl)==0) {printf("Transfering 2LPT results to output area\n"); fflush(stdout); }
	for( auto index=R[7].begin(); index!=R[7].end(); ++index ) {
	    auto pos = index.position();
	    float x = R[7](pos) * fftNormalize;
	    float y = R[8](pos) * fftNormalize;
	    float z = R[9](pos) * fftNormalize;
	    output(pos).dr[0] += x;
	    output(pos).dr[1] += y;
	    output(pos).dr[2] += z;
	    output(pos).v[0] += f2_a * x * velFactor;
	    output(pos).v[1] += f2_a * y * velFactor;
	    output(pos).v[2] += f2_a * z * velFactor;
	    }
	}
    csmFinish(csm);

    return nLocal;
    }


/*
** Just like pkdGenerateIC(), this function will generate the particle ICs,
** the difference being that this will use the CLASS interface.
** Note that no back-scaling is performed; the CLASS transfer function
** is realized at the given scale factor `a` as is.
** Currently, only 1LPT is implemented.
*/
#ifdef __cplusplus
extern "C"
#endif
int pkdGenerateClassICm(PKD pkd, MDLFFT fft, int iSeed, int bFixed, float fPhase, int nGrid,
    double dBoxSize, struct csmVariables *cosmo, double a, double *noiseMean, double *noiseCSQ) {
    MDL mdl = pkd->mdl;
    float kx, ky, kz;
    int iNyquist = nGrid / 2;
    float iLbox = 2*M_PI / dBoxSize;
    float k2, amp;
    float dOmega;

    basicParticle *p;
    int nLocal;
    CSM csm;

    csmInitialize(&csm);
    csm->val = *cosmo;
    assert(csm->val.classData.bClass);
    csmClassGslInitialize(csm);

    // Construct the arrays. "data" points to the same block for all threads!
    GridInfo G(pkd->mdl,fft);
    complex_array_t K[10];
    real_array_t R[10];
    auto data = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkdParticleBase(pkd)));
    for(auto i=0; i<10; ++i) {
	G.setupArray(data,K[i]);
	G.setupArray(data,R[i]);
	data += fft->rgrid->nLocal;
	}

    /* Particles will overlap K[0] through K[5] eventually */
    nLocal = blitz::product(R[7].shape());
    blitz::Array<basicParticle,3> output(
        reinterpret_cast<basicParticle *>(mdlSetArray(pkd->mdl,0,0,pkdParticleBase(pkd))),
        blitz::shape(G.n1r(),G.n2(),G.nz()), blitz::neverDeleteData,
        RegularArray(G.sz()));

    /* Generate white noise realization -> K[6] */
    if (mdlSelf(mdl)==0) {printf("Generating random noise\n"); fflush(stdout);}
    pkdGenerateNoise(pkd,iSeed,bFixed,fPhase,fft,K[6],noiseMean,noiseCSQ);

    if (mdlSelf(mdl)==0) {printf("Imprinting power\n"); fflush(stdout);}

    /* Particle positions */
    for( auto index=K[6].begin(); index!=K[6].end(); ++index ) {
	auto pos = index.position();
	kz = fwrap(pos[2],fft->rgrid->n3) * iLbox;
	ky = fwrap(pos[1],fft->rgrid->n2) * iLbox;
	kx = fwrap(pos[0],fft->rgrid->n1) * iLbox;
	k2 = kx*kx + ky*ky + kz*kz;
	if (k2>0) {
	    amp = csmDelta_m(csm, a, sqrt(k2));
	    K[7](pos) = *index * amp * kx/k2 * I;
	    K[8](pos) = *index * amp * ky/k2 * I;
	    K[9](pos) = *index * amp * kz/k2 * I;
	}
	else{
	    K[7](pos) = .0;
	    K[8](pos) = .0;
	    K[9](pos) = .0;
	}
    }

    if (mdlSelf(mdl)==0) {printf("Generating x displacements\n"); fflush(stdout);}
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[7].dataFirst());
    if (mdlSelf(mdl)==0) {printf("Generating y displacements\n"); fflush(stdout);}
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[8].dataFirst());
    if (mdlSelf(mdl)==0) {printf("Generating z displacements\n"); fflush(stdout);}
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[9].dataFirst());

    /* Move the 1LPT positions/velocities to the particle area */
    if (mdlSelf(mdl)==0) {printf("Transfering 1LPT positions to output area\n"); fflush(stdout);}
    for( auto index=R[7].begin(); index!=R[7].end(); ++index ) {
	auto pos = index.position();
	output(pos).dr[0] = R[7](pos);
	output(pos).dr[1] = R[8](pos);
	output(pos).dr[2] = R[9](pos);
    }

    /* Particle velocities */
    for( auto index=K[6].begin(); index!=K[6].end(); ++index ) {
	auto pos = index.position();
	kz = fwrap(pos[2],fft->rgrid->n3) * iLbox;
	ky = fwrap(pos[1],fft->rgrid->n2) * iLbox;
	kx = fwrap(pos[0],fft->rgrid->n1) * iLbox;
	k2 = kx*kx + ky*ky + kz*kz;
	if (k2>0) {
	    amp = csmTheta_m(csm, a, sqrt(k2));
	    K[7](pos) = *index * amp * kx/k2 * -I;
	    K[8](pos) = *index * amp * ky/k2 * -I;
	    K[9](pos) = *index * amp * kz/k2 * -I;
	}
	else{
	    K[7](pos) = .0;
	    K[8](pos) = .0;
	    K[9](pos) = .0;
        }
    }
    if (mdlSelf(mdl)==0) {printf("Generating x velocities\n"); fflush(stdout);}
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[7].dataFirst());
    if (mdlSelf(mdl)==0) {printf("Generating y velocitites\n"); fflush(stdout);}
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[8].dataFirst());
    if (mdlSelf(mdl)==0) {printf("Generating z velocities\n"); fflush(stdout);}
    mdlIFFT(mdl, fft, (FFTW3(complex)*)K[9].dataFirst());

    /* Move the 1LPT positions/velocities to the particle area */
    if (mdlSelf(mdl)==0) {printf("Transfering 1LPT velocities to output area\n"); fflush(stdout); }
    for( auto index=R[7].begin(); index!=R[7].end(); ++index ) {
	auto pos = index.position();
	output(pos).v[0] = R[7](pos);
	output(pos).v[1] = R[8](pos);
	output(pos).v[2] = R[9](pos);
	}

    csmFinish(csm);

    return nLocal;
}

void pkdGenerateLinGrid(PKD pkd, MDLFFT fft, double a, double a_next, double Lbox, int iSeed,
    int bFixed, float fPhase, int bRho){
    /* If bRho == 0, we generate the \delta field,
    ** otherwise we generate the \delta\rho field.
    ** If a == a_next, we generate the field at this a.
    ** Otherwize, we generate the weighted average of the
    ** field over the interval [a, a_next]. Note that
    ** averaging is only implemented for bRho == 1.
    */
    if (!bRho && a != a_next){
        fprintf(stderr,
            "WARNING: In pkdGenerateLinGrid(): Averaging of \\delta fields not implemented\n");
        abort();
    }
    /* For the sake of performance we precompute the (possibly averaged)
    ** field at the |k|'s at which the perturbations are tabulated.
    */
    size_t size = pkd->param.csm->val.classData.perturbations.size_k;
    double *logk = (double*)malloc(sizeof(double)*size);
    double *field = (double*)malloc(sizeof(double)*size);
    size_t i;
    double k;
    gsl_interp_accel *acc;
    gsl_spline *spline;
    for (i = 0; i < size; i++){
        k = pkd->param.csm->val.classData.perturbations.k[i];
        logk[i] = log(k);
        if (bRho)
            /* Use the \delta\rho field, possibly averaged */
            field[i] = csmDeltaRho_lin(pkd->param.csm, a, a_next, k);
        else
            /* Use the \delta field */
            field[i] = csmDelta_lin(pkd->param.csm, a, k);
        /* The cubic splining is much better without the analytic zeta */
        field[i] /= csmZeta(pkd->param.csm, k);
    }
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(spline, logk, field, size);
    /* Generate grid */
    double noiseMean, noiseCSQ;
    int idx;
    double ix, iy, iz, i2;
    double iLbox = 2*M_PI/Lbox;
    int iNuquist = fft->rgrid->n3 / 2;

    GridInfo G(pkd->mdl,fft);
    complex_array_t K;
    G.setupArray((FFTW3(real) *)mdlSetArray(pkd->mdl,0,0,pkd->pLite),K);
    pkdGenerateNoise(pkd, iSeed, bFixed, fPhase, fft, K, &noiseMean, &noiseCSQ);
    for( auto index=K.begin(); index!=K.end(); ++index ) {
	auto pos = index.position();
	iz = fwrap(pos[2],fft->rgrid->n3); // Range: (-iNyquist,iNyquist]
	iy = fwrap(pos[1],fft->rgrid->n2);
	ix = fwrap(pos[0],fft->rgrid->n1);
        i2 = ix*ix + iy*iy + iz*iz;
        if (i2>0){
            k = sqrt((double)i2)*iLbox;
            *index *= csmZeta(pkd->param.csm, k)*gsl_spline_eval(spline, log(k), acc);
        }
        else
            *index = 0.0;
    }
    /* Cleanup */
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
    free(logk);
    free(field);
}

void pltMoveIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inMoveIC *in = reinterpret_cast<struct inMoveIC *>(vin);
    int i;

    mdlassert(pst->mdl,nIn == sizeof(struct inMoveIC));
    assert(pstOnNode(pst)); /* We pass around pointers! */
    if (pstNotCore(pst)) {
	struct inMoveIC icUp;
	
	icUp.pBase = in->pBase;
	icUp.nMove = pst->nUpper * in->nMove / pst->nLeaves;
	in->nMove -= icUp.nMove;
	icUp.iStart = in->iStart + in->nMove;
	icUp.fMass = in->fMass;
	icUp.fSoft = in->fSoft;
	icUp.nGrid = in->nGrid;
	icUp.nInflateFactor = in->nInflateFactor;

	int rID = mdlReqService(pst->mdl,pst->idUpper,PLT_MOVEIC,&icUp,nIn);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	pltMoveIC(pst->pstLower,in,nIn,vout,pnOut);
	}
    else {
	PKD pkd = plcl->pkd;
	assert(in->nInflateFactor>0);
	if (in->nMove <= (pkd->nStore/in->nInflateFactor)) {}
	else {
	    printf("nMove=%" PRIu64 " nStore=%d nInflateFactor=%d\n", in->nMove, pkd->nStore, in->nInflateFactor);
	}
	assert(in->nMove <= (pkd->nStore/in->nInflateFactor));
	double inGrid = 1.0 / in->nGrid;
	for(i=in->nMove-1; i>=0; --i) {
	    PARTICLE *p = pkdParticle(pkd,i);
	    vel_t *pVel = pkdVel(pkd,p);
#ifdef INTEGER_POSITION
	    // If we have no particle order convert directly to Integerized positions.
	    // We do this to save space as an "Integer" particle is small.
	    if (pkd->bNoParticleOrder) {
		integerParticle *b = ((integerParticle *)in->pBase) + in->iStart + i;
		integerParticle temp;
		memcpy(&temp,b,sizeof(temp));
		pVel[2] = temp.v[2];
		pVel[1] = temp.v[1];
		pVel[0] = temp.v[0];
		pkdSetPosRaw(pkd,p,2,temp.r[2]);
		pkdSetPosRaw(pkd,p,1,temp.r[1]);
		pkdSetPosRaw(pkd,p,0,temp.r[0]);
		}
	    else
#endif
	    {
		expandParticle *b = ((expandParticle *)in->pBase) + in->iStart + i;
		expandParticle temp;
		memcpy(&temp,b,sizeof(temp));
		pVel[2] = temp.v[2];
		pVel[1] = temp.v[1];
		pVel[0] = temp.v[0];
		pkdSetPos(pkd,p,2,temp.dr[2] + (temp.iz+0.5) * inGrid - 0.5);
		pkdSetPos(pkd,p,1,temp.dr[1] + (temp.iy+0.5) * inGrid - 0.5);
		pkdSetPos(pkd,p,0,temp.dr[0] + (temp.ix+0.5) * inGrid - 0.5);
		if (!pkd->bNoParticleOrder)
		    p->iOrder = temp.ix + in->nGrid*(temp.iy + 1ul*in->nGrid*temp.iz);
		}
	    pkdSetClass(pkd,in->fMass,in->fSoft,FIO_SPECIES_DARK,p);
	    p->bMarked = 1;
	    p->uRung = 0;
	    if (pkd->bNoParticleOrder) ((UPARTICLE *)p)->iGroup = 0;
	    else p->uNewRung = 0;
	    float *pPot = pkdPot(pkd,p);
	    if (pPot) pPot = 0;
	    }
	pkd->nLocal = pkd->nActive = in->nMove;
	}
    if (pnOut) *pnOut = 0;
    }

void pstMoveIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inGenerateIC *in = reinterpret_cast<struct inGenerateIC *>(vin);

    if (pstOffNode(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MOVEIC,vin,nIn);
	pstMoveIC(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	MDLFFT fft = pkd->fft;
	int myProc = mdlProc(pst->mdl);
	int iProc;
	uint64_t iUnderBeg=0, iOverBeg=0;
	uint64_t iUnderEnd, iOverEnd;
	uint64_t iBeg, iEnd;

	int *scount = reinterpret_cast<int *>(malloc(sizeof(int)*mdlProcs(pst->mdl))); assert(scount!=NULL);
	int *sdisps = reinterpret_cast<int *>(malloc(sizeof(int)*mdlProcs(pst->mdl))); assert(sdisps!=NULL);
	int *rcount = reinterpret_cast<int *>(malloc(sizeof(int)*mdlProcs(pst->mdl))); assert(rcount!=NULL);
	int *rdisps = reinterpret_cast<int *>(malloc(sizeof(int)*mdlProcs(pst->mdl))); assert(rdisps!=NULL);

	assert(pstAmNode(pst));
	assert(fft != NULL);

	assert(in->nInflateFactor>0);
	uint64_t nPerNode = (uint64_t)mdlCores(pst->mdl) * (pkd->nStore / in->nInflateFactor);
	uint64_t nLocal = (int64_t)fft->rgrid->rn[myProc] * in->nGrid*in->nGrid;

	/* Calculate how many slots are free (under) and how many need to be sent (over) before my rank */
	iUnderBeg = iOverBeg = 0;
	for(iProc=0; iProc<myProc; ++iProc) {
	    uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
	    if (nOnNode>nPerNode) iOverBeg += nOnNode - nPerNode;
	    else iUnderBeg += nPerNode - nOnNode;
	    }
	size_t nSize = sizeof(expandParticle);
#ifdef INTEGER_POSITION
	if (pkd->bNoParticleOrder) nSize = sizeof(integerParticle);
#endif
	char *pBase = (char *)pkdParticleBase(pkd);
	char *pRecv = pBase + nSize*nLocal;
	char *eBase;
	if (nLocal > nPerNode) {      /* Too much here: send extra particles to other nodes */
	    eBase = pBase + nSize*nPerNode;
	    iUnderBeg = 0;
	    iOverEnd = iOverBeg + nLocal - nPerNode;
	    for(iProc=0; iProc<mdlProcs(pst->mdl); ++iProc) {
		rcount[iProc] = rdisps[iProc] = 0; // We cannot receive anything
		uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
		if (nOnNode<nPerNode) {
		    iUnderEnd = iUnderBeg + nPerNode - nOnNode;
		    /* The transfer condition */
		    if (iUnderEnd>iOverBeg && iUnderBeg<iOverEnd) {
			iBeg = iOverBeg>iUnderBeg ? iOverBeg : iUnderBeg;
			iEnd = iOverEnd<iUnderEnd ? iOverEnd : iUnderEnd;
			scount[iProc] = (iEnd-iBeg);
			sdisps[iProc] = (iBeg-iOverBeg);
			nLocal -= iEnd-iBeg;
			}
		    else scount[iProc] = sdisps[iProc] = 0;
		    iUnderBeg = iUnderEnd;
		    }
		else scount[iProc] = sdisps[iProc] = 0;
		}
	    assert(nLocal == nPerNode);
	    }
	else if (nLocal < nPerNode) { /* We have room: *maybe* receive particles from other nodes */
	    eBase = pBase + nSize*nLocal;
	    iOverBeg = 0;
	    iUnderEnd = iUnderBeg + nPerNode - nLocal;
	    for(iProc=0; iProc<mdlProcs(pst->mdl); ++iProc) {
		scount[iProc] = sdisps[iProc] = 0; // We have nothing to send
		uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
		if (nOnNode>nPerNode) {
		    iOverEnd = iOverBeg + nOnNode - nPerNode;
		    if (iOverEnd>iUnderBeg && iOverBeg<iUnderEnd) {
			iBeg = iOverBeg>iUnderBeg ? iOverBeg : iUnderBeg;
			iEnd = iOverEnd<iUnderEnd ? iOverEnd : iUnderEnd;
			rcount[iProc] = (iEnd-iBeg);
			rdisps[iProc] = (iBeg-iUnderBeg);
			nLocal += iEnd-iBeg;
			}
		    else rcount[iProc] = rdisps[iProc] = 0;
		    iOverBeg = iOverEnd;
		    }
		else rcount[iProc] = rdisps[iProc] = 0;
		}
	    assert(nLocal <= nPerNode);
	    }
	else {
	    for(iProc=0; iProc<mdlProcs(pst->mdl); ++iProc) {
		rcount[iProc] = rdisps[iProc] = 0; // We cannot receive anything
		scount[iProc] = sdisps[iProc] = 0; // We have nothing to send
		}
	    }
	mdlAlltoallv(pst->mdl, nSize,
	    pBase + nSize*nPerNode, scount, sdisps,
	    pRecv,            rcount, rdisps);
	free(scount);
	free(sdisps);
	free(rcount);
	free(rdisps);
	mdlFFTNodeFinish(pst->mdl,fft);
	pkd->fft = NULL;

	/* We need to relocate the particles */
	struct inMoveIC move;
	uint64_t nTotal;
	nTotal = in->nGrid; /* Careful: 32 bit integer cubed => 64 bit integer */
	nTotal *= in->nGrid;
	nTotal *= in->nGrid;
	move.pBase = (overlayedParticle *)pkdParticleBase(pkd);
	move.iStart = 0;
	move.nMove = nLocal;
	move.fMass = in->dBoxMass;
	move.fSoft = 1.0 / (50.0*in->nGrid);
	move.nGrid = in->nGrid;
	move.nInflateFactor = in->nInflateFactor;
	pltMoveIC(pst,&move,sizeof(move),NULL,0);
	}
    }

/* NOTE: only called when on-node -- pointers are passed around. */
void pltGenerateIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inGenerateICthread *tin = reinterpret_cast<struct inGenerateICthread *>(vin);
    struct inGenerateIC *in = tin->ic;
    struct outGenerateIC *out = reinterpret_cast<struct outGenerateIC *>(vout), outUp;
    mdlassert(pst->mdl,nIn == sizeof(struct inGenerateICthread));
    mdlassert(pst->mdl,vout != NULL);
    assert(pstOnNode(pst)); /* We pass around pointers! */

    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PLT_GENERATEIC,vin,nIn);
	pltGenerateIC(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUp,pnOut);
	out->N += outUp.N;
	out->noiseMean += outUp.noiseMean;
	out->noiseCSQ += outUp.noiseCSQ;
	}
    else {
	if (in->bClass)
	    out->N = pkdGenerateClassICm(plcl->pkd,tin->fft,in->iSeed, in->bFixed,in->fPhase,
	        in->nGrid, in->dBoxSize,&in->cosmo,in->dExpansion,&out->noiseMean,&out->noiseCSQ);
	else
	    out->N = pkdGenerateIC(plcl->pkd,tin->fft,in->iSeed,in->bFixed,in->fPhase,
	        in->nGrid,in->b2LPT,in->dBoxSize, &in->cosmo,in->dExpansion,in->nTf,
	        in->k, in->tf,&out->noiseMean,&out->noiseCSQ);
	out->dExpansion = in->dExpansion;
	}

    if (pnOut) *pnOut = sizeof(struct outGenerateIC);
    }

void pstGenerateIC(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inGenerateIC *in = reinterpret_cast<struct inGenerateIC *>(vin);
    struct outGenerateIC *out = reinterpret_cast<struct outGenerateIC *>(vout), outUp;
    int64_t i;

    mdlassert(pst->mdl,nIn == sizeof(struct inGenerateIC));
    mdlassert(pst->mdl,vout != NULL);

    if (pstAmNode(pst)) {
	struct inGenerateICthread tin;
	MDLFFT fft = mdlFFTNodeInitialize(pst->mdl,in->nGrid,in->nGrid,in->nGrid,0,0);
	tin.ic = reinterpret_cast<struct inGenerateIC *>(vin);
	tin.fft = fft;
	pltGenerateIC(pst,&tin,sizeof(tin),vout,pnOut);

	int myProc = mdlProc(pst->mdl);
	uint64_t nLocal = (int64_t)fft->rgrid->rn[myProc] * in->nGrid*in->nGrid;

	/* Expand the particles by adding an iOrder */
	assert(sizeof(expandParticle) >= sizeof(basicParticle));
	overlayedParticle  * pbBase = (overlayedParticle *)pkdParticleBase(pkd);
	int iz = fft->rgrid->rs[myProc] + fft->rgrid->rn[myProc];
	int iy=0, ix=0;
	float inGrid = 1.0 / in->nGrid;
	for(i=nLocal-1; i>=0; --i) {
	    basicParticle  *b = &pbBase->b + i;
	    basicParticle temp;
	    memcpy(&temp,b,sizeof(temp));
	    if (ix>0) --ix;
	    else {
		ix = in->nGrid-1;
		if (iy>0) --iy;
		else {
		    iy = in->nGrid-1;
		    --iz;
		    assert(iz>=0);
		    }
		}
#ifdef INTEGER_POSITION
	    // If we have no particle order convert directly to Integerized positions.
	    // We do this to save space as an "Integer" particle is small.
	    if (pkd->bNoParticleOrder) {
		integerParticle *p = &pbBase->i + i;
		p->v[2] = temp.v[2];
		p->v[1] = temp.v[1];
		p->v[0] = temp.v[0];
		p->r[2] = pkdDblToPos(pkd,temp.dr[2] + (iz+0.5) * inGrid - 0.5);
		p->r[1] = pkdDblToPos(pkd,temp.dr[1] + (iy+0.5) * inGrid - 0.5);
		p->r[0] = pkdDblToPos(pkd,temp.dr[0] + (ix+0.5) * inGrid - 0.5);
		}
	    else
#endif
	    {
		expandParticle *p = &pbBase->e + i;
		p->v[2] = temp.v[2];
		p->v[1] = temp.v[1];
		p->v[0] = temp.v[0];
		p->dr[2] = temp.dr[2];
		p->dr[1] = temp.dr[1];
		p->dr[0] = temp.dr[0];
		p->ix = ix;
		p->iy = iy;
		p->iz = iz;
		}
	    }
	assert(ix==0 && iy==0 && iz==fft->rgrid->rs[myProc]);
	/* Now we need to move excess particles between nodes so nStore is obeyed. */
	pkd->fft = fft; /* This is freed in pstMoveIC() */
	}
    else if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_GENERATEIC,in,nIn);
	pstGenerateIC(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,&outUp,pnOut);
	out->N += outUp.N;
	out->noiseMean += outUp.noiseMean;
	out->noiseCSQ += outUp.noiseCSQ;
	}
    if (pnOut) *pnOut = sizeof(struct outGenerateIC);
    }
