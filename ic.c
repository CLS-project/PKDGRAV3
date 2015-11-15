#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "pkd.h"
#include "ic.h"
#include "RngStream.h"
#define RE 0
#define IM 1

typedef struct {
    double omegam;
    double omegav;
    } ddplus_ctx;

static double ddplus(double a,void *ctx) {
    ddplus_ctx *dc = ctx;
    double eta;
    if ( a == 0.0 ) return 0.0;
    eta = sqrt(dc->omegam/a + dc->omegav*a*a + 1.0 - dc->omegam - dc->omegav);
    return 2.5/(eta*eta*eta);
    }
static double dplus(double a,double omegam,double omegav) {
    double eta;
    ddplus_ctx ddctx;
    gsl_function F;
    double result, error;
    gsl_integration_workspace *W;

    ddctx.omegam = omegam;
    ddctx.omegav = omegav;
    eta = sqrt(omegam/a + omegav*a*a + 1.0 - omegam - omegav);
    W = gsl_integration_workspace_alloc (1000);
    F.function = &ddplus;
    F.params = &ddctx;
    gsl_integration_qag(&F, 0.0, a, 0.0, 1e-8, 1000,
        GSL_INTEG_GAUSS61, W, &result, &error); 
    gsl_integration_workspace_free(W);
    return eta/a * result;
    }

typedef struct {
    double m_scale, m_lasta;

    double m_f_baryon,m_k_equality,m_sound_horizon,
	m_beta_c,m_alpha_c,m_beta_node,m_alpha_b,m_beta_b,
	m_k_silk;
    double m_omega, m_omegab, m_omegav, m_h0, m_sigma8, m_ns;
    double m_pnorm;

    } PowerEHU97;

static double TF_pressureless(double q, double a, double b) {
    double v;
    v = log(M_E+1.8*b*q);
    return v/(v + (14.2/a + 386/(1.+69.9*pow(q,1.08)))*pow(q,2));
    }

static double TFmdm_onek_mpc(PowerEHU97 *power,double k) {

    //  Calculate transfer function
    //
    //  Input: 
    //!  k -- wavenumber in Mpc^{-1}  
    //!        omhh -- The density of CDM and baryons, in units of critical dens,
    //!                multiplied by the square of the Hubble constant, in units
    //!                of 100 km/s/Mpc
    //!        m_f_baryon -- The fraction of baryons to CDM
    //! 
    //!  Output:
    //!  tf_full -- The full fitting formula, eq. (16), for the matter
    //!             transfer function. 
    //!  tf_baryon -- The baryonic piece of the full fitting formula, eq. 21.
    //!  tf_cdm -- The CDM piece of the full fitting formula, eq. 17.
    //!
    //!       Tcmb -- The temperature of the CMB in Kelvin, 2.728(4) is COBE and is
    //!         the default reached by inputing Tcmb=0 -- reset on output. */
    //!
    //! theta_cmb,      /* Tcmb in units of 2.7 K */ 
    //! z_equality,     /* Redshift of matter-radiation equality, really 1+z */
    //! k_equality,     /* Scale of equality, in Mpc^-1 */
    //! z_drag,         /* Redshift of drag epoch */
    //! R_drag,         /* Photon-baryon ratio at drag epoch */
    //! R_equality,     /* Photon-baryon ratio at equality epoch */
    //! sound_horizon,  /* Sound horizon at drag epoch, in Mpc */
    //! k_silk,         /* Silk damping scale, in Mpc^-1 */
    //! alpha_c,        /* CDM suppression */
    //! beta_c,         /* CDM log shift */
    //! alpha_b,        /* Baryon suppression */
    //! beta_b,         /* Baryon envelope shift */


    double tf_full,tf_cdm,q,ks,tf_baryon,s_tilde;

    assert(k>0.0);

    //!  Auxiliary Variables
           
    q = k/13.41/power->m_k_equality;
    ks = k*power->m_sound_horizon;

    //!  Main Variables

    tf_cdm = 1./(1.+pow(ks/5.4,4.));
    tf_cdm = tf_cdm*TF_pressureless(q,1.,power->m_beta_c) +
        (1.-tf_cdm)*TF_pressureless(q,power->m_alpha_c,power->m_beta_c);

    s_tilde = power->m_sound_horizon/pow(1.+pow(power->m_beta_node/ks,3.),1./3.) ;
    tf_baryon = TF_pressureless(q,1.,1.)/(1.+pow(ks/5.2,2.));
    tf_baryon = tf_baryon + power->m_alpha_b/(1.+pow(power->m_beta_b/ks,3))
        *exp(-pow(k/power->m_k_silk,1.4));
    tf_baryon = tf_baryon *(sin(k*s_tilde)/(k*s_tilde));
    tf_full = power->m_f_baryon*tf_baryon + (1-power->m_f_baryon)*tf_cdm;

    return tf_full;
    }

double PowerEHU97pmat(PowerEHU97 *power,double ak, double a) {
    if ( ak <= 0.0 ) return 0.0;
    double t = TFmdm_onek_mpc(power,ak);
    double p,tpois;


    //  Apply transfer function to primordial power spectrum.
    //  Primordial spectrum of psi (or entropy, in the isocurvature case):
    p = power->m_pnorm*pow(ak,power->m_ns-4.0);
    //  Apply transfer function to get spectrum of phi at a=1.
    p = p * t*t;
    //  Convert to density fluctuation power spectrum.  Note that k^2 is
    //  corrected for an open universe.  Scale to a using linear theory.
    tpois = -(2.0/3.0) / power->m_omega*(pow(ak*2.99793e5/power->m_h0*100,2.0) - 4.*(power->m_omega-1.0));
    if ( power->m_scale==0.0 || power->m_lasta!=a) {
        power->m_scale = dplus(a,power->m_omega,power->m_omegav)/dplus(1.0,power->m_omega,power->m_omegav);
        power->m_lasta = a;
        }

    return p*tpois*tpois*power->m_scale*power->m_scale;
    }

static double dsig8_gsl(double ak, void * params) {
    PowerEHU97 *power = params;
    double x,w;
    if (ak <= 0.0) return 0.0;
    // Window function for spherical tophat of radius 8 Mpc/h.
    x = ak*800.0/(power->m_h0/100);
    w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
    return ak*ak* PowerEHU97pmat(power,ak,1.0)*w*w;
}

static void calculateNormalization(PowerEHU97 *power) {
    power->m_pnorm = 1.0;
    assert(power->m_sigma8>0.0);

    gsl_function F;
    double result, error;
    gsl_integration_workspace *W = gsl_integration_workspace_alloc (1000);
    F.function = &dsig8_gsl;
    F.params = power;
    gsl_error_handler_t *oldhand = gsl_set_error_handler_off();
    gsl_integration_qag(&F, 0.0, 10.0, 0.0, 1e-6, 1000,
        GSL_INTEG_GAUSS61, W, &result, &error); 
    gsl_set_error_handler(oldhand);
    power->m_pnorm = power->m_sigma8*power->m_sigma8/((4.0*M_PI)*result);
    gsl_integration_workspace_free(W);
    }

void PowerEHU97Setup(PowerEHU97 *power,double h0,double omega,double omegab,double Tcmb,double sigma8,double ns) {
    double omhh = omega * h0 * h0;
    double y,obhh;
    double theta_cmb,z_equality,z_drag,R_drag,R_equality;

    power->m_scale = 0.0;
    power->m_lasta = 0.0;
    power->m_omega = omega;
    power->m_omegab = omegab;
    power->m_omegav = 1.0 - omega;
    power->m_h0 = h0;
    power->m_ns = ns;
    power->m_sigma8 = sigma8;

    power->m_f_baryon = omegab / omega;
    if (power->m_f_baryon <= 0.0) power->m_f_baryon = 1.e-5;
    assert(omhh>0.0);

    //! Auxiliary variables
    obhh = omhh * power->m_f_baryon;
    theta_cmb = Tcmb / 2.7;

    //! Main variables
    z_equality = 2.50e4 * omhh * pow(theta_cmb,-4.) - 1.0;
    power->m_k_equality = 0.0746 * omhh * pow(theta_cmb,-2.); 

    z_drag = 0.313 * pow(omhh,-0.419) * (1.+0.607 * pow(omhh,0.674));
    z_drag = 1e0 + z_drag * pow(obhh,0.238 * pow(omhh,0.223));
    z_drag = 1291e0 * pow(omhh,0.251)/
        (1e0 + 0.659 * pow(omhh,0.828)) * z_drag;

    R_drag = 31.5 * obhh * pow(theta_cmb,-4.)*1000e0/(1e0 + z_drag) ;
    R_equality = 31.5*obhh * pow(theta_cmb,-4.)*1000e0/(1e0 + z_equality) ;

    power->m_sound_horizon = 2./3./power->m_k_equality*sqrt(6./R_equality)*
        log(( sqrt(1.+R_drag)+sqrt(R_drag+R_equality) )
            /(1.+sqrt(R_equality)));

    power->m_k_silk = 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1e0 + pow(10.4*omhh,-0.95));

    power->m_alpha_c = (pow(46.9*omhh,0.670)*(1e0+pow(32.1*omhh,-0.532)));
    power->m_alpha_c = pow(power->m_alpha_c,-power->m_f_baryon);
    power->m_alpha_c = power->m_alpha_c*pow(pow(12.0*omhh,0.424)*(1e0 +
            pow(45.0*omhh,-0.582)),-pow(power->m_f_baryon,3.));


    power->m_beta_c = 0.944/(1+pow(458.*omhh,-0.708));
    power->m_beta_c = 1.+power->m_beta_c*(pow(1.-power->m_f_baryon,pow(0.395*omhh,-0.0266)) - 1e0);
    power->m_beta_c = 1./power->m_beta_c;

    y = (1e0+z_equality)/(1e0+z_drag);
    power->m_alpha_b = y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)
            /(sqrt(1.+y)-1.)));
    power->m_alpha_b = 2.07*power->m_k_equality*power->m_sound_horizon*
        pow(1.+R_drag,-0.75)*power->m_alpha_b;


    power->m_beta_b = 0.5+power->m_f_baryon+(3.-2.*power->m_f_baryon)*
        sqrt(pow(17.2*omhh,2.+1e0));

    power->m_beta_node = 8.41*pow(omhh,0.435);

    calculateNormalization(power);

    }

typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double *tk, *tf;
    double spectral;
    double normalization;
    int nTf;
    } powerParameters;

static double power(powerParameters *P,double k, double a) {
    double T = gsl_spline_eval(P->spline,log(k),P->acc);
    return pow(k,P->spectral) * P->normalization * T * T / P->tf[0] / P->tf[0];
    }

typedef struct {
    powerParameters *P;
    double r;
    } varianceParameters;

static double variance_integrand(double ak, void * params) {
    varianceParameters *vprm = params;
    double x, w;
    // Window function for spherical tophat of given radius
    x = ak * vprm->r;
    w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
    return power(vprm->P,ak,1.0)*ak*ak*w*w*4.0*M_PI;
    }

static double variance(powerParameters *P,double dRadius,double dSigma8) {
    varianceParameters vprm;
    gsl_function F;
    double result, error;
    gsl_integration_workspace *W = gsl_integration_workspace_alloc (1000);
    vprm.P = P;
    vprm.r = dRadius; /* 8 Mpc/h for example */
    F.function = &variance_integrand;
    F.params = &vprm;
    gsl_integration_qag(&F, exp(P->tk[0]), exp(P->tk[P->nTf-1]),
	0.0, 1e-6, 1000, GSL_INTEG_GAUSS61, W, &result, &error);
    gsl_integration_workspace_free(W);
    return result;
    }

static void pairg( RngStream g, FFTW3(real) *y1, FFTW3(real) *y2 ) {
    double x1, x2, w;
    // Instead of this:
    //   y1 = sqrt( - 2 ln(x1) ) cos( 2 pi x2 )
    //   y2 = sqrt( - 2 ln(x1) ) sin( 2 pi x2 )
    // we use this:
    do {
	x1 = 2.0 * RngStream_RandU01(g) - 1.0;
	x2 = 2.0 * RngStream_RandU01(g) - 1.0;
	w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 || w == 0.0 ); // Loop ~ 21.5% of the time
    w = sqrt( (-2.0 * log( w ) ) / w );
    *y1 = x1 * w;
    *y2 = x2 * w;
    }

static void mrandg( RngStream g, int n, FFTW3(real) *y ) {
    int i;
    FFTW3(real) y2;
    int n1 = n & ~1;
    
    for( i=0; i<n1; i += 2 ) pairg(g,y+i,y+i+1);
    if ( n&1 ) pairg(g,y+n1, &y2 );
    }

static double clockNow() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec+(tv.tv_usec*1e-6));
    }


static int wrap(int v,int h,int m) {
    return v - (v >= h ? m : 0);
    }

static int getk(int v,int h,int m) {
    return v <= h ? v : m - v;
    }

void pkdGenerateNoise(PKD pkd,unsigned long seed,MDLFFT fft,float complex *ic) {
    MDL mdl = pkd->mdl;
    const int nGrid = fft->kgrid->n3;
    const int iNyquist = nGrid / 2;
    RngStream g;
    unsigned long fullKey[6];
    int j,k,jj,kk;
    float v1,v2;

    mdlGridCoord kfirst, klast, kindex;
    mdlGridCoordFirstLast(mdl,fft->kgrid,&kfirst,&klast);

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

    j = k = nGrid; /* Start with invalid values so we advance the RNG correctly. */
    for( kindex=kfirst; !mdlGridCoordCompare(&kindex,&klast); mdlGridCoordIncrement(&kindex) ) {
	if (j!=kindex.z || k!=kindex.y) {
	    assert(kindex.x==0);
	    j = kindex.z; /* Remember: z and y indexes are permuted in k-space */
	    k = kindex.y;

	    jj = j<=iNyquist ? j*2 : (nGrid-j)*2 + 1;
	    kk = k<=iNyquist ? k*2 : (nGrid-k)*2 + 1;

	    RngStream_ResetStartStream (g);

	    /* We need to generate the correct complex conjugates */
	    if (k && j > iNyquist) {
		int jjc = j<iNyquist ? j*2 + 1 : (nGrid-j)*2;
		int kkc = k<iNyquist ? k*2 + 1 : (nGrid-k)*2;
		RngStream_AdvanceState (g, 0, (1L<<40)*(jjc) + (1L<<20)*(kkc) );
		pairg(g,&v1,&v2);
		ic[kindex.i+iNyquist] = v1 - I * v2;
		pairg(g,&v1,&v2);
		ic[kindex.i] = v1 - I * v2;
		RngStream_ResetStartStream (g);
		RngStream_AdvanceState (g, 0, (1L<<40)*jj + (1L<<20)*kk );
		pairg(g,&v1,&v2);
		pairg(g,&v1,&v2);
		}
	    else {
		RngStream_AdvanceState (g, 0, (1L<<40)*jj + (1L<<20)*kk );
		pairg(g,&v1,&v2);
		ic[kindex.i+iNyquist] = v1 + I * v2;
		pairg(g,&v1,&v2);
		ic[kindex.i] = v1 + I * v2;
		}
	    }
	else if (kindex.i!=iNyquist) {
	    pairg(g,&v1,&v2);
	    ic[kindex.i] = v1 + I * v2;
	    }
	if (kindex.x%iNyquist + kindex.y%iNyquist + kindex.z%iNyquist == 0) {
	    ic[kindex.i] = creal(ic[kindex.i]);
	    }
	}
    RngStream_DeleteStream(&g);
    if (kfirst.x+kfirst.y+kfirst.z == 0)
	ic[0] = 0;
    }

/* Approximation */
static double Growth(double Om0, double OL0, double a, double *Om, double *OL) {
    double Hsq,OK0;
    OK0 = 1 - Om0 - OL0;
    Hsq = Om0 / (a*a*a) + OK0 / (a*a) + OL0;
    *Om = Om0 / (a*a*a*Hsq);
    *OL = OL0 / Hsq;
    return 2.5 * a * *Om / (pow(*Om,4./7.) - *OL + (1.+0.5* *Om)*(1.+ *OL/70.));
    }

int pkdGenerateIC(PKD pkd,int iSeed,int nGrid,int b2LPT,double dBoxSize,
    double dOmega0,double dLambda0,double dSigma8,double dSpectral,double h,
    double a,int nTf, double *tk, double *tf) {
    MDL mdl = pkd->mdl;
    double twopi = 2.0 * 4.0 * atan(1.0);
    double itwopi = 1.0 / twopi;
    float inGrid = 1.0 / nGrid;
    float fftNormalize = inGrid*inGrid*inGrid;
    int i,j,k,idx;
    int ix, iy, iz;
    int iNyquist = nGrid / 2;
    double iLbox = twopi / dBoxSize;
    double iLbox32 = pow(iLbox,3.0) / 2.0;
    double ak, ak2, xfac, yfac, zfac, amp;
    double dOmega, dLambda, D0, Da;
    double velFactor;
    float f1, f2;
    basicParticle *p;
    int nLocal;

    MDLFFT fft;
    mdlGridCoord kfirst, klast, kindex;
    mdlGridCoord rfirst, rlast, rindex;
    gridptr ic[10];

    powerParameters P;

    D0 = Growth(dOmega0,dLambda0,1,&dOmega,&dLambda);
    Da = Growth(dOmega0,dLambda0,a,&dOmega,&dLambda);
    dSigma8 *= Da/D0;

    P.normalization = 1e4;
    P.spectral = dSpectral;
    P.nTf = nTf;
    P.tk = tk;
    P.tf = tf;
    P.acc = gsl_interp_accel_alloc();
    P.spline = gsl_spline_alloc (gsl_interp_cspline, nTf);
    gsl_spline_init(P.spline, P.tk, P.tf, P.nTf);
    P.normalization *= dSigma8*dSigma8 / variance(&P,8.0,dSigma8);
    f1 = pow(dOmega,5.0/9.0);
    f2 = 2.0 * pow(dOmega,6.0/11.0);

    velFactor = (sqrt(dOmega0/(a*a*a) + (1-dOmega0-dLambda0)/(a*a) + dLambda0)) * sqrt(8.0/3.0*M_PI);
    velFactor *= a*a; /* Comoving */

    mdlFFTInitialize(mdl,&fft,nGrid,nGrid,nGrid,0,0);
    mdlGridCoordFirstLast(mdl,fft->kgrid,&kfirst,&klast);
    mdlGridCoordFirstLast(mdl,fft->rgrid,&rfirst,&rlast);

    /* The mdlSetArray will use the values from thread 0 */
    ic[0].r = (FFTW3(real)*)pkdParticleBase(pkd);
    ic[1].r = ic[0].r + fft->rgrid->nLocal;
    ic[2].r = ic[1].r + fft->rgrid->nLocal;
    ic[3].r = ic[2].r + fft->rgrid->nLocal;
    ic[4].r = ic[3].r + fft->rgrid->nLocal;
    ic[5].r = ic[4].r + fft->rgrid->nLocal;
    ic[6].r = ic[5].r + fft->rgrid->nLocal;
    ic[7].r = ic[6].r + fft->rgrid->nLocal;
    ic[8].r = (FFTW3(real)*)pkd->pLite;;
    ic[9].r = ic[8].r + fft->rgrid->nLocal;

    ic[0].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[0].r);
    ic[1].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[1].r);
    ic[2].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[2].r);
    ic[3].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[3].r);
    ic[4].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[4].r);
    ic[5].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[5].r);
    ic[6].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[6].r);
    ic[7].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[7].r);
    ic[8].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[8].r);
    ic[9].r = mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),ic[9].r);

    /* Particles will overlap ic[0] through ic[5] eventually */
    nLocal = rlast.i / fft->rgrid->a1 * fft->rgrid->n1;
    p = mdlSetArray(pkd->mdl,nLocal,sizeof(basicParticle),pkdParticleBase(pkd));

    /* Generate white noise realization -> ic[5] */
    if (mdlSelf(mdl)==0) {printf("Generating random noise\n"); fflush(stdout); }
    pkdGenerateNoise(pkd,iSeed,fft,ic[5].k);

    if (mdlSelf(mdl)==0) {printf("Imprinting power\n"); fflush(stdout); }
    for( kindex=kfirst; !mdlGridCoordCompare(&kindex,&klast); mdlGridCoordIncrement(&kindex) ) {
	/* Range: [0,iNyquist] */
	iy = getk(kindex.z,iNyquist,fft->rgrid->n3);
	iz = getk(kindex.y,iNyquist,fft->rgrid->n2);
	ix = getk(kindex.x,iNyquist,fft->rgrid->n1);
	idx = kindex.i;

	ak2 = ix*ix + iy*iy + iz*iz;
	if (ak2>0) {
	    ak = sqrt(ak2) * iLbox;
	    amp = sqrt(power(&P,ak,a) * iLbox32) * itwopi / ak2;
	    }
	else amp = 0.0;
	/*
	** Note: noise for Nyquist is complex conjugate, so:
	** (a + ib) * +x * -i = (b - ia) * x
	** (a - ib) * -x * -i = (b + ia) * x [conjugate]
	*/
	xfac = amp * ix;
	yfac = amp * iy;
	zfac = amp * iz;

	ic[7].k[idx] = ic[5].k[idx] * xfac;
	ic[8].k[idx] = ic[5].k[idx] * yfac;
	ic[9].k[idx] = ic[5].k[idx] * zfac;

	if (b2LPT) {
	    ic[0].k[idx] = ic[7].k[idx] * ix * twopi * -I; /* xx */
	    ic[1].k[idx] = ic[8].k[idx] * iy * twopi * -I; /* yy */
	    ic[2].k[idx] = ic[9].k[idx] * iz * twopi * -I; /* zz */
	    ic[3].k[idx] = ic[7].k[idx] * iy * twopi * -I; /* xy */
	    ic[4].k[idx] = ic[8].k[idx] * iz * twopi * -I; /* yz */
	    ic[5].k[idx] = ic[9].k[idx] * ix * twopi * -I; /* zx */
	    }
	}

    if (mdlSelf(mdl)==0) {printf("Generating x velocities\n"); fflush(stdout); }
    mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[7].k );
    if (mdlSelf(mdl)==0) {printf("Generating y velocities\n"); fflush(stdout); }
    mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[8].k );
    if (mdlSelf(mdl)==0) {printf("Generating z velocities\n"); fflush(stdout); }
    mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[9].k );

    if (b2LPT) {
	if (mdlSelf(mdl)==0) {printf("Generating xx term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[0].k );
	if (mdlSelf(mdl)==0) {printf("Generating yy term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[1].k );
	if (mdlSelf(mdl)==0) {printf("Generating zz term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[2].k );
	if (mdlSelf(mdl)==0) {printf("Generating xy term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[3].k );
	if (mdlSelf(mdl)==0) {printf("Generating yz term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[4].k );
	if (mdlSelf(mdl)==0) {printf("Generating xz term\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[5].k );

	/* Calculate the source term */
	if (mdlSelf(mdl)==0) {printf("Generating source term\n"); fflush(stdout); }
	for( rindex=rfirst; !mdlGridCoordCompare(&rindex,&rlast); mdlGridCoordIncrement(&rindex) ) {
	    int i = rindex.i;
	    ic[6].r[i] = ic[0].r[i]*ic[1].r[i] + ic[0].r[i]*ic[2].r[i] + ic[1].r[i]*ic[2].r[i]
		- ic[3].r[i]*ic[3].r[i] - ic[4].r[i]*ic[4].r[i] - ic[5].r[i]*ic[5].r[i];
	    }
	mdlFFT(mdl, fft, ic[6].r );
	}

    /* Move the 1LPT positions/velocities to the particle area */
    idx = 0;
    for( rindex=rfirst; !mdlGridCoordCompare(&rindex,&rlast); mdlGridCoordIncrement(&rindex) ) {
	float x = ic[7].r[rindex.i];
	float y = ic[8].r[rindex.i];
	float z = ic[9].r[rindex.i];
	p[idx].r[0] = (rindex.x+0.5) * inGrid + x - 0.5;
	p[idx].r[1] = (rindex.y+0.5) * inGrid + y - 0.5;
	p[idx].r[2] = (rindex.z+0.5) * inGrid + z - 0.5;
	p[idx].v[0] = f1 * x * velFactor;
	p[idx].v[1] = f1 * y * velFactor;
	p[idx].v[2] = f1 * z * velFactor;
	++idx;
	}
    assert(idx == nLocal);

    if (b2LPT) {
	for(kindex=kfirst; !mdlGridCoordCompare(&kindex,&klast); mdlGridCoordIncrement(&kindex)) {
	    iy = wrap(kindex.z,iNyquist,fft->rgrid->n3);
	    iz = wrap(kindex.y,iNyquist,fft->rgrid->n2);
	    ix = wrap(kindex.x,iNyquist,fft->rgrid->n1);
	    idx = kindex.i;
	    ak2 = ix*ix + iy*iy + iz*iz;
	    ak = sqrt(ak2);

	    if (ak2>0.0) {
		ic[7].k[idx] = -3.0/ak2 * ix * ic[6].k[idx] * -I;
		ic[8].k[idx] = -3.0/ak2 * iy * ic[6].k[idx] * -I;
		ic[9].k[idx] = -3.0/ak2 * iz * ic[6].k[idx] * -I;
		}
	    else {
		ic[7].k[idx] = 0.0;
		ic[8].k[idx] = 0.0;
		ic[9].k[idx] = 0.0;
		}
	    }

	if (mdlSelf(mdl)==0) {printf("Generating x2 velocities\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[7].k );
	if (mdlSelf(mdl)==0) {printf("Generating y2 velocities\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[8].k );
	if (mdlSelf(mdl)==0) {printf("Generating z2 velocities\n"); fflush(stdout); }
	mdlIFFT(mdl, fft, (FFTW3(complex)*)ic[9].k );

	/* Add the 2LPT positions/velocities corrections to the particle area */
	idx = 0;
	for( rindex=rfirst; !mdlGridCoordCompare(&rindex,&rlast); mdlGridCoordIncrement(&rindex) ) {
	    float x = ic[7].r[rindex.i] * fftNormalize;
	    float y = ic[8].r[rindex.i] * fftNormalize;
	    float z = ic[9].r[rindex.i] * fftNormalize;

	    p[idx].r[0] += x;
	    p[idx].r[1] += y;
	    p[idx].r[2] += z;
	    p[idx].v[0] += f2 * x * velFactor;
	    p[idx].v[1] += f2 * y * velFactor;
	    p[idx].v[2] += f2 * z * velFactor;
	    ++idx;
	    }
	}
    gsl_spline_free(P.spline);
    gsl_interp_accel_free(P.acc);

    return nLocal;
    }
