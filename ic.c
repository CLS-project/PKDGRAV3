#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
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


void pkdGenerateNoise(PKD pkd,MDLFFT fft,FFTW3(real) *ic,
    unsigned long seed,int sy,int ey,int sz,int ez) {
    RngStream g;
    unsigned long fullKey[6];
    int i,j,k,idx;

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
    RngStream_IncreasedPrecis (g, 1);
    RngStream_SetSeed(g,fullKey);


    double b = clockNow();

    for(k=sz; k<ez; ++k) {
	RngStream_ResetStartStream (g);
	RngStream_AdvanceState (g, 0, (1L<<40)*k );
	for(j=sy; j<ey; ++j) {
	    RngStream_ResetStartStream (g);
	    RngStream_AdvanceState (g, 0, (1L<<40)*k + (1L<<20)*j );
	    for(i=0; i<fft->rgrid->n1; i += 2) {
		idx = mdlFFTrIdx(pkd->mdl,fft,i,j,k);
		pairg(g,&ic[idx],&ic[idx+1]);
		}
	    }
	}
    RngStream_DeleteStream(&g);

    double a = clockNow();
    printf("Wallclock %f seconds\n", a-b);


    }

static int wrap(int v,int h,int m) {
    return v - (v >= h ? m : 0);
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



void pkdGenerateIC(PKD pkd,int iSeed,double dBoxSize,double dOmega0,double dLambda0,double a,
    MDLFFT fft,int iBegYr,int iEndYr,int iBegZk,int iEndZk,gridptr dic[],gridpos *pos) {
    MDL mdl = pkd->mdl;
    double twopi = 2.0 * 4.0 * atan(1.0);
    double itwopi = 1.0 / twopi;
    int i,j,k,sy,ey,sz,ez,idx;
    int ix, iy, iz;
    int nyx = fft->rgrid->n1 / 2;
    int nyy = fft->rgrid->n2 / 2;
    int nyz = fft->rgrid->n3 / 2;
    double iLbox = twopi / dBoxSize;
    double iLbox32 = pow(iLbox,1.5);
    double ak, ak2, xfac, yfac, zfac;
    double dOmega, dLambda, D0, Da;

    D0 = Growth(dOmega0,dLambda0,1,&dOmega,&dLambda);
    Da = Growth(dOmega0,dLambda0,a,&dOmega,&dLambda);

    printf("a=%g D0=%g Da=%g\n", a, D0, Da );
    printf("%g %g\n", Da / D0, dplus(a,dOmega0,dLambda0) / dplus(1.0,dOmega0,dLambda0) );


    sy = fft->kgrid->sSlab;
    ey = sy + fft->kgrid->nSlab;

    sz = fft->rgrid->sSlab;
    ez = sz + fft->rgrid->nSlab;

    /* Generate white noise realization -> dic[5] */
    if (mdlSelf(mdl)==0) printf("Generating noise in r-space\n");
    pkdGenerateNoise(pkd,fft,dic[5].r,iSeed,sy,ey,iBegZk,iEndZk);
    if (mdlSelf(mdl)==0) printf("Transforming to k-space\n");
    mdlFFT( mdl, fft, dic[5].r );

    for(j=sy; j<ey; ++j) {
	iy = wrap(j,nyy,fft->rgrid->n2);
	for(k=iBegZk; k<iEndZk; ++k) {
	    iz = wrap(j,nyy,fft->rgrid->n2);
	    for(i=0; i<=fft->kgrid->a1; ++i) {
		ix = wrap(j,nyy,fft->rgrid->n2);
		idx = mdlFFTkIdx(mdl,fft,i,j,k);
		ak2 = ix*ix + iy*iy + iz*iz;
		ak = sqrt(ak2);
		double amp = sqrt(1.0*(ak*iLbox) * iLbox32) * itwopi;

		amp /= ak2;
		xfac = amp * ix; // THESE CAN BE NEGATIVE!
		yfac = amp * iy;
		zfac = amp * iz;
		dic[0].k[idx][RE] = dic[5].k[idx][RE] * xfac;
		dic[0].k[idx][IM] = dic[5].k[idx][IM] * xfac;
		dic[1].k[idx][RE] = dic[5].k[idx][RE] * yfac;
		dic[1].k[idx][IM] = dic[5].k[idx][IM] * yfac;
		dic[2].k[idx][RE] = dic[5].k[idx][RE] * zfac;
		dic[2].k[idx][IM] = dic[5].k[idx][IM] * zfac;
		dic[3].k[idx][RE] = dic[0].k[idx][RE];
		dic[3].k[idx][IM] = dic[0].k[idx][IM];
		dic[4].k[idx][RE] = dic[1].k[idx][RE];
		dic[4].k[idx][IM] = dic[1].k[idx][IM];
		dic[5].k[idx][RE] = dic[2].k[idx][RE];
		dic[5].k[idx][IM] = dic[2].k[idx][IM];
		}
	    }
	}
    if (iBegZk==0) {
	for(i=0; i<6; ++i) dic[i].k[0][RE] = dic[i].k[0][IM] = 0.0;
	}


    if (mdlSelf(mdl)==0) printf("FFT 1\n");
    mdlIFFT(mdl, fft, dic[3].k );
    if (mdlSelf(mdl)==0) printf("FFT 2\n");
    mdlIFFT(mdl, fft, dic[4].k );
    if (mdlSelf(mdl)==0) printf("FFT 3\n");
    mdlIFFT(mdl, fft, dic[5].k );
    if (mdlSelf(mdl)==0) printf("FFT done\n");

    idx = 0;
    for(k=sz; k<ez; ++k) {
	for(j=iBegYr; j<iEndYr; ++j) {
	    for(i=0; i<=fft->rgrid->n1; ++i) {
		pos[idx].x = dic[3].r[mdlFFTrIdx(mdl, fft,  i,  j,  k)];
		pos[idx].y = dic[4].r[mdlFFTrIdx(mdl, fft,  i,  j,  k)];
		pos[idx].z = dic[5].r[mdlFFTrIdx(mdl, fft,  i,  j,  k)];
		++idx;
		}
	    }
	}

    }
