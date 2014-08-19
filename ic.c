#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "pkd.h"
#include "ic.h"

typedef uint32_t rand_t;
typedef uint64_t randl_t;
static const rand_t rand_max = 0xffffffff;
static const randl_t randl_max = 0xffffffffffffffff;
#define mt_table_size 624

typedef struct {
    rand_t mt[mt_table_size];
    int mi;
    } rand_state_t;

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


void generate(rand_state_t *state) {
    rand_t y;
    int i;
    for( i=0; i<mt_table_size; i++) {
	y = ((state->mt[i]>>31)&1) + (state->mt[(i+1)%mt_table_size] & 0x7fffffff);
	state->mt[i] = state->mt[(i + 397) % mt_table_size] ^ (y>>1);
	if ( y & 1)
	    state->mt[i] ^= 0x9908b0df;
        }
    state->mi = 0;
    }


static void initializeRNG(rand_state_t *state, int key_length, const rand_t *key) {
    int i, j, k;

    state->mt[0] = 19650218UL;
    for( i=1; i<mt_table_size; i++ )
	state->mt[i] = (1812433253 * (state->mt[i-1] ^ (state->mt[i-1]>>30)) + i) & 0xffffffff;
    state->mi = mt_table_size;

    i=1; j=0;
    k = (mt_table_size>key_length ? mt_table_size : key_length);
    for (; k; k--) {
	state->mt[i] = (state->mt[i] ^ ((state->mt[i-1] ^ (state->mt[i-1] >> 30)) * 1664525UL))
	    + key[j] + j; /* non linear */
	state->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
	i++; j++;
	if (i>=mt_table_size) { state->mt[0] = state->mt[mt_table_size-1]; i=1; }
	if (j>=key_length) j=0;
        }
    for (k=mt_table_size-1; k; k--) {
	state->mt[i] = (state->mt[i] ^ ((state->mt[i-1] ^ (state->mt[i-1] >> 30)) * 1566083941UL))
	    - i; /* non linear */
	state->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
	i++;
	if (i>=mt_table_size) { state->mt[0] = state->mt[mt_table_size-1]; i=1; }
        }
    
        state->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    }

rand_t rand32(rand_state_t *state) {
    rand_t y;
    
    if ( state->mi == mt_table_size ) generate(state);
    y = state->mt[state->mi++];
    y = y ^ (y>>11);
    y = y ^ ((y<<7) & 0x9d2c5680);
    y = y ^ ((y<<15) & 0xefc60000);
    y = y ^ (y>>18);
    return y;
    }

// Random number between 0 and 2^64 - 1
static randl_t randl(rand_state_t *state) {
    randl_t y1, y2;
    y1 = rand32(state);
    y2 = rand32(state);
    return y1 | (y2<<32);
    }

static double randd(rand_state_t *state) {
    return randl(state)/(randl_max+1.0);
    }

static void pairg( rand_state_t *state, double *y1, double *y2 ) {
    double x1, x2, w;
    // Instead of this:
    //   y1 = sqrt( - 2 ln(x1) ) cos( 2 pi x2 )
    //   y2 = sqrt( - 2 ln(x1) ) sin( 2 pi x2 )
    // we use this:
    do {
	x1 = 2.0 * randd(state) - 1.0;
	x2 = 2.0 * randd(state) - 1.0;
	w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 ); // Loop ~ 21.5% of the time
    w = sqrt( (-2.0 * log( w ) ) / w );
    *y1 = x1 * w;
    *y2 = x2 * w;
    }

static void mrandg( rand_state_t *state, int n, double *y ) {
    int i;
    double y2;
    int n1 = n & ~1;
    
    for( i=0; i<n1; i += 2 ) pairg(state,y+i,y+i+1);
    if ( n&1 ) pairg(state,y+n1, &y2 );
    }

#define k_loop_begin(i,j,k,idx,ex,sy,ey,sz,ez,YINIT,ZINIT)	\
    do { for(j=(sy); j<(ey); ++j) {				\
    YINIT							\
    for(k=0; k<(ez); ++k) {					\
    ZINIT							\
    for(i=0; i<(ex); ++i) {					\
    idx = mdlFFTkIdx(fft,i,j,k);				\
    do

#define k_loop_end while(0); }}} } while(0)

void pkdGenerateNoise(PKD pkd,MDLFFT fft,fftw_complex *ic,
    rand_t seed,int sy,int ey,int sz,int ez) {
    rand_state_t state;
    rand_t fullKey[3];
    int i,j,k,idx;


    fullKey[0] = seed;

//    k_loop_begin(i,j,k,idx,fft->kgrid->a1,sy,ey,sz,ez,
//	{fullKey[1] = j;},
//	{fullKey[2] = k; initializeRNG(&state,3,fullKey);}) {
//	pairg(&state,ic[idx]+0,ic[idx]+1);
//	} k_loop_end;

    for(j=sy; j<ey; ++j) {
	fullKey[1] = j;
	for(k=sz; k<ez; ++k) {
	    fullKey[2] = k;
	    initializeRNG(&state,3,fullKey);
	    for(i=0; i<=fft->kgrid->a1; ++i) {
		idx = mdlFFTkIdx(fft,i,j,k);
		pairg(&state,ic[idx]+0,ic[idx]+1);
		}
	    }
	}
    }

static int wrap(int v,int h,int m) {
    return v - (v >= h ? m : 0);
    }


#define RE 0
#define IM 0
void pkdGenerateIC(PKD pkd,MDLFFT fft,int iBegYr,int iEndYr,int iBegZk,int iEndZk,gridptr dic[]) {
    double twopi = 2.0 * 4.0 * atan(1.0);
    double itwopi = 1.0 / twopi;
    int i,j,k,sy,ey,idx;
    int ix, iy, iz;
    int nyx = fft->rgrid->n1 / 2;
    int nyy = fft->rgrid->n2 / 2;
    int nyz = fft->rgrid->n3 / 2;
    double ak, ak2, xfac, yfac, zfac;

    sy = fft->kgrid->s;
    ey = sy + fft->kgrid->n;

    /* Generate white noise realization -> dic[5] */
    pkdGenerateNoise(pkd,fft,dic[5].k,12345,sy,ey,iBegZk,iEndZk);

    /* Nyquist modes need to be real */
    dic[5].k[mdlFFTkIdx(fft,  0,  0,  0)][IM] = 0.0;
    dic[5].k[mdlFFTkIdx(fft,nyx,  0,  0)][IM] = 0.0;
    dic[5].k[mdlFFTkIdx(fft,  0,nyy,  0)][IM] = 0.0;
    dic[5].k[mdlFFTkIdx(fft,nyx,nyy,  0)][IM] = 0.0;
    dic[5].k[mdlFFTkIdx(fft,  0,  0,nyz)][IM] = 0.0;
    dic[5].k[mdlFFTkIdx(fft,nyx,  0,nyz)][IM] = 0.0;
    dic[5].k[mdlFFTkIdx(fft,  0,nyy,nyz)][IM] = 0.0;
    dic[5].k[mdlFFTkIdx(fft,nyx,nyy,nyz)][IM] = 0.0;

    for(j=sy; j<ey; ++j) {
	iy = wrap(j,nyy,fft->rgrid->n2);
	for(k=iBegZk; k<iEndZk; ++k) {
	    iz = wrap(j,nyy,fft->rgrid->n2);
	    for(i=0; i<=fft->kgrid->a1; ++i) {
		ix = wrap(j,nyy,fft->rgrid->n2);
		idx = mdlFFTkIdx(fft,i,j,k);
		ak2 = ix*ix + iy*iy + iz*iz;
		ak = sqrt(ak2);
		double amp = 1.0;

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


    if (mdlSelf(pkd->mdl)==0) printf("FFT\n");
    mdlIFFT( pkd->mdl, fft, dic[3].k );
    mdlIFFT( pkd->mdl, fft, dic[4].k );
    mdlIFFT( pkd->mdl, fft, dic[5].k );
    if (mdlSelf(pkd->mdl)==0) printf("FFT\n");


    }
