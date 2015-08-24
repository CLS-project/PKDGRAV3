#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include "ewald.h"
#include "pkd.h"
#include "meval.h"
#include "qeval.h"
#include "moments.h"
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
#ifdef __AVX__
#include "avx_mathfun.h"
#else
#include "sse_mathfun.h"
#endif
#endif

#ifdef __INTEL_COMPILER
#ifndef USE_SVML
/*#define USE_SVML*/
#endif
#endif

static int evalEwald(struct EwaldVariables *ew,double *ax, double *ay, double *az, double *fPot,
    double x, double y, double z,
    double g0, double g1, double g2, double g3,double g4, double g5) {
    const MOMC * restrict mom = &ew->mom;
    double onethird = 1.0/3.0;
    double xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
    double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
    double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;

    xx = 0.5*x*x;
    xxx = onethird*xx*x;
    xxy = xx*y;
    xxz = xx*z;
    yy = 0.5*y*y;
    yyy = onethird*yy*y;
    xyy = yy*x;
    yyz = yy*z;
    zz = 0.5*z*z;
    zzz = onethird*zz*z;
    xzz = zz*x;
    yzz = zz*y;
    xy = x*y;
    xyz = xy*z;
    xz = x*z;
    yz = y*z;
    Q2mirx = mom->xx*x + mom->xy*y + mom->xz*z;
    Q2miry = mom->xy*x + mom->yy*y + mom->yz*z;
    Q2mirz = mom->xz*x + mom->yz*y + mom->zz*z;
    Q3mirx = mom->xxx*xx + mom->xxy*xy + mom->xxz*xz + mom->xyy*yy + mom->xyz*yz + mom->xzz*zz;
    Q3miry = mom->xxy*xx + mom->xyy*xy + mom->xyz*xz + mom->yyy*yy + mom->yyz*yz + mom->yzz*zz;
    Q3mirz = mom->xxz*xx + mom->xyz*xy + mom->xzz*xz + mom->yyz*yy + mom->yzz*yz + mom->zzz*zz;
    Q4mirx = mom->xxxx*xxx + mom->xxxy*xxy + mom->xxxz*xxz + mom->xxyy*xyy + mom->xxyz*xyz +
	mom->xxzz*xzz + mom->xyyy*yyy + mom->xyyz*yyz + mom->xyzz*yzz + mom->xzzz*zzz;
    Q4miry = mom->xxxy*xxx + mom->xxyy*xxy + mom->xxyz*xxz + mom->xyyy*xyy + mom->xyyz*xyz +
	mom->xyzz*xzz + mom->yyyy*yyy + mom->yyyz*yyz + mom->yyzz*yzz + mom->yzzz*zzz;
    Q4mirz = mom->xxxz*xxx + mom->xxyz*xxy + mom->xxzz*xxz + mom->xyyz*xyy + mom->xyzz*xyz +
	mom->xzzz*xzz + mom->yyyz*yyy + mom->yyzz*yyz + mom->yzzz*yzz + mom->zzzz*zzz;
    Q4x = ew->Q4xx*x + ew->Q4xy*y + ew->Q4xz*z;
    Q4y = ew->Q4xy*x + ew->Q4yy*y + ew->Q4yz*z;
    Q4z = ew->Q4xz*x + ew->Q4yz*y + ew->Q4zz*z;
    Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ew->Q3x*x + ew->Q3y*y + ew->Q3z*z) + ew->Q4;
    Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
    Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
    Qta = g1*mom->m - g2*ew->Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
    *fPot -= g0*mom->m - g1*ew->Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
    *ax += g2*(Q2mirx - ew->Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
    *ay += g2*(Q2miry - ew->Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
    *az += g2*(Q2mirz - ew->Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;

    return 447;
    }

/* Once CUDA has completed, we need to acumulate */
void pkdAccumulateCUDA(PKD pkd,workEwald *we,double *pax,double *pay,double *paz,double *pot) {
    int i;

    for(i=0; i<we->nP; ++i) {
	workParticle *wp = we->ppWorkPart[i];
	int pi = we->piWorkPart[i];
	PINFOOUT *out = &wp->pInfoOut[pi];
	out->a[0] += pax[i];
	out->a[1] += pay[i];
	out->a[2] += paz[i];
	out->fPot += pot[i];
	pkdParticleWorkDone(wp);
	}
   }


#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
static const struct CONSTS {
    vdouble onequarter,onethird,half,one,two,three,five,seven,nine;
    vdouble log2e,C1,C2;
    vdouble threshold0,threshold1,threshold2;
#if defined(__AVX__)
    vdouble p0ab,p0cd,p1ab,p1cd,p2ab,p2cd,p3ab,p3cd,p4ab,p4cd,p5ab,p5cd;
    vdouble q0ab,q0cd,q1ab,q1cd,q2ab,q2cd,q3ab,q3cd,q4ab,q4cd,q5ab,q5cd;
#else
    vdouble p0a,p0b,p0c,p0d, p1a,p1b,p1c,p1d, p2a,p2b,p2c,p2d;
    vdouble p3a,p3b,p3c,p3d, p4a,p4b,p4c,p4d, p5a,p5b,p5c,p5d;
    vdouble q0a,q0b,q0c,q0d, q1a,q1b,q1c,q1d, q2a,q2b,q2c,q2d;
    vdouble q3a,q3b,q3c,q3d, q4a,q4b,q4c,q4d, q5a,q5b,q5c,q5d;
#endif
    vdouble p0,p1,p2,q0,q1,q2,q3;
    } consts = {
        {SIMD_DCONST(0.25)},
        {SIMD_DCONST(1.0/3.0)},
        {SIMD_DCONST(0.5)},
        {SIMD_DCONST(1.0)},
        {SIMD_DCONST(2.0)},
        {SIMD_DCONST(3.0)},
        {SIMD_DCONST(5.0)},
        {SIMD_DCONST(7.0)},
        {SIMD_DCONST(9.0)},

#ifndef USE_SVML
	{SIMD_DCONST(1.4426950408889634073599)},     /* 1/log(2) */
	{SIMD_DCONST(6.93145751953125E-1)},
	{SIMD_DCONST(1.42860682030941723212E-6)},

        {SIMD_DCONST(0.65)},
        {SIMD_DCONST(2.2)},
        {SIMD_DCONST(6.0)},

#ifdef __AVX__
#define ERF_CONSTS(d,c,b,a) {{a,b,a,b}},{{c,d,c,d}}
#elif defined(slower__SSE4_1__)
#define ERF_CONSTS(d,c,b,a) {a,a},{b,b},{c,c},{d,d}
#else
#define ERF_CONSTS(d,c,b,a) {a,a},{b-a,b-a},{c-b,c-b},{d-c,d-c}
#endif
	/* erf/erfc polynomial tables */
	/*         [6.0,)                  [2.2,6.0)               [0.65,2.2)              [0,0.64)   */
	ERF_CONSTS(0.00000000000000000e+0, 2.25716982919217555e-2, 7.06940843763253131e-3, 0.00000000000000000e+0),
	ERF_CONSTS(8.08040729052301677e+0, 1.57289620742838702e-1, 7.14193832506776067e-2, 6.49254556481904354e-5),
	ERF_CONSTS(4.77209965874436377e+1, 5.81528574177741135e-1, 3.31899559578213215e-1, 1.20339380863079457e-3),
	ERF_CONSTS(3.84683103716117320e+1, 1.26739901455873222e+0, 8.78115804155881782e-1, 4.03259488531795274e-2),
	ERF_CONSTS(8.80253746105525775e+0, 1.62356584489366647e+0, 1.33154163936765307e+0, 1.35894887627277916e-1),
	ERF_CONSTS(5.64189583547756078e-1, 9.99921140009714409e-1, 9.99999992049799098e-1, 1.12837916709551256e+0),

	ERF_CONSTS(0.00000000000000000e+0, 4.00072964526861362e-2, 1.25304936549413393e-2, 0.00000000000000000e+0),
	ERF_CONSTS(0.00000000000000000e+0, 2.78788439273628983e-1, 1.26579413030177940e-1, 0.00000000000000000e+0),
	ERF_CONSTS(3.73997570145040850e+1, 1.05074004614827206e+0, 5.94651311286481502e-1, 3.64915280629351082e-4),
	ERF_CONSTS(1.12123870801026015e+2, 2.38574194785344389e+0, 1.61876655543871376e+0, 8.49717371168693357e-3),
	ERF_CONSTS(7.54843505665954743e+1, 3.37367334657284535e+0, 2.65383972869775752e+0, 8.69936222615385890e-2),
	ERF_CONSTS(1.61020914205869003e+1, 2.75143870676376208e+0, 2.45992070144245533e+0, 4.53767041780002545e-1),

	/* exp polynomial table */
	{SIMD_DCONST(1.26177193074810590878E-4)},
	{SIMD_DCONST(3.02994407707441961300E-2)},
	{SIMD_DCONST(9.99999999999999999910E-1)},
	{SIMD_DCONST(3.00198505138664455042E-6)},
	{SIMD_DCONST(2.52448340349684104192E-3)},
	{SIMD_DCONST(2.27265548208155028766E-1)},
	{SIMD_DCONST(2.00000000000000000009E0)},
#endif
    };

static const struct ICONSTS {
    vint64 isignmask;
    vint64 fixmasks[4];
    vint64 keepmask[4];
    vint i0x7f;
    vint i1023;
    } iconsts = {
        {SIMD_CONST(0x8000000000000000)},
	{{0,0,0,0},{0,-1,-1,-1},{0,0,-1,-1},{0,0,0,-1}},
	{{0,0,0,0},{-1,0,0,0},{-1,-1,0,0},{-1,-1,-1,0}},
        {SIMD_CONST(0x7f)},
        {SIMD_CONST(1023)},
    };

#ifdef USE_SVML
extern __m256d __svml_erf4(__m256d a);
extern __m256d __svml_exp4(__m256d a);
extern __m256d __svml_invsqrt4(__m256d a);
/*extern __m256 __svml_sincosf8(__m256,__m256 *);*/
#define verf __svml_erf4
#define vexp __svml_exp4
#define vrsqrt __svml_invsqrt4
/*#define vsincos(v,s,c) s = __svml_sincosf8(v,&c)*/
#ifdef __AVX__
#define vsincos(v,s,c) sincos256_ps(v,&s,&c)
#else
#define vsincos(v,s,c) sincos_ps(v,&s,&c);
#endif
#else
#ifdef __AVX__
#define vsincos(v,s,c) sincos256_ps(v,&s,&c)
#else
#define vsincos(v,s,c) sincos_ps(v,&s,&c);
#endif
#define vrsqrt SIMD_DRSQRT_EXACT

#ifdef __AVX__
#define SET_PREFACTOR(q) q = _mm256_blendv_pd(_mm256_permutevar_pd(consts.q##ab.p,pred0i),_mm256_permutevar_pd(consts.q##cd.p,pred2i),pred1)
#elif defined(slower__SSE4_1__)
#define SET_PREFACTOR(q) q = _mm_blendv_pd(_mm_blendv_pd(consts.q##a.p,consts.q##b.p,pred0),_mm_blendv_pd(consts.q##c.p,consts.q##d.p,pred2),pred1);
#else
#define SET_PREFACTOR(q) q = SIMD_DADD(SIMD_DADD(consts.q##a.p,SIMD_DAND(pred0,consts.q##b.p)),SIMD_DADD(SIMD_DAND(pred1,consts.q##c.p),SIMD_DAND(pred2,consts.q##d.p)))
#endif

/* Cost: AVX: 30 ops, SSE: 33 */
v_df vexp(v_df x) {
    v_df d,xx,pow2n,Pexp,Qexp;
    v_i4 n1,n2;
    v_i n;

    d = SIMD_DMUL(consts.log2e.p,x);
    n1 = MM_FCN(cvttpd,epi32)(d);
    d = MM_FCN(cvtepi32,pd)(n1);
    n1  = _mm_add_epi32(n1, iconsts.i1023.ph);
#ifdef __AVX__
    n2 = _mm_shuffle_epi32(n1,_MM_SHUFFLE(3, 3, 2, 2));
    n1 = _mm_shuffle_epi32(n1,_MM_SHUFFLE(1, 1, 0, 0));
    n2 = _mm_slli_epi64(n2, 52);
    n1 = _mm_slli_epi64(n1, 52);
    n = _mm256_castsi128_si256(n1);
    n = _mm256_insertf128_si256(n, n2, 0x1);
    pow2n = _mm256_castsi256_pd(n);
#else
    n1  = _mm_shuffle_epi32(n1, _MM_SHUFFLE(2, 1, 2, 0));
    n1  = _mm_slli_epi64(n1, 52);
    pow2n = MM_FCN(castsi128,pd)(n1);
#endif
    x = SIMD_DSUB(x,SIMD_DMUL(d,consts.C1.p));
    x = SIMD_DSUB(x,SIMD_DMUL(d,consts.C2.p));
    xx = SIMD_DMUL(x,x);

    Pexp = SIMD_DMUL(SIMD_DMADD(SIMD_DMADD(consts.p0.p,xx,consts.p1.p),xx,consts.p2.p),x);
    Qexp = SIMD_DSUB(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(consts.q0.p,xx,consts.q1.p),xx,consts.q2.p),xx,consts.q3.p),Pexp);
    x = SIMD_DDIV(Pexp,Qexp);
    x = SIMD_DMUL(SIMD_DMADD(x,consts.two.p,consts.one.p),pow2n);
    return x;
    }

/*
** This is an optimized version of erf that is accurate for all ranges.
** The exp() value is passed in because it it needed outside as well.
** Cost: AVX: 105, SSE: 172
*/
v_df verf(v_df v,v_df iv,v_df ex2,v_df *r_erf,v_df *r_erfc) {
    v_df p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,q5,Perf,Qerf,v2;
    v_df pred0,  pred1,  pred2;
    v_i  pred0i, pred2i;
    v_df t,nt;

    pred0 = SIMD_DCMP_GE(v,consts.threshold0.p);
    pred0i = SIMD_D2I(pred0);
    pred1 = SIMD_DCMP_GE(v,consts.threshold1.p);
    pred2 = SIMD_DCMP_GE(v,consts.threshold2.p);
    pred2i = SIMD_D2I(pred2);

    SET_PREFACTOR(p0);
    SET_PREFACTOR(p1);
    SET_PREFACTOR(p2);
    SET_PREFACTOR(p3);
    SET_PREFACTOR(p4);
    SET_PREFACTOR(p5);

    /* Gives x^2, x, x, 1/x^2 */
    v2 = SIMD_DSELECT(v,iv,pred2); /* x unless x>6.0, then 1/x */
    v2 = SIMD_DSELECT(SIMD_DMUL(v2,v2),v2,SIMD_DXOR(pred0,pred2));

    SET_PREFACTOR(q0);
    SET_PREFACTOR(q1);
    SET_PREFACTOR(q2);
    SET_PREFACTOR(q3);
    SET_PREFACTOR(q4);
    SET_PREFACTOR(q5);

    Perf = SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(p0,v2,p1),v2,p2),v2,p3),v2,p4),v2,p5);
    Qerf = SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(q0,v2,q1),v2,q2),v2,q3),v2,q4),v2,q5),v2,consts.one.p);

    t = SIMD_DMUL(SIMD_DMUL(SIMD_DSELECT(consts.one.p,iv,pred2),SIMD_DSELECT(v,ex2,pred0)),SIMD_DDIV(Perf,Qerf));
    nt = SIMD_DSUB(consts.one.p,t);
    *r_erf = SIMD_DSELECT(t,nt,pred0);
    *r_erfc = SIMD_DSELECT(nt,t,pred0);

    return t;
    }
#endif

int evalEwaldSIMD( ewaldSIMD *ew, v_df *ax, v_df *ay, v_df *az, v_df *dPot, v_df x, v_df y, v_df z, v_df r2, v_df doerfc ) {
    v_df dir,dir2,a,g0,g1,g2,g3,g4,g5,alphan;
    v_df xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
    v_df Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
    v_df Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;
    v_df rerf,rerfc,ex2,t,tx,ty,tz,tpot;
    v_df alpha2x2 = SIMD_DMUL(consts.two.p,ew->ewp.alpha2.p);

    dir = vrsqrt(r2);
    dir2 = SIMD_DMUL(dir,dir);
    ex2 = vexp(SIMD_DMUL(SIMD_DXOR(r2,iconsts.isignmask.pd),ew->ewp.alpha2.p));
    a = SIMD_DMUL(ex2,SIMD_DMUL(ew->ewp.ka.p,dir2));
    verf(SIMD_DMUL(ew->ewp.alpha.p,SIMD_DMUL(r2,dir)),SIMD_DMUL(ew->ewp.ialpha.p,dir),ex2,&rerf,&rerfc);
    g0 = SIMD_DMUL(dir,SIMD_DSELECT(SIMD_DXOR(rerf,iconsts.isignmask.pd),rerfc,doerfc));
    g1 = SIMD_DMADD(g0,dir2,a);
    alphan = alpha2x2;
    g2 = SIMD_DMADD(SIMD_DMUL(consts.three.p,dir2),g1,SIMD_DMUL(alphan,a));
    alphan = SIMD_DMUL(alphan,alpha2x2);
    g3 = SIMD_DMADD(SIMD_DMUL(consts.five.p,dir2),g2,SIMD_DMUL(alphan,a));
    alphan = SIMD_DMUL(alphan,alpha2x2);
    g4 = SIMD_DMADD(SIMD_DMUL(consts.seven.p,dir2),g3,SIMD_DMUL(alphan,a));
    alphan = SIMD_DMUL(alphan,alpha2x2);
    g5 = SIMD_DMADD(SIMD_DMUL(consts.nine.p,dir2),g4,SIMD_DMUL(alphan,a));
    xx = SIMD_DMUL(x,SIMD_DMUL(x,consts.half.p));
    xxx = SIMD_DMUL(xx,SIMD_DMUL(x,consts.onethird.p));
    xxy = SIMD_DMUL(xx,y);
    xxz = SIMD_DMUL(xx,z);
    yy = SIMD_DMUL(y,SIMD_DMUL(y,consts.half.p));
    yyy = SIMD_DMUL(yy,SIMD_DMUL(y,consts.onethird.p));
    xyy = SIMD_DMUL(yy,x);
    yyz = SIMD_DMUL(yy,z);
    zz = SIMD_DMUL(z,SIMD_DMUL(z,consts.half.p));
    zzz = SIMD_DMUL(zz,SIMD_DMUL(z,consts.onethird.p));
    xzz = SIMD_DMUL(zz,x);
    yzz = SIMD_DMUL(zz,y);
    xy = SIMD_DMUL(x,y);
    xyz = SIMD_DMUL(xy,z);
    xz = SIMD_DMUL(x,z);
    yz = SIMD_DMUL(y,z);
    Q2mirx = SIMD_DMUL(ew->ewm.xx.p,x);
    Q2mirx = SIMD_DMADD(ew->ewm.xy.p,y,Q2mirx);
    Q2mirx = SIMD_DMADD(ew->ewm.xz.p,z,Q2mirx);
    Q2miry = SIMD_DMUL(ew->ewm.xy.p,x);
    Q2miry = SIMD_DMADD(ew->ewm.yy.p,y,Q2miry);
    Q2miry = SIMD_DMADD(ew->ewm.yz.p,z,Q2miry);
    Q2mirz = SIMD_DMUL(ew->ewm.xz.p,x);
    Q2mirz = SIMD_DMADD(ew->ewm.yz.p,y,Q2mirz);
    Q2mirz = SIMD_DMADD(ew->ewm.zz.p,z,Q2mirz);
    Q3mirx = SIMD_DMUL(ew->ewm.xxx.p,xx);
    Q3mirx = SIMD_DMADD(ew->ewm.xxy.p,xy,Q3mirx);
    Q3mirx = SIMD_DMADD(ew->ewm.xxz.p,xz,Q3mirx);
    Q3mirx = SIMD_DMADD(ew->ewm.xyy.p,yy,Q3mirx);
    Q3mirx = SIMD_DMADD(ew->ewm.xyz.p,yz,Q3mirx);
    Q3mirx = SIMD_DMADD(ew->ewm.xzz.p,zz,Q3mirx);
    Q3miry = SIMD_DMUL(ew->ewm.xxy.p,xx);
    Q3miry = SIMD_DMADD(ew->ewm.xyy.p,xy,Q3miry);
    Q3miry = SIMD_DMADD(ew->ewm.xyz.p,xz,Q3miry);
    Q3miry = SIMD_DMADD(ew->ewm.yyy.p,yy,Q3miry);
    Q3miry = SIMD_DMADD(ew->ewm.yyz.p,yz,Q3miry);
    Q3miry = SIMD_DMADD(ew->ewm.yzz.p,zz,Q3miry);
    Q3mirz = SIMD_DMUL(ew->ewm.xxz.p,xx);
    Q3mirz = SIMD_DMADD(ew->ewm.xyz.p,xy,Q3mirz);
    Q3mirz = SIMD_DMADD(ew->ewm.xzz.p,xz,Q3mirz);
    Q3mirz = SIMD_DMADD(ew->ewm.yyz.p,yy,Q3mirz);
    Q3mirz = SIMD_DMADD(ew->ewm.yzz.p,yz,Q3mirz);
    Q3mirz = SIMD_DMADD(ew->ewm.zzz.p,zz,Q3mirz);
    Q4mirx = SIMD_DMUL(ew->ewm.xxxx.p,xxx);
    Q4mirx = SIMD_DMADD(ew->ewm.xxxy.p,xxy,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xxxz.p,xxz,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xxyy.p,xyy,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xxyz.p,xyz,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xxzz.p,xzz,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xyyy.p,yyy,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xyyz.p,yyz,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xyzz.p,yzz,Q4mirx);
    Q4mirx = SIMD_DMADD(ew->ewm.xzzz.p,zzz,Q4mirx);
    Q4miry = SIMD_DMUL(ew->ewm.xxxy.p,xxx);
    Q4miry = SIMD_DMADD(ew->ewm.xxyy.p,xxy,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.xxyz.p,xxz,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.xyyy.p,xyy,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.xyyz.p,xyz,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.xyzz.p,xzz,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.yyyy.p,yyy,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.yyyz.p,yyz,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.yyzz.p,yzz,Q4miry);
    Q4miry = SIMD_DMADD(ew->ewm.yzzz.p,zzz,Q4miry);
    Q4mirz = SIMD_DMUL(ew->ewm.xxxz.p,xxx);
    Q4mirz = SIMD_DMADD(ew->ewm.xxyz.p,xxy,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.xxzz.p,xxz,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.xyyz.p,xyy,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.xyzz.p,xyz,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.xzzz.p,xzz,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.yyyz.p,yyy,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.yyzz.p,yyz,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.yzzz.p,yzz,Q4mirz);
    Q4mirz = SIMD_DMADD(ew->ewm.zzzz.p,zzz,Q4mirz);
    Q4x = SIMD_DMUL(ew->ewp.Q4xx.p,x);
    Q4x = SIMD_DMADD(ew->ewp.Q4xy.p,y,Q4x);
    Q4x = SIMD_DMADD(ew->ewp.Q4xz.p,z,Q4x);
    Q4y = SIMD_DMUL(ew->ewp.Q4xy.p,x);
    Q4y = SIMD_DMADD(ew->ewp.Q4yy.p,y,Q4y);
    Q4y = SIMD_DMADD(ew->ewp.Q4yz.p,z,Q4y);
    Q4z = SIMD_DMUL(ew->ewp.Q4xz.p,x);
    Q4z = SIMD_DMADD(ew->ewp.Q4yz.p,y,Q4z);
    Q4z = SIMD_DMADD(ew->ewp.Q4zz.p,z,Q4z);
    Q2mir = SIMD_DMUL(Q2mirx,x);
    Q2mir = SIMD_DMADD(Q2miry,y,Q2mir);
    Q2mir = SIMD_DMADD(Q2mirz,z,Q2mir);
    Q2mir = SIMD_DMUL(Q2mir,consts.half.p);
    Q2mir = SIMD_DNMADD(ew->ewp.Q3x.p,x,Q2mir);
    Q2mir = SIMD_DNMADD(ew->ewp.Q3y.p,y,Q2mir);
    Q2mir = SIMD_DNMADD(ew->ewp.Q3z.p,z,Q2mir);
    Q2mir = SIMD_DADD(Q2mir,ew->ewp.Q4.p);
    Q3mir = SIMD_DMUL(Q3mirx,x);
    Q3mir = SIMD_DMADD(Q3miry,y,Q3mir);
    Q3mir = SIMD_DMADD(Q3mirz,z,Q3mir);
    Q3mir = SIMD_DMUL(Q3mir,consts.onethird.p);
    t =             SIMD_DMUL(Q4x,x);
    t = SIMD_DMADD(Q4y,y,t);
    t = SIMD_DMADD(Q4z,z,t);
    Q3mir = SIMD_DNMADD(t,consts.half.p,Q3mir);
    Q4mir = SIMD_DMUL(Q4mirx,x);
    Q4mir = SIMD_DMADD(Q4miry,y,Q4mir);
    Q4mir = SIMD_DMADD(Q4mirz,z,Q4mir);
    Q4mir = SIMD_DMUL(Q4mir,consts.onequarter.p);
    Qta = SIMD_DMUL(g1,ew->ewm.m.p);
    Qta = SIMD_DNMADD(g2,ew->ewp.Q2.p,Qta);
    Qta = SIMD_DMADD(g3,Q2mir,Qta);
    Qta = SIMD_DMADD(g4,Q3mir,Qta);
    Qta = SIMD_DMADD(g5,Q4mir,Qta);
    tpot = SIMD_DMUL(g0,ew->ewm.m.p);
    tpot = SIMD_DNMADD(g1,ew->ewp.Q2.p,tpot);
    tpot = SIMD_DMADD(g2,Q2mir,tpot);
    tpot = SIMD_DMADD(g3,Q3mir,tpot);
    tpot = SIMD_DMADD(g4,Q4mir,tpot);
    tx = SIMD_DMUL(g2,SIMD_DSUB(Q2mirx,ew->ewp.Q3x.p));
    tx = SIMD_DMADD(g3,SIMD_DSUB(Q3mirx,Q4x),tx);
    tx = SIMD_DMADD(g4,Q4mirx,tx);
    tx = SIMD_DNMADD(x,Qta,tx);
    ty = SIMD_DMUL(g2,SIMD_DSUB(Q2miry,ew->ewp.Q3y.p));
    ty = SIMD_DMADD(g3,SIMD_DSUB(Q3miry,Q4y),ty);
    ty = SIMD_DMADD(g4,Q4miry,ty);
    ty = SIMD_DNMADD(y,Qta,ty);
    tz = SIMD_DMUL(g2,SIMD_DSUB(Q2mirz,ew->ewp.Q3z.p));
    tz = SIMD_DMADD(g3,SIMD_DSUB(Q3mirz,Q4z),tz);
    tz = SIMD_DMADD(g4,Q4mirz,tz);
    tz = SIMD_DNMADD(z,Qta,tz);

    *dPot = tpot;
    *ax = tx;
    *ay = ty;
    *az = tz;
    return 447 * SIMD_WIDTH;
    }
#endif

int pkdParticleEwald(PKD pkd,PARTICLE *p, float *pa, float *pPot) {
    struct EwaldVariables *ew = &pkd->ew;
    EwaldTable *ewt = &pkd->ewt;
    const MOMC * restrict mom = &ew->mom;
    double L,Pot,ax,ay,az,dx,dy,dz,x,y,z,r2;
#ifdef USE_SIMD_EWALD
    v_df dPot,dax,day,daz;
    v_sf fPot,fax,fay,faz,fx,fy,fz;
    vdouble px, py, pz, pr2, doerfc;
    int nSIMD = 0;
#endif
    int i,ix,iy,iz;
    int bInHole,bInHolex,bInHolexy;
    int nFlop = 0;
    int nLoop = 0;

    assert(pkd->oAcceleration); /* Validate memory model */

    L = ew->Lbox;
    dx = pkdPos(p->r,0) - ew->r[0];
    dy = pkdPos(p->r,1) - ew->r[1];
    dz = pkdPos(p->r,2) - ew->r[2];

    ax = ay = az = 0.0;
    Pot = mom->m*ew->k1;
#ifdef USE_SIMD_EWALD
    dPot = SIMD_DSPLAT(0.0);
    dax = SIMD_DSPLAT(0.0);
    day = SIMD_DSPLAT(0.0);
    daz = SIMD_DSPLAT(0.0);
#endif
    for (ix=-ew->nEwReps;ix<=ew->nEwReps;++ix) {
	bInHolex = (abs(ix) <= ew->nReps);
	x = dx + ix*L;
	for (iy=-ew->nEwReps;iy<=ew->nEwReps;++iy) {
	    bInHolexy = (bInHolex && abs(iy) <= ew->nReps);
	    y = dy + iy*L;
	    for (iz=-ew->nEwReps;iz<=ew->nEwReps;++iz) {
		bInHole = (bInHolexy && abs(iz) <= ew->nReps);
		z = dz + iz*L;
		r2 = x*x + y*y + z*z;
		if (r2 > ew->fEwCut2 && !bInHole) continue;
		if (r2 < ew->fInner2) {
		    double g0,g1,g2,g3,g4,g5,alphan;
		    /*
		     * For small r, series expand about
		     * the origin to avoid errors caused
		     * by cancellation of large terms.
		     */
		    alphan = ew->ka;
		    r2 *= ew->alpha2;
		    g0 = alphan*((1.0/3.0)*r2 - 1.0);
		    alphan *= 2*ew->alpha2;
		    g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
		    alphan *= 2*ew->alpha2;
		    g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
		    alphan *= 2*ew->alpha2;
		    g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
		    alphan *= 2*ew->alpha2;
		    g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
		    alphan *= 2*ew->alpha2;
		    g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));

		    nFlop += evalEwald(ew,&ax,&ay,&az,&Pot,x,y,z,g0,g1,g2,g3,g4,g5);
		    }
		else {
#if defined(USE_SIMD_EWALD)
		    px.d[nSIMD] = x;
		    py.d[nSIMD] = y;
		    pz.d[nSIMD] = z;
		    pr2.d[nSIMD] = r2;
		    doerfc.i[nSIMD] = bInHole ? 0 : UINT64_MAX;
		    if (++nSIMD == SIMD_DWIDTH) {
			v_df tax, tay, taz, tpot;
			nFlop += evalEwaldSIMD(&pkd->es,&tax,&tay,&taz,&tpot,px.p,py.p,pz.p,pr2.p,doerfc.p);
 			dax = SIMD_DADD(dax,tax);
			day = SIMD_DADD(day,tay);
			daz = SIMD_DADD(daz,taz);
			dPot= SIMD_DSUB(dPot,tpot);
			nSIMD = 0;
			}
#else
		    double dir,dir2,a;
		    double g0,g1,g2,g3,g4,g5,alphan;
		    dir = 1/sqrt(r2);
		    dir2 = dir*dir;
		    a = exp(-r2*ew->alpha2);
		    a *= ew->ka*dir2;


		    if (bInHole) {
			g0 = -erf(ew->alpha*r2*dir);
			}
		    else {
			g0 = erfc(ew->alpha*r2*dir);
			}
		    g0 *= dir;
		    g1 = g0*dir2 + a;
		    alphan = 2*ew->alpha2;
		    g2 = 3*g1*dir2 + alphan*a;
		    alphan *= 2*ew->alpha2;
		    g3 = 5*g2*dir2 + alphan*a;
		    alphan *= 2*ew->alpha2;
		    g4 = 7*g3*dir2 + alphan*a;
		    alphan *= 2*ew->alpha2;
		    g5 = 9*g4*dir2 + alphan*a;
		    nFlop += evalEwald(ew,&ax,&ay,&az,&Pot,x,y,z,g0,g1,g2,g3,g4,g5);
#endif
		    }
		++nLoop;
		}
	    }
	}
#if defined(USE_SIMD_EWALD)
    /* Finish remaining SIMD operations if necessary */
    if (nSIMD) { /* nSIMD can be 0, 1, 2 or 3 */
	v_df t, tax, tay, taz, tpot;
	evalEwaldSIMD(&pkd->es,&tax,&tay,&taz,&tpot,px.p,py.p,pz.p,pr2.p,doerfc.p);
	t = iconsts.keepmask[nSIMD].pd;
	tax = SIMD_DAND(tax,t);
	tay = SIMD_DAND(tay,t);
	taz = SIMD_DAND(taz,t);
	tpot = SIMD_DAND(tpot,t);
	dax = SIMD_DADD(dax,tax);
	day = SIMD_DADD(day,tay);
	daz = SIMD_DADD(daz,taz);
	dPot= SIMD_DSUB(dPot,tpot);
	nSIMD = 0;
	}
#endif
#ifdef USE_SIMD_EWALD
    /* h-loop is done in float precision */
    fax = SIMD_D2F(dax);
    fay = SIMD_D2F(day);
    faz = SIMD_D2F(daz);
    fPot = SIMD_D2F(dPot);

    fx = SIMD_SPLAT(dx);
    fy = SIMD_SPLAT(dy);
    fz = SIMD_SPLAT(dz);

    nLoop = (ew->nEwhLoop+SIMD_MASK) >> SIMD_BITS;
    i = 0;
    do {
	v_sf hdotx,s,c,t;
	hdotx = SIMD_MADD(ewt->hz.p[i],fz,SIMD_MADD(ewt->hy.p[i],fy,SIMD_MUL(ewt->hx.p[i],fx)));
	vsincos(hdotx,s,c);
	fPot = SIMD_ADD(fPot,SIMD_ADD(SIMD_MUL(ewt->hSfac.p[i],s),SIMD_MUL(ewt->hCfac.p[i],c)));
	s = SIMD_MUL(ewt->hCfac.p[i],s);
	c = SIMD_MUL(ewt->hSfac.p[i],c);
	t = SIMD_SUB(s,c);
	fax = SIMD_ADD(fax,SIMD_MUL(ewt->hx.p[i],t));
	fay = SIMD_ADD(fay,SIMD_MUL(ewt->hy.p[i],t));
	faz = SIMD_ADD(faz,SIMD_MUL(ewt->hz.p[i],t));
	} while(++i < nLoop);
    nFlop += nLoop*58*SIMD_WIDTH;

    ax += SIMD_HADD(fax);
    ay += SIMD_HADD(fay);
    az += SIMD_HADD(faz);
    Pot += SIMD_HADD(fPot);
#else
    /*
    ** Scoring for the h-loop (+,*)
    ** 	Without trig = (10,14)
    **	    Trig est.	 = 2*(6,11)  same as 1/sqrt scoring.
    **		Total        = (22,36)
    **					 = 58
    */
    for (i=0;i<ew->nEwhLoop;++i) {
	double hdotx,s,c,t;
	hdotx = ewt->hx.f[i]*dx + ewt->hy.f[i]*dy + ewt->hz.f[i]*dz;
	c = cos(hdotx);
	s = sin(hdotx);
	Pot += ewt->hCfac.f[i]*c + ewt->hSfac.f[i]*s;
	t = ewt->hCfac.f[i]*s - ewt->hSfac.f[i]*c;
	ax += ewt->hx.f[i]*t;
	ay += ewt->hy.f[i]*t;
	az += ewt->hz.f[i]*t;
	}
    nFlop += ew->nEwhLoop*58;
#endif
    pa[0] += ax;
    pa[1] += ay;
    pa[2] += az;
    *pPot += Pot;
    return(nFlop);
    }

void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut) {
    struct EwaldVariables *ew = &pkd->ew;
    EwaldTable *ewt = &pkd->ewt;
    const MOMC * restrict mom = &ew->mom;
    int i,hReps,hx,hy,hz,h2;
    double k4,L;
    double gam[6],mfacc,mfacs;
    double ax,ay,az;
    const int iOrder = 4;

    L = pkd->fPeriod[0];
    ew->Lbox = L;
    /*
    ** Create SIMD versions of the moments.
    */
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    pkd->es.ewm.m.p = SIMD_DSPLAT(mom->m);
    pkd->es.ewm.xx.p = SIMD_DSPLAT(mom->xx);
    pkd->es.ewm.yy.p = SIMD_DSPLAT(mom->yy);
    pkd->es.ewm.xy.p = SIMD_DSPLAT(mom->xy);
    pkd->es.ewm.xz.p = SIMD_DSPLAT(mom->xz);
    pkd->es.ewm.yz.p = SIMD_DSPLAT(mom->yz);
    pkd->es.ewm.xxx.p = SIMD_DSPLAT(mom->xxx);
    pkd->es.ewm.xyy.p = SIMD_DSPLAT(mom->xyy);
    pkd->es.ewm.xxy.p = SIMD_DSPLAT(mom->xxy);
    pkd->es.ewm.yyy.p = SIMD_DSPLAT(mom->yyy);
    pkd->es.ewm.xxz.p = SIMD_DSPLAT(mom->xxz);
    pkd->es.ewm.yyz.p = SIMD_DSPLAT(mom->yyz);
    pkd->es.ewm.xyz.p = SIMD_DSPLAT(mom->xyz);
    pkd->es.ewm.xxxx.p = SIMD_DSPLAT(mom->xxxx);
    pkd->es.ewm.xyyy.p = SIMD_DSPLAT(mom->xyyy);
    pkd->es.ewm.xxxy.p = SIMD_DSPLAT(mom->xxxy);
    pkd->es.ewm.yyyy.p = SIMD_DSPLAT(mom->yyyy);
    pkd->es.ewm.xxxz.p = SIMD_DSPLAT(mom->xxxz);
    pkd->es.ewm.yyyz.p = SIMD_DSPLAT(mom->yyyz);
    pkd->es.ewm.xxyy.p = SIMD_DSPLAT(mom->xxyy);
    pkd->es.ewm.xxyz.p = SIMD_DSPLAT(mom->xxyz);
    pkd->es.ewm.xyyz.p = SIMD_DSPLAT(mom->xyyz);
    pkd->es.ewm.zz.p = SIMD_DSPLAT(mom->zz);
    pkd->es.ewm.xzz.p = SIMD_DSPLAT(mom->xzz);
    pkd->es.ewm.yzz.p = SIMD_DSPLAT(mom->yzz);
    pkd->es.ewm.zzz.p = SIMD_DSPLAT(mom->zzz);
    pkd->es.ewm.xxzz.p = SIMD_DSPLAT(mom->xxzz);
    pkd->es.ewm.xyzz.p = SIMD_DSPLAT(mom->xyzz);
    pkd->es.ewm.xzzz.p = SIMD_DSPLAT(mom->xzzz);
    pkd->es.ewm.yyzz.p = SIMD_DSPLAT(mom->yyzz);
    pkd->es.ewm.yzzz.p = SIMD_DSPLAT(mom->yzzz);
    pkd->es.ewm.zzzz.p = SIMD_DSPLAT(mom->zzzz);
#endif

    /*
    ** Set up traces of the complete multipole moments.
    */
    ew->Q4xx = 0.5*(mom->xxxx + mom->xxyy + mom->xxzz);
    ew->Q4xy = 0.5*(mom->xxxy + mom->xyyy + mom->xyzz);
    ew->Q4xz = 0.5*(mom->xxxz + mom->xyyz + mom->xzzz);
    ew->Q4yy = 0.5*(mom->xxyy + mom->yyyy + mom->yyzz);
    ew->Q4yz = 0.5*(mom->xxyz + mom->yyyz + mom->yzzz);
    ew->Q4zz = 0.5*(mom->xxzz + mom->yyzz + mom->zzzz);
    ew->Q4 = 0.25*(ew->Q4xx + ew->Q4yy + ew->Q4zz);
    ew->Q3x = 0.5*(mom->xxx + mom->xyy + mom->xzz);
    ew->Q3y = 0.5*(mom->xxy + mom->yyy + mom->yzz);
    ew->Q3z = 0.5*(mom->xxz + mom->yyz + mom->zzz);
    ew->Q2 = 0.5*(mom->xx + mom->yy + mom->zz);
    ew->nReps = nReps;
    ew->nEwReps = d2i(ceil(fEwCut));
    ew->nEwReps = ew->nEwReps > nReps ? ew->nEwReps : nReps;
    ew->fEwCut2 = fEwCut*fEwCut*L*L;
    ew->fInner2 = 1.2e-3*L*L;
    ew->alpha = 2.0/L;
    ew->ialpha = 0.5 * L;
    ew->alpha2 = ew->alpha*ew->alpha;
    ew->k1 = M_PI/(ew->alpha2*L*L*L);
    ew->ka = 2.0*ew->alpha/sqrt(M_PI);
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    pkd->es.ewp.Q4xx.p = SIMD_DSPLAT(ew->Q4xx);
    pkd->es.ewp.Q4xy.p = SIMD_DSPLAT(ew->Q4xy);
    pkd->es.ewp.Q4xz.p = SIMD_DSPLAT(ew->Q4xz);
    pkd->es.ewp.Q4yy.p = SIMD_DSPLAT(ew->Q4yy);
    pkd->es.ewp.Q4yz.p = SIMD_DSPLAT(ew->Q4yz);
    pkd->es.ewp.Q4zz.p = SIMD_DSPLAT(ew->Q4zz);
    pkd->es.ewp.Q4.p = SIMD_DSPLAT(ew->Q4);
    pkd->es.ewp.Q3x.p = SIMD_DSPLAT(ew->Q3x);
    pkd->es.ewp.Q3y.p = SIMD_DSPLAT(ew->Q3y);
    pkd->es.ewp.Q3z.p = SIMD_DSPLAT(ew->Q3z);
    pkd->es.ewp.Q2.p = SIMD_DSPLAT(ew->Q2);
    pkd->es.ewp.fEwCut2.p = SIMD_DSPLAT(ew->fEwCut2);
    pkd->es.ewp.fInner2.p = SIMD_DSPLAT(ew->fInner2);
    pkd->es.ewp.alpha.p = SIMD_DSPLAT(ew->alpha);
    pkd->es.ewp.ialpha.p = SIMD_DSPLAT(ew->ialpha);
    pkd->es.ewp.alpha2.p = SIMD_DSPLAT(ew->alpha2);
    pkd->es.ewp.k1.p = SIMD_DSPLAT(ew->k1);
    pkd->es.ewp.ka.p = SIMD_DSPLAT(ew->ka);
#endif


    /*
    ** Now setup stuff for the h-loop.
    */
    hReps = d2i(ceil(fhCut));
    k4 = M_PI*M_PI/(ew->alpha*ew->alpha*L*L);

    i = (int)pow(1+2*hReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE__)
    i = (i + SIMD_MASK) & ~SIMD_MASK;
#endif
    if ( i>ew->nMaxEwhLoop ) {
	ew->nMaxEwhLoop = i;
	ewt->hx.f = SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hx.f));
	assert(ewt->hx.f != NULL);
	ewt->hy.f = SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hy.f));
	assert(ewt->hy.f != NULL);
	ewt->hz.f = SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hz.f));
	assert(ewt->hz.f != NULL);
	ewt->hCfac.f = SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hCfac.f));
	assert(ewt->hCfac.f != NULL);
	ewt->hSfac.f = SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hSfac.f));
	assert(ewt->hSfac.f != NULL);
	}
    ew->nEwhLoop = i;
    i = (int)pow(1+2*ew->nEwReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    i = (i + SIMD_MASK) & ~SIMD_MASK;
#endif
    i = 0;
    for (hx=-hReps;hx<=hReps;++hx) {
	for (hy=-hReps;hy<=hReps;++hy) {
	    for (hz=-hReps;hz<=hReps;++hz) {
		h2 = hx*hx + hy*hy + hz*hz;
		if (h2 == 0) continue;
		if (h2 > fhCut*fhCut) continue;
		assert (i < ew->nMaxEwhLoop);
		gam[0] = exp(-k4*h2)/(M_PI*h2*L);
		gam[1] = 2*M_PI/L*gam[0];
		gam[2] = -2*M_PI/L*gam[1];
		gam[3] = 2*M_PI/L*gam[2];
		gam[4] = -2*M_PI/L*gam[3];
		gam[5] = 2*M_PI/L*gam[4];
		gam[1] = 0.0;
		gam[3] = 0.0;
		gam[5] = 0.0;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		mfacc = 0.0;
		QEVAL(iOrder,ew->mom,gam,hx,hy,hz,ax,ay,az,mfacc);
		gam[0] = exp(-k4*h2)/(M_PI*h2*L);
		gam[1] = 2*M_PI/L*gam[0];
		gam[2] = -2*M_PI/L*gam[1];
		gam[3] = 2*M_PI/L*gam[2];
		gam[4] = -2*M_PI/L*gam[3];
		gam[5] = 2*M_PI/L*gam[4];
		gam[0] = 0.0;
		gam[2] = 0.0;
		gam[4] = 0.0;
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
		mfacs = 0.0;
		QEVAL(iOrder,ew->mom,gam,hx,hy,hz,ax,ay,az,mfacs);
		ewt->hx.f[i] = 2*M_PI/L*hx;
		ewt->hy.f[i] = 2*M_PI/L*hy;
		ewt->hz.f[i] = 2*M_PI/L*hz;
		ewt->hCfac.f[i] = mfacc;
		ewt->hSfac.f[i] = mfacs;
		++i;
		}
	    }
	}
    ew->nEwhLoop = i;
    while(i<ew->nMaxEwhLoop) {
	ewt->hx.f[i] = 0;
	ewt->hy.f[i] = 0;
	ewt->hz.f[i] = 0;
	ewt->hCfac.f[i] = 0;
	ewt->hSfac.f[i] = 0;
	++i;
	}
#ifdef USE_CUDA
    if (mdlCore(pkd->mdl)==0) {
	cudaEwaldInit(ew,ewt);
	}
    mdlThreadBarrier(pkd->mdl);
#endif


    }
