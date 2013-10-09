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
    double x, double y, double z, const MOMC * restrict mom,
    double g0, double g1, double g2, double g3,double g4, double g5) {
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

static int evalMainBox(struct EwaldVariables *ew,double *ax, double *ay, double *az, double *fPot,
    double x, double y, double z, const MOMC * restrict mom) {
    double r2,dir,dir2,a,alphan,L;
    double g0,g1,g2,g3,g4,g5;

    r2 = x*x + y*y + z*z;
    if (r2 < ew->fInner2) {
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
	}
    else {
	dir = 1/sqrt(r2);
	dir2 = dir*dir;
	a = exp(-r2*ew->alpha2);
	a *= ew->ka*dir2;
	/*
	** We can get away with not using erfc because it is vanishing
	** if the following are true:
	** 1. The center of mass is near/at the centre of the box,
	** 2. We use 1 replica for RMS error < ~1e-4 (erfc 3 = 2.2e-5)
	** 3. We use 2 replicas for smaller RMS  (erfc 5 = 1.5e-12)
	** Worst case with COM near the edge:
	** 1 replica:  erfc 2 = 4.6e-3 (this would be wrong)
	** 2 replicas: erfc 4 = 1.5e-8
	** We currently switch from 1 replica to 2 replicates at at
	** theta of 0.5 which corresponds to an RMS error of ~1e-4.
	*/
	g0 = -erf(ew->alpha*r2*dir);
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
	}
    return 58 + evalEwald(ew,ax,ay,az,fPot,x,y,z,mom,g0,g1,g2,g3,g4,g5);
    }

#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
static const struct CONSTS {
    vdouble onequarter,onethird,half,one,two,three,five,seven,nine;
#ifndef USE_SVML
    vdouble log2e,C1,C2;
    vdouble threshold0,threshold1,threshold2;
#ifdef __AVX__
#ifdef ERF_FULL_RANGE
    vdouble p0bc,p0a,p1bc,p1a,p2bc,p2a,p3bc,p3a,p4bc,p4a,p5bc,p5a;
    vdouble q0bc,q0a,q1bc,q1a,q2bc,q2a,q3bc,q3a,q4bc,q4a,q5bc,q5a;
#else
    vdouble p0bc,p1bc,p2bc,p3bc,p4bc,p5bc;
    vdouble q0bc,q1bc,q2bc,q3bc,q4bc,q5bc;
#endif
#else
    vdouble p0b,p0c,p1b,p1c,p2b,p2c,p3b,p3c,p4b,p4c,p5b,p5c;
    vdouble q0b,q0c,q1b,q1c,q2b,q2c,q3b,q3c,q4b,q4c,q5b,q5c;
#endif
    vdouble p0,p1,p2,q0,q1,q2,q3;
#endif
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
#ifdef ERF_FULL_RANGE
#define ERF_CONSTS(c,b,a) {{c,b,c,b}},{a,a,a,a}}
#else
#define ERF_CONSTS(c,b,a) {{c,b,c,b}}
#endif
#else
#ifdef ERF_FULL_RANGE
#define ERF_CONSTS(c,b,a) {{c,c}},{{b-c,b-c}},{{b-c-a,b-c-a}}
#else
#define ERF_CONSTS(c,b,a) {{c,c}},{{b-c,b-c}}
#endif
#endif
	ERF_CONSTS(2.25716982919217555e-2,7.06940843763253131e-3,0),
	ERF_CONSTS(1.57289620742838702e-1,7.14193832506776067e-2,6.49254556481904354e-5),
	ERF_CONSTS(5.81528574177741135e-1,3.31899559578213215e-1,1.20339380863079457e-3),
	ERF_CONSTS(1.26739901455873222e+0,8.78115804155881782e-1,4.03259488531795274e-2),
	ERF_CONSTS(1.62356584489366647e+0,1.33154163936765307,1.35894887627277916e-1),
	ERF_CONSTS(9.99921140009714409e-1,9.99999992049799098e-1,1.12837916709551256e+0),

	ERF_CONSTS(4.00072964526861362e-2,1.25304936549413393e-2,0),
	ERF_CONSTS(2.78788439273628983e-1,1.26579413030177940e-1,0),
	ERF_CONSTS(1.05074004614827206,5.94651311286481502e-1,3.64915280629351082e-4),
	ERF_CONSTS(2.38574194785344389,1.61876655543871376,8.49717371168693357e-3),
	ERF_CONSTS(3.37367334657284535,2.65383972869775752,8.69936222615385890e-2),
	ERF_CONSTS(2.75143870676376208,2.45992070144245533,4.53767041780002545e-1),

	SIMD_DCONST(1.26177193074810590878E-4),
	SIMD_DCONST(3.02994407707441961300E-2),
	SIMD_DCONST(9.99999999999999999910E-1),
	SIMD_DCONST(3.00198505138664455042E-6),
	SIMD_DCONST(2.52448340349684104192E-3),
	SIMD_DCONST(2.27265548208155028766E-1),
	SIMD_DCONST(2.00000000000000000009E0),
#endif
    };

static const struct ICONSTS {
    vint64 isignmask;
    vint64 fixmasks[4];
    vint i0x7f;
    } iconsts = {
        {SIMD_CONST(0x8000000000000000)},
	{{0,0,0,0},{0,-1,-1,-1},{0,0,-1,-1},{0,0,0,-1}},
        {SIMD_CONST(0x7f)},
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
#ifdef ERF_FULL_RANGE
#ifdef __AVX__
#define SET_PREFACTOR(q) q = _mm256_blendv_pd(_mm256_permutevar_pd(consts.q##bc.p,pred1i),consts.q##a.p,pred0)
#else
#endif
#else
#ifdef __AVX__
#define SET_PREFACTOR(q) q = _mm256_permutevar_pd(consts.q##bc.p,pred1i)
#else
#define SET_PREFACTOR(q) q = SIMD_DADD(consts.q##b.p,SIMD_DANDNOT(pred1,consts.q##c.p))
#endif
#endif


v_df vexp(v_df x) {
    v_df d,xx,pow2n,Pexp,Qexp;
    __m128i n;

    d = SIMD_DFLOOR(SIMD_DMADD(consts.log2e.p,x,consts.half.p));
    n = MM_FCN(cvttpd,epi32)(d);
    n = _mm_add_epi32(n, iconsts.i0x7f.ph );
    n = _mm_slli_epi32(n, 23);
    pow2n = MM_FCN(cvtps,pd)(_mm_castsi128_ps(n));
    x = SIMD_DSUB(x,SIMD_DMUL(d,consts.C1.p));
    x = SIMD_DSUB(x,SIMD_DMUL(d,consts.C2.p));
    xx = SIMD_DMUL(x,x);

    Pexp = SIMD_DMUL(SIMD_DMADD(SIMD_DMADD(consts.p0.p,xx,consts.p1.p),xx,consts.p2.p),x);
    Qexp = SIMD_DSUB(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(consts.q0.p,xx,consts.q1.p),xx,consts.q2.p),xx,consts.q3.p),Pexp);
//    x = SIMD_DMUL(Pexp,SIMD_DRE_GOOD(Qexp));
    x = SIMD_DDIV(Pexp,Qexp);
    x = SIMD_DMUL(SIMD_DMADD(x,consts.two.p,consts.one.p),pow2n);
    return x;
    }


/*
** This is an optimized version of erf that is accurate (without ERF_FULL_RANGE)
** in the range [0.65,9.0]. This is ideal for ewald replicas.
*/
v_df verf(v_df v) {
    v_df p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,q5,Perf,Qerf,Pexp,Qexp,v2;
    v_i pred1i;
    v_i4 n;
    v_df pred1,t,t2,pow2n;
#ifdef ERF_FULL_RANGE
    v_df pred0;
#endif

#ifdef ERF_FULL_RANGE
    pred0 = SIMD_DCMP_LT(v,consts.threshold0.p);
#endif
    pred1 = SIMD_DCMP_LT(v,consts.threshold1.p);
    pred1i = SIMD_D2I(pred1);

    SET_PREFACTOR(p0);
    SET_PREFACTOR(p1);
    SET_PREFACTOR(p2);
    SET_PREFACTOR(p3);
    SET_PREFACTOR(p4);
    SET_PREFACTOR(p5);
#ifdef ERF_FULL_RANGE
    v = SIMD_DMIN(v,consts.threshold2.p);
    v2 = SIMD_DMUL(v,SIMD_DOR(SIMD_DAND(pred0,v),SIMD_DANDNOT(pred0,consts.one.p)));
#else
    v2 = v;
#endif
    SET_PREFACTOR(q0);
    SET_PREFACTOR(q1);
    SET_PREFACTOR(q2);
    SET_PREFACTOR(q3);
    SET_PREFACTOR(q4);
    SET_PREFACTOR(q5);

    Perf = SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(p0,v2,p1),v2,p2),v2,p3),v2,p4),v2,p5);
    Qerf = SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(q0,v2,q1),v2,q2),v2,q3),v2,q4),v2,q5),v2,consts.one.p);

    t = SIMD_DMUL(SIMD_DXOR(v,iconsts.isignmask.pd),v);
    Pexp = SIMD_DFLOOR(SIMD_DMADD(consts.log2e.p,t,consts.half.p));
    t = SIMD_DSUB(t,SIMD_DMUL(Pexp,consts.C1.p));
    t = SIMD_DSUB(t,SIMD_DMUL(Pexp,consts.C2.p));
    n = MM_FCN(cvttpd,epi32)(Pexp);
    n = _mm_add_epi32(n, iconsts.i0x7f.ph );
    n = _mm_slli_epi32(n, 23);
    pow2n = MM_FCN(cvtps,pd)(_mm_castsi128_ps(n));
    t2 = SIMD_DMUL(t,t);

    Pexp = SIMD_DMUL(SIMD_DMADD(SIMD_DMADD(consts.p0.p,t2,consts.p1.p),t2,consts.p2.p),t);
    Qexp = SIMD_DSUB(SIMD_DMADD(SIMD_DMADD(SIMD_DMADD(consts.q0.p,t2,consts.q1.p),t2,consts.q2.p),t2,consts.q3.p),Pexp);

//    t = SIMD_DMUL(Perf,SIMD_DRE_EXACT(SIMD_DMUL(Qerf,Qexp)));
    t = SIMD_DDIV(Perf,SIMD_DMUL(Qerf,Qexp));
    v2 = SIMD_DNMSUB(SIMD_DMADD(consts.two.p,Pexp,Qexp),SIMD_DMUL(pow2n,t),consts.one.p);

#ifdef ERF_FULL_RANGE
    t2 = SIMD_DMUL(v,SIMD_DMUL(t,Qexp));
    t = SIMD_DOR(SIMD_DAND(pred0,t2),SIMD_DANDNOT(pred0,v2));
    return t;
#else
    return v2;
#endif
    }
#endif

int pkdParticleEwaldSIMD(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
    PARTICLE *p, float *pa, float *pPot) {
    struct EwaldVariables *ew = &pkd->ew;
    v_df ax,ay,az,L,dPot,dx,dy,dz,alpha2x2;
    v_df tx,ty,tz,tpot,t;
    v_sf fPot,fax,fay,faz,fx,fy,fz;
    double px,py,pz,a1,a2,a3,p1;
    int i, nLoop;
    const v_df *Lx = ew->ewt.Lx.p;
    const v_df *Ly = ew->ewt.Ly.p;
    const v_df *Lz = ew->ewt.Lz.p;

    assert(pkd->oAcceleration); /* Validate memory model */
    assert(pkd->oPotential); /* Validate memory model */

    if (!pkdIsDstActive(p,uRungLo,uRungHi)) return 0;

    alpha2x2 = SIMD_DMUL(consts.two.p,ew->ewp.alpha2.p);
    L = SIMD_DSPLAT(pkd->fPeriod[0]);
    dPot = SIMD_DSPLAT(ew->ewm.m.d[0]*pkd->ew.k1);
    ax = SIMD_DSPLAT(0.0);
    ay = SIMD_DSPLAT(0.0);
    az = SIMD_DSPLAT(0.0);
    px = p->r[0] - pkdTopNode(pkd,ROOT)->r[0];
    py = p->r[1] - pkdTopNode(pkd,ROOT)->r[1];
    pz = p->r[2] - pkdTopNode(pkd,ROOT)->r[2];
    dx = SIMD_DSPLAT(px);
    dy = SIMD_DSPLAT(py);
    dz = SIMD_DSPLAT(pz);
    fx = SIMD_SPLAT(px);
    fy = SIMD_SPLAT(py);
    fz = SIMD_SPLAT(pz);

    nLoop = (pkd->ew.nEwLoopInner+SIMD_DMASK) >> SIMD_DBITS;
    i = 0;
    do {
	v_df x,y,z,r2,dir,dir2,a,g0,g1,g2,g3,g4,g5,alphan;
	v_df xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
	v_df Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
	v_df Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;

	x = SIMD_DADD(dx,Lx[i]);
	y = SIMD_DADD(dy,Ly[i]);
	z = SIMD_DADD(dz,Lz[i]);
	r2 = SIMD_DMADD(z,z,SIMD_DMADD(y,y,SIMD_DMUL(x,x)));
	dir = vrsqrt(r2);
	dir2 = SIMD_DMUL(dir,dir);
	a = SIMD_DMUL(vexp(SIMD_DMUL(SIMD_DXOR(r2,iconsts.isignmask.pd),ew->ewp.alpha2.p)),SIMD_DMUL(ew->ewp.ka.p,dir2));
	g0 = SIMD_DMUL(dir,SIMD_DXOR(verf(SIMD_DMUL(ew->ewp.alpha.p,SIMD_DMUL(r2,dir))),iconsts.isignmask.pd));
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
	Q2mirx = SIMD_DADD(Q2mirx,SIMD_DMUL(ew->ewm.xy.p,y));
	Q2mirx = SIMD_DADD(Q2mirx,SIMD_DMUL(ew->ewm.xz.p,z));
	Q2miry = SIMD_DMUL(ew->ewm.xy.p,x);
	Q2miry = SIMD_DADD(Q2miry,SIMD_DMUL(ew->ewm.yy.p,y));
	Q2miry = SIMD_DADD(Q2miry,SIMD_DMUL(ew->ewm.yz.p,z));
	Q2mirz = SIMD_DMUL(ew->ewm.xz.p,x);
	Q2mirz = SIMD_DADD(Q2mirz,SIMD_DMUL(ew->ewm.yz.p,y));
	Q2mirz = SIMD_DADD(Q2mirz,SIMD_DMUL(ew->ewm.zz.p,z));
	Q3mirx = SIMD_DMUL(ew->ewm.xxx.p,xx);
	Q3mirx = SIMD_DADD(Q3mirx,SIMD_DMUL(ew->ewm.xxy.p,xy));
	Q3mirx = SIMD_DADD(Q3mirx,SIMD_DMUL(ew->ewm.xxz.p,xz));
	Q3mirx = SIMD_DADD(Q3mirx,SIMD_DMUL(ew->ewm.xyy.p,yy));
	Q3mirx = SIMD_DADD(Q3mirx,SIMD_DMUL(ew->ewm.xyz.p,yz));
	Q3mirx = SIMD_DADD(Q3mirx,SIMD_DMUL(ew->ewm.xzz.p,zz));
	Q3miry = SIMD_DMUL(ew->ewm.xxy.p,xx);
	Q3miry = SIMD_DADD(Q3miry,SIMD_DMUL(ew->ewm.xyy.p,xy));
	Q3miry = SIMD_DADD(Q3miry,SIMD_DMUL(ew->ewm.xyz.p,xz));
	Q3miry = SIMD_DADD(Q3miry,SIMD_DMUL(ew->ewm.yyy.p,yy));
	Q3miry = SIMD_DADD(Q3miry,SIMD_DMUL(ew->ewm.yyz.p,yz));
	Q3miry = SIMD_DADD(Q3miry,SIMD_DMUL(ew->ewm.yzz.p,zz));
	Q3mirz = SIMD_DMUL(ew->ewm.xxz.p,xx);
	Q3mirz = SIMD_DADD(Q3mirz,SIMD_DMUL(ew->ewm.xyz.p,xy));
	Q3mirz = SIMD_DADD(Q3mirz,SIMD_DMUL(ew->ewm.xzz.p,xz));
	Q3mirz = SIMD_DADD(Q3mirz,SIMD_DMUL(ew->ewm.yyz.p,yy));
	Q3mirz = SIMD_DADD(Q3mirz,SIMD_DMUL(ew->ewm.yzz.p,yz));
	Q3mirz = SIMD_DADD(Q3mirz,SIMD_DMUL(ew->ewm.zzz.p,zz));
	Q4mirx = SIMD_DMUL(ew->ewm.xxxx.p,xxx);
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xxxy.p,xxy));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xxxz.p,xxz));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xxyy.p,xyy));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xxyz.p,xyz));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xxzz.p,xzz));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xyyy.p,yyy));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xyyz.p,yyz));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xyzz.p,yzz));
	Q4mirx = SIMD_DADD(Q4mirx,SIMD_DMUL(ew->ewm.xzzz.p,zzz));
	Q4miry = SIMD_DMUL(ew->ewm.xxxy.p,xxx);
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.xxyy.p,xxy));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.xxyz.p,xxz));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.xyyy.p,xyy));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.xyyz.p,xyz));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.xyzz.p,xzz));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.yyyy.p,yyy));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.yyyz.p,yyz));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.yyzz.p,yzz));
	Q4miry = SIMD_DADD(Q4miry,SIMD_DMUL(ew->ewm.yzzz.p,zzz));
	Q4mirz = SIMD_DMUL(ew->ewm.xxxz.p,xxx);
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.xxyz.p,xxy));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.xxzz.p,xxz));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.xyyz.p,xyy));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.xyzz.p,xyz));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.xzzz.p,xzz));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.yyyz.p,yyy));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.yyzz.p,yyz));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.yzzz.p,yzz));
	Q4mirz = SIMD_DADD(Q4mirz,SIMD_DMUL(ew->ewm.zzzz.p,zzz));
	Q4x = SIMD_DMUL(ew->ewp.Q4xx.p,x);
	Q4x = SIMD_DADD(Q4x,SIMD_DMUL(ew->ewp.Q4xy.p,y));
	Q4x = SIMD_DADD(Q4x,SIMD_DMUL(ew->ewp.Q4xz.p,z));
	Q4y = SIMD_DMUL(ew->ewp.Q4xy.p,x);
	Q4y = SIMD_DADD(Q4y,SIMD_DMUL(ew->ewp.Q4yy.p,y));
	Q4y = SIMD_DADD(Q4y,SIMD_DMUL(ew->ewp.Q4yz.p,z));
	Q4z = SIMD_DMUL(ew->ewp.Q4xz.p,x);
	Q4z = SIMD_DADD(Q4z,SIMD_DMUL(ew->ewp.Q4yz.p,y));
	Q4z = SIMD_DADD(Q4z,SIMD_DMUL(ew->ewp.Q4zz.p,z));
	Q2mir =                  SIMD_DMUL(Q2mirx,x);
	Q2mir = SIMD_DADD(Q2mir,SIMD_DMUL(Q2miry,y));
	Q2mir = SIMD_DADD(Q2mir,SIMD_DMUL(Q2mirz,z));
	Q2mir = SIMD_DMUL(Q2mir,consts.half.p);
	Q2mir = SIMD_DSUB(Q2mir,SIMD_DMUL(ew->ewp.Q3x.p,x));
	Q2mir = SIMD_DSUB(Q2mir,SIMD_DMUL(ew->ewp.Q3y.p,y));
	Q2mir = SIMD_DSUB(Q2mir,SIMD_DMUL(ew->ewp.Q3z.p,z));
	Q2mir = SIMD_DADD(Q2mir,ew->ewp.Q4.p);
	Q3mir =                 SIMD_DMUL(Q3mirx,x);
	Q3mir = SIMD_DADD(Q3mir,SIMD_DMUL(Q3miry,y));
	Q3mir = SIMD_DADD(Q3mir,SIMD_DMUL(Q3mirz,z));
	Q3mir = SIMD_DMUL(Q3mir,consts.onethird.p);
	t =             SIMD_DMUL(Q4x,x);
	t = SIMD_DADD(t,SIMD_DMUL(Q4y,y));
	t = SIMD_DADD(t,SIMD_DMUL(Q4z,z));
	t = SIMD_DMUL(t,consts.half.p);
	Q3mir = SIMD_DSUB(Q3mir,t);
	Q4mir =                 SIMD_DMUL(Q4mirx,x);
	Q4mir = SIMD_DADD(Q4mir,SIMD_DMUL(Q4miry,y));
	Q4mir = SIMD_DADD(Q4mir,SIMD_DMUL(Q4mirz,z));
	Q4mir = SIMD_DMUL(Q4mir,consts.onequarter.p);
	Qta = SIMD_DMUL(g1,ew->ewm.m.p);
	Qta = SIMD_DSUB(Qta,SIMD_DMUL(g2,ew->ewp.Q2.p));
	Qta = SIMD_DADD(Qta,SIMD_DMUL(g3,Q2mir));
	Qta = SIMD_DADD(Qta,SIMD_DMUL(g4,Q3mir));
	Qta = SIMD_DADD(Qta,SIMD_DMUL(g5,Q4mir));
	tpot = SIMD_DMUL(g0,ew->ewm.m.p);
	tpot = SIMD_DSUB(tpot,SIMD_DMUL(g1,ew->ewp.Q2.p));
	tpot = SIMD_DADD(tpot,SIMD_DMUL(g2,Q2mir));
	tpot = SIMD_DADD(tpot,SIMD_DMUL(g3,Q3mir));
	tpot = SIMD_DADD(tpot,SIMD_DMUL(g4,Q4mir));
	tx = SIMD_DMUL(g2,SIMD_DSUB(Q2mirx,ew->ewp.Q3x.p));
	tx = SIMD_DADD(tx,SIMD_DMUL(g3,SIMD_DSUB(Q3mirx,Q4x)));
	tx = SIMD_DADD(tx,SIMD_DMUL(g4,Q4mirx));
	tx = SIMD_DSUB(tx,SIMD_DMUL(x,Qta));
	ty = SIMD_DMUL(g2,SIMD_DSUB(Q2miry,ew->ewp.Q3y.p));
	ty = SIMD_DADD(ty,SIMD_DMUL(g3,SIMD_DSUB(Q3miry,Q4y)));
	ty = SIMD_DADD(ty,SIMD_DMUL(g4,Q4miry));
	ty = SIMD_DSUB(ty,SIMD_DMUL(y,Qta));
	tz = SIMD_DMUL(g2,SIMD_DSUB(Q2mirz,ew->ewp.Q3z.p));
	tz = SIMD_DADD(tz,SIMD_DMUL(g3,SIMD_DSUB(Q3mirz,Q4z)));
	tz = SIMD_DADD(tz,SIMD_DMUL(g4,Q4mirz));
	tz = SIMD_DSUB(tz,SIMD_DMUL(z,Qta));
	dPot = SIMD_DSUB(dPot,tpot);
	ax = SIMD_DADD(ax,tx);
	ay = SIMD_DADD(ay,ty);
	az = SIMD_DADD(az,tz);
	} while(++i < nLoop);

    /* Correct the unused terms */
    nLoop = pkd->ew.nEwLoopInner;
    t = iconsts.fixmasks[nLoop&SIMD_DMASK].pd;
    tx = SIMD_DAND(tx,t);
    ty = SIMD_DAND(ty,t);
    tz = SIMD_DAND(tz,t);
    tpot = SIMD_DAND(tpot,t);
    ax = SIMD_DSUB(ax,tx);
    ay = SIMD_DSUB(ay,ty);
    az = SIMD_DSUB(az,tz);
    dPot = SIMD_DADD(dPot,tpot);
    fax = SIMD_D2F(ax);
    fay = SIMD_D2F(ay);
    faz = SIMD_D2F(az);
    fPot = SIMD_D2F(dPot);

    /* Add the contribution from the main box */
    a1 = a2 = a3 = p1 = 0;
    evalMainBox(&pkd->ew,&a1,&a2,&a3,&p1,px,py,pz,&pkd->momRoot);
    fax = SIMD_ADD(fax,SIMD_LOADS(a1));
    fay = SIMD_ADD(fay,SIMD_LOADS(a2));
    faz = SIMD_ADD(faz,SIMD_LOADS(a3));
    fPot= SIMD_ADD(fPot,SIMD_LOADS(p1));
    nLoop = (pkd->ew.nEwhLoop+SIMD_MASK) >> SIMD_BITS;
    i = 0;
    do {
	v_sf hdotx,s,c,t;
	hdotx = SIMD_MADD(pkd->ew.ewt.hz.p[i],fz,SIMD_MADD(pkd->ew.ewt.hy.p[i],fy,SIMD_MUL(pkd->ew.ewt.hx.p[i],fx)));
	vsincos(hdotx,s,c);
	fPot = SIMD_ADD(fPot,SIMD_ADD(SIMD_MUL(pkd->ew.ewt.hSfac.p[i],s),SIMD_MUL(pkd->ew.ewt.hCfac.p[i],c)));
	s = SIMD_MUL(pkd->ew.ewt.hCfac.p[i],s);
	c = SIMD_MUL(pkd->ew.ewt.hSfac.p[i],c);
	t = SIMD_SUB(s,c);
	fax = SIMD_ADD(fax,SIMD_MUL(pkd->ew.ewt.hx.p[i],t));
	fay = SIMD_ADD(fay,SIMD_MUL(pkd->ew.ewt.hy.p[i],t));
	faz = SIMD_ADD(faz,SIMD_MUL(pkd->ew.ewt.hz.p[i],t));
	} while(++i < nLoop);

    pa[0] += SIMD_HADD(fax);
    pa[1] += SIMD_HADD(fay);
    pa[2] += SIMD_HADD(faz);
    *pPot += SIMD_HADD(fPot);
    return pkd->ew.nEwLoopInner*447 + pkd->ew.nEwhLoop*58;
    }
#endif

int pkdParticleEwald(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
    PARTICLE *p, float *pa, float *pPot) {
    MOMC mom = pkd->momRoot;
    double fPot,ax,ay,az;
    double dx,dy,dz,x,y,z,r2,dir,dir2,a,alphan,L;
    double xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
    double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
    double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;
    double g0,g1,g2,g3,g4,g5;
    double onethird = 1.0/3.0;
    int i,ix,iy,iz;
#ifdef USE_FORMALLY_CORRECT_EWALD
    int bInHole,bInHolex,bInHolexy;
#endif
    int nFlop;
    int nLoop = 0;
    float hdotx,s,c,t;

    assert(pkd->oAcceleration); /* Validate memory model */
    assert(pkd->oPotential); /* Validate memory model */

    if (!pkdIsDstActive(p,uRungLo,uRungHi)) return 0;

    L = pkd->fPeriod[0];
    fPot = mom.m*pkd->ew.k1;
    ax = 0.0;
    ay = 0.0;
    az = 0.0;
    dx = p->r[0] - pkdTopNode(pkd,ROOT)->r[0];
    dy = p->r[1] - pkdTopNode(pkd,ROOT)->r[1];
    dz = p->r[2] - pkdTopNode(pkd,ROOT)->r[2];

#ifndef USE_FORMALLY_CORRECT_EWALD
    nLoop = pkd->ew.nEwLoopInner;
    r2 = dx*dx + dy*dy + dz*dz;
    /* We will handle this at the end if necessary */
    if (r2 >= pkd->ew.fInner2) ++nLoop;
    for( i=0; i<nLoop; ++i) {
	x = dx + pkd->ew.ewt.Lx.d[i];
	y = dy + pkd->ew.ewt.Ly.d[i];
	z = dz + pkd->ew.ewt.Lz.d[i];
	r2 = x*x + y*y + z*z;
	dir = 1/sqrt(r2);
	dir2 = dir*dir;
	a = exp(-r2*pkd->ew.alpha2);
	a *= pkd->ew.ka*dir2;
	g0 = -erf(pkd->ew.alpha*r2*dir);
	g0 *= dir;
	g1 = g0*dir2 + a;
	alphan = 2*pkd->ew.alpha2;
	g2 = 3*g1*dir2 + alphan*a;
	alphan *= 2*pkd->ew.alpha2;
	g3 = 5*g2*dir2 + alphan*a;
	alphan *= 2*pkd->ew.alpha2;
	g4 = 7*g3*dir2 + alphan*a;
	alphan *= 2*pkd->ew.alpha2;
	g5 = 9*g4*dir2 + alphan*a;

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
	Q2mirx = mom.xx*x + mom.xy*y + mom.xz*z;
	Q2miry = mom.xy*x + mom.yy*y + mom.yz*z;
	Q2mirz = mom.xz*x + mom.yz*y + mom.zz*z;
	Q3mirx = mom.xxx*xx + mom.xxy*xy + mom.xxz*xz + mom.xyy*yy + mom.xyz*yz + mom.xzz*zz;
	Q3miry = mom.xxy*xx + mom.xyy*xy + mom.xyz*xz + mom.yyy*yy + mom.yyz*yz + mom.yzz*zz;
	Q3mirz = mom.xxz*xx + mom.xyz*xy + mom.xzz*xz + mom.yyz*yy + mom.yzz*yz + mom.zzz*zz;
	Q4mirx = mom.xxxx*xxx + mom.xxxy*xxy + mom.xxxz*xxz + mom.xxyy*xyy + mom.xxyz*xyz +
	    mom.xxzz*xzz + mom.xyyy*yyy + mom.xyyz*yyz + mom.xyzz*yzz + mom.xzzz*zzz;
	Q4miry = mom.xxxy*xxx + mom.xxyy*xxy + mom.xxyz*xxz + mom.xyyy*xyy + mom.xyyz*xyz +
	    mom.xyzz*xzz + mom.yyyy*yyy + mom.yyyz*yyz + mom.yyzz*yzz + mom.yzzz*zzz;
	Q4mirz = mom.xxxz*xxx + mom.xxyz*xxy + mom.xxzz*xxz + mom.xyyz*xyy + mom.xyzz*xyz +
	    mom.xzzz*xzz + mom.yyyz*yyy + mom.yyzz*yyz + mom.yzzz*yzz + mom.zzzz*zzz;
	Q4x = pkd->ew.Q4xx*x + pkd->ew.Q4xy*y + pkd->ew.Q4xz*z;
	Q4y = pkd->ew.Q4xy*x + pkd->ew.Q4yy*y + pkd->ew.Q4yz*z;
	Q4z = pkd->ew.Q4xz*x + pkd->ew.Q4yz*y + pkd->ew.Q4zz*z;
	Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (pkd->ew.Q3x*x + pkd->ew.Q3y*y + pkd->ew.Q3z*z) + pkd->ew.Q4;
	Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
	Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
	Qta = g1*mom.m - g2*pkd->ew.Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
	fPot -= g0*mom.m - g1*pkd->ew.Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
	ax += (g2*(Q2mirx - pkd->ew.Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta);
	ay += (g2*(Q2miry - pkd->ew.Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta);
	az += (g2*(Q2mirz - pkd->ew.Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta);
	}
    if ( nLoop == pkd->ew.nEwLoopInner) {
	alphan = pkd->ew.ka;
	r2 = dx*dx + dy*dy + dz*dz;
	r2 *= pkd->ew.alpha2;
	g0 = alphan*((1.0/3.0)*r2 - 1.0);
	alphan *= 2*pkd->ew.alpha2;
	g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
	alphan *= 2*pkd->ew.alpha2;
	g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
	alphan *= 2*pkd->ew.alpha2;
	g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
	alphan *= 2*pkd->ew.alpha2;
	g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
	alphan *= 2*pkd->ew.alpha2;
	g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
	evalEwald(&pkd->ew,&ax,&ay,&az,&fPot,dx,dy,dz,&mom,g0,g1,g2,g3,g4,g5);
	}
#else

    for (ix=-pkd->ew.nEwReps;ix<=pkd->ew.nEwReps;++ix) {
#ifdef USE_FORMALLY_CORRECT_EWALD
	bInHolex = (abs(ix) <= pkd->ew.nReps);
#endif
	x = dx + ix*L;
	for (iy=-pkd->ew.nEwReps;iy<=pkd->ew.nEwReps;++iy) {
#ifdef USE_FORMALLY_CORRECT_EWALD
	    bInHolexy = (bInHolex && abs(iy) <= pkd->ew.nReps);
#endif
	    y = dy + iy*L;
	    for (iz=-pkd->ew.nEwReps;iz<=pkd->ew.nEwReps;++iz) {
#ifdef USE_FORMALLY_CORRECT_EWALD
		bInHole = (bInHolexy && abs(iz) <= pkd->ew.nReps);
#endif
		/*
		** Scoring for Ewald inner stuff = (+,*)
		**		Visible ops 		= (104,161)
		**		sqrt, 1/sqrt est. 	= (6,11)
		**     division            = (6,11)  same as sqrt.
		**		exp est.			= (6,11)  same as sqrt.
		**		erf/erfc est.		= (12,22) twice a sqrt.
		**		Total			= (128,205) = 333
		**     Old scoring				    = 447
		*/
		z = dz + iz*L;
		r2 = x*x + y*y + z*z;
#ifdef USE_FORMALLY_CORRECT_EWALD
		if (r2 > pkd->ew.fEwCut2 && !bInHole) continue;
#endif
		if (r2 < pkd->ew.fInner2) {
		    /*
		     * For small r, series expand about
		     * the origin to avoid errors caused
		     * by cancellation of large terms.
		     */
		    alphan = pkd->ew.ka;
		    r2 *= pkd->ew.alpha2;
		    g0 = alphan*((1.0/3.0)*r2 - 1.0);
		    alphan *= 2*pkd->ew.alpha2;
		    g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
		    alphan *= 2*pkd->ew.alpha2;
		    g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
		    }
		else {
		    dir = 1/sqrt(r2);
		    dir2 = dir*dir;
		    a = exp(-r2*pkd->ew.alpha2);
		    a *= pkd->ew.ka*dir2;
		    /*
		    ** We can get away with not using erfc because it is vanishing
		    ** if the following are true:
		    ** 1. The center of mass is near/at the centre of the box,
		    ** 2. We use 1 replica for RMS error < ~1e-4 (erfc 3 = 2.2e-5)
		    ** 3. We use 2 replicas for smaller RMS  (erfc 5 = 1.5e-12)
		    ** Worst case with COM near the edge:
		    ** 1 replica:  erfc 2 = 4.6e-3 (this would be wrong)
		    ** 2 replicas: erfc 4 = 1.5e-8
		    ** We currently switch from 1 replica to 2 replicates at at
		    ** theta of 0.5 which corresponds to an RMS error of ~1e-4.
		    */
#ifdef USE_FORMALLY_CORRECT_EWALD
		    if (bInHole) {
			g0 = -erf(pkd->ew.alpha*r2*dir);
			}
		    else {
			g0 = erfc(pkd->ew.alpha*r2*dir);
			}
#else
		    g0 = -erf(pkd->ew.alpha*r2*dir);
#endif
		    g0 *= dir;
		    g1 = g0*dir2 + a;
		    alphan = 2*pkd->ew.alpha2;
		    g2 = 3*g1*dir2 + alphan*a;
		    alphan *= 2*pkd->ew.alpha2;
		    g3 = 5*g2*dir2 + alphan*a;
		    alphan *= 2*pkd->ew.alpha2;
		    g4 = 7*g3*dir2 + alphan*a;
		    alphan *= 2*pkd->ew.alpha2;
		    g5 = 9*g4*dir2 + alphan*a;
		    }
		evalEwald(&pkd->ew,&ax,&ay,&az,&fPot,x,y,z,&mom,g0,g1,g2,g3,g4,g5);
		++nLoop;
		}
	    }
	}
#endif

    /*
    ** Scoring for the h-loop (+,*)
    ** 	Without trig = (10,14)
    **	    Trig est.	 = 2*(6,11)  same as 1/sqrt scoring.
    **		Total        = (22,36)
    **					 = 58
    */
    for (i=0;i<pkd->ew.nEwhLoop;++i) {
	hdotx = pkd->ew.ewt.hx.f[i]*dx + pkd->ew.ewt.hy.f[i]*dy + pkd->ew.ewt.hz.f[i]*dz;
	c = cosf(hdotx);
	s = sinf(hdotx);
	fPot += pkd->ew.ewt.hCfac.f[i]*c + pkd->ew.ewt.hSfac.f[i]*s;
	t = pkd->ew.ewt.hCfac.f[i]*s - pkd->ew.ewt.hSfac.f[i]*c;
	ax += pkd->ew.ewt.hx.f[i]*t;
	ay += pkd->ew.ewt.hy.f[i]*t;
	az += pkd->ew.ewt.hz.f[i]*t;
	}
    *pPot += fPot;
    pa[0] += ax;
    pa[1] += ay;
    pa[2] += az;
    nFlop = nLoop*447 + pkd->ew.nEwhLoop*58;
    return(nFlop);
    }

void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut) {
    MOMC mom = pkd->momRoot;
    int i,ix,iy,iz,hReps,hx,hy,hz,h2;
    double k4,L;
    double gam[6],mfacc,mfacs;
    double ax,ay,az,dx,dy,dz;
    const int iOrder = 4;

    L = pkd->fPeriod[0];

    /*
    ** We have made certain assumptions about the COM in order to optimize
    ** the erf/erfc parts. In particular, the erf function is accurate only
    ** in the range [0.65,6.0]. The range for the replicas are:
    **   1 replica:  [1.0,3.0] - COM
    **   2 replicas: [3.0,5.0] - COM
    **   3 replicas: [5.0,7.0] - COM
    **   4 replicas: [7.0,9.0] - COM
    */
    if (nReps>0) {
	assert(pkdTopNode(pkd,ROOT)->r[0]/L * pkd->ew.alpha < 0.35);
	assert(pkdTopNode(pkd,ROOT)->r[1]/L * pkd->ew.alpha < 0.35);
	assert(pkdTopNode(pkd,ROOT)->r[1]/L * pkd->ew.alpha < 0.35);
	assert(nReps<=4);
	}

    /*
    ** Create SIMD versions of the moments.
    */
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    pkd->ew.ewm.m.p = SIMD_DSPLAT(mom.m);
    pkd->ew.ewm.xx.p = SIMD_DSPLAT(mom.xx);
    pkd->ew.ewm.yy.p = SIMD_DSPLAT(mom.yy);
    pkd->ew.ewm.xy.p = SIMD_DSPLAT(mom.xy);
    pkd->ew.ewm.xz.p = SIMD_DSPLAT(mom.xz);
    pkd->ew.ewm.yz.p = SIMD_DSPLAT(mom.yz);
    pkd->ew.ewm.xxx.p = SIMD_DSPLAT(mom.xxx);
    pkd->ew.ewm.xyy.p = SIMD_DSPLAT(mom.xyy);
    pkd->ew.ewm.xxy.p = SIMD_DSPLAT(mom.xxy);
    pkd->ew.ewm.yyy.p = SIMD_DSPLAT(mom.yyy);
    pkd->ew.ewm.xxz.p = SIMD_DSPLAT(mom.xxz);
    pkd->ew.ewm.yyz.p = SIMD_DSPLAT(mom.yyz);
    pkd->ew.ewm.xyz.p = SIMD_DSPLAT(mom.xyz);
    pkd->ew.ewm.xxxx.p = SIMD_DSPLAT(mom.xxxx);
    pkd->ew.ewm.xyyy.p = SIMD_DSPLAT(mom.xyyy);
    pkd->ew.ewm.xxxy.p = SIMD_DSPLAT(mom.xxxy);
    pkd->ew.ewm.yyyy.p = SIMD_DSPLAT(mom.yyyy);
    pkd->ew.ewm.xxxz.p = SIMD_DSPLAT(mom.xxxz);
    pkd->ew.ewm.yyyz.p = SIMD_DSPLAT(mom.yyyz);
    pkd->ew.ewm.xxyy.p = SIMD_DSPLAT(mom.xxyy);
    pkd->ew.ewm.xxyz.p = SIMD_DSPLAT(mom.xxyz);
    pkd->ew.ewm.xyyz.p = SIMD_DSPLAT(mom.xyyz);
    pkd->ew.ewm.zz.p = SIMD_DSPLAT(mom.zz);
    pkd->ew.ewm.xzz.p = SIMD_DSPLAT(mom.xzz);
    pkd->ew.ewm.yzz.p = SIMD_DSPLAT(mom.yzz);
    pkd->ew.ewm.zzz.p = SIMD_DSPLAT(mom.zzz);
    pkd->ew.ewm.xxzz.p = SIMD_DSPLAT(mom.xxzz);
    pkd->ew.ewm.xyzz.p = SIMD_DSPLAT(mom.xyzz);
    pkd->ew.ewm.xzzz.p = SIMD_DSPLAT(mom.xzzz);
    pkd->ew.ewm.yyzz.p = SIMD_DSPLAT(mom.yyzz);
    pkd->ew.ewm.yzzz.p = SIMD_DSPLAT(mom.yzzz);
    pkd->ew.ewm.zzzz.p = SIMD_DSPLAT(mom.zzzz);
#endif

    /*
    ** Set up traces of the complete multipole moments.
    */
    pkd->ew.Q4xx = 0.5*(mom.xxxx + mom.xxyy + mom.xxzz);
    pkd->ew.Q4xy = 0.5*(mom.xxxy + mom.xyyy + mom.xyzz);
    pkd->ew.Q4xz = 0.5*(mom.xxxz + mom.xyyz + mom.xzzz);
    pkd->ew.Q4yy = 0.5*(mom.xxyy + mom.yyyy + mom.yyzz);
    pkd->ew.Q4yz = 0.5*(mom.xxyz + mom.yyyz + mom.yzzz);
    pkd->ew.Q4zz = 0.5*(mom.xxzz + mom.yyzz + mom.zzzz);
    pkd->ew.Q4 = 0.25*(pkd->ew.Q4xx + pkd->ew.Q4yy + pkd->ew.Q4zz);
    pkd->ew.Q3x = 0.5*(mom.xxx + mom.xyy + mom.xzz);
    pkd->ew.Q3y = 0.5*(mom.xxy + mom.yyy + mom.yzz);
    pkd->ew.Q3z = 0.5*(mom.xxz + mom.yyz + mom.zzz);
    pkd->ew.Q2 = 0.5*(mom.xx + mom.yy + mom.zz);
    pkd->ew.nReps = nReps;
#ifdef USE_FORMALLY_CORRECT_EWALD
    pkd->ew.nEwReps = d2i(ceil(fEwCut));
    pkd->ew.nEwReps = pkd->ew.nEwReps > nReps ? pkd->ew.nEwReps : nReps;
#else
    pkd->ew.nEwReps = nReps;
#endif
    pkd->ew.fEwCut2 = fEwCut*fEwCut*L*L;
    pkd->ew.fInner2 = 1.2e-3*L*L;
    pkd->ew.alpha = 2.0/L;
    pkd->ew.alpha2 = pkd->ew.alpha*pkd->ew.alpha;
    pkd->ew.k1 = M_PI/(pkd->ew.alpha2*L*L*L);
    pkd->ew.ka = 2.0*pkd->ew.alpha/sqrt(M_PI);

#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    pkd->ew.ewp.Q4xx.p = SIMD_DSPLAT(pkd->ew.Q4xx);
    pkd->ew.ewp.Q4xy.p = SIMD_DSPLAT(pkd->ew.Q4xy);
    pkd->ew.ewp.Q4xz.p = SIMD_DSPLAT(pkd->ew.Q4xz);
    pkd->ew.ewp.Q4yy.p = SIMD_DSPLAT(pkd->ew.Q4yy);
    pkd->ew.ewp.Q4yz.p = SIMD_DSPLAT(pkd->ew.Q4yz);
    pkd->ew.ewp.Q4zz.p = SIMD_DSPLAT(pkd->ew.Q4zz);
    pkd->ew.ewp.Q4.p = SIMD_DSPLAT(pkd->ew.Q4);
    pkd->ew.ewp.Q3x.p = SIMD_DSPLAT(pkd->ew.Q3x);
    pkd->ew.ewp.Q3y.p = SIMD_DSPLAT(pkd->ew.Q3y);
    pkd->ew.ewp.Q3z.p = SIMD_DSPLAT(pkd->ew.Q3z);
    pkd->ew.ewp.Q2.p = SIMD_DSPLAT(pkd->ew.Q2);
    pkd->ew.ewp.fEwCut2.p = SIMD_DSPLAT(pkd->ew.fEwCut2);
    pkd->ew.ewp.fInner2.p = SIMD_DSPLAT(pkd->ew.fInner2);
    pkd->ew.ewp.alpha.p = SIMD_DSPLAT(pkd->ew.alpha);
    pkd->ew.ewp.alpha2.p = SIMD_DSPLAT(pkd->ew.alpha2);
    pkd->ew.ewp.k1.p = SIMD_DSPLAT(pkd->ew.k1);
    pkd->ew.ewp.ka.p = SIMD_DSPLAT(pkd->ew.ka);
#endif


    /*
    ** Now setup stuff for the h-loop.
    */
    hReps = d2i(ceil(fhCut));
    k4 = M_PI*M_PI/(pkd->ew.alpha*pkd->ew.alpha*L*L);

    i = (int)pow(1+2*hReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE__)
    i = (i + SIMD_MASK) & ~SIMD_MASK;
#endif
    if ( i>pkd->ew.nMaxEwhLoop ) {
	pkd->ew.nMaxEwhLoop = i;
	pkd->ew.ewt.hx.f = SIMD_malloc(pkd->ew.nMaxEwhLoop*sizeof(pkd->ew.ewt.hx.f));
	assert(pkd->ew.ewt.hx.f != NULL);
	pkd->ew.ewt.hy.f = SIMD_malloc(pkd->ew.nMaxEwhLoop*sizeof(pkd->ew.ewt.hy.f));
	assert(pkd->ew.ewt.hy.f != NULL);
	pkd->ew.ewt.hz.f = SIMD_malloc(pkd->ew.nMaxEwhLoop*sizeof(pkd->ew.ewt.hz.f));
	assert(pkd->ew.ewt.hz.f != NULL);
	pkd->ew.ewt.hCfac.f = SIMD_malloc(pkd->ew.nMaxEwhLoop*sizeof(pkd->ew.ewt.hCfac.f));
	assert(pkd->ew.ewt.hCfac.f != NULL);
	pkd->ew.ewt.hSfac.f = SIMD_malloc(pkd->ew.nMaxEwhLoop*sizeof(pkd->ew.ewt.hSfac.f));
	assert(pkd->ew.ewt.hSfac.f != NULL);
	}
    pkd->ew.nEwhLoop = i;
    i = (int)pow(1+2*pkd->ew.nEwReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    i = (i + SIMD_MASK) & ~SIMD_MASK;
#endif
    if ( i>pkd->ew.nMaxEwLoopInner ) {
	pkd->ew.nMaxEwLoopInner = i;
	pkd->ew.ewt.Lx.d = SIMD_malloc(pkd->ew.nMaxEwLoopInner*sizeof(pkd->ew.ewt.Lx.d));
	assert(pkd->ew.ewt.Lx.d!=NULL);
	pkd->ew.ewt.Ly.d = SIMD_malloc(pkd->ew.nMaxEwLoopInner*sizeof(pkd->ew.ewt.Ly.d));
	assert(pkd->ew.ewt.Ly.d!=NULL);
	pkd->ew.ewt.Lz.d = SIMD_malloc(pkd->ew.nMaxEwLoopInner*sizeof(pkd->ew.ewt.Lz.d));
	assert(pkd->ew.ewt.Lz.d!=NULL);
	i = 0;
	for (ix=-pkd->ew.nEwReps;ix<=pkd->ew.nEwReps;++ix) {
	    dx = ix*L;
	    for (iy=-pkd->ew.nEwReps;iy<=pkd->ew.nEwReps;++iy) {
		dy = iy*L;
		for (iz=-pkd->ew.nEwReps;iz<=pkd->ew.nEwReps;++iz) {
		    if (ix||iy||iz) { /* We treat 0,0,0 specially */
			dz = iz*L;
			pkd->ew.ewt.Lx.d[i] = dx;
			pkd->ew.ewt.Ly.d[i] = dy;
			pkd->ew.ewt.Lz.d[i] = dz;
			++i;
			}
		    }
		}
	    }
	pkd->ew.nEwLoopInner = i;
	}
    while ( i<pkd->ew.nMaxEwLoopInner ) {
	pkd->ew.ewt.Lx.d[i] = 0;
	pkd->ew.ewt.Ly.d[i] = 0;
	pkd->ew.ewt.Lz.d[i] = 0;
	++i;
	}
    i = 0;
    for (hx=-hReps;hx<=hReps;++hx) {
	for (hy=-hReps;hy<=hReps;++hy) {
	    for (hz=-hReps;hz<=hReps;++hz) {
		h2 = hx*hx + hy*hy + hz*hz;
		if (h2 == 0) continue;
		if (h2 > fhCut*fhCut) continue;
		assert (i < pkd->ew.nMaxEwhLoop);
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
		QEVAL(iOrder,pkd->momRoot,gam,hx,hy,hz,ax,ay,az,mfacc);
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
		QEVAL(iOrder,pkd->momRoot,gam,hx,hy,hz,ax,ay,az,mfacs);
		pkd->ew.ewt.hx.f[i] = 2*M_PI/L*hx;
		pkd->ew.ewt.hy.f[i] = 2*M_PI/L*hy;
		pkd->ew.ewt.hz.f[i] = 2*M_PI/L*hz;
		pkd->ew.ewt.hCfac.f[i] = mfacc;
		pkd->ew.ewt.hSfac.f[i] = mfacs;
		++i;
		}
	    }
	}
    pkd->ew.nEwhLoop = i;
    while(i<pkd->ew.nMaxEwhLoop) {
	pkd->ew.ewt.hx.f[i] = 0;
	pkd->ew.ewt.hy.f[i] = 0;
	pkd->ew.ewt.hz.f[i] = 0;
	pkd->ew.ewt.hCfac.f[i] = 0;
	pkd->ew.ewt.hSfac.f[i] = 0;
	++i;
	}
    }
