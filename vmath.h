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

#ifndef VMATH_H
#define VMATH_H
#include <math.h>
#include "simd.h"

const float vmath_DP1 = 2*0.78515625;
const float vmath_DP2 = 2*2.4187564849853515625e-4;
const float vmath_DP3 = 2*3.77489497744594108e-8;

static inline
fvec sinf(const fvec &xx) {
    fvec sin_sign = xx & fvec::sign_mask(); // Save sign of argument
    fvec x = abs(xx);
    fvec y = (x * (1.27323954473516f*0.5f) );
    i32v j = cvt_i32v(y);
    y = cvt_fvec(j);
    j = j & 3;
    fvec fj = cvt_fvec(j);
    fmask mask = fj >= 2.0f;
    sin_sign = mask_xor(mask,sin_sign,fvec::sign_mask());
    fj = mask_sub(mask,fj,fvec(2.0f));
    fmask poly = fj >= 1.0f;
    x = ((x - y * vmath_DP1) - y * vmath_DP2) - y * vmath_DP3;
    fvec z = x * x;
    fvec y1 = (((
	  2.443315711809948E-005f * z
        - 1.388731625493765E-003f) * z
        + 4.166664568298827E-002f) * z
        - 0.5f) * z + 1.0f;
    fvec y2 = ((
	 -1.9515295891E-4f * z
	+ 8.3321608736E-3f) * z
	- 1.6666654611E-1f) * z * x + x;
    y = mask_mov(y2,poly,y1);
    return y ^ sin_sign;
    }

static inline
fvec cosf(const fvec &xx) {
    fvec cos_sign; // assume positive
    fvec x = abs(xx);
    fvec y = (x * (1.27323954473516f*0.5f) );
    i32v j = cvt_i32v(y);
    y = cvt_fvec(j);
    j = j & 3;
    fvec fj = cvt_fvec(j);
    fmask mask = fj >= 2.0f;
    fj = mask_sub(mask,fj,fvec(2.0f));
    mask ^= (fj >= 1.0f);
    cos_sign = maskz_mov(mask,fvec::sign_mask());
    fmask poly = fj >= 1.0f;
    x = ((x - y * vmath_DP1) - y * vmath_DP2) - y * vmath_DP3;
    fvec z = x * x;
    fvec y1 = (((
	  2.443315711809948E-005f * z
	- 1.388731625493765E-003f) * z
	+ 4.166664568298827E-002f) * z
	- 0.5f) * z + 1.0f;
    fvec y2 = ((
	 -1.9515295891E-4f * z
	+ 8.3321608736E-3f) * z
	- 1.6666654611E-1f) * z * x + x;
    y = mask_mov(y1,poly,y2);
    return y ^ cos_sign;
    }

static inline
void sincosf(const fvec &xx, fvec &sin, fvec &cos) {
    fvec cos_sign;
    fvec sin_sign = xx & fvec::sign_mask(); // Save sign of argument
    fvec x = abs(xx);
    fvec y = (x * (1.27323954473516f*0.5f) );
    i32v j = cvt_i32v(y);
    y = cvt_fvec(j);
    j = j & 3;
    fvec fj = cvt_fvec(j);
    fmask mask = fj >= 2.0f;
    sin_sign = mask_xor(mask,sin_sign,fvec::sign_mask());
    fj = mask_sub(mask,fj,fvec(2.0f));
    mask ^= (fj >= 1.0f);
    cos_sign = maskz_mov(mask,fvec::sign_mask());
    fmask poly = fj >= 1.0f;
    x = ((x - y * vmath_DP1) - y * vmath_DP2) - y * vmath_DP3;
    fvec z = x * x;
    fvec y1 = (((
	  2.443315711809948E-005f * z
        - 1.388731625493765E-003f) * z
        + 4.166664568298827E-002f) * z
        - 0.5f) * z + 1.0f;
    fvec y2 = ((
	 -1.9515295891E-4f * z
	+ 8.3321608736E-3f) * z
	- 1.6666654611E-1f) * z * x + x;
    sin = mask_mov(y2,poly,y1) ^ sin_sign;
    cos = mask_mov(y1,poly,y2) ^ cos_sign;
    }

/* Cost: AVX: 30 ops, SSE: 33 */
static inline
dvec exp(const dvec &x) {
    dvec d,xx,pow2n,Pexp,Qexp;

    const dvec C1 = 6.93145751953125E-1;
    const dvec C2 = 1.42860682030941723212E-6;
    /* exp polynomial table */
    const dvec p0 = 1.26177193074810590878E-4;
    const dvec p1 = 3.02994407707441961300E-2;
    const dvec p2 = 9.99999999999999999910E-1;
    const dvec q0 = 3.00198505138664455042E-6;
    const dvec q1= 2.52448340349684104192E-3;
    const dvec q2 = 2.27265548208155028766E-1;
    const dvec q3 = 2.00000000000000000009E0;

    d = M_LOG2E * x;
#if defined(__AVX512F__)
    __m256i n1 = _mm512_cvtpd_epi32(d);
    d = _mm512_cvtepi32_pd(n1);
    n1  = _mm256_add_epi32(n1, _mm256_set1_epi32(1023));
    __m512i n = _mm512_permutexvar_epi32(
	_mm512_setr_epi32(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7),
	_mm512_castsi256_si512(n1));
    pow2n = _mm512_castsi512_pd(_mm512_slli_epi64(n,52));
#elif defined(__AVX__)
   __m128i n1 = _mm256_cvttpd_epi32(d);
    d = _mm256_cvtepi32_pd(n1);
    n1  = _mm_add_epi32(n1, _mm_set1_epi32(1023));
    __m128i n2 = _mm_shuffle_epi32(n1,_MM_SHUFFLE(3, 3, 2, 2));
    n1 = _mm_shuffle_epi32(n1,_MM_SHUFFLE(1, 1, 0, 0));
    n2 = _mm_slli_epi64(n2, 52);
    n1 = _mm_slli_epi64(n1, 52);
    __m256i n = _mm256_castsi128_si256(n1);
    n = _mm256_insertf128_si256(n, n2, 0x1);
    pow2n = _mm256_castsi256_pd(n);
#else
    __m128i n1 = _mm_cvttpd_epi32(d);
    d = _mm_cvtepi32_pd(n1);
    n1  = _mm_add_epi32(n1, _mm_set1_epi32(1023));
    n1  = _mm_shuffle_epi32(n1, _MM_SHUFFLE(2, 1, 2, 0));
    n1  = _mm_slli_epi64(n1, 52);
    pow2n = _mm_castsi128_pd(n1);
#endif
    dvec y = x - d*C1 - d*C2;
    xx = y * y;

    Pexp = ((p0*xx + p1)*xx + p2)*y;
    Qexp = (((q0*xx + q1)*xx + q2)*xx + q3) - Pexp;
    dvec r = Pexp / Qexp;
    r = (r*2.0f + 1.0f)*pow2n;

    return r;
    }

/*
** This is an optimized version of erf that is accurate for all ranges.
** The exp() value is passed in because it it needed outside as well.
** Cost: AVX: 105, SSE: 172
*/
static inline
dvec verf(const dvec &v,const dvec &iv,const dvec &ex2,dvec &r_erf,dvec &r_erfc) {
    dvec Perf,Qerf,v2;
    dvec t,nt;
    const dvec threshold0 = 0.65;
    const dvec threshold1 = 2.2;
    const dvec threshold2 = 6.0;

    dmask pred0 = v >= threshold0;
    dmask pred1 = v >= threshold1;
    dmask pred2 = v >= threshold2;

const struct CONSTS {
#if defined(__AVX512F__)
    dvec::array_t p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,q5;
    i64v::array_t one,two;
#elif defined(__AVX2__)
    dvec::array_t p0,p1,p2,p3,p4,p5,q0,q1,q2,q3,q4,q5;
    i32v::array_t init,two;
#elif defined(__AVX__)
    dvec::array_t p0ab,p0cd,p1ab,p1cd,p2ab,p2cd,p3ab,p3cd,p4ab,p4cd,p5ab,p5cd;
    dvec::array_t q0ab,q0cd,q1ab,q1cd,q2ab,q2cd,q3ab,q3cd,q4ab,q4cd,q5ab,q5cd;
#else
    dvec::array_t p0a,p0b,p0c,p0d, p1a,p1b,p1c,p1d, p2a,p2b,p2c,p2d;
    dvec::array_t p3a,p3b,p3c,p3d, p4a,p4b,p4c,p4d, p5a,p5b,p5c,p5d;
    dvec::array_t q0a,q0b,q0c,q0d, q1a,q1b,q1c,q1d, q2a,q2b,q2c,q2d;
    dvec::array_t q3a,q3b,q3c,q3d, q4a,q4b,q4c,q4d, q5a,q5b,q5c,q5d;
#endif
    } consts = {
#if defined(__AVX512F__)
#define ERF_CONSTS(d,c,b,a) {a,b,c,d,a,b,c,d}
#elif defined(__AVX2__)
#define ERF_CONSTS(d,c,b,a) {a,b,c,d}
#elif defined(__AVX__)
#define ERF_CONSTS(d,c,b,a) {a,b,a,b},{c,d,c,d}
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
#if defined(__AVX512F__)
	{1,1,1,1,1,1,1,1},{2,2,2,2,2,2,2,2},
#elif defined(__AVX2__)
	{0,1,0,1,0,1,0,1},{2,2,2,2,2,2,2,2}
#endif
    };
#if defined(__AVX512F__)
#define SET_PREFACTOR(q) dvec q = _mm512_permutexvar_pd(idx,dvec(consts.q))
#elif defined(__AVX2__)
#define SET_PREFACTOR(q) dvec q = _mm256_castsi256_pd(_mm256_permutevar8x32_epi32(_mm256_castpd_si256(dvec(consts.q)),idx));
#elif defined(__AVX__)
#define SET_PREFACTOR(q) dvec q = _mm256_blendv_pd(_mm256_permutevar_pd(dvec(consts.q##ab),_mm256_castpd_si256(pred0)),_mm256_permutevar_pd(dvec(consts.q##cd),_mm256_castpd_si256(pred2)),pred1)
#else
#define SET_PREFACTOR(q) dvec q = dvec(consts.q##a) + (pred0&dvec(consts.q##b)) + (pred1&dvec(consts.q##c)) + (pred2&dvec(consts.q##d))
#endif

#if defined(__AVX512F__)
i64v idx =
	_mm512_mask_add_epi64(
	    _mm512_maskz_mov_epi64(pred0,i64v(consts.one)),
	    pred1,i64v(consts.two),
	    _mm512_maskz_mov_epi64(pred2,i64v(consts.one)));
#elif defined(__AVX2__)
    i32v idx = i32v(consts.init)
	+ i32v(_mm256_and_si256(i32v(consts.two),_mm256_castpd_si256(pred0)))
	+ i32v(_mm256_and_si256(i32v(consts.two),_mm256_castpd_si256(pred1)))
	+ i32v(_mm256_and_si256(i32v(consts.two),_mm256_castpd_si256(pred2)));
#endif


    SET_PREFACTOR(p0);
    SET_PREFACTOR(p1);
    SET_PREFACTOR(p2);
    SET_PREFACTOR(p3);
    SET_PREFACTOR(p4);
    SET_PREFACTOR(p5);

    SET_PREFACTOR(q0);
    SET_PREFACTOR(q1);
    SET_PREFACTOR(q2);
    SET_PREFACTOR(q3);
    SET_PREFACTOR(q4);
    SET_PREFACTOR(q5);

    /* Gives x^2, x, x, 1/x^2 */
    v2 = mask_mov(v,pred2,iv); /* x unless x>6.0, then 1/x */
    v2 = mask_mov(v2*v2,pred0^pred2,v2);

    Perf = (((((p0*v2 + p1)*v2 + p2)*v2 + p3)*v2 + p4)*v2 + p5);
    Qerf = (((((q0*v2 + q1)*v2 + q2)*v2 + q3)*v2 + q4)*v2 + q5)*v2 + 1.0;
    t = mask_mov(dvec(1.0),pred2,iv) * mask_mov(v,pred0,ex2) * (Perf/Qerf);
    nt = 1.0 - t;
    r_erf = mask_mov(t,pred0,nt);
    r_erfc = mask_mov(nt,pred0,t);

    return t;
    }
#endif/*VMATH_H*/
