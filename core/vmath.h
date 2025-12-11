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
    mask ^= (fmask)(fj >= 1.0f);
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
    mask ^= (fmask)(fj >= 1.0f);
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
    dvec xx,pow2n,Pexp,Qexp;

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

    dvec d = M_LOG2E * x;
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

static inline dvec log_poly(const dvec &m) {
    // m in [1,2), t = m - 1 in [0,1]
    const dvec t  = m - 1.0;
    dvec tpow = 1.0;

    // Q(t) = a0 + a1 t + ... + a12 t^12
    // log(1+t) ≈ t * Q(t)
    const dvec a0  = 9.99999999940960782e-01;
    const dvec a1  = -4.99999988894528569e-01;
    const dvec a2  = 3.33332814506182407e-01;
    const dvec a3  = -2.49989386749519288e-01;
    const dvec a4  = 1.99880786360059820e-01;
    const dvec a5  = -1.65837785714768954e-01;
    const dvec a6  = 1.39020881042502709e-01;
    const dvec a7  = -1.12595369219835192e-01;
    const dvec a8  = 8.20632177197394008e-02;
    const dvec a9  = -4.91258177129472937e-02;
    const dvec a10 = 2.16911029784248484e-02;
    const dvec a11 = -6.09830024373608868e-03;
    const dvec a12 = 8.05026589618270359e-04;

    dvec q = a0;
    tpow = t;           q += a1  * tpow;
    tpow = tpow * t;    q += a2  * tpow;
    tpow = tpow * t;    q += a3  * tpow;
    tpow = tpow * t;    q += a4  * tpow;
    tpow = tpow * t;    q += a5  * tpow;
    tpow = tpow * t;    q += a6  * tpow;
    tpow = tpow * t;    q += a7  * tpow;
    tpow = tpow * t;    q += a8  * tpow;
    tpow = tpow * t;    q += a9  * tpow;
    tpow = tpow * t;    q += a10 * tpow;
    tpow = tpow * t;    q += a11 * tpow;
    tpow = tpow * t;    q += a12 * tpow;

    // log(m) ≈ t * Q(t)
    return t * q;
}

#if defined(__AVX512F__)
static inline dvec log(const dvec &x) {
    const __m512i exp_mask = _mm512_set1_epi64(0x7FF0000000000000ULL);
    const __m512i man_mask = _mm512_set1_epi64(0x000FFFFFFFFFFFFFULL);
    const __m512i one_exp  = _mm512_set1_epi64(1023ULL << 52);

    __m512i xi = _mm512_castpd_si512(x);

    // Extract exponent (unbiased)
    __m512i e = _mm512_srli_epi64(_mm512_and_si512(xi, exp_mask), 52);
    e = _mm512_sub_epi64(e, _mm512_set1_epi64(1023));

    // Normalized mantissa m ∈ [1,2)
    __m512i mi = _mm512_or_si512(_mm512_and_si512(xi, man_mask), one_exp);
    dvec m = _mm512_castsi512_pd(mi);

    dvec logm = log_poly(m);

    // e → double
    dvec ed = _mm512_cvtepi64_pd(e);

    const dvec ln2 = 0.6931471805599453094172321215;
    return logm + ed * ln2;
}
#elif defined(__AVX__)
static inline dvec log(const dvec &x) {
    const __m256d exp_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FF0000000000000ULL));
    const __m256d man_mask = _mm256_castsi256_pd(_mm256_set1_epi64x(0x000FFFFFFFFFFFFFULL));
    const __m256d one_exp  = _mm256_castsi256_pd(_mm256_set1_epi64x(1023ULL << 52));

    // Extract exponent (unbiased)
    __m256i e64 = _mm256_castpd_si256(_mm256_and_pd(x, exp_mask));
#if defined(__AVX2__) // AVX2 has 64-bit integer shift and subtract
    e64 = _mm256_srli_epi64(e64, 52);
    e64 = _mm256_sub_epi64(e64, _mm256_set1_epi64x(1023));
#endif
    // Convert exponent lanes → doubles (via two SSE halves)
    // We want _mm256_cvtepi64_pd which only appears in AVX512
    __m128i e_lo = _mm256_castsi256_si128(e64);
    __m128i e_hi = _mm256_extractf128_si256(e64, 1);
#if defined(__AVX2__)
    // e → 32-bit integers for conversion (lanes 0 and 2)
    __m128i e32_lo = _mm_shuffle_epi32(e_lo, _MM_SHUFFLE(3,3,2,0));
    __m128i e32_hi = _mm_shuffle_epi32(e_hi, _MM_SHUFFLE(3,3,2,0));
#else
    // e → 32-bit integers for conversion (lanes 1 and 3)
    __m128i e32_lo = _mm_shuffle_epi32(e_lo, _MM_SHUFFLE(0,0,3,1));
    __m128i e32_hi = _mm_shuffle_epi32(e_hi, _MM_SHUFFLE(0,0,3,1));
    e32_lo = _mm_srli_epi32(e32_lo, 20);
    e32_hi = _mm_srli_epi32(e32_hi, 20);
    e32_lo = _mm_sub_epi32(e32_lo, _mm_set1_epi32(1023));
    e32_hi = _mm_sub_epi32(e32_hi, _mm_set1_epi32(1023));
#endif

    __m128d ed_lo = _mm_cvtepi32_pd(e32_lo);
    __m128d ed_hi = _mm_cvtepi32_pd(e32_hi);

    dvec ed = _mm256_castpd128_pd256(ed_lo);
    ed = _mm256_insertf128_pd(ed, ed_hi, 1);

    // Normalized mantissa m ∈ [1,2)
    dvec m = _mm256_or_pd(_mm256_and_pd(x, man_mask), one_exp);

    dvec logm = log_poly(m);

    const dvec ln2 = 0.6931471805599453094172321215;
    return logm + ed * ln2;
}
#elif defined(__SSE2__) || defined(__ARM_NEON)
static inline dvec log(const dvec &x) {
    const __m128i exp_mask = _mm_set1_epi64x(0x7FF0000000000000ULL);
    const __m128i man_mask = _mm_set1_epi64x(0x000FFFFFFFFFFFFFULL);
    const __m128i one_exp  = _mm_set1_epi64x(1023ULL << 52);

    __m128i xi = _mm_castpd_si128(x);

    // Extract exponent (unbiased)
    __m128i e64 = _mm_srli_epi64(_mm_and_si128(xi, exp_mask), 52);
    e64 = _mm_sub_epi64(e64, _mm_set1_epi64x(1023));

    // Pack lower 32-bit parts for conversion
    __m128i e32 = _mm_shuffle_epi32(e64, _MM_SHUFFLE(2,2,2,0));
    dvec ed = _mm_cvtepi32_pd(e32);

    // Normalized mantissa m in [1,2)
    __m128i mi = _mm_or_si128(_mm_and_si128(xi, man_mask), one_exp);
    dvec m = _mm_castsi128_pd(mi);

    dvec logm = log_poly(m);

    const dvec ln2 = 0.6931471805599453094172321215;
    return logm + ed * ln2;
}

#endif

static inline fvec exp( fvec xx ) {
    static const float LOG2EF = 1.44269504088896341f;
    static const fvec C1 =   0.693359375f;
    static const fvec C2 =  -2.12194440e-4f;

    fvec z = floor( LOG2EF * xx + 0.5f ); /* floor() truncates toward -infinity. */
    const fvec z_min = -126.0f;
    const fvec z_max =  127.0f;
    z = max(z, z_min);
    z = min(z, z_max);
    i32v n = cvt_i32v(z);   // now safely in valid exponent range

    auto x = xx - z * C1 - z * C2; // Cody–Waite reduction

    z = x * x;
    z =
        ((((( 1.9875691500E-4  * x
              + 1.3981999507E-3) * x
            + 8.3334519073E-3) * x
           + 4.1665795894E-2) * x
          + 1.6666665459E-1) * x
         + 5.0000001201E-1) * z
        + x
        + 1.0;

    // exponent bits: (n + 127) << 23
    const i32v bias = 127;
    i32v e = (n + i32v(127)) << 23;

    // reinterpret as float -> scale = 2^n
    fvec scale = cast_fvec(e);
    z *= scale;
    z = mask_mov(z, xx >  88.72283905206835, std::numeric_limits<float>::infinity());
    z = mask_mov(z, xx < -103.278929903431851103, 0.0f);
    return z;
}

static inline fvec log(const fvec &x_in) {
    // Save original input for edge-case handling (x <= 0)
    const i32v exp_mask  = 0x7F800000;     // exponent bits
    const i32v mans_mask = 0x807FFFFF;     // mantissa bits + sign bit
    const i32v half_exp  = 0x3F000000;     // exponent = 126 → mantissa in [0.5,1)
    const fvec zero      = 0.0f;

    // Subnormal handling
    const fvec flt_min  = std::numeric_limits<float>::min();  // 1.175494e-38
    const fvec scale24  = 16777216.0f;
    auto pos_mask       = x_in > 0.0f;
    auto sub_mask       = pos_mask & (x_in < flt_min);
    fvec x_scaled       = x_in * scale24;
    auto x_norm         = mask_mov(x_in, sub_mask, x_scaled);
    //auto sub_adj        = maskz_mov(cast_i32v(sub_mask),i32v(24));
    auto sub_adj        = maskz_mov(sub_mask,i32v(24));

    // Now treat x_norm as the input for frexp()
    i32v xi = cast_i32v(x_norm);

    // 1. Extract exponent: (exp_bits >> 23) - 126
    i32v e32 = ((xi & exp_mask) >> 23) - 126 - sub_adj;

    // 2. Build mantissa in [0.5,1):
    //    Replace exponent with 126, keep sign & mantissa.
    //
    i32v mi = (xi & mans_mask) | half_exp;
    fvec m = cast_fvec(mi);

    fvec SQRTHF_v = 0.707106781186547524f;
    auto mask = m < SQRTHF_v;

    i32v one_i  = 1;
    //i32v mask_i = cast_i32v(mask);
    //e32 -= mask_i & one_i;
    e32 -= maskz_mov(mask,one_i);

    fvec e = cvt_fvec(e32);

    fvec x = mask_mov(m - 1.0f, mask, 2*m - 1.0f);

    auto z = x * x;

    auto y =
        (((((((( 7.0376836292E-2 * x
                 - 1.1514610310E-1) * x
               + 1.1676998740E-1) * x
              - 1.2420140846E-1) * x
             + 1.4249322787E-1) * x
            - 1.6668057665E-1) * x
           + 2.0000714765E-1) * x
          - 2.4999993993E-1) * x
         + 3.3333331174E-1) * x * z;

    auto nonzero_mask = e != 0.0f;

    y += maskz_mov(nonzero_mask,-2.12194440e-4f * e);
    y -= 0.5 * z;
    z = x + y;
    z += maskz_mov(nonzero_mask, 0.693359375f * e);

    // Handle edge cases
    const fvec nanv    = std::numeric_limits<float>::quiet_NaN();
    const fvec pos_inf = std::numeric_limits<float>::infinity();
    const fvec neg_inf = -std::numeric_limits<float>::infinity();
    z = mask_mov(z, x_in == 0.0f,    neg_inf);  // x == 0 → −inf
    z = mask_mov(z, x_in < 0.0f,     nanv);     // x < 0 → NaN (includes −inf)
    z = mask_mov(z, x_in == pos_inf, pos_inf);  // +inf → +inf
    z = mask_mov(z, x_in != x_in,    nanv);     // NaN → NaN

    return z;
}

static inline dvec pow(const dvec &x, const dvec &y) {
    return exp( y * log(x) );
}

static inline fvec pow(const fvec &x, const fvec &y) {
    return exp( y * log(x) );
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
    } consts alignas(dvec) = {
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

extern "C" {
#if defined(HAVE_LIBVECM_POW)
    __m128d _ZGVbN2vv_pow(__m128d x, __m128d y);
#endif
#if defined(HAVE_LIBVECM256_POW)
    __m256d _ZGVdN4vv_pow(__m256d x, __m256d y);
#endif
#if defined(HAVE_LIBVECM512_POW)
    __m512d _ZGVeN8vv_pow(__m512d x, __m512d y);
#endif
}

#endif/*VMATH_H*/
