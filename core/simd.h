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

#ifndef SIMD_H
#define SIMD_H

#include "pkd_config.h"
#include <cstdint>
#include <math.h>

#ifdef USE_SIMD
    #if defined(__SSE__)
        #include <xmmintrin.h>
        #ifdef __SSE2__
            #include <emmintrin.h>
        #endif
        #ifdef __SSE3__
            #include <pmmintrin.h>
        #endif
        #if defined(__SSE4_1__)
            #include <smmintrin.h>
        #endif
        #ifdef __AVX__
            #include <immintrin.h>
        #endif
        #ifdef __FMA4__
            #include <x86intrin.h>
        #endif
    #elif defined(__ARM_NEON)
        #include "sse2neon.h"
        #define __SSE__
        #define __SSE2__
        #define __SSE3__
        #define __SSE4_1__
    #endif/*__SSE__*/
#endif/*USE_SIMD*/

#if __GNUC__ > 5
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

/**********************************************************************\
* SIMD Vector class template
\**********************************************************************/
template<typename vtype,typename ftype>
struct vec {
public:
    using vector_t = vtype;
    using scalar_t = ftype;
    static constexpr int width() { return sizeof(vector_t)/sizeof(scalar_t); }
    static constexpr int mask()  { return width()-1; }
    using array_t  = scalar_t[width()];
private:
    vector_t ymm;
public:
    vec &load1(scalar_t f);
    vec &load(const scalar_t *f);
    vec &zero();

    vec() = default;
#ifdef USE_SIMD
    vec(scalar_t const &d);
#endif
    vec(vector_t const &d) { ymm = d; }
    explicit vec(scalar_t const *d) { load(d); }
    operator vector_t() const { return ymm; }
    scalar_t operator [] (std::uint32_t idx) const;
    const vec &store(scalar_t *f) const;
    static const vec sign_mask();
    vec<vector_t,scalar_t> operator-() const;
};

template<typename vtype>
class mmask {
public:
    typedef vtype vector_t;
private:
    vector_t ymm;
public:
    mmask() {}
    mmask(vector_t const &d) { ymm = d; }
    operator vector_t() const { return ymm; }
};

#if defined(__AVX512F__) && defined(USE_SIMD)
typedef vec<__m512i,std::int32_t> i32v;
typedef vec<__m512,float> fvec;
typedef mmask<__mmask16> fmask;
typedef vec<__m512i,std::int64_t> i64v;
typedef vec<__m512d,double> dvec;
typedef mmask<__mmask8> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v(_mm512_cvtps_epi32(a)); }
//inline i64v cvt_i64v(const fvec &a) { return i64v(_mm512_cvtps_epi64(a)); }
inline fvec cvt_fvec(const i32v &a) { return fvec(_mm512_cvtepi32_ps(a)); }
inline fvec cvt_fvec(const dvec &a) { return fvec(_mm512_castps256_ps512(_mm512_cvtpd_ps(a))); }
inline fvec cast_fvec(const i32v &a) { return fvec(_mm512_castsi512_ps(a)); }
inline dvec cast_dvec(const i64v &a) { return dvec(_mm512_castsi512_pd(a)); }

//inline fvec cvt_fvec(const i64v &a) { return fvec(_mm512_cvtepi64_ps(a)); }
#elif defined(__AVX__) && defined(USE_SIMD)
typedef vec<__m256i,std::int32_t> i32v;
typedef vec<__m256,float> fvec;
typedef vec<__m256,float> fmask;
typedef vec<__m256i,std::int64_t> i64v;
typedef vec<__m256d,double> dvec;
typedef vec<__m256d,double> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v(_mm256_cvtps_epi32(a)); }
inline fvec cvt_fvec(const i32v &a) { return fvec(_mm256_cvtepi32_ps(a)); }
inline fvec cvt_fvec(const dvec &a) { return fvec(_mm256_castps128_ps256(_mm256_cvtpd_ps(a))); }
inline fvec cast_fvec(const i32v &a) { return fvec(_mm256_castsi256_ps(a)); }
inline dvec cast_dvec(const i64v &a) { return dvec(_mm256_castsi256_pd(a)); }
#elif defined(__SSE__) && defined(USE_SIMD)
typedef vec<__m128i,std::int32_t> i32v;
typedef vec<__m128,float> fvec;
typedef vec<__m128,float> fmask;
#if defined(__SSE2__)
typedef vec<__m128i,std::int64_t> i64v;
typedef vec<__m128d,double> dvec;
typedef vec<__m128d,double> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v(_mm_cvtps_epi32(a)); }
inline fvec cvt_fvec(const i32v &a) { return fvec(_mm_cvtepi32_ps(a)); }
inline fvec cvt_fvec(const dvec &a) { return fvec(_mm_cvtpd_ps(a)); }
inline fvec cast_fvec(const i32v &a) { return fvec(_mm_castsi128_ps(a)); }
inline dvec cast_dvec(const i64v &a) { return dvec(_mm_castsi128_pd(a)); }
#endif/*__SSE2__*/
#else/*__AVX512F__,__AVX__,__SSE2__*/
typedef vec<std::int32_t,std::int32_t> i32v;
typedef vec<float,float> fvec;
typedef vec<double,double> dvec;
typedef vec<std::int64_t,std::int64_t> i64v;
typedef mmask<bool> fmask;
typedef mmask<bool> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v((std::int32_t)a); }
inline fvec cvt_fvec(const i32v &a) { return fvec((float)a); }
inline fvec cvt_dvec(const dvec &a) { return fvec((float)a); }
#endif/*__AVX512F__,__AVX__,__SSE2__*/

#if !defined(__CUDACC__) && !defined(__METAL_VERSION__)

/**********************************************************************\
* Generic Operators
\**********************************************************************/

template<typename v,typename ftype> inline vec<v,ftype> abs(vec<v,ftype> const &a) { return a & ~vec<v,ftype>::sign_mask(); }
template<typename v,typename ftype> inline vec<v,ftype> &operator+=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a + b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator+=(vec<v,ftype> &a,ftype const &b) { return a = a + b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator-=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a - b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator-=(vec<v,ftype> &a,ftype const &b) { return a = a - b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator*=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a * b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator*=(vec<v,ftype> &a,ftype const &b) { return a = a * b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator/=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a / b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator/=(vec<v,ftype> &a,ftype const &b) { return a = a / b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator&=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a & b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator&=(vec<v,ftype> &a,ftype const &b) { return a = a & b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator|=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a | b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator|=(vec<v,ftype> &a,ftype const &b) { return a = a | b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator^=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a ^ b; }
template<typename v,typename ftype> inline vec<v,ftype> &operator^=(vec<v,ftype> &a,ftype const &b) { return a = a ^ b; }
template<typename v,typename ftype> inline ftype vec<v,ftype>::operator [] (std::uint32_t idx) const {
    ftype d[sizeof(vector_t)/sizeof(scalar_t)];
    store(d);
    return d[idx];
}
template<typename v> inline mmask<v> &operator&=(mmask<v> &a,mmask<v> const &b) { return a = a & b; }
template<typename v> inline mmask<v> &operator|=(mmask<v> &a,mmask<v> const &b) { return a = a | b; }
template<typename v> inline mmask<v> &operator^=(mmask<v> &a,mmask<v> const &b) { return a = a ^ b; }

#if defined(__AVX512F__) && defined(USE_SIMD)

/**********************************************************************\
* AVX512 single precision
\**********************************************************************/

template<> inline vec<__m512,float>::vec(const float &d) { ymm = _mm512_set1_ps(d); }
template<> inline vec<__m512,float> &vec<__m512,float>::zero() { ymm = _mm512_setzero_ps(); return *this; }
template<> inline vec<__m512,float> &vec<__m512,float>::load1(float f) { ymm = _mm512_setr_ps(f,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512,float> &vec<__m512,float>::load(const float *pf) { ymm = _mm512_loadu_ps(pf); return *this; }
template<> inline const vec<__m512,float> &vec<__m512,float>::store(float *pf) const { _mm512_storeu_ps(pf,ymm); return *this; }
template<> inline const vec<__m512,float> vec<__m512,float>::sign_mask() { return _mm512_castsi512_ps(_mm512_set1_epi32(0x80000000)); }
inline vec<__m512,float> min(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_min_ps(a,b); }
inline vec<__m512,float> max(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_max_ps(a,b); }
//inline vec<__m512,float> abs(vec<__m512,float> const &a) { return _mm512_abs_ps(a); }
inline vec<__m512,float> operator*(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_mul_ps(a,b); }
inline vec<__m512,float> operator/(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_div_ps(a,b); }
inline vec<__m512,float> operator+(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_add_ps(a,b); }
inline vec<__m512,float> operator-(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_sub_ps(a,b); }
#if defined(__AVX512ER__)
inline vec<__m512,float> rsqrt(vec<__m512,float> const &a) { return _mm512_rsqrt28_ps(a); }
#else
inline vec<__m512,float> rsqrt(vec<__m512,float> const &a2) {
    vec<__m512,float> ia =  _mm512_rsqrt14_ps(a2); // ia = (approx) 1.0 / sqrt(a2)
    ia *= 1.5f - 0.5f*a2*ia*ia; // Do one Newton step to double the precision
    return ia;
}
#endif
inline mmask<__mmask16> operator==(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_cmp_ps_mask(a,b,_CMP_EQ_OQ); }
inline mmask<__mmask16> operator!=(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_cmp_ps_mask(a,b,_CMP_NEQ_OQ); }
inline mmask<__mmask16> operator>(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_cmp_ps_mask(a,b,_CMP_GT_OQ); }
inline mmask<__mmask16> operator<(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_cmp_ps_mask(a,b,_CMP_LT_OQ); }
inline mmask<__mmask16> operator>=(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_cmp_ps_mask(a,b,_CMP_GE_OQ); }
inline mmask<__mmask16> operator<=(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_cmp_ps_mask(a,b,_CMP_LE_OQ); }
inline vec<__m512,float> sqrt(vec<__m512,float> const &r2) { return _mm512_sqrt_ps(r2); }
#ifdef __AVX512DQ__
inline vec<__m512,float> operator&(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_and_ps(a,b); }
inline vec<__m512,float> operator|(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_or_ps(a,b); }
inline vec<__m512,float> operator^(vec<__m512,float> const &a,vec<__m512,float> const &b) { return _mm512_xor_ps(a,b); }
#else
inline vec<__m512,float> operator&(vec<__m512,float> const &a,vec<__m512,float> const &b)
{ return _mm512_castsi512_ps(_mm512_and_epi32( _mm512_castps_si512(a),_mm512_castps_si512(b))); }
inline vec<__m512,float> operator|(vec<__m512,float> const &a,vec<__m512,float> const &b)
{ return _mm512_castsi512_ps(_mm512_or_epi32( _mm512_castps_si512(a),_mm512_castps_si512(b))); }
inline vec<__m512,float> operator^(vec<__m512,float> const &a,vec<__m512,float> const &b)
{ return _mm512_castsi512_ps(_mm512_xor_epi32( _mm512_castps_si512(a),_mm512_castps_si512(b))); }
#endif
template<> inline vec<__m512,float> vec<__m512,float>::operator-() const {return ymm ^ sign_mask(); }
inline vec<__m512,float> operator~(vec<__m512,float> const &a)
{ return _mm512_castsi512_ps(_mm512_xor_epi32( _mm512_castps_si512(a),_mm512_set1_epi32(0xffffffff))); }
inline float hadd(__m256 p) {
    __m128 sum;
    __m128 upper = _mm256_extractf128_ps(p, 1);
    __m128 lower = _mm256_castps256_ps128(p);
    sum = _mm_add_ps(lower,upper);
    __m128 shuf = _mm_movehdup_ps(sum);
    sum = _mm_add_ps(sum, shuf);
    shuf = _mm_movehl_ps(sum, sum);
    sum = _mm_add_ps(sum, shuf);
    return _mm_cvtss_f32(sum);
}
inline float hadd(vec<__m512,float> const &a) {
    __m256 lower = _mm512_castps512_ps256(a);
    __m256 upper = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(a),1));
    return hadd(_mm256_add_ps(lower,upper));
}

inline vec<__m512,float> mask_xor(mmask<__mmask16> const &k,vec<__m512,float> const &a,vec<__m512,float> const &b)
{ return _mm512_castsi512_ps(_mm512_mask_xor_epi32( _mm512_castps_si512(a),k,_mm512_castps_si512(a),_mm512_castps_si512(b))); }
inline vec<__m512,float> mask_sub(mmask<__mmask16> const &k,vec<__m512,float> const &a,vec<__m512,float> const &b)
{ return _mm512_mask_sub_ps(a,k,a,b); }
inline vec<__m512,float> mask_mov(vec<__m512,float> const &src,mmask<__mmask16> const &k,vec<__m512,float> const &a)
{ return _mm512_mask_mov_ps(src,k,a); }
inline vec<__m512,float> maskz_mov(mmask<__mmask16> const &k,vec<__m512,float> const &a)
{ return _mm512_maskz_mov_ps(k,a); }

/**********************************************************************\
* AVX512 32-bit integer
\**********************************************************************/
template<> inline vec<__m512i,std::int32_t>::vec(const std::int32_t &d) { ymm = _mm512_set1_epi32(d); }
template<> inline vec<__m512i,std::int32_t> &vec<__m512i,std::int32_t>::zero() { ymm = _mm512_setzero_si512(); return *this; }
template<> inline vec<__m512i,std::int32_t> &vec<__m512i,std::int32_t>::load1(std::int32_t f) { ymm = _mm512_setr_epi32(f,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512i,std::int32_t> &vec<__m512i,std::int32_t>::load(const std::int32_t *pf) { ymm = _mm512_loadu_si512((__m512i *)pf); return *this; }
template<> inline const vec<__m512i,std::int32_t> &vec<__m512i,std::int32_t>::store(std::int32_t *pf) const { _mm512_storeu_si512((__m512i *)pf,ymm); return *this; }
template<> inline const vec<__m512i,std::int32_t> vec<__m512i,std::int32_t>::sign_mask() { return _mm512_set1_epi32(0x80000000); }
template<> inline vec<__m512i,std::int32_t> vec<__m512i,std::int32_t>::operator-() const {
    return _mm512_castps_si512(_mm512_xor_ps(_mm512_castsi512_ps(ymm),_mm512_castsi512_ps(sign_mask())));
}
inline vec<__m512i,std::int32_t> operator&(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_and_epi32(a,b); }
inline vec<__m512i,std::int32_t> operator|(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_or_epi32(a,b); }
inline vec<__m512i,std::int32_t> operator^(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_xor_epi32(a,b); }

inline mmask<__mmask16> operator==(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ); }
inline mmask<__mmask16> operator!=(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE); }
inline mmask<__mmask16> operator>(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GT); }
inline mmask<__mmask16> operator<(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LT); }
inline mmask<__mmask16> operator>=(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GE); }
inline mmask<__mmask16> operator<=(vec<__m512i,std::int32_t> const &a,vec<__m512i,std::int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LE); }

inline vec<__m512i,std::int32_t> mask_mov(vec<__m512i,std::int32_t> const &src,mmask<__mmask16> const &k,vec<__m512i,std::int32_t> const &a)
{ return _mm512_mask_mov_epi32(src,k,a); }

/**********************************************************************\
* AVX512 64-bit integer
\**********************************************************************/
template<> inline vec<__m512i,std::int64_t>::vec(const std::int64_t &d) { ymm = _mm512_set1_epi64(d); }
template<> inline vec<__m512i,std::int64_t> &vec<__m512i,std::int64_t>::zero() { ymm = _mm512_setzero_si512(); return *this; }
template<> inline vec<__m512i,std::int64_t> &vec<__m512i,std::int64_t>::load1(std::int64_t f) { ymm = _mm512_setr_epi64(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512i,std::int64_t> &vec<__m512i,std::int64_t>::load(const std::int64_t *pf) { ymm = _mm512_loadu_si512((__m512i *)pf); return *this; }
template<> inline const vec<__m512i,std::int64_t> &vec<__m512i,std::int64_t>::store(std::int64_t *pf) const { _mm512_storeu_si512((__m512i *)pf,ymm); return *this; }
template<> inline const vec<__m512i,std::int64_t> vec<__m512i,std::int64_t>::sign_mask() { return _mm512_set1_epi64(0x8000000000000000); }
template<> inline vec<__m512i,std::int64_t> vec<__m512i,std::int64_t>::operator-() const {
    return _mm512_castps_si512(_mm512_xor_ps(_mm512_castsi512_ps(ymm),_mm512_castsi512_ps(sign_mask())));
}
inline vec<__m512i,std::int64_t> operator&(vec<__m512i,std::int64_t> const &a,vec<__m512i,std::int64_t> const &b) { return _mm512_and_epi64(a,b); }
inline vec<__m512i,std::int64_t> operator|(vec<__m512i,std::int64_t> const &a,vec<__m512i,std::int64_t> const &b) { return _mm512_or_epi64(a,b); }
inline vec<__m512i,std::int64_t> operator^(vec<__m512i,std::int64_t> const &a,vec<__m512i,std::int64_t> const &b) { return _mm512_xor_epi64(a,b); }

/**********************************************************************\
* AVX512 single precision mask
\**********************************************************************/

inline mmask<__mmask16> operator&(mmask<__mmask16> const &a,mmask<__mmask16> const &b) { return _mm512_kand(a,b); }
inline mmask<__mmask16> operator|(mmask<__mmask16> const &a,mmask<__mmask16> const &b) { return _mm512_kor(a,b); }
inline mmask<__mmask16> operator^(mmask<__mmask16> const &a,mmask<__mmask16> const &b) { return _mm512_kxor(a,b); }
inline mmask<__mmask16> operator~(mmask<__mmask16> const &a) { return _mm512_knot(a); }
inline int testz(mmask<__mmask16> const &a) { return _mm512_kortestz(a,a); }

/**********************************************************************\
* AVX512 double precision
\**********************************************************************/

template<> inline vec<__m512d,double>::vec(const double &d) { ymm = _mm512_set1_pd(d); }
template<> inline vec<__m512d,double> &vec<__m512d,double>::zero() { ymm = _mm512_setzero_pd(); return *this; }
template<> inline vec<__m512d,double> &vec<__m512d,double>::load1(double f) { ymm = _mm512_setr_pd(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512d,double> &vec<__m512d,double>::load(const double *pf) { ymm = _mm512_loadu_pd(pf); return *this; }
template<> inline const vec<__m512d,double> &vec<__m512d,double>::store(double *pf) const { _mm512_storeu_pd(pf,ymm); return *this; }
template<> inline const vec<__m512d,double> vec<__m512d,double>::sign_mask() { return _mm512_castsi512_pd(_mm512_set1_epi64(0x8000000000000000)); }
inline vec<__m512d,double> min(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_min_pd(a,b); }
inline vec<__m512d,double> max(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_max_pd(a,b); }
//inline vec<__m512d,double> abs(vec<__m512d,double> const &a) { return _mm512_abs_pd(a); }
inline vec<__m512d,double> operator*(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_mul_pd(a,b); }
inline vec<__m512d,double> operator/(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_div_pd(a,b); }
inline vec<__m512d,double> operator+(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_add_pd(a,b); }
inline vec<__m512d,double> operator-(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_sub_pd(a,b); }
#if defined(__AVX512ER__)
inline vec<__m512d,double> rsqrt(vec<__m512d,double> const &a2) {
    vec<__m512d,double> ia = _mm512_rsqrt28_pd(a2); // ia = (approx) 1.0 / sqrt(a2)
    ia *= 1.5 - 0.5*a2*ia*ia; // Do one Newton step to double the precision
    return ia;
}
#else
inline vec<__m512d,double> rsqrt(vec<__m512d,double> const &a2) {
    vec<__m512d,double> ia =  _mm512_rsqrt14_pd(a2); // ia = (approx) 1.0 / sqrt(a2)
    ia *= 1.5 - 0.5*a2*ia*ia; // Do one Newton step to double the precision
    ia *= 1.5 - 0.5*a2*ia*ia; // Do one Newton step to double the precision
    return ia;
}
#endif

inline mmask<__mmask8> operator==(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_cmp_pd_mask(a,b,_CMP_EQ_OQ); }
inline mmask<__mmask8> operator!=(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_cmp_pd_mask(a,b,_CMP_NEQ_OQ); }
inline mmask<__mmask8> operator>(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_cmp_pd_mask(a,b,_CMP_GT_OQ); }
inline mmask<__mmask8> operator<(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_cmp_pd_mask(a,b,_CMP_LT_OQ); }
inline mmask<__mmask8> operator>=(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_cmp_pd_mask(a,b,_CMP_GE_OQ); }
inline mmask<__mmask8> operator<=(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_cmp_pd_mask(a,b,_CMP_LE_OQ); }
inline vec<__m512d,double> sqrt(vec<__m512d,double> const &r2) { return _mm512_sqrt_pd(r2); }
#ifdef __AVX512DQ__
inline vec<__m512d,double> operator&(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_and_pd(a,b); }
inline vec<__m512d,double> operator|(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_or_pd(a,b); }
inline vec<__m512d,double> operator^(vec<__m512d,double> const &a,vec<__m512d,double> const &b) { return _mm512_xor_pd(a,b); }
#else
inline vec<__m512d,double> operator&(vec<__m512d,double> const &a,vec<__m512d,double> const &b)
{ return _mm512_castsi512_pd(_mm512_and_epi32( _mm512_castpd_si512(a),_mm512_castpd_si512(b))); }
inline vec<__m512d,double> operator|(vec<__m512d,double> const &a,vec<__m512d,double> const &b)
{ return _mm512_castsi512_pd(_mm512_or_epi32( _mm512_castpd_si512(a),_mm512_castpd_si512(b))); }
inline vec<__m512d,double> operator^(vec<__m512d,double> const &a,vec<__m512d,double> const &b)
{ return _mm512_castsi512_pd(_mm512_xor_epi32( _mm512_castpd_si512(a),_mm512_castpd_si512(b))); }
#endif
template<> inline vec<__m512d,double> vec<__m512d,double>::operator-() const {return ymm ^ sign_mask(); }

inline vec<__m512d,double> mask_mov(vec<__m512d,double> const &src,mmask<__mmask8> const &k,vec<__m512d,double> const &a)
{ return _mm512_mask_mov_pd(src,k,a); }

/**********************************************************************\
* AVX512 double precision mask
\**********************************************************************/

inline mmask<__mmask8> operator&(mmask<__mmask8> const &a,mmask<__mmask8> const &b) { return _mm512_kand(a,b); }
inline mmask<__mmask8> operator|(mmask<__mmask8> const &a,mmask<__mmask8> const &b) { return _mm512_kor(a,b); }
inline mmask<__mmask8> operator^(mmask<__mmask8> const &a,mmask<__mmask8> const &b) { return _mm512_kxor(a,b); }
inline mmask<__mmask8> operator~(mmask<__mmask8> const &a) { return _mm512_knot(a); }
inline int testz(mmask<__mmask8> const &a) { return _mm512_kortestz(a,a); }
inline int movemask(mmask<__mmask8> const &k) { return (int)(k); }

#else/*__AVX512F__*/

/**********************************************************************\
* Mask Operators
\**********************************************************************/

//inline mmask<__m256> operator&(mmask<__m256> const &a,mmask<__m256> const &b) { return _mm512_kand(a,b); }
//inline mmask<__m256> operator|(mmask<__m256> const &a,mmask<__m256> const &b) { return _mm512_kor(a,b); }
//inline mmask<__m256> operator^(mmask<__m256> const &a,mmask<__m256> const &b) { return _mm512_kxor(a,b); }

template<typename v,typename ftype,typename mtype> inline
vec<v,ftype> mask_xor(const mtype &k,vec<v,ftype> &a,vec<v,ftype> const &b)
{ return a ^ (b&k); }
template<typename v,typename ftype,typename mtype> inline
vec<v,ftype> mask_sub(const mtype &k,vec<v,ftype> &a,vec<v,ftype> const &b)
{ return a - (b&k); }

#if defined(__AVX__) && defined(USE_SIMD)
/**********************************************************************\
* AVX single precision
\**********************************************************************/

template<> inline vec<__m256,float>::vec(const float &d) { ymm = _mm256_set1_ps(d); }
template<> inline vec<__m256,float> &vec<__m256,float>::zero() { ymm = _mm256_setzero_ps(); return *this; }
template<> inline vec<__m256,float> &vec<__m256,float>::load1(float f) { ymm = _mm256_setr_ps(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m256,float> &vec<__m256,float>::load(const float *pf) { ymm = _mm256_loadu_ps(pf); return *this; }
template<> inline const vec<__m256,float> &vec<__m256,float>::store(float *pf) const { _mm256_storeu_ps(pf,ymm); return *this; }
template<> inline const vec<__m256,float> vec<__m256,float>::sign_mask() { return _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)); }
template<> inline vec<__m256,float> vec<__m256,float>::operator-() const { return _mm256_xor_ps(ymm,sign_mask()); }
inline vec<__m256,float> min(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_min_ps(a,b); }
inline vec<__m256,float> max(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_max_ps(a,b); }
inline vec<__m256,float> sqrt(vec<__m256,float> const &r2) { return _mm256_sqrt_ps(r2); }
inline vec<__m256,float> operator*(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_mul_ps(a,b); }
inline vec<__m256,float> operator/(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_div_ps(a,b); }
inline vec<__m256,float> operator+(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_add_ps(a,b); }
inline vec<__m256,float> operator-(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_sub_ps(a,b); }
inline vec<__m256,float> rsqrt(vec<__m256,float> const &a2) {
    vec<__m256,float> ia =  _mm256_rsqrt_ps(a2); // ia = (approx) 1.0 / sqrt(a2)
    ia *= 1.5f - 0.5f*a2*ia*ia; // Do one Newton step to double the precision
    return ia;
}
inline vec<__m256,float> operator==(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_EQ_OQ); }
inline vec<__m256,float> operator!=(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_NEQ_OQ); }
inline vec<__m256,float> operator>(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_GT_OQ); }
inline vec<__m256,float> operator<(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_LT_OQ); }
inline vec<__m256,float> operator>=(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_GE_OQ); }
inline vec<__m256,float> operator<=(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_LE_OQ); }
inline vec<__m256,float> operator&(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_and_ps(a,b); }
inline vec<__m256,float> operator|(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_or_ps(a,b); }
inline vec<__m256,float> operator^(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_xor_ps(a,b); }
inline vec<__m256,float> operator~(vec<__m256,float> const &a) { return _mm256_xor_ps(a,_mm256_castsi256_ps(_mm256_set1_epi32(0xffffffff))); }
inline float hadd(vec<__m256,float> const &p) {
    __m128 sum;
    __m128 upper = _mm256_extractf128_ps(p, 1);
    __m128 lower = _mm256_castps256_ps128(p);
    sum = _mm_add_ps(lower,upper);
    __m128 shuf = _mm_movehdup_ps(sum);
    sum = _mm_add_ps(sum, shuf);
    shuf = _mm_movehl_ps(sum, sum);
    sum = _mm_add_ps(sum, shuf);
    return _mm_cvtss_f32(sum);
}
inline vec<__m256,float> mask_mov(vec<__m256,float> const &src,vec<__m256,float> const &p,vec<__m256,float> const &a)
{ return _mm256_blendv_ps(src,a,p); }
inline vec<__m256,float> maskz_mov(vec<__m256,float> const &p,vec<__m256,float> const &a) { return a & p; }
inline int testz(vec<__m256,float> const &a) { return !_mm256_movemask_ps(a); }

/**********************************************************************\
* AVX 32-bit integer
\**********************************************************************/
template<> inline vec<__m256i,std::int32_t>::vec(const std::int32_t &d) { ymm = _mm256_set1_epi32(d); }
template<> inline vec<__m256i,std::int32_t> &vec<__m256i,std::int32_t>::zero() { ymm = _mm256_setzero_si256(); return *this; }
template<> inline vec<__m256i,std::int32_t> &vec<__m256i,std::int32_t>::load1(std::int32_t f) { ymm = _mm256_setr_epi32(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m256i,std::int32_t> &vec<__m256i,std::int32_t>::load(const std::int32_t *pf) { ymm = _mm256_loadu_si256((__m256i *)pf); return *this; }
template<> inline const vec<__m256i,std::int32_t> &vec<__m256i,std::int32_t>::store(std::int32_t *pf) const { _mm256_storeu_si256((__m256i *)pf,ymm); return *this; }
template<> inline const vec<__m256i,std::int32_t> vec<__m256i,std::int32_t>::sign_mask() { return _mm256_set1_epi32(0x80000000); }
#ifdef __AVX2__
template<> inline vec<__m256i,std::int32_t> vec<__m256i,std::int32_t>::operator-() const { return _mm256_xor_si256(ymm,sign_mask()); }
inline vec<__m256i,std::int32_t> operator+(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_add_epi32(a,b); }
inline vec<__m256i,std::int32_t> operator-(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_sub_epi32(a,b); }
inline vec<__m256i,std::int32_t> operator&(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_and_si256(a,b); }
inline vec<__m256i,std::int32_t> mask_mov(vec<__m256i,std::int32_t> const &src,vec<__m256i,std::int32_t> const &p,vec<__m256i,std::int32_t> const &a)
{ return _mm256_blendv_epi8(src,a,p); }
inline vec<__m256i,std::int32_t> mask_mov(vec<__m256i,std::int32_t> const &src,vec<__m256,float> const &p,vec<__m256i,std::int32_t> const &a)
{ return _mm256_blendv_epi8(src,a,_mm256_castps_si256(p)); }
#else
template<> inline vec<__m256i,std::int32_t> vec<__m256i,std::int32_t>::operator-() const {
    return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(ymm),_mm256_castsi256_ps(sign_mask())));
}
inline vec<__m256i,std::int32_t> operator&(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b)
{ return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(a),_mm256_castsi256_ps(b))); }
inline vec<__m256i,std::int32_t> mask_mov(vec<__m256i,std::int32_t> const &src,vec<__m256i,std::int32_t> const &p,vec<__m256i,std::int32_t> const &a)
{ return _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(src),_mm256_castsi256_ps(a),_mm256_castsi256_ps(p))); }
inline vec<__m256i,std::int32_t> mask_mov(vec<__m256i,std::int32_t> const &src,vec<__m256,float> const &p,vec<__m256i,std::int32_t> const &a)
{ return _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(src),_mm256_castsi256_ps(a),p)); }
#endif
#ifdef __AVX2__
inline vec<__m256i,std::int32_t> operator==(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_cmpeq_epi32(a,b); }
//inline vec<__m256i,std::int32_t> operator!=(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_cmpneq_epi32(a,b); }
inline vec<__m256i,std::int32_t> operator>(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_cmpgt_epi32(a,b); }
//inline vec<__m256i,std::int32_t> operator<(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_cmplt_epi32(a,b); }
//inline vec<__m256i,std::int32_t> operator>=(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_cmpge_epi32(a,b); }
//inline vec<__m256i,std::int32_t> operator<=(vec<__m256i,std::int32_t> const &a,vec<__m256i,std::int32_t> const &b) { return _mm256_cmple_epi32(a,b); }
#else
#endif

/**********************************************************************\
* AVX double precision
\**********************************************************************/

template<> inline vec<__m256d,double>::vec(const double &d) { ymm = _mm256_set1_pd(d); }
template<> inline vec<__m256d,double> &vec<__m256d,double>::zero() { ymm = _mm256_setzero_pd(); return *this; }
template<> inline vec<__m256d,double> &vec<__m256d,double>::load1(double f) { ymm = _mm256_setr_pd(f,0,0,0); return *this; }
template<> inline vec<__m256d,double> &vec<__m256d,double>::load(const double *pf) { ymm = _mm256_loadu_pd(pf); return *this; }
template<> inline const vec<__m256d,double> &vec<__m256d,double>::store(double *pf) const { _mm256_storeu_pd(pf,ymm); return *this; }
template<> inline const vec<__m256d,double> vec<__m256d,double>::sign_mask() { return _mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000)); }
template<> inline vec<__m256d,double> vec<__m256d,double>::operator-() const { return _mm256_xor_pd(ymm,sign_mask()); }
inline vec<__m256d,double> min(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_min_pd(a,b); }
inline vec<__m256d,double> max(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_max_pd(a,b); }
inline vec<__m256d,double> sqrt(vec<__m256d,double> const &r2) { return _mm256_sqrt_pd(r2); }
inline vec<__m256d,double> operator*(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_mul_pd(a,b); }
inline vec<__m256d,double> operator/(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_div_pd(a,b); }
inline vec<__m256d,double> operator+(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_add_pd(a,b); }
inline vec<__m256d,double> operator-(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_sub_pd(a,b); }
inline vec<__m256d,double> rsqrt(vec<__m256d,double> const &a2) {
    vec<__m256d,double> ia =  _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(a2))); // ia = (approx) 1.0 / sqrt(a2)
    ia *= 1.5 - 0.5*a2*ia*ia; // Do one Newton step to double the precision
    ia *= 1.5 - 0.5*a2*ia*ia; // Do one Newton step to double the precision
    return ia;
}
inline vec<__m256d,double> operator==(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_EQ_OQ); }
inline vec<__m256d,double> operator!=(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_NEQ_OQ); }
inline vec<__m256d,double> operator>(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_GT_OQ); }
inline vec<__m256d,double> operator<(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_LT_OQ); }
inline vec<__m256d,double> operator>=(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_GE_OQ); }
inline vec<__m256d,double> operator<=(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_LE_OQ); }
inline vec<__m256d,double> operator&(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_and_pd(a,b); }
inline vec<__m256d,double> operator|(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_or_pd(a,b); }
inline vec<__m256d,double> operator^(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_xor_pd(a,b); }
inline vec<__m256d,double> operator~(vec<__m256d,double> const &a) { return _mm256_xor_pd(a,_mm256_castsi256_pd(_mm256_set1_epi32(0xffffffff))); }

inline int movemask(vec<__m256d,double> const &r2) { return _mm256_movemask_pd(r2); }
inline vec<__m256d,double> mask_mov(vec<__m256d,double> const &src,vec<__m256d,double> const &k,vec<__m256d,double> const &a)
{ return _mm256_blendv_pd(src,a,k); }
inline int testz(vec<__m256d,double> const &a) { return !_mm256_movemask_pd(a); }

/**********************************************************************\
* AVX 64-bit integer
\**********************************************************************/
template<> inline vec<__m256i,std::int64_t> &vec<__m256i,std::int64_t>::zero() { ymm = _mm256_setzero_si256(); return *this; }
template<> inline vec<__m256i,std::int64_t> &vec<__m256i,std::int64_t>::load1(std::int64_t f) { ymm = _mm256_setr_epi64x(f,0,0,0); return *this; }
template<> inline vec<__m256i,std::int64_t> &vec<__m256i,std::int64_t>::load(const std::int64_t *pf) { ymm = _mm256_loadu_si256((__m256i *)pf); return *this; }

#elif defined(__SSE__) && defined(USE_SIMD)

/**********************************************************************\
* SSE single precision
\**********************************************************************/

template<> inline vec<__m128,float>::vec(const float &d) { ymm = _mm_set1_ps(d); }
template<> inline vec<__m128,float> &vec<__m128,float>::zero() { ymm = _mm_setzero_ps(); return *this; }
template<> inline vec<__m128,float> &vec<__m128,float>::load1(float f) { ymm = _mm_setr_ps(f,0,0,0); return *this; }
template<> inline vec<__m128,float> &vec<__m128,float>::load(const float *pf) { ymm = _mm_loadu_ps(pf); return *this; }
template<> inline const vec<__m128,float> &vec<__m128,float>::store(float *pf) const { _mm_storeu_ps(pf,ymm); return *this; }
template<> inline const vec<__m128,float> vec<__m128,float>::sign_mask() { return _mm_castsi128_ps(_mm_set1_epi32(0x80000000)); }
template<> inline vec<__m128,float> vec<__m128,float>::operator-() const { return _mm_xor_ps(ymm,sign_mask()); }
inline vec<__m128,float> min(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_min_ps(a,b); }
inline vec<__m128,float> max(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_max_ps(a,b); }
inline vec<__m128,float> operator*(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_mul_ps(a,b); }
inline vec<__m128,float> operator/(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_div_ps(a,b); }
inline vec<__m128,float> operator+(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_add_ps(a,b); }
inline vec<__m128,float> operator-(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_sub_ps(a,b); }
inline vec<__m128,float> rsqrt(vec<__m128,float> const &a2) {
    vec<__m128,float> ia =  _mm_rsqrt_ps(a2); // ia = (approx) 1.0 / sqrt(a2)
    ia *= 1.5f - 0.5f*a2*ia*ia; // Do one Newton step to double the precision
    return ia;
}
inline vec<__m128,float> operator==(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpeq_ps(a,b); }
inline vec<__m128,float> operator!=(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpneq_ps(a,b); }
inline vec<__m128,float> operator>(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpgt_ps(a,b); }
inline vec<__m128,float> operator<(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmplt_ps(a,b); }
inline vec<__m128,float> operator>=(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpge_ps(a,b); }
inline vec<__m128,float> operator<=(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmple_ps(a,b); }
inline vec<__m128,float> operator&(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_and_ps(a,b); }
inline vec<__m128,float> operator|(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_or_ps(a,b); }
inline vec<__m128,float> operator^(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_xor_ps(a,b); }
inline vec<__m128,float> operator~(vec<__m128,float> const &a) { return _mm_xor_ps(a,_mm_castsi128_ps(_mm_set1_epi32(0xffffffff))); }
inline float hadd(vec<__m128,float> const &p) {
    __m128 sum = p;
#ifdef __SSE3__
    __m128 shuf = _mm_movehdup_ps(sum);
    sum = _mm_add_ps(sum, shuf);
    shuf = _mm_movehl_ps(sum, sum);
#else
    __m128 shuf = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(0,1,2,3));
    sum = _mm_add_ps(sum, shuf);
    shuf = _mm_shuffle_ps(sum, sum, _MM_SHUFFLE(2,3,0,1));
#endif
    sum = _mm_add_ps(sum, shuf);
    return _mm_cvtss_f32(sum);
}
inline int movemask(vec<__m128,float> const &r2) { return _mm_movemask_ps(r2); }
inline vec<__m128,float> sqrt(vec<__m128,float> const &r2) { return _mm_sqrt_ps(r2); }
inline int testz(vec<__m128,float> const &a) { return !_mm_movemask_ps(a); }
inline vec<__m128,float> maskz_mov(vec<__m128,float> const &p,vec<__m128,float> const &a) { return a & p; }
inline vec<__m128,float> mask_mov(vec<__m128,float> const &src,vec<__m128,float> const &p,vec<__m128,float> const &a) {
#ifdef __SSE4_1__
    return _mm_blendv_ps(src,a,p);
#else
    return _mm_or_ps(_mm_and_ps(p,a),_mm_andnot_ps(p,src));
#endif
}

/**********************************************************************\
* SSE 32-bit integer
\**********************************************************************/
template<> inline vec<__m128i,std::int32_t>::vec(const std::int32_t &d) { ymm = _mm_set1_epi32(d); }
template<> inline vec<__m128i,std::int32_t> &vec<__m128i,std::int32_t>::zero() { ymm = _mm_setzero_si128(); return *this; }
template<> inline vec<__m128i,std::int32_t> &vec<__m128i,std::int32_t>::load1(std::int32_t f) { ymm = _mm_setr_epi32(f,0,0,0); return *this; }
template<> inline vec<__m128i,std::int32_t> &vec<__m128i,std::int32_t>::load(const std::int32_t *pf) { ymm = _mm_loadu_si128((__m128i *)pf); return *this; }
template<> inline vec<__m128i,std::int64_t> &vec<__m128i,std::int64_t>::load(const std::int64_t *pf) { ymm = _mm_loadu_si128((__m128i *)pf); return *this; }
template<> inline const vec<__m128i,std::int32_t> &vec<__m128i,std::int32_t>::store(std::int32_t *pf) const { _mm_storeu_si128((__m128i *)pf,ymm); return *this; }
template<> inline const vec<__m128i,std::int32_t> vec<__m128i,std::int32_t>::sign_mask() { return _mm_set1_epi32(0x80000000); }
#ifdef __SSE2__
template<> inline vec<__m128i,std::int32_t> vec<__m128i,std::int32_t>::operator-() const { return _mm_xor_si128(ymm,sign_mask()); }
inline vec<__m128i,std::int32_t> operator&(vec<__m128i,std::int32_t> const &a,vec<__m128i,std::int32_t> const &b) { return _mm_and_si128(a,b); }
#else
template<> inline vec<__m128i,std::int32_t> vec<__m128i,std::int32_t>::operator-() {
    return _mm_castps_si128(_mm_xor_ps(_mm_castsi128_ps(ymm),_mm_castsi128_ps(sign_mask())));
}
inline vec<__m128i,std::int32_t> operator&(vec<__m128i,std::int32_t> const &a,vec<__m128i,std::int32_t> const &b)
{ return _mm_castps_si128(_mm_and_ps(_mm_castsi128_ps(a),_mm_castsi128_ps(b))); }
#endif

inline vec<__m128i,std::int32_t> mask_mov(vec<__m128i,std::int32_t> const &src,vec<__m128i,std::int32_t> const &p,vec<__m128i,std::int32_t> const &a) {
#ifdef __SSE4_1__
    return _mm_blendv_epi8(src,a,p);
#else
    return _mm_castps_si128(_mm_or_ps(_mm_and_ps(_mm_castsi128_ps(p),_mm_castsi128_ps(a)),_mm_andnot_ps(_mm_castsi128_ps(p),_mm_castsi128_ps(src))));
#endif
}
inline vec<__m128i,std::int32_t> mask_mov(vec<__m128i,std::int32_t> const &src,vec<__m128,float> const &p,vec<__m128i,std::int32_t> const &a) {
#ifdef __SSE4_1__
    return _mm_blendv_epi8(src,a,_mm_castps_si128(p));
#else
    return _mm_castps_si128(_mm_or_ps(_mm_and_ps(p,_mm_castsi128_ps(a)),_mm_andnot_ps(p,_mm_castsi128_ps(src))));
#endif
}

#if defined(__SSE2__)
/**********************************************************************\
* SSE double precision
\**********************************************************************/

template<> inline vec<__m128d,double>::vec(const double &d) { ymm = _mm_set1_pd(d); }
template<> inline vec<__m128d,double> &vec<__m128d,double>::zero() { ymm = _mm_setzero_pd(); return *this; }
template<> inline vec<__m128d,double> &vec<__m128d,double>::load1(double f) { ymm = _mm_setr_pd(f,0); return *this; }
template<> inline vec<__m128d,double> &vec<__m128d,double>::load(const double *pf) { ymm = _mm_loadu_pd(pf); return *this; }
template<> inline const vec<__m128d,double> &vec<__m128d,double>::store(double *pf) const { _mm_storeu_pd(pf,ymm); return *this; }
template<> inline const vec<__m128d,double> vec<__m128d,double>::sign_mask() { return _mm_castsi128_pd(_mm_set1_epi64x(0x8000000000000000)); }
template<> inline vec<__m128d,double> vec<__m128d,double>::operator-() const { return _mm_xor_pd(ymm,sign_mask()); }
inline vec<__m128d,double> min(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_min_pd(a,b); }
inline vec<__m128d,double> max(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_max_pd(a,b); }
inline vec<__m128d,double> operator*(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_mul_pd(a,b); }
inline vec<__m128d,double> operator/(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_div_pd(a,b); }
inline vec<__m128d,double> operator+(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_add_pd(a,b); }
inline vec<__m128d,double> operator-(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_sub_pd(a,b); }
inline vec<__m128d,double> operator==(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_cmpeq_pd(a,b); }
inline vec<__m128d,double> operator!=(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_cmpneq_pd(a,b); }
inline vec<__m128d,double> operator>(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_cmpgt_pd(a,b); }
inline vec<__m128d,double> operator<(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_cmplt_pd(a,b); }
inline vec<__m128d,double> operator>=(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_cmpge_pd(a,b); }
inline vec<__m128d,double> operator<=(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_cmple_pd(a,b); }
inline vec<__m128d,double> operator&(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_and_pd(a,b); }
inline vec<__m128d,double> operator|(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_or_pd(a,b); }
inline vec<__m128d,double> operator^(vec<__m128d,double> const &a,vec<__m128d,double> const &b) { return _mm_xor_pd(a,b); }
inline vec<__m128d,double> operator~(vec<__m128d,double> const &a) { return _mm_xor_pd(a,_mm_castsi128_pd(_mm_set1_epi32(0xffffffff))); }
#if defined(__SSE3__)
inline double hadd(vec<__m128d,double> const &a) {
    return _mm_cvtsd_f64(_mm_hadd_pd(a,a));
}
#endif
inline int movemask(vec<__m128d,double> const &r2) { return _mm_movemask_pd(r2); }
inline int testz(vec<__m128d,double> const &a) { return !_mm_movemask_pd(a); }
inline vec<__m128d,double> sqrt(vec<__m128d,double> const &r2) { return _mm_sqrt_pd(r2); }
inline vec<__m128d,double> rsqrt(vec<__m128d,double> const &r2) {
    vec<__m128d,double> r = _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(r2))); /* Approximation */
    r *= 1.5 - 0.5*r*r*r2; /* Newton step correction */
    return r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
}
#endif/*__SSE2__*/

inline vec<__m128d,double> mask_mov(vec<__m128d,double> const &src,vec<__m128d,double> const &p,vec<__m128d,double> const &a) {
#ifdef __SSE4_1__
    return _mm_blendv_pd(src,a,p);
#else
    return _mm_or_pd(_mm_and_pd(p,a),_mm_andnot_pd(p,src));
#endif
}

#else/*#elif defined(__SSE__) && defined(USE_SIMD)*/
/**********************************************************************\
* single precision
\**********************************************************************/

//template<> inline vec<float,float>::vec(const float &d) { ymm = d; }
template<> inline vec<float,float> &vec<float,float>::zero() { ymm = 0; return *this; }
template<> inline vec<float,float> &vec<float,float>::load1(float f) { ymm = f; return *this; }
template<> inline vec<float,float> &vec<float,float>::load(const float *pf) { ymm = *pf; return *this; }
template<> inline const vec<float,float> &vec<float,float>::store(float *pf) const { *pf = ymm; return *this; }
template<> inline const vec<float,float> vec<float,float>::sign_mask() { return 0x80000000; }
template<> inline vec<float,float> vec<float,float>::operator-() const { return -ymm; }
inline vec<float,float> min(vec<float,float> const &a,vec<float,float> const &b) { return a<b?a:b; }
inline vec<float,float> max(vec<float,float> const &a,vec<float,float> const &b) { return a>b?a:b; }
//inline vec<float,float> operator*(vec<float,float> const &a,vec<float,float> const &b) { return a*b; }
//inline vec<float,float> operator/(vec<float,float> const &a,vec<float,float> const &b) { return a/b; }
//inline vec<float,float> operator+(vec<float,float> const &a,vec<float,float> const &b) { return a+b; }
//inline vec<float,float> operator-(vec<float,float> const &a,vec<float,float> const &b) { return a-b; }

//LET FLOAT DO THE JOB
//inline mmask<bool> operator==(vec<float,float> const &a,vec<float,float> const &b) { return a==b; }
//inline mmask<bool> operator!=(vec<float,float> const &a,vec<float,float> const &b) { return a!=b; }
//inline mmask<bool> operator>(vec<float,float> const &a,vec<float,float> const &b) { return a>b; }
//inline mmask<bool> operator<(vec<float,float> const &a,vec<float,float> const &b) { return a<b; }
//inline mmask<bool> operator>=(vec<float,float> const &a,vec<float,float> const &b) { return a>=b; }
//inline mmask<bool> operator<=(vec<float,float> const &a,vec<float,float> const &b) { return a<=b; }
//inline vec<float,float> operator&(vec<float,float> const &a,vec<float,float> const &b) { return a&b; }
//inline vec<float,float> operator|(vec<float,float> const &a,vec<float,float> const &b) { return a|b; }
//inline vec<float,float> operator^(vec<float,float> const &a,vec<float,float> const &b) { return a^b; }
//inline vec<float,float> operator~(vec<float,float> const &a) { return _mm_xor_ps(a,_mm_castsi128_ps(_mm_set1_epi32(0xffffffff))); }
inline float hadd(vec<float,float> const &a) { return a; }
inline vec<float,float> sqrt(vec<float,float> const &r2) { return sqrtf(r2); }
inline vec<float,float> rsqrt(vec<float,float> const &r2) { return 1.0f / sqrtf(r2); }
inline int testz(mmask<bool> const &a) { return a==0; }
inline vec<float,float> maskz_mov(mmask<bool> const &p,vec<float,float> const &a) {
    return p ? a : vec<float,float>(0.0f);
}
inline vec<float,float> maskz_mov(bool const &p,vec<float,float> const &a) {
    return p ? a : vec<float,float>(0.0f);
}
inline vec<float,float> mask_mov(vec<float,float> const &src,mmask<bool> const &p,vec<float,float> const &a) {
    return p ? a : src;
}
inline vec<float,float> mask_mov(vec<float,float> const &src,bool const &p,vec<float,float> const &a) {
    return p ? a : src;
}

/**********************************************************************\
* double precision
\**********************************************************************/

//template<> inline vec<double,double>::vec(const double &d) { ymm = d; }
template<> inline vec<double,double> &vec<double,double>::zero() { ymm = 0; return *this; }
template<> inline vec<double,double> &vec<double,double>::load1(double f) { ymm = f; return *this; }
template<> inline vec<double,double> &vec<double,double>::load(const double *pf) { ymm = *pf; return *this; }
template<> inline const vec<double,double> &vec<double,double>::store(double *pf) const { *pf = ymm; return *this; }

inline vec<double,double> min(vec<double,double> const &a,vec<double,double> const &b) { return a<b?a:b; }
inline vec<double,double> max(vec<double,double> const &a,vec<double,double> const &b) { return a>b?a:b; }
inline int movemask(mmask<bool> const &k) { return (int)(k); }

/**********************************************************************\
* 32-bit integer
\**********************************************************************/
//template<> inline vec<std::int32_t,std::int32_t>::vec(const std::int32_t &d) { ymm = d; }
template<> inline vec<std::int32_t,std::int32_t> &vec<std::int32_t,std::int32_t>::zero() { ymm = 0; return *this; }
template<> inline vec<std::int32_t,std::int32_t> &vec<std::int32_t,std::int32_t>::load1(std::int32_t f) { ymm = f; return *this; }
template<> inline vec<std::int32_t,std::int32_t> &vec<std::int32_t,std::int32_t>::load(const std::int32_t *pf) { ymm = *pf; return *this; }
template<> inline const vec<std::int32_t,std::int32_t> &vec<std::int32_t,std::int32_t>::store(std::int32_t *pf) const { *pf = ymm; return *this; }
inline vec<std::int32_t,std::int32_t> maskz_mov(mmask<bool> const &p,vec<std::int32_t,std::int32_t> const &a) {
    return p ? a : vec<std::int32_t,std::int32_t>(0);
}
inline vec<std::int32_t,std::int32_t> mask_mov(vec<std::int32_t,std::int32_t> const &src,mmask<bool> const &p,vec<std::int32_t,std::int32_t> const &a) {
    return p ? a : src;
}

#endif/*#elif defined(__SSE__)*/
#endif /*__AVX512F__*/

#if __GNUC__ > 4
    #pragma GCC diagnostic pop
#endif

#endif/*__CUDACC__*/

#endif/*SIMD_H*/
