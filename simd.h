#ifndef SIMD_H
#define SIMD_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_SIMD
#include <stdint.h>

#ifdef __AVX__
#define SIMD_BITS 3
#else
#define SIMD_BITS 2
#endif

#define SIMD_WIDTH (1<<SIMD_BITS)
#define SIMD_MASK (SIMD_WIDTH-1)

#ifdef HAVE_ANSIDECL_H
#include <ansidecl.h>
#else
#define ATTRIBUTE_ALIGNED_ALIGNOF(m)
#endif

#ifdef _MSC_VER
#define __SSE__
#define __SSE2__
#endif

#if defined(__SSE__)
#include <xmmintrin.h>
#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __SSE3__
#include <pmmintrin.h>
#endif
#ifdef __AVX__
#include <immintrin.h>
#endif

#elif defined(__ALTIVEC__)
#include <altivec.h>
#include <math.h> /* for sqrtf() */
#else
#error 'SIMD selected, but no known SIMD supported'
#endif

#if defined(HAVE_POSIX_MEMALIGN) || defined(HAVE_MEMALIGN)
#include <stdlib.h>
#endif

#if defined(__SSE__)
#ifdef __AVX__
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256 v4sf;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256 v4bool;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256i v4i;
#define MM_FCN(f) _mm256_##f
#define MM_CMP(f,F,a,b) _mm256_cmp_ps(a,b,_CMP_##F##_OQ)
#else
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128 v4sf;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128 v4bool;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128i v4i;
#define MM_FCN(f) _mm_##f
#define MM_CMP(f,F,a,b) _mm_cmp##f##_ps(a,b)
#endif
#else
typedef vector float v4sf;
typedef vector bool int v4bool;
typedef vector bool int v4i;
#endif

#if SIMD_WIDTH==4
#define SIMD_CONST(c) {c,c,c,c}
#elif SIMD_WIDTH==8
#define SIMD_CONST(c) {c,c,c,c,c,c,c,c}
#else
#error Invalid SIMD_WIDTH
#endif

typedef union {
    float f[SIMD_WIDTH];
    v4sf p;
    } v4;
typedef union {
    int32_t  i[SIMD_WIDTH];
    uint32_t u[SIMD_WIDTH];
    v4i      p;
    v4sf     pf;
    } i4;

#if defined(HAVE_POSIX_MEMALIGN)
static inline void * SIMD_malloc( size_t newSize ) {
    void *np;

    if ( posix_memalign( &np, sizeof(v4sf), newSize ) == 0 )
	return np;
    else
	return NULL;
    }

static inline void SIMD_free(void *p) {
    free(p);
    }
#else
#define SIMD_malloc(a) malloc(a)
#define SIMD_free(a) free(a)
#endif

static inline v4sf SIMD_SPLAT(float f) {
#ifdef __SSE__
    return MM_FCN(set1_ps)(f);
#else
    int i;
    v4 r;
    for(i=0; i<SIMD_WIDTH; i++) r.f[i] = f;
    return r.p;
    /*return (v4sf)(f,f,f,f);*/
#endif
    }

static inline v4i SIMD_SPLATI32(int i) {
#ifdef __SSE__
    return MM_FCN(set1_epi32)(i);
#else
    typedef union {
	int i[SIMD_WIDTH];
	v4i p;
	} i4;
    int j;
    i4 r;
    for(j=0; j<SIMD_WIDTH; j++) r.i[j] = i;
    return r.p;
    /*return (v4sf)(f,f,f,f);*/
#endif
    }

static inline v4sf SIMD_LOADS(float f) {
#ifdef __AVX__
    return _mm256_set_ps (0.0,0.0,0.0,0.0,0.0,0.0,0.0,f);
#elif defined(__SSE__)
    return MM_FCN(set_ss)(f);
#else
    int i;
    v4 r;
    r.f[0] = f;
    for(i=1; i<SIMD_WIDTH; i++) r.f[i] = 0;
    return r.p;
    /*return (v4sf)(f,0,0,0);*/
#endif
    }

#if defined(__SSE__)
#define SIMD_MUL(a,b) MM_FCN(mul_ps)(a,b)
#define SIMD_ADD(a,b) MM_FCN(add_ps)(a,b)
#define SIMD_SUB(a,b) MM_FCN(sub_ps)(a,b)
#ifdef __AVX__
#define SIMD_MADD(a,b,c) MM_FCN(fmadd_ps)(a,b,c)
#define SIMD_NMSUB(a,b,c) MM_FCN(fnmadd_ps)(a,b,c)
#else
#define SIMD_MADD(a,b,c) MM_FCN(add_ps)(MM_FCN(mul_ps)(a,b),c)
#define SIMD_NMSUB(a,b,c) MM_FCN(sub_ps)(c,MM_FCN(mul_ps)(a,b))
#endif
#define SIMD_DIV(a,b) MM_FCN(div_ps)(a,b)
#define SIMD_RSQRT(a) MM_FCN(rsqrt_ps)(a)
#define SIMD_RE(a) MM_FCN(rcp_ps)(a);
static inline v4sf SIMD_RE_EXACT(v4sf a) {
    static const v4 one = {SIMD_CONST(1.0)};
    v4sf r = SIMD_RE(a);
    return MM_FCN(add_ps)(MM_FCN(mul_ps)(r,MM_FCN(sub_ps)(one.p,MM_FCN(mul_ps)(r,a))),r);
    }
#define SIMD_MAX(a,b) MM_FCN(max_ps)(a,b)
#define SIMD_MIN(a,b) MM_FCN(min_ps)(a,b)
#define SIMD_CMP_EQ(a,b) MM_CMP(eq,EQ,a,b)
#define SIMD_CMP_NE(a,b) MM_CMP(ne,NE,a,b)
#define SIMD_CMP_LE(a,b) MM_CMP(le,LE,a,b)
#define SIMD_CMP_LT(a,b) MM_CMP(lt,LT,a,b)
#define SIMD_CMP_GE(a,b) MM_CMP(ge,GE,a,b)
#define SIMD_CMP_GT(a,b) MM_CMP(gt,GT,a,b)
#define SIMD_AND(a,b) MM_FCN(and_ps)(a,b)
#define SIMD_ANDNOT(a,b) MM_FCN(andnot_ps)(a,b)
#define SIMD_OR(a,b) MM_FCN(or_ps)(a,b)
#define SIMD_XOR(a,b) MM_FCN(xor_ps)(a,b)
#define SIMD_ALL_ZERO(a) MM_FCN(movemask_ps)(a)
#ifdef __AVX__
#define SIMD_I2F(a) MM_FCN(castsi256_ps)(a)
#define SIMD_F2I(a) MM_FCN(castps_si256)(a)

typedef union {
    ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256i p8;
    ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m128i p4[2];
    } avx_punner;
static inline v4i SIMD_CMP_EQ_EPI32(v4i a,v4i b) {
    register avx_punner x,y;
    x.p8 = a;
    y.p8 = b;
    x.p4[0] = _mm_cmpeq_epi32(x.p4[0],y.p4[0]);
    x.p4[1] = _mm_cmpeq_epi32(x.p4[1],y.p4[1]);
    return x.p8;
    }
static inline v4i SIMD_CMP_GT_EPI32(v4i a,v4i b) {
    typedef union {
	ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256i p8;
	ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m128i p4[2];
	} punner;
    register avx_punner x,y;
    x.p8 = a;
    y.p8 = b;
    x.p4[0] = _mm_cmpgt_epi32(x.p4[0],y.p4[0]);
    x.p4[1] = _mm_cmpgt_epi32(x.p4[1],y.p4[1]);
    return x.p8;
    }
#else
#define SIMD_I2F(a) MM_FCN(castsi128_ps)(a)
#define SIMD_F2I(a) MM_FCN(castps_si128)(a)
#define SIMD_CMP_EQ_EPI32(a,b) MM_FCN(cmpeq_epi32)(a,b)
#define SIMD_CMP_GT_EPI32(a,b) MM_FCN(cmpgt_epi32)(a,b)
#define SIMD_AND_EPI32(a,b) MM_FCN(and_si128)(a,b)
#define SIMD_ANDNOT_EPI32(a,b) MM_FCN(andnot_si128)(a,b)
#define SIMD_OR_EPI32(a,b) MM_FCN(or_si128)(a,b)
#endif
#else
static v4sf   simd_zero  = {SIMD_CONST(0)};
static v4bool simd_false = {SIMD_CONST(0)};
#define SIMD_MUL(a,b) vec_madd(a,b,simd_zero)
#define SIMD_ADD(a,b) vec_add(a,b)
#define SIMD_SUB(a,b) vec_sub(a,b)
#define SIMD_MADD(a,b,c) vec_madd(a,b,c)
#define SIMD_NMSUB(a,b,c) vec_nmsub(a,b,c)
#define SIMD_RSQRT(a) vec_rsqrte(a)
#define SIMD_RE(a) vec_re(a)
static inline v4sf SIMD_RE_EXACT( v4sf a) {
    static const v4sf one = {SIMD_CONST(1.0)};
    v4sf r = SIMD_RE(a);
    return vec_madd(r,vec_nmsub(r,a,one),r);
    }
#define SIMD_MAX(a,b) vec_max(a,b)
#define SIMD_MIN(a,b) vec_min(a,b)
#define SIMD_CMP_EQ(a,b) vec_cmpeq(a,b)
#define SIMD_CMP_NE(a,b) vec_cmpne(a,b)
#define SIMD_CMP_LE(a,b) vec_cmple(a,b)
#define SIMD_CMP_LT(a,b) vec_cmplt(a,b)
#define SIMD_CMP_GE(a,b) vec_cmpge(a,b)
#define SIMD_CMP_GT(a,b) vec_cmpgt(a,b)
#define SIMD_AND(a,b) vec_and(a,b)
#define SIMD_ANDNOT(a,b) vec_andc(b,a)
#define SIMD_OR(a,b) vec_or(a,b)
#define SIMD_XOR(a,b) vec_xor(a,b)
#define SIMD_ALL_ZERO(a) (!vec_all_eq(a,simd_false))
#define SIMD_F_TO_I(a) (a)
#define SIMD_CMP_EQ_EPI32(a,b) vec_cmpeq(a,b)
#define SIMD_CMP_GT_EPI32(a,b) vec_cmpgt(a,b)
#endif

static inline v4sf SIMD_RSQRT_EXACT(v4sf B) {
    static const v4 threehalves = {SIMD_CONST(1.5)};
    static const v4 half= {SIMD_CONST(0.5)};
    v4sf r, t1, t2;
    r = SIMD_RSQRT(B);          /* dir = 1/sqrt(d2) */
    t2 = SIMD_MUL(B,half.p);      /* 0.5*d2 */
    t1 = SIMD_MUL(r,r);         /* dir*dir */
    t1 = SIMD_NMSUB(t1,t2,threehalves.p);  /* 1.5 - 0.5*d2*dir*dir */
    r = SIMD_MUL(r,t1);     /* dir*(1.5 - 0.5*d2*dir*dir) */
    return r;
    }

#ifdef __SSE3__
static inline float SIMD_HADD(v4sf p) {
    register v4 r;
    r.p = p;
    r.p = MM_FCN(hadd_ps)(r.p,r.p);
    r.p = MM_FCN(hadd_ps)(r.p,r.p);
#ifdef __AVX__
    r.f[0] += r.f[4];
#endif
    return r.f[0];
    }
#else
static inline float SIMD_HADD(v4sf p) {
    v4 r;
    float f;
    int i;
    r.p = p;
    f = r.f[0];
    for(i=1; i<SIMD_WIDTH; i++) f += r.f[i];
    return f;
    }
#endif


#else
#define SIMD_malloc(a) malloc(a)
#define SIMD_free(a) free(a)
#endif

#endif
