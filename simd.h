#ifndef SIMD_H
#define SIMD_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif
#if !defined(_STDINT_H) && !defined(_STDINT_H_) && !defined(_STDINT_H_INCLUDED) && !defined(_STDINT) && !defined(__STDINT_H_)
#include <stdint.h>
#endif
#include <math.h>

/* Define if SIMD optimizations should be used. */
/* Visual Studio definition */
#if defined(_MSC_VER) && defined(__AVX__)
#define USE_SIMD 1
#ifndef __SSE__
#define __SSE__
#endif
#ifndef __SSE2__
#define __SSE2__
#endif
#ifdef __AVX2__
#define __FMA__
#endif
#elif defined(_M_X64)
#define __SSE__
#define __SSE2__
#define __AVX__
#ifdef __AVX2__
#define __FMA__
#endif
#define USE_SIMD 1
#elif _M_IX86_FP >= 1
#define USE_SIMD 1
#define __SSE__
#if _M_IX86_FP >= 2
#define __SSE2__
#endif
#endif

#ifdef USE_SIMD
#if defined(__AVX512F__)
#define SIMD_BITS 4
#elif defined(__AVX__)
#define SIMD_BITS 3
#elif defined(__SSE__) || defined(__ALTIVEC__)
#define SIMD_BITS 2
#endif
#endif
#ifdef SIMD_BITS
#define SIMD_DBITS (SIMD_BITS-1)
#else
#define SIMD_BITS 0
#define SIMD_DBITS 0
#endif
#define SIMD_WIDTH (1<<SIMD_BITS)
#define SIMD_DWIDTH (1<<SIMD_DBITS)
#define SIMD_MASK (SIMD_WIDTH-1)
#define SIMD_DMASK (SIMD_DWIDTH-1)

#if SIMD_WIDTH==1
#define SIMD_CONST(c) {c}
#define SIMD_DCONST(c) {c}
#elif SIMD_WIDTH==4
#define SIMD_CONST(c) {c,c,c,c}
#define SIMD_DCONST(c) {c,c}
#elif SIMD_WIDTH==8
#define SIMD_CONST(c) {c,c,c,c,c,c,c,c}
#define SIMD_DCONST(c) {c,c,c,c}
#elif SIMD_WIDTH==16
#define SIMD_CONST(c) {c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c}
#define SIMD_DCONST(c) {c,c,c,c,c,c,c,c}
#else
#error Invalid SIMD_WIDTH
#endif

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
#elif defined(__ALTIVEC__)/*not defined(__SSE__)*/
#include <altivec.h>
#include <math.h> /* for sqrtf() */
#endif/*__SSE__,__ALTIVEC__*/
#endif/*USE_SIMD*/

#ifdef USE_SIMD
#ifdef HAVE_ANSIDECL_H
#include <ansidecl.h>
#else/*!HAVE_ANSIDECL_H*/
#define ATTRIBUTE_ALIGNED_ALIGNOF(m)
#endif/*HAVE_ANSIDECL_H*/
#if defined(__SSE__)
#if defined(__AVX512F__)
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m512) __m512  v_sf;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m512) __m512d v_df;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m512) __m512  v_bool;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m512) __m512i v_i;
#define MM_FCN(f,p) _mm512_##f##_##p
#define MM_CMP(F,f,p,a,b) _mm512_cmp_##p(a,b,_CMP_##F##_OQ)
#elif defined(__AVX__)
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256  v_sf;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256d v_df;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256  v_bool;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256i v_i;
#define MM_FCN(f,p) _mm256_##f##_##p
#define MM_CMP(F,f,p,a,b) _mm256_cmp_##p(a,b,_CMP_##F##_OQ)
#else/*must be __SSE__*/
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128  v_sf;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128d v_df;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128  v_bool;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128i v_i;
#define MM_FCN(f,p) _mm_##f##_##p
#define MM_CMP(F,f,p,a,b) _mm_cmp##f##_##p(a,b)
#endif/*__AVX512F*/
#else/*__SSE__ must be altivec */
typedef vector float v_sf;
typedef vector bool int v_bool;
typedef vector bool int v_i;
#endif/*__SSE__*/
#else/*USE_SIMD*/
typedef float v_sf;
typedef double v_df;
typedef int32_t v_i;
#endif/*USE_SIMD*/

typedef union {
    float f[SIMD_WIDTH];
    v_sf p;
    float i[SIMD_WIDTH];
    } vfloat;

typedef union {
    double d[SIMD_DWIDTH];
    v_df p;
    uint64_t i[SIMD_DWIDTH];
    } vdouble;

typedef union {
    int32_t  i[SIMD_WIDTH];
    uint32_t u[SIMD_WIDTH];
    v_i      p;
    v_sf     pf;
    v_df     pd;
    } vint;

typedef union {
    uint64_t u[SIMD_WIDTH];
    v_i      pi;
    v_df     pd;
    } vint64;

#ifdef USE_SIMD

#if defined(HAVE_POSIX_MEMALIGN) || defined(HAVE_MEMALIGN)
#include <stdlib.h>
#endif
#ifdef HAVE_LIBMEMKIND
#include <hbwmalloc.h>
static inline void * SIMD_malloc( size_t newSize ) {
    void *np;

    if ( hbw_posix_memalign( &np, sizeof(v_sf), newSize ) == 0 )
	return np;
    else
	return NULL;
    }

static inline void SIMD_free(void *p) {
    hbw_free(p);
    }
#elif defined(HAVE_POSIX_MEMALIGN)
static inline void * SIMD_malloc( size_t newSize ) {
    void *np;

    if ( posix_memalign( &np, 2*sizeof(v_sf), newSize ) == 0 )
	return np;
    else
	return NULL;
    }

static inline void SIMD_free(void *p) {
    free(p);
    }
#else
static inline void * SIMD_malloc(size_t newSize) {
    return _mm_malloc(newSize, 2*sizeof(vfloat));
}

static inline void SIMD_free(void *p) {
    _mm_free(p);
}
#endif

#ifndef NO_C_MACROS

static inline v_sf SIMD_SPLAT(float f) {
#ifdef __SSE__
    return MM_FCN(set1,ps)(f);
#else
    int i;
    vfloat r;
    for(i=0; i<SIMD_WIDTH; i++) r.f[i] = f;
    return r.p;
    /*return (v_sf)(f,f,f,f);*/
#endif
    }

static inline v_i SIMD_SPLATI32(int i) {
#ifdef __SSE__
    return MM_FCN(set1,epi32)(i);
#else
    int j;
    vint r;
    for(j=0; j<SIMD_WIDTH; j++) r.i[j] = i;
    return r.p;
    /*return (v_sf)(f,f,f,f);*/
#endif
    }

static inline v_sf SIMD_LOADS(float f) {
#if defined(__AVX512F__)
    return _mm512_set_ps (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,f);
#elif defined(__AVX__)
    return _mm256_set_ps (0.0,0.0,0.0,0.0,0.0,0.0,0.0,f);
#elif defined(__SSE__)
    return MM_FCN(set,ss)(f);
#else
    int i;
    vfloat r;
    r.f[0] = f;
    for(i=1; i<SIMD_WIDTH; i++) r.f[i] = 0;
    return r.p;
    /*return (v_sf)(f,0,0,0);*/
#endif
    }

#if defined(__SSE__)
#define SIMD_MUL(a,b) MM_FCN(mul,ps)(a,b)
#define SIMD_ADD(a,b) MM_FCN(add,ps)(a,b)
#define SIMD_SUB(a,b) MM_FCN(sub,ps)(a,b)
#if defined(__FMA4__)
#define SIMD_MADD(a,b,c) MM_FCN(macc,ps)(a,b,c)    /*  (a*b) + c */
#define SIMD_MSUB(a,b,c) MM_FCN(msub,ps)(a,b,c)    /*  (a*b) - c */
#define SIMD_NMADD(a,b,c) MM_FCN(nmacc,ps)(a,b,c)  /* -(a*b) + c */
#define SIMD_NMSUB(a,b,c) MM_FCN(nmsub,ps)(a,b,c)  /* -(a*b) - c */
#elif defined(__FMA__)
#define SIMD_MADD(a,b,c) MM_FCN(fmadd,ps)(a,b,c)   /*  (a*b) + c */
#define SIMD_MSUB(a,b,c) MM_FCN(fmsub,ps)(a,b,c)   /*  (a*b) - c */
#define SIMD_NMADD(a,b,c) MM_FCN(fnmadd,ps)(a,b,c) /* -(a*b) + c */
#define SIMD_NMSUB(a,b,c) MM_FCN(fnmsub,ps)(a,b,c) /* -(a*b) - c */
#else
static const vfloat   simd_minus1  = {SIMD_CONST(-1.0)};
#define SIMD_MADD(a,b,c) MM_FCN(add,ps)(MM_FCN(mul,ps)(a,b),c)
#define SIMD_MSUB(a,b,c) MM_FCN(sub,ps)(MM_FCN(mul,ps)(a,b),c)
#define SIMD_NMADD(a,b,c) MM_FCN(sub,ps)(c,MM_FCN(mul,ps)(a,b))
#define SIMD_NMSUB(a,b,c) MM_FCN(mul,ps)(MM_FCN(add,ps)(c,MM_FCN(mul,ps)(a,b)),simd_minus1)
#endif/*defined(__FMA4__)*/
#define SIMD_DIV(a,b) MM_FCN(div,ps)(a,b)
#define SIMD_SQRT(a) MM_FCN(sqrt,ps)(a)
#ifdef __AVX512F__
#define SIMD_RSQRT(a) _mm512_rsqrt14_ps(a)
#define SIMD_RE(a) _mm512_rcp14_ps(a)
#else
#define SIMD_RSQRT(a) MM_FCN(rsqrt,ps)(a)
#define SIMD_RE(a) MM_FCN(rcp,ps)(a)
#endif
static inline v_sf SIMD_RE_EXACT(v_sf a) {
    static const vfloat one = {SIMD_CONST(1.0)};
    v_sf r = SIMD_RE(a);
    return MM_FCN(add,ps)(MM_FCN(mul,ps)(r,MM_FCN(sub,ps)(one.p,MM_FCN(mul,ps)(r,a))),r);
    }
#define SIMD_MAX(a,b) MM_FCN(max,ps)(a,b)
#define SIMD_MIN(a,b) MM_FCN(min,ps)(a,b)
#define SIMD_CMP_EQ(a,b) MM_CMP(EQ,eq,ps,a,b)
#define SIMD_CMP_NE(a,b) MM_CMP(NE,ne,ps,a,b)
#define SIMD_CMP_LE(a,b) MM_CMP(LE,le,ps,a,b)
#define SIMD_CMP_LT(a,b) MM_CMP(LT,lt,ps,a,b)
#define SIMD_CMP_GE(a,b) MM_CMP(GE,ge,ps,a,b)
#define SIMD_CMP_GT(a,b) MM_CMP(GT,gt,ps,a,b)
#define SIMD_AND(a,b) MM_FCN(and,ps)(a,b)
#define SIMD_ANDNOT(a,b) MM_FCN(andnot,ps)(a,b)
#define SIMD_OR(a,b) MM_FCN(or,ps)(a,b)
#define SIMD_XOR(a,b) MM_FCN(xor,ps)(a,b)
#define SIMD_ALL_ZERO(a) MM_FCN(movemask,ps)(a)
#ifdef __AVX__
#define SIMD_I2F(a) MM_FCN(castsi256,ps)(a)
#define SIMD_F2I(a) MM_FCN(castps,si256)(a)
#ifdef __AVX2__
#define SIMD_CMP_EQ_EPI32(a,b) MM_FCN(cmpeq,epi32)(a,b)
#define SIMD_CMP_GT_EPI32(a,b) MM_FCN(cmpgt,epi32)(a,b)
#else/*__AVX2__*/
typedef union {
    ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256i p8;
    ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m128i p4[2];
    } avx_punner;
static inline v_i SIMD_CMP_EQ_EPI32(v_i a,v_i b) {
    avx_punner x,y;
    x.p8 = a;
    y.p8 = b;
    x.p4[0] = _mm_cmpeq_epi32(x.p4[0],y.p4[0]);
    x.p4[1] = _mm_cmpeq_epi32(x.p4[1],y.p4[1]);
    return x.p8;
    }
static inline v_i SIMD_CMP_GT_EPI32(v_i a,v_i b) {
    typedef union {
	ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256i p8;
	ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m128i p4[2];
	} punner;
    avx_punner x,y;
    x.p8 = a;
    y.p8 = b;
    x.p4[0] = _mm_cmpgt_epi32(x.p4[0],y.p4[0]);
    x.p4[1] = _mm_cmpgt_epi32(x.p4[1],y.p4[1]);
    return x.p8;
    }
#endif
#else/*__AVX__*/
#define SIMD_I2F(a) MM_FCN(castsi128,ps)(a)
#define SIMD_F2I(a) MM_FCN(castps,si128)(a)
#define SIMD_CMP_EQ_EPI32(a,b) MM_FCN(cmpeq,epi32)(a,b)
#define SIMD_CMP_GT_EPI32(a,b) MM_FCN(cmpgt,epi32)(a,b)
#define SIMD_AND_EPI32(a,b) MM_FCN(and,si128)(a,b)
#define SIMD_ANDNOT_EPI32(a,b) MM_FCN(andnot,si128)(a,b)
#define SIMD_OR_EPI32(a,b) MM_FCN(or,si128)(a,b)
#endif
#else/*__SSE__*/
static const v_sf   simd_zero    = {SIMD_CONST(0)};
static const v_sf   simd_minus1  = {SIMD_CONST(-1)};
static const v_bool simd_false   = {SIMD_CONST(0)};
#define SIMD_MUL(a,b) vec_madd(a,b,simd_zero)
#define SIMD_ADD(a,b) vec_add(a,b)
#define SIMD_SUB(a,b) vec_sub(a,b)
#define SIMD_MADD(a,b,c) vec_madd(a,b,c)                                  /*  (a*b) + c */
#define SIMD_MSUB(a,b,c) vec_madd(vec_nmsub(a,b,c),simd_minus1,simd_zero) /*  (a*b) - c */
#define SIMD_NMADD(a,b,c) vec_nmsub(a,b,c)                                /* -(a*b) + c */
#define SIMD_NMSUB(a,b,c) vec_madd(vec_madd(a,b,c),simd_minus1,simd_zero) /* -(a*b) - c */
#define SIMD_RSQRT(a) vec_rsqrte(a)
#define SIMD_RE(a) vec_re(a)
static inline v_sf SIMD_RE_EXACT( v_sf a) {
    static const v_sf one = {SIMD_CONST(1.0)};
    v_sf r = SIMD_RE(a);
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
#endif/*__SSE__*/

static inline v_sf SIMD_RSQRT_EXACT(v_sf B) {
    static const vfloat threehalves = {SIMD_CONST(1.5)};
    static const vfloat half= {SIMD_CONST(0.5)};
    v_sf r, t1, t2;
    r = SIMD_RSQRT(B);          /* dir = 1/sqrt(d2) */
    t2 = SIMD_MUL(B,half.p);      /* 0.5*d2 */
    t1 = SIMD_MUL(r,r);         /* dir*dir */
    t1 = SIMD_NMADD(t1,t2,threehalves.p);  /* 1.5 - 0.5*d2*dir*dir */
    r = SIMD_MUL(r,t1);     /* dir*(1.5 - 0.5*d2*dir*dir) */
    return r;
    }

#ifndef __AVX512F__
#ifdef __SSE__
/**
 * latencies and throughputs:
 * _extractf128_ps : latency: 3, throughput: 1
 * _castps256_ps128: no latency or throughput
 * all other instructions have latency and throughput of 1
 */
static inline float SIMD_HADD(v_sf p) {
    __m128 sum;
    #ifdef __AVX__
    __m128 upper = _mm256_extractf128_ps(p, 1);
    __m128 lower = _mm256_castps256_ps128(p);
    sum = _mm_add_ps(lower,upper);
    #else
    sum = p;
    #endif

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
#else/*__SSE__*/
static inline float SIMD_HADD(v_sf p) {
    vfloat r;
    float f;
    int i;
    r.p = p;
    f = r.f[0];
    for(i=1; i<SIMD_WIDTH; i++) f += r.f[i];
    return f;
    }
#endif/*__SSE__*/
#endif/*__AVX512F__*/

/* With SSE2 and beyond we have double support */
#ifdef __SSE2__

static inline v_df SIMD_DSPLAT(double f) {
    return MM_FCN(set1,pd)(f);
    }
static inline v_df SIMD_DLOADS(double f) {
#if defined(__AVX512F__)
    return _mm512_set_pd (0.0,0.0,0.0,0.0,0.0,0.0,0.0,f);
#elif defined(__AVX__)
    return _mm256_set_pd (0.0,0.0,0.0,f);
#else
    return MM_FCN(set,sd)(f);
#endif
    }

#define SIMD_DMUL(a,b) MM_FCN(mul,pd)(a,b)
#define SIMD_DADD(a,b) MM_FCN(add,pd)(a,b)
#define SIMD_DSUB(a,b) MM_FCN(sub,pd)(a,b)
#if defined(__FMA4__)
#define SIMD_DMADD(a,b,c) MM_FCN(macc,pd)(a,b,c)    /*  (a*b) + c */
#define SIMD_DMSUB(a,b,c) MM_FCN(msub,pd)(a,b,c)    /*  (a*b) - c */
#define SIMD_DNMADD(a,b,c) MM_FCN(nmacc,pd)(a,b,c)  /* -(a*b) + c */
#define SIMD_DNMSUB(a,b,c) MM_FCN(nmsub,pd)(a,b,c)  /* -(a*b) - c */
#elif defined(__FMA__)
#define SIMD_DMADD(a,b,c) MM_FCN(fmadd,pd)(a,b,c)   /*  (a*b) + c */
#define SIMD_DMSUB(a,b,c) MM_FCN(fmsub,pd)(a,b,c)   /*  (a*b) - c */
#define SIMD_DNMADD(a,b,c) MM_FCN(fnmadd,pd)(a,b,c) /* -(a*b) + c */
#define SIMD_DNMSUB(a,b,c) MM_FCN(fnmsub,pd)(a,b,c) /* -(a*b) - c */
#else
#define SIMD_DMADD(a,b,c) MM_FCN(add,pd)(MM_FCN(mul,pd)(a,b),c)
#define SIMD_DNMADD(a,b,c) MM_FCN(sub,pd)(c,MM_FCN(mul,pd)(a,b))
#endif
#define SIMD_DDIV(a,b) MM_FCN(div,pd)(a,b)
#define SIMD_DSQRT(a) MM_FCN(sqrt,pd)(a)
#define SIMD_DMAX(a,b) MM_FCN(max,pd)(a,b)
#define SIMD_DMIN(a,b) MM_FCN(min,pd)(a,b)
#define SIMD_DCMP_EQ(a,b) MM_CMP(EQ,eq,pd,a,b)
#define SIMD_DCMP_NE(a,b) MM_CMP(NE,ne,pd,a,b)
#define SIMD_DCMP_LE(a,b) MM_CMP(LE,le,pd,a,b)
#define SIMD_DCMP_LT(a,b) MM_CMP(LT,lt,pd,a,b)
#define SIMD_DCMP_GE(a,b) MM_CMP(GE,ge,pd,a,b)
#define SIMD_DCMP_GT(a,b) MM_CMP(GT,gt,pd,a,b)
#define SIMD_DAND(a,b) MM_FCN(and,pd)(a,b)
#define SIMD_DANDNOT(a,b) MM_FCN(andnot,pd)(a,b)
#define SIMD_DOR(a,b) MM_FCN(or,pd)(a,b)
#define SIMD_DXOR(a,b) MM_FCN(xor,pd)(a,b)
#define SIMD_DKXOR(a,b) MM_FCN(xor,pd)(a,b)
#define SIMD_DALL_ZERO(a) MM_FCN(movemask,pd)(a)

/* p==false then select a -- p==true, select b */
#ifdef __AVX__
#define SIMD_SELECT(a,b,p) MM_FCN(blendv,ps)(a,b,p)
#define SIMD_DSELECT(a,b,p) MM_FCN(blendv,pd)(a,b,p)
#else
#define SIMD_SELECT(a,b,p) SIMD_OR(SIMD_AND(p,b),SIMD_ANDNOT(p,a))
#define SIMD_DSELECT(a,b,p) SIMD_DOR(SIMD_DAND(p,b),SIMD_DANDNOT(p,a))
#endif

#if defined(__AVX512F__)
#define SIMD_I2D(a) MM_FCN(castsi512,pd)(a)
#define SIMD_D2I(a) MM_FCN(castpd,si512)(a)
#define SIMD_D2F(a) MM_FCN(castps256,ps512)(MM_FCN(cvtpd,ps)(a))
#elif defined(__AVX__)
#define SIMD_I2D(a) MM_FCN(castsi256,pd)(a)
#define SIMD_D2I(a) MM_FCN(castpd,si256)(a)
#define SIMD_D2F(a) MM_FCN(castps128,ps256)(MM_FCN(cvtpd,ps)(a))
#else
#define SIMD_I2D(a) MM_FCN(castsi128,pd)(a)
#define SIMD_D2I(a) MM_FCN(castpd,si128)(a)
#define SIMD_D2F(a) MM_FCN(cvtpd,ps)(a)
#endif

#if !defined(__AVX512F__)
static inline v_df SIMD_DRSQRT(v_df B) {
    static const vdouble half = {SIMD_DCONST(0.5)};
    static const vdouble three = {SIMD_DCONST(3.0)};
    v_df r;
    r = MM_FCN(cvtps,pd)(_mm_rsqrt_ps(MM_FCN(cvtpd,ps)(B))); 
    /* vvv ** Up to two more interations for full precision are needed ** vvv*/
    r = SIMD_DMUL(SIMD_DMUL(half.p,r),SIMD_DSUB(three.p,SIMD_DMUL(SIMD_DMUL(B,r),r)));
    return r;
    }
static inline v_df SIMD_DRSQRT_EXACT(v_df B) {
    static const vdouble half = {SIMD_DCONST(0.5)};
    static const vdouble three = {SIMD_DCONST(3.0)};
    v_df r = SIMD_DRSQRT(B);
    r = SIMD_DMUL(SIMD_DMUL(half.p,r),SIMD_DSUB(three.p,SIMD_DMUL(SIMD_DMUL(B,r),r)));
    /*r = SIMD_DMUL(SIMD_DMUL(half.p,r),SIMD_DSUB(three.p,SIMD_DMUL(SIMD_DMUL(B,r),r)));*/
    return r;
    }
static inline v_df SIMD_DRE(v_df a) {
    v_df r = MM_FCN(cvtps,pd)(_mm_rcp_ps(MM_FCN(cvtpd,ps)(a))); 
    /* vvv ** Up to two more interations for full precision are needed ** vvv*/
    r = SIMD_DSUB(SIMD_DADD(r,r),SIMD_DMUL(SIMD_DMUL(a,r),r));
    return r;
    }
static inline v_df SIMD_DRE_GOOD(v_df a) {
    v_df r = SIMD_DRE(a);
    r = SIMD_DSUB(SIMD_DADD(r,r),SIMD_DMUL(SIMD_DMUL(a,r),r));
    return r;
    }
static inline v_df SIMD_DRE_EXACT(v_df a) {
    v_df r = SIMD_DRE_GOOD(a);
    r = SIMD_DSUB(SIMD_DADD(r,r),SIMD_DMUL(SIMD_DMUL(a,r),r));
    return r;
    }
#endif/*__AVX512F__*/

#endif/*__SSE2*/
#endif/*NO_C_MACROS*/
#else/*USE_SIMD*/
#define SIMD_malloc malloc
#define SIMD_free free
#endif/*USE_SIMD*/

#if defined(__cplusplus)

#if __GNUC__ > 5
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

/**********************************************************************\
* SIMD Vector class template
\**********************************************************************/
template<typename vtype,typename ftype>
struct vec {
    vtype ymm;
public:
    vec() {}
#ifdef USE_SIMD
    vec(ftype const &d);
#endif
    vec(vtype const &d) { ymm = d; }
    operator vtype() const { return ymm; }
    ftype operator [] (uint32_t idx) const;
    vec & zero();
    vec & load1(ftype f);
    vec & load(ftype *f);
    const vec & store(ftype *f) const;
    static int width() { return sizeof(vtype)/sizeof(ftype); }
    static int mask()  { return sizeof(vtype)/sizeof(ftype)-1; }
    static const vec sign_mask();
    vec<vtype,ftype> operator-() const;
    };

template<typename vtype>
class mmask {
    vtype ymm;
public:
    mmask() {}
    mmask(vtype const &d) { ymm = d; }
    operator vtype() const { return ymm; }
    };

#if defined(__AVX512F__) && defined(USE_SIMD)
typedef vec<__m512i,int32_t> i32v;
typedef vec<__m512,float> fvec;
typedef mmask<__mmask16> fmask;
typedef vec<__m512i,int64_t> i64v;
typedef vec<__m512d,double> dvec;
typedef mmask<__mmask8> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v(_mm512_cvtps_epi32(a)); }
//inline i64v cvt_i64v(const fvec &a) { return i64v(_mm512_cvtps_epi64(a)); }
inline fvec cvt_fvec(const i32v &a) { return fvec(_mm512_cvtepi32_ps(a)); }
inline fvec cvt_fvec(const dvec &a) { return fvec(_mm512_castps256_ps512(_mm512_cvtpd_ps(a))); }
//#define SIMD_D2F(a) MM_FCN(castps256,ps512)(MM_FCN(cvtpd,ps)(a))

//inline fvec cvt_fvec(const i64v &a) { return fvec(_mm512_cvtepi64_ps(a)); }
#elif defined(__AVX__) && defined(USE_SIMD)
typedef vec<__m256i,int32_t> i32v;
typedef vec<__m256,float> fvec;
typedef vec<__m256,float> fmask;
typedef vec<__m256i,int64_t> i64v;
typedef vec<__m256d,double> dvec;
typedef vec<__m256d,double> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v(_mm256_cvtps_epi32(a)); }
inline fvec cvt_fvec(const i32v &a) { return fvec(_mm256_cvtepi32_ps(a)); }
inline fvec cvt_fvec(const dvec &a) { return fvec(_mm256_castps128_ps256(_mm256_cvtpd_ps(a))); }
#elif defined(__SSE__) && defined(USE_SIMD)
typedef vec<__m128i,int32_t> i32v;
typedef vec<__m128,float> fvec;
typedef vec<__m128,float> fmask;
#if defined(__SSE2__)
typedef vec<__m128i,int64_t> i64v;
typedef vec<__m128d,double> dvec;
typedef vec<__m128d,double> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v(_mm_cvtps_epi32(a)); }
inline fvec cvt_fvec(const i32v &a) { return fvec(_mm_cvtepi32_ps(a)); }
inline fvec cvt_fvec(const dvec &a) { return fvec(_mm_cvtpd_ps(a)); }
#endif/*__SSE2__*/
#else/*__AVX512F__,__AVX__,__SSE2__*/
typedef vec<int32_t,int32_t> i32v;
typedef vec<float,float> fvec;
typedef vec<double,double> dvec;
typedef mmask<bool> fmask;
typedef mmask<bool> dmask;
inline i32v cvt_i32v(const fvec &a) { return i32v((int32_t)a); }
inline fvec cvt_fvec(const i32v &a) { return fvec((float)a); }
inline fvec cvt_fvec(const dvec &a) { return fvec((float)a); }
#endif/*__AVX512F__,__AVX__,__SSE2__*/

/**********************************************************************\
* Generic Operators
\**********************************************************************/

template<typename v,typename ftype> inline vec<v,ftype> abs(vec<v,ftype> const &a) { return a & ~vec<v,ftype>::sign_mask(); }
template<typename v,typename ftype> inline vec<v,ftype> & operator+=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a + b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator+=(vec<v,ftype> &a,ftype const &b) { return a = a + b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator-=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a - b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator-=(vec<v,ftype> &a,ftype const &b) { return a = a - b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator*=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a * b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator*=(vec<v,ftype> &a,ftype const &b) { return a = a * b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator/=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a / b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator/=(vec<v,ftype> &a,ftype const &b) { return a = a / b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator&=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a & b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator&=(vec<v,ftype> &a,ftype const &b) { return a = a & b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator|=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a | b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator|=(vec<v,ftype> &a,ftype const &b) { return a = a | b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator^=(vec<v,ftype> &a,vec<v,ftype> const &b) { return a = a ^ b; }
template<typename v,typename ftype> inline vec<v,ftype> & operator^=(vec<v,ftype> &a,ftype const &b) { return a = a ^ b; }
template<typename v,typename ftype> inline ftype vec<v,ftype>::operator [] (uint32_t idx) const {
//	ftype d[width()];
	ftype d[SIMD_WIDTH];
	store(d);
    return d[idx];
    }
template<typename v> inline mmask<v> & operator&=(mmask<v> &a,mmask<v> const &b) { return a = a & b; }
template<typename v> inline mmask<v> & operator|=(mmask<v> &a,mmask<v> const &b) { return a = a | b; }
template<typename v> inline mmask<v> & operator^=(mmask<v> &a,mmask<v> const &b) { return a = a ^ b; }


#if defined(__AVX512F__) && defined(USE_SIMD)

/**********************************************************************\
* AVX512 single precision
\**********************************************************************/

template<> inline vec<__m512,float>::vec(const float &d) { ymm = _mm512_set1_ps(d); }
template<> inline vec<__m512,float> & vec<__m512,float>::zero() { ymm = _mm512_setzero_ps(); return *this; }
template<> inline vec<__m512,float> & vec<__m512,float>::load1(float f) { ymm = _mm512_setr_ps(f,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512,float> & vec<__m512,float>::load(float *pf) { ymm = _mm512_loadu_ps(pf); return *this; }
template<> inline const vec<__m512,float> & vec<__m512,float>::store(float *pf) const { _mm512_storeu_ps(pf,ymm); return *this; }
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
    __m256 low  = _mm512_castps512_ps256(a);
    __m256 high = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(a),1));
    return hadd(low) + hadd(high);
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
template<> inline vec<__m512i,int32_t>::vec(const int32_t &d) { ymm = _mm512_set1_epi32(d); }
template<> inline vec<__m512i,int32_t> & vec<__m512i,int32_t>::zero() { ymm = _mm512_setzero_si512(); return *this; }
template<> inline vec<__m512i,int32_t> & vec<__m512i,int32_t>::load1(int32_t f) { ymm = _mm512_setr_epi32(f,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512i,int32_t> & vec<__m512i,int32_t>::load(int32_t *pf) { ymm = _mm512_loadu_si512((__m512i*)pf); return *this; }
template<> inline const vec<__m512i,int32_t> & vec<__m512i,int32_t>::store(int32_t *pf) const { _mm512_storeu_si512((__m512i*)pf,ymm); return *this; }
template<> inline const vec<__m512i,int32_t> vec<__m512i,int32_t>::sign_mask() { return _mm512_set1_epi32(0x80000000); }
template<> inline vec<__m512i,int32_t> vec<__m512i,int32_t>::operator-() const {
    return _mm512_castps_si512(_mm512_xor_ps(_mm512_castsi512_ps(ymm),_mm512_castsi512_ps(sign_mask()))); }
inline vec<__m512i,int32_t> operator&(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_and_epi32(a,b); }
inline vec<__m512i,int32_t> operator|(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_or_epi32(a,b); }
inline vec<__m512i,int32_t> operator^(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_xor_epi32(a,b); }

inline mmask<__mmask16> operator==(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_EQ); }
inline mmask<__mmask16> operator!=(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_NE); }
inline mmask<__mmask16> operator>(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GT); }
inline mmask<__mmask16> operator<(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LT); }
inline mmask<__mmask16> operator>=(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_GE); }
inline mmask<__mmask16> operator<=(vec<__m512i,int32_t> const &a,vec<__m512i,int32_t> const &b) { return _mm512_cmp_epi32_mask(a,b,_MM_CMPINT_LE); }


/**********************************************************************\
* AVX512 64-bit integer
\**********************************************************************/
template<> inline vec<__m512i,int64_t>::vec(const int64_t &d) { ymm = _mm512_set1_epi64(d); }
template<> inline vec<__m512i,int64_t> & vec<__m512i,int64_t>::zero() { ymm = _mm512_setzero_si512(); return *this; }
template<> inline vec<__m512i,int64_t> & vec<__m512i,int64_t>::load1(int64_t f) { ymm = _mm512_setr_epi64(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512i,int64_t> & vec<__m512i,int64_t>::load(int64_t *pf) { ymm = _mm512_loadu_si512((__m512i*)pf); return *this; }
template<> inline const vec<__m512i,int64_t> & vec<__m512i,int64_t>::store(int64_t *pf) const { _mm512_storeu_si512((__m512i*)pf,ymm); return *this; }
template<> inline const vec<__m512i,int64_t> vec<__m512i,int64_t>::sign_mask() { return _mm512_set1_epi64(0x8000000000000000); }
template<> inline vec<__m512i,int64_t> vec<__m512i,int64_t>::operator-() const {
    return _mm512_castps_si512(_mm512_xor_ps(_mm512_castsi512_ps(ymm),_mm512_castsi512_ps(sign_mask()))); }
inline vec<__m512i,int64_t> operator&(vec<__m512i,int64_t> const &a,vec<__m512i,int64_t> const &b) { return _mm512_and_epi64(a,b); }
inline vec<__m512i,int64_t> operator|(vec<__m512i,int64_t> const &a,vec<__m512i,int64_t> const &b) { return _mm512_or_epi64(a,b); }
inline vec<__m512i,int64_t> operator^(vec<__m512i,int64_t> const &a,vec<__m512i,int64_t> const &b) { return _mm512_xor_epi64(a,b); }

/**********************************************************************\
* AVX512 single precision mask
\**********************************************************************/

inline mmask<__mmask16> operator&(mmask<__mmask16> const &a,mmask<__mmask16> const &b) { return _mm512_kand(a,b); }
inline mmask<__mmask16> operator|(mmask<__mmask16> const &a,mmask<__mmask16> const &b) { return _mm512_kor(a,b); }
inline mmask<__mmask16> operator^(mmask<__mmask16> const &a,mmask<__mmask16> const &b) { return _mm512_kxor(a,b); }
inline int testz(mmask<__mmask16> const &a) { return _mm512_kortestz(a,a); }

/**********************************************************************\
* AVX512 double precision
\**********************************************************************/

template<> inline vec<__m512d,double>::vec(const double &d) { ymm = _mm512_set1_pd(d); }
template<> inline vec<__m512d,double> & vec<__m512d,double>::zero() { ymm = _mm512_setzero_pd(); return *this; }
template<> inline vec<__m512d,double> & vec<__m512d,double>::load1(double f) { ymm = _mm512_setr_pd(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m512d,double> & vec<__m512d,double>::load(double *pf) { ymm = _mm512_loadu_pd(pf); return *this; }
template<> inline const vec<__m512d,double> & vec<__m512d,double>::store(double *pf) const { _mm512_storeu_pd(pf,ymm); return *this; }
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
template<> inline vec<__m256,float> & vec<__m256,float>::zero() { ymm = _mm256_setzero_ps(); return *this; }
template<> inline vec<__m256,float> & vec<__m256,float>::load1(float f) { ymm = _mm256_setr_ps(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m256,float> & vec<__m256,float>::load(float *pf) { ymm = _mm256_loadu_ps(pf); return *this; }
template<> inline const vec<__m256,float> & vec<__m256,float>::store(float *pf) const { _mm256_storeu_ps(pf,ymm); return *this; }
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
inline vec<__m256,float> blend(vec<__m256,float> const &a,vec<__m256,float> const &b,vec<__m256,float> const &p) { return _mm256_blendv_ps(a,b,p); }
inline vec<__m256,float> mask_mov(vec<__m256,float> const &src,vec<__m256,float> const &p,vec<__m256,float> const &a)
    { return _mm256_blendv_ps(src,a,p); }
inline vec<__m256,float> maskz_mov(vec<__m256,float> const &p,vec<__m256,float> const &a) { return a & p; }
inline int testz(vec<__m256,float> const &a) { return !_mm256_movemask_ps(a); }

/**********************************************************************\
* AVX 32-bit integer
\**********************************************************************/
template<> inline vec<__m256i,int32_t>::vec(const int32_t &d) { ymm = _mm256_set1_epi32(d); }
template<> inline vec<__m256i,int32_t> & vec<__m256i,int32_t>::zero() { ymm = _mm256_setzero_si256(); return *this; }
template<> inline vec<__m256i,int32_t> & vec<__m256i,int32_t>::load1(int32_t f) { ymm = _mm256_setr_epi32(f,0,0,0,0,0,0,0); return *this; }
template<> inline vec<__m256i,int32_t> & vec<__m256i,int32_t>::load(int32_t *pf) { ymm = _mm256_loadu_si256((__m256i*)pf); return *this; }
template<> inline const vec<__m256i,int32_t> & vec<__m256i,int32_t>::store(int32_t *pf) const { _mm256_storeu_si256((__m256i*)pf,ymm); return *this; }
template<> inline const vec<__m256i,int32_t> vec<__m256i,int32_t>::sign_mask() { return _mm256_set1_epi32(0x80000000); }
#ifdef __AVX2__
template<> inline vec<__m256i,int32_t> vec<__m256i,int32_t>::operator-() const { return _mm256_xor_si256(ymm,sign_mask()); }
inline vec<__m256i,int32_t> operator+(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_add_epi32(a,b); }
inline vec<__m256i,int32_t> operator-(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_sub_epi32(a,b); }
inline vec<__m256i,int32_t> operator&(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_and_si256(a,b); }
#else
template<> inline vec<__m256i,int32_t> vec<__m256i,int32_t>::operator-() const {
    return _mm256_castps_si256(_mm256_xor_ps(_mm256_castsi256_ps(ymm),_mm256_castsi256_ps(sign_mask()))); }
inline vec<__m256i,int32_t> operator&(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b)
    { return _mm256_castps_si256(_mm256_and_ps(_mm256_castsi256_ps(a),_mm256_castsi256_ps(b))); }
#endif

#ifdef __AVX2__
inline vec<__m256i,int32_t> operator==(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_cmpeq_epi32(a,b); }
//inline vec<__m256i,int32_t> operator!=(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_cmpneq_epi32(a,b); }
inline vec<__m256i,int32_t> operator>(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_cmpgt_epi32(a,b); }
//inline vec<__m256i,int32_t> operator<(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_cmplt_epi32(a,b); }
//inline vec<__m256i,int32_t> operator>=(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_cmpge_epi32(a,b); }
//inline vec<__m256i,int32_t> operator<=(vec<__m256i,int32_t> const &a,vec<__m256i,int32_t> const &b) { return _mm256_cmple_epi32(a,b); }
#else
#endif


/**********************************************************************\
* AVX double precision
\**********************************************************************/

template<> inline vec<__m256d,double>::vec(const double &d) { ymm = _mm256_set1_pd(d); }
template<> inline vec<__m256d,double> & vec<__m256d,double>::zero() { ymm = _mm256_setzero_pd(); return *this; }
template<> inline vec<__m256d,double> & vec<__m256d,double>::load1(double f) { ymm = _mm256_setr_pd(f,0,0,0); return *this; }
template<> inline vec<__m256d,double> & vec<__m256d,double>::load(double *pf) { ymm = _mm256_loadu_pd(pf); return *this; }
template<> inline const vec<__m256d,double> & vec<__m256d,double>::store(double *pf) const { _mm256_storeu_pd(pf,ymm); return *this; }
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

inline vec<__m256d,double> blend(vec<__m256d,double> const &a,vec<__m256d,double> const &b,vec<__m256d,double> const &p) { return _mm256_blendv_pd(a,b,p); }
inline int movemask(vec<__m256d,double> const &r2) { return _mm256_movemask_pd(r2); }

inline vec<__m256d,double> mask_mov(vec<__m256d,double> const &src,vec<__m256d,double> const &k,vec<__m256d,double> const &a)
    { return _mm256_blendv_pd(src,a,k); }
inline int testz(vec<__m256d,double> const &a) { return !_mm256_movemask_pd(a); }

#elif defined(__SSE__) && defined(USE_SIMD)

/**********************************************************************\
* SSE single precision
\**********************************************************************/

template<> inline vec<__m128,float>::vec(const float &d) { ymm = _mm_set1_ps(d); }
template<> inline vec<__m128,float> & vec<__m128,float>::zero() { ymm = _mm_setzero_ps(); return *this; }
template<> inline vec<__m128,float> & vec<__m128,float>::load1(float f) { ymm = _mm_setr_ps(f,0,0,0); return *this; }
template<> inline vec<__m128,float> & vec<__m128,float>::load(float *pf) { ymm = _mm_loadu_ps(pf); return *this; }
template<> inline const vec<__m128,float> & vec<__m128,float>::store(float *pf) const { _mm_storeu_ps(pf,ymm); return *this; }
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
inline vec<__m128,float> blend(vec<__m128,float> const &a,vec<__m128,float> const &b,vec<__m128,float> const &p) {
#ifdef __SSE4_1__
    return _mm_blendv_ps(a,b,p);
#else
    return _mm_or_ps(_mm_and_ps(p,b),_mm_andnot_ps(p,a));
#endif
    }
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
template<> inline vec<__m128i,int32_t>::vec(const int32_t &d) { ymm = _mm_set1_epi32(d); }
template<> inline vec<__m128i,int32_t> & vec<__m128i,int32_t>::zero() { ymm = _mm_setzero_si128(); return *this; }
template<> inline vec<__m128i,int32_t> & vec<__m128i,int32_t>::load1(int32_t f) { ymm = _mm_setr_epi32(f,0,0,0); return *this; }
template<> inline vec<__m128i,int32_t> & vec<__m128i,int32_t>::load(int32_t *pf) { ymm = _mm_loadu_si128((__m128i*)pf); return *this; }
template<> inline const vec<__m128i,int32_t> & vec<__m128i,int32_t>::store(int32_t *pf) const { _mm_storeu_si128((__m128i*)pf,ymm); return *this; }
template<> inline const vec<__m128i,int32_t> vec<__m128i,int32_t>::sign_mask() { return _mm_set1_epi32(0x80000000); }
#ifdef __SSE2__
template<> inline vec<__m128i,int32_t> vec<__m128i,int32_t>::operator-() const { return _mm_xor_si128(ymm,sign_mask()); }
inline vec<__m128i,int32_t> operator&(vec<__m128i,int32_t> const &a,vec<__m128i,int32_t> const &b) { return _mm_and_si128(a,b); }
#else
template<> inline vec<__m128i,int32_t> vec<__m128i,int32_t>::operator-() {
    return _mm_castps_si128(_mm_xor_ps(_mm_castsi128_ps(ymm),_mm_castsi128_ps(sign_mask()))); }
inline vec<__m128i,int32_t> operator&(vec<__m128i,int32_t> const &a,vec<__m128i,int32_t> const &b)
    { return _mm_castps_si128(_mm_and_ps(_mm_castsi128_ps(a),_mm_castsi128_ps(b))); }
#endif

#if defined(__SSE2__)
/**********************************************************************\
* SSE double precision
\**********************************************************************/

template<> inline vec<__m128d,double>::vec(const double &d) { ymm = _mm_set1_pd(d); }
template<> inline vec<__m128d,double> & vec<__m128d,double>::zero() { ymm = _mm_setzero_pd(); return *this; }
template<> inline vec<__m128d,double> & vec<__m128d,double>::load1(double f) { ymm = _mm_setr_pd(f,0); return *this; }
template<> inline vec<__m128d,double> & vec<__m128d,double>::load(double *pf) { ymm = _mm_loadu_pd(pf); return *this; }
template<> inline const vec<__m128d,double> & vec<__m128d,double>::store(double *pf) const { _mm_storeu_pd(pf,ymm); return *this; }
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
inline vec<__m128d,double> blend(vec<__m128d,double> const &a,vec<__m128d,double> const &b,vec<__m128d,double> const &p) {
#ifdef __SSE4_1__
    return _mm_blendv_pd(a,b,p);
#else
    return _mm_or_pd(_mm_and_pd(p,b),_mm_andnot_pd(p,a));
#endif
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
template<> inline vec<float,float> & vec<float,float>::zero() { ymm = 0; return *this; }
template<> inline vec<float,float> & vec<float,float>::load1(float f) { ymm = f; return *this; }
template<> inline vec<float,float> & vec<float,float>::load(float *pf) { ymm = *pf; return *this; }
template<> inline const vec<float,float> & vec<float,float>::store(float *pf) const { *pf = ymm; return *this; }
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
* 32-bit integer
\**********************************************************************/
//template<> inline vec<int32_t,int32_t>::vec(const int32_t &d) { ymm = d; }
template<> inline vec<int32_t,int32_t> & vec<int32_t,int32_t>::zero() { ymm = 0; return *this; }
template<> inline vec<int32_t,int32_t> & vec<int32_t,int32_t>::load1(int32_t f) { ymm = f; return *this; }
template<> inline vec<int32_t,int32_t> & vec<int32_t,int32_t>::load(int32_t *pf) { ymm = *pf; return *this; }
template<> inline const vec<int32_t,int32_t> & vec<int32_t,int32_t>::store(int32_t *pf) const { *pf = ymm; return *this; }
inline vec<int32_t,int32_t> maskz_mov(mmask<bool> const &p,vec<int32_t,int32_t> const &a) {
    return p ? a : vec<int32_t,int32_t>(0); 
    }
inline vec<int32_t,int32_t> mask_mov(vec<int32_t,int32_t> const &src,mmask<bool> const &p,vec<int32_t,int32_t> const &a) {
	return p ? a : src;
    }

#endif/*#elif defined(__SSE__)*/
#endif /*__AVX512F__*/
#if __GNUC__ > 4
#pragma GCC diagnostic pop
#endif

#endif/*defined(__cplusplus)*/


#endif/*SIMD_H*/
