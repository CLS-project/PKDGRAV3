#ifndef SIMD_H
#define SIMD_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_SIMD
#include <stdint.h>

#if defined(__AVX__)
#define SIMD_BITS 3
#elif defined(__SSE__) || defined(__ALTIVEC__)
#define SIMD_BITS 2
#else
#define SIMD_BITS 0
#endif
#define SIMD_DBITS (SIMD_BITS-1)


#define SIMD_WIDTH (1<<SIMD_BITS)
#define SIMD_DWIDTH (1<<SIMD_DBITS)
#define SIMD_MASK (SIMD_WIDTH-1)
#define SIMD_DMASK (SIMD_DWIDTH-1)

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
#ifdef __FMA4__
#include <x86intrin.h>
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
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256  v_sf;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256d v_df;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256  v_bool;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m256) __m256i v_i;
#define MM_FCN(f,p) _mm256_##f##_##p
#define MM_CMP(F,f,p,a,b) _mm256_cmp_##p(a,b,_CMP_##F##_OQ)
#else
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128  v_sf;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128d v_df;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128  v_bool;
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128i v_i;
#define MM_FCN(f,p) _mm_##f##_##p
#define MM_CMP(F,f,p,a,b) _mm_cmp##f##_##p(a,b)
#endif
typedef ATTRIBUTE_ALIGNED_ALIGNOF(__m128) __m128i v_i4;
#else
typedef vector float v_sf;
typedef vector bool int v_bool;
typedef vector bool int v_i;
#endif

#if SIMD_WIDTH==4
#define SIMD_CONST(c) {c,c,c,c}
#define SIMD_DCONST(c) {c,c}
#elif SIMD_WIDTH==8
#define SIMD_CONST(c) {c,c,c,c,c,c,c,c}
#define SIMD_DCONST(c) {c,c,c,c}
#else
#error Invalid SIMD_WIDTH
#endif

typedef union {
    float f[SIMD_WIDTH];
    v_sf p;
    } vfloat;

typedef union {
    double d[SIMD_DWIDTH];
    v_df p;
    } vdouble;

typedef union {
    int32_t  i[SIMD_WIDTH];
    uint32_t u[SIMD_WIDTH];
    v_i      p;
    v_i4     ph;
    v_sf     pf;
    v_df     pd;
    } vint;

typedef union {
    uint64_t u[SIMD_WIDTH];
    v_i      pi;
    v_df     pd;
    } vint64;

#if defined(HAVE_POSIX_MEMALIGN)
static inline void * SIMD_malloc( size_t newSize ) {
    void *np;

    if ( posix_memalign( &np, sizeof(v_sf), newSize ) == 0 )
	return np;
    else
	return NULL;
    }

static inline void SIMD_free(void *p) {
    free(p);
    }
#else
#define SIMD_malloc malloc
#define SIMD_free free
#endif

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
#ifdef __AVX__
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
#define SIMD_MADD(a,b,c) MM_FCN(macc,ps)(a,b,c)
#define SIMD_MSUB(a,b,c) MM_FCN(msub,ps)(a,b,c)
#define SIMD_NMADD(a,b,c) MM_FCN(nmacc,ps)(a,b,c)
#elif defined(__FMA__)
#define SIMD_MADD(a,b,c) MM_FCN(fmadd,ps)(a,b,c)  /*  (a*b) + c */
#define SIMD_MSUB(a,b,c) MM_FCN(fmsub,ps)(a,b,c)  /*  (a*b) - c */
#define SIMD_NMADD(a,b,c) MM_FCN(fnmadd,ps)(a,b,c) /* -(a*b) + c*/
#else
#define SIMD_MADD(a,b,c) MM_FCN(add,ps)(MM_FCN(mul,ps)(a,b),c)
#define SIMD_MSUB(a,b,c) MM_FCN(sub,ps)(MM_FCN(mul,ps)(a,b),c)
#define SIMD_NMADD(a,b,c) MM_FCN(sub,ps)(c,MM_FCN(mul,ps)(a,b))
#endif
#define SIMD_DIV(a,b) MM_FCN(div,ps)(a,b)
#define SIMD_SQRT(a) MM_FCN(sqrt,ps)(a)
#define SIMD_RSQRT(a) MM_FCN(rsqrt,ps)(a)
#define SIMD_RE(a) MM_FCN(rcp,ps)(a);
static inline v_sf SIMD_RE_EXACT(v_sf a) {
    static const vfloat one = {SIMD_CONST(1.0)};
    v_sf r = SIMD_RE(a);
    return MM_FCN(add,ps)(MM_FCN(mul,ps)(r,MM_FCN(sub,ps)(one.p,MM_FCN(mul,ps)(r,a))),r);
    }
#define SIMD_MAX(a,b) MM_FCN(max,ps)(a,b)
#define SIMD_MIN(a,b) MM_FCN(min,ps)(a,b)
#ifdef __AVX__
#define SIMD_FLOOR(a) MM_FCN(floor,ps)(a)
#else
#define SIMD_FLOOR(a) _mm_cvtepi32_ps(_mm_cvtps_epi32(a))
#endif
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
#else
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
#else
#define SIMD_I2F(a) MM_FCN(castsi128,ps)(a)
#define SIMD_F2I(a) MM_FCN(castps,si128)(a)
#define SIMD_CMP_EQ_EPI32(a,b) MM_FCN(cmpeq,epi32)(a,b)
#define SIMD_CMP_GT_EPI32(a,b) MM_FCN(cmpgt,epi32)(a,b)
#define SIMD_AND_EPI32(a,b) MM_FCN(and,si128)(a,b)
#define SIMD_ANDNOT_EPI32(a,b) MM_FCN(andnot,si128)(a,b)
#define SIMD_OR_EPI32(a,b) MM_FCN(or,si128)(a,b)
#endif
#else
static v_sf   simd_zero  = {SIMD_CONST(0)};
static v_bool simd_false = {SIMD_CONST(0)};
#define SIMD_MUL(a,b) vec_madd(a,b,simd_zero)
#define SIMD_ADD(a,b) vec_add(a,b)
#define SIMD_SUB(a,b) vec_sub(a,b)
#define SIMD_MADD(a,b,c) vec_madd(a,b,c)
#define SIMD_NMADD(a,b,c) vec_nmsub(a,b,c)
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
#endif

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

#ifdef __SSE3__
static inline float SIMD_HADD(v_sf p) {
    vfloat r;
    r.p = p;
    r.p = MM_FCN(hadd,ps)(r.p,r.p);
    r.p = MM_FCN(hadd,ps)(r.p,r.p);
#ifdef __AVX__
    r.f[0] += r.f[4];
#endif
    return r.f[0];
    }
#else
static inline float SIMD_HADD(v_sf p) {
    vfloat r;
    float f;
    int i;
    r.p = p;
    f = r.f[0];
    for(i=1; i<SIMD_WIDTH; i++) f += r.f[i];
    return f;
    }
#endif

/* With SSE2 and beyond we have double support */
#ifdef __SSE2__
static inline v_df SIMD_DSPLAT(double f) {
    return MM_FCN(set1,pd)(f);
    }
static inline v_df SIMD_DLOADS(double f) {
#ifdef __AVX__
    return _mm256_set_pd (0.0,0.0,0.0,f);
#else
    return MM_FCN(set,sd)(f);
#endif
    }

#define SIMD_DMUL(a,b) MM_FCN(mul,pd)(a,b)
#define SIMD_DADD(a,b) MM_FCN(add,pd)(a,b)
#define SIMD_DSUB(a,b) MM_FCN(sub,pd)(a,b)
#ifdef __FMA4__
#define SIMD_DMADD(a,b,c) MM_FCN(macc,pd)(a,b,c)
#define SIMD_DNMADD(a,b,c) MM_FCN(nmacc,pd)(a,b,c)
#else
#define SIMD_DMADD(a,b,c) MM_FCN(add,pd)(MM_FCN(mul,pd)(a,b),c)
#define SIMD_DNMADD(a,b,c) MM_FCN(sub,pd)(c,MM_FCN(mul,pd)(a,b))
#endif
#define SIMD_DDIV(a,b) MM_FCN(div,pd)(a,b)
#define SIMD_DSQRT(a) MM_FCN(sqrt,pd)(a)
#define SIMD_DMAX(a,b) MM_FCN(max,pd)(a,b)
#define SIMD_DMIN(a,b) MM_FCN(min,pd)(a,b)
#ifdef __AVX__
#define SIMD_DFLOOR(a) MM_FCN(floor,pd)(a)
#else
#define SIMD_DFLOOR(a)  _mm_cvtepi32_pd(_mm_cvtpd_epi32(a))
#endif
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
#define SIMD_DALL_ZERO(a) MM_FCN(movemask,pd)(a)

#ifdef __AVX__
#define SIMD_I2D(a) MM_FCN(castsi256,pd)(a)
#define SIMD_D2I(a) MM_FCN(castpd,si256)(a)
#define SIMD_D2F(a) MM_FCN(castps128,ps256)(MM_FCN(cvtpd,ps)(a))
#else
#define SIMD_I2D(a) MM_FCN(castsi128,pd)(a)
#define SIMD_D2I(a) MM_FCN(castpd,si128)(a)
#define SIMD_D2F(a) MM_FCN(cvtpd,ps)(a)
#endif
static inline v_df SIMD_DRSQRT(v_df B) {
    static const vdouble one = {SIMD_DCONST(1.0)};
    static const vdouble half = {SIMD_DCONST(0.5)};
    static const vdouble three = {SIMD_DCONST(3.0)};
    v_df r;
    r = MM_FCN(cvtps,pd)(_mm_rsqrt_ps(MM_FCN(cvtpd,ps)(B))); 
    /* vvv ** Up to two more interations for full precision are needed ** vvv*/
    r = SIMD_DMUL(SIMD_DMUL(half.p,r),SIMD_DSUB(three.p,SIMD_DMUL(SIMD_DMUL(B,r),r)));
    return r;
    }
static inline v_df SIMD_DRSQRT_EXACT(v_df B) {
    static const vdouble one = {SIMD_DCONST(1.0)};
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
#endif


#else
#define SIMD_malloc malloc
#define SIMD_free free
#endif

#endif
