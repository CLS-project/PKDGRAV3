#ifndef SIMD_H
#define SIMD_H
#ifdef HAVE_CONFIG_H
#include "config.h"
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
#elif defined(__ALTIVEC__)
#include <altivec.h>
#include <math.h> /* for sqrtf() */
#else
#error 'SIMD selected, but no known SIMD supported'
#endif

#ifdef HAVE_POSIX_MEMALIGN
#include <stdlib.h>
#endif

#if defined(__SSE__)
typedef __m128 v4sf;
typedef __m128 v4bool;
typedef int v4i __attribute__ ((vector_size(16)));
#else
typedef vector float v4sf;
typedef vector bool int v4bool;
typedef vector bool int v4i;
#endif

#if defined(HAVE_POSIX_MEMALIGN)
static inline void * SIMD_malloc( size_t newSize )
{
    void *np;

    if ( posix_memalign( &np, sizeof(v4sf), newSize ) == 0 )
        return np;
    else
        return NULL;
}

static inline void SIMD_free(void *p)
{
    free(p);
}
#else
#define SIMD_malloc(a) malloc(a)
#define SIMD_free(a) free(a)
#endif

static inline v4sf SIMD_SPLAT(float f) {
#ifdef __SSE__
    return _mm_set1_ps(f);
#else
    v4 r;
    r.f[0] = r.f[1] = r.f[2] = r.f[3] = f;
    return r.p;
    /*return (v4sf)(f,f,f,f);*/
#endif
}

static inline v4sf SIMD_LOADS(float f) {
#ifdef __SSE__
    return _mm_set_ss(f);
#else
    v4 r;
    r.f[0] = f;
    r.f[1] = r.f[2] = r.f[3] = 0;
    return r.p;
    /*return (v4sf)(f,0,0,0);*/
#endif
}

static inline v4sf SIMD_LOAD(float a, float b, float c, float d ) {
#ifdef __SSE__
    return _mm_set_ps(a,b,c,d);
#else
    v4 r;
    r.f[0] = a;
    r.f[1] = b;
    r.f[2] = c;
    r.f[3] = d;
    return r.p;
    /*return (v4sf)(a,b,c,d);*/
#endif
}


#if defined(__SSE__)
#define SIMD_MUL(a,b) _mm_mul_ps(a,b)
#define SIMD_ADD(a,b) _mm_add_ps(a,b)
#define SIMD_SUB(a,b) _mm_sub_ps(a,b)
#define SIMD_MADD(a,b,c) _mm_add_ps(_mm_mul_ps(a,b),c)
#define SIMD_NMSUB(a,b,c) _mm_sub_ps(c,_mm_mul_ps(a,b))
#define SIMD_DIV(a,b) _mm_div_ps(a,b)
#define SIMD_RSQRT(a) _mm_rsqrt_ps(a)
#define SIMD_RE(a) _mm_rcp_ps(a);
static inline v4sf SIMD_RE_EXACT(v4sf a) {
    static const v4sf one = {1.0,1.0,1.0,1.0};
    v4sf r = SIMD_RE(a);
    return _mm_add_ps(_mm_mul_ps(r,_mm_sub_ps(one,_mm_mul_ps(r,a))),r);
}
#define SIMD_MAX(a,b) _mm_max_ps(a,b)
#define SIMD_MIN(a,b) _mm_min_ps(a,b)
#define SIMD_CMP_LE(a,b) _mm_cmple_ps(a,b)
#define SIMD_CMP_LT(a,b) _mm_cmplt_ps(a,b)
#define SIMD_CMP_GE(a,b) _mm_cmpge_ps(a,b)
#define SIMD_CMP_GT(a,b) _mm_cmpgt_ps(a,b)
#define SIMD_AND(a,b) _mm_and_ps(a,b)
#define SIMD_ANDNOT(a,b) _mm_andnot_ps(a,b)
#define SIMD_OR(a,b) _mm_or_ps(a,b)
#define SIMD_XOR(a,b) _mm_xor_ps(a,b)
#define SIMD_ALL_ZERO(a) _mm_movemask_ps(a)
#else
#define SIMD_MUL(a,b) vec_madd(a,b,(v4sf)(0,0,0,0))
#define SIMD_ADD(a,b) vec_add(a,b)
#define SIMD_SUB(a,b) vec_sub(a,b)
#define SIMD_MADD(a,b,c) vec_madd(a,b,c)
#define SIMD_NMSUB(a,b,c) vec_nmsub(a,b,c)
#define SIMD_RSQRT(a) vec_rsqrte(a)
#define SIMD_RE(a) vec_re(a)
static inline v4sf SIMD_RE_EXACT( v4sf a) {
    static const v4sf one = {1.0,1.0,1.0,1.0};
    v4sf r = SIMD_RE(a);
    return vec_madd(r,vec_nmsub(r,a,one),r);
}
#define SIMD_MAX(a,b) vec_max(a,b)
#define SIMD_MIN(a,b) vec_min(a,b)
#define SIMD_CMP_LE(a,b) vec_cmple(a,b)
#define SIMD_CMP_LT(a,b) vec_cmplt(a,b)
#define SIMD_CMP_GE(a,b) vec_cmpge(a,b)
#define SIMD_CMP_GT(a,b) vec_cmpgt(a,b)
#define SIMD_AND(a,b) vec_and(a,b)
#define SIMD_ANDNOT(a,b) vec_andc(b,a)
#define SIMD_OR(a,b) vec_or(a,b)
#define SIMD_XOR(a,b) vec_xor(a,b)
#define SIMD_ALL_ZERO(a) (!vec_all_eq(a,(v4bool)(0,0,0,0)))
#endif

static inline v4sf SIMD_RSQRT_EXACT(v4sf B) {
    static const v4sf one = {1.0,1.0,1.0,1.0};
    static const v4sf half= {0.5,0.5,0.5,0.5};
    v4sf r, t1, t2;
    r = SIMD_RSQRT(B);          /* dir = 1/sqrt(d2) */
    t1 = SIMD_MUL(r,r);         /* dir*dir */
    t2 = SIMD_MUL(r,half);      /* 0.5*dir */
    t1 = SIMD_NMSUB(B,t1,one);  /* 1 - dir*dir*a */
    r = SIMD_MADD(t1,t2,r);     /* 0.5*dir*(1-dir*dir*a)+dir */
    return r;
}

#ifdef __SSE3__
static inline float SIMD_HADD(v4sf p) {
    register v4sf r = p;
    float f;

    r = _mm_hadd_ps(r,r);
    r = _mm_hadd_ps(r,r);
    _mm_store_ss(&f,r);
    return f;
}
#else
static inline float SIMD_HADD(v4sf p) {
    union {
	v4sf p;
	float f[4];
    } r;
    r.p = p;
    return = r.f[0] + r.f[1] + r.f[2] + r.f[3];
}
#endif


#else
#define SIMD_malloc(a) malloc(a)
#define SIMD_free(a) free(a)
#endif

#endif
