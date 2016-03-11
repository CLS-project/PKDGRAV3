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
#if defined(__SSE4_1__)
#include <smmintrin.h>
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
static inline void * SIMD_malloc(size_t newSize) {
	return _mm_malloc(newSize, sizeof(vfloat));
}

static inline void SIMD_free(void *p) {
	_mm_free(p);
}
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

/* p==false then select a -- p==true, select b */
#ifdef __AVX__
#define SIMD_SELECT(a,b,p) MM_FCN(blendv,ps)(a,b,p)
#define SIMD_DSELECT(a,b,p) MM_FCN(blendv,pd)(a,b,p)
#else
#define SIMD_SELECT(a,b,p) SIMD_OR(SIMD_AND(p,b),SIMD_ANDNOT(p,a))
#define SIMD_DSELECT(a,b,p) SIMD_DOR(SIMD_DAND(p,b),SIMD_DANDNOT(p,a))
#endif

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
#endif


#else
#define SIMD_malloc malloc
#define SIMD_free free
#endif


#ifdef __cplusplus
/**********************************************************************\
* SIMD Vector class template
\**********************************************************************/
template<typename vtype,typename ftype>
class vec {
    vtype ymm;
public:
    vec() {}
    vec(const ftype &d);
    vec(vtype const &d) { ymm = d; }
    operator vtype() const { return ymm; }
    ftype operator [] (uint32_t idx) const;
    vec & zero();
    vec & load1(ftype f);
    vec & load(ftype *f);
    const vec & store(ftype *f) const;
    static int width() { return sizeof(vtype)/sizeof(ftype); }
    static int mask()  { width()-1; }
    };

#if defined(__AVX__)
/**********************************************************************\
* AVX single precision
\**********************************************************************/

inline vec<__m256,float>::vec(const float &d) { ymm = _mm256_set1_ps(d); }
inline vec<__m256,float> & vec<__m256,float>::zero() { ymm = _mm256_setzero_ps(); return *this; }
inline vec<__m256,float> & vec<__m256,float>::load1(float f) { ymm = _mm256_setr_ps(f,0,0,0,0,0,0,0); return *this; }
inline vec<__m256,float> & vec<__m256,float>::load(float *pf) { ymm = _mm256_loadu_ps(pf); return *this; }
inline const vec<__m256,float> & vec<__m256,float>::store(float *pf) const { _mm256_storeu_ps(pf,ymm); return *this; }
inline vec<__m256,float> operator-(vec<__m256,float> const &a) {
    return _mm256_xor_ps(a,_mm256_castsi256_ps(_mm256_set1_epi32(0x80000000)));
    }

inline vec<__m256,float> min(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_min_ps(a,b); }
inline vec<__m256,float> max(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_max_ps(a,b); }
inline vec<__m256,float> cmp(vec<__m256,float> const &a,vec<__m256,float> const &b, const int imm8) { return _mm256_cmp_ps(a,b,imm8); }
inline vec<__m256,float> operator*(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_mul_ps(a,b); }
inline vec<__m256,float> operator/(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_div_ps(a,b); }
inline vec<__m256,float> operator+(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_add_ps(a,b); }
inline vec<__m256,float> operator-(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_sub_ps(a,b); }
inline vec<__m256,float> operator==(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_EQ_OQ); }
inline vec<__m256,float> operator!=(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_NEQ_OQ); }
inline vec<__m256,float> operator>(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_GT_OQ); }
inline vec<__m256,float> operator<(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_LT_OQ); }
inline vec<__m256,float> operator>=(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_GE_OQ); }
inline vec<__m256,float> operator<=(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_cmp_ps(a,b,_CMP_LE_OQ); }
inline vec<__m256,float> operator&(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_and_ps(a,b); }
inline vec<__m256,float> operator|(vec<__m256,float> const &a,vec<__m256,float> const &b) { return _mm256_or_ps(a,b); }
inline float hadd(vec<__m256,float> const &a) {
    __m256 t1 = _mm256_hadd_ps(a,a);
    __m256 t2 = _mm256_hadd_ps(t1,t1);
    __m128 t3 = _mm256_extractf128_ps(t2,1);
    return _mm_cvtss_f32(_mm_add_ss(_mm256_castps256_ps128(t2),t3));
    }
inline int movemask(vec<__m256,float> const &r2) { return _mm256_movemask_ps(r2); }
inline vec<__m256,float> sqrt(vec<__m256,float> const &r2) { return _mm256_sqrt_ps(r2); }
inline vec<__m256,float> rsqrt(vec<__m256,float> const &r2) {
    vec<__m256,float> r = _mm256_rsqrt_ps(r2); /* Approximation */
    return r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
    }

/**********************************************************************\
* AVX double precision
\**********************************************************************/

inline vec<__m256d,double>::vec(const double &d) { ymm = _mm256_set1_pd(d); }
inline vec<__m256d,double> & vec<__m256d,double>::zero() { ymm = _mm256_setzero_pd(); return *this; }
inline vec<__m256d,double> & vec<__m256d,double>::load1(double f) { ymm = _mm256_setr_pd(f,0,0,0); return *this; }
inline vec<__m256d,double> & vec<__m256d,double>::load(double *pf) { ymm = _mm256_loadu_pd(pf); return *this; }
inline const vec<__m256d,double> & vec<__m256d,double>::store(double *pf) const { _mm256_storeu_pd(pf,ymm); return *this; }
inline vec<__m256d,double>  operator-(vec<__m256d,double> const &a) {
    return _mm256_xor_pd(a,_mm256_castsi256_pd(_mm256_set1_epi64x(0x8000000000000000)));
    }
inline vec<__m256d,double> min(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_min_pd(a,b); }
inline vec<__m256d,double> max(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_max_pd(a,b); }
inline vec<__m256d,double> cmp(vec<__m256d,double> const &a,vec<__m256d,double> const &b, const int imm8) { return _mm256_cmp_pd(a,b,imm8); }
inline vec<__m256d,double> operator*(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_mul_pd(a,b); }
inline vec<__m256d,double> operator/(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_div_pd(a,b); }
inline vec<__m256d,double> operator+(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_add_pd(a,b); }
inline vec<__m256d,double> operator-(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_sub_pd(a,b); }
inline vec<__m256d,double> operator==(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_EQ_OQ); }
inline vec<__m256d,double> operator!=(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_NEQ_OQ); }
inline vec<__m256d,double> operator>(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_GT_OQ); }
inline vec<__m256d,double> operator<(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_LT_OQ); }
inline vec<__m256d,double> operator>=(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_GE_OQ); }
inline vec<__m256d,double> operator<=(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_cmp_pd(a,b,_CMP_LE_OQ); }
inline vec<__m256d,double> operator&(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_and_pd(a,b); }
inline vec<__m256d,double> operator|(vec<__m256d,double> const &a,vec<__m256d,double> const &b) { return _mm256_or_pd(a,b); }

inline int movemask(vec<__m256d,double> const &r2) { return _mm256_movemask_pd(r2); }
inline vec<__m256d,double> sqrt(vec<__m256d,double> const &r2) { return _mm256_sqrt_pd(r2); }
inline vec<__m256d,double> rsqrt(vec<__m256d,double> const &r2) {
    vec<__m256d,double> r = _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(r2))); /* Approximation */
    r = r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
    return r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
    }
#endif

#if defined(__SSE__)

/**********************************************************************\
* SSE single precision
\**********************************************************************/

inline vec<__m128,float>::vec(const float &d) { ymm = _mm_set1_ps(d); }
inline vec<__m128,float> & vec<__m128,float>::zero() { ymm = _mm_setzero_ps(); return *this; }
inline vec<__m128,float> & vec<__m128,float>::load1(float f) { ymm = _mm_setr_ps(f,0,0,0); return *this; }
inline vec<__m128,float> & vec<__m128,float>::load(float *pf) { ymm = _mm_loadu_ps(pf); return *this; }
inline const vec<__m128,float> & vec<__m128,float>::store(float *pf) const { _mm_storeu_ps(pf,ymm); return *this; }
inline vec<__m128,float> operator-(vec<__m128,float> const &a) {
    return _mm_xor_ps(a,_mm_castsi128_ps(_mm_set1_epi32(0x80000000)));
    }
inline vec<__m128,float> min(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_min_ps(a,b); }
inline vec<__m128,float> max(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_max_ps(a,b); }
inline vec<__m128,float> operator*(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_mul_ps(a,b); }
inline vec<__m128,float> operator/(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_div_ps(a,b); }
inline vec<__m128,float> operator+(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_add_ps(a,b); }
inline vec<__m128,float> operator-(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_sub_ps(a,b); }
inline vec<__m128,float> operator==(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpeq_ps(a,b); }
inline vec<__m128,float> operator!=(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpneq_ps(a,b); }
inline vec<__m128,float> operator>(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpgt_ps(a,b); }
inline vec<__m128,float> operator<(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmplt_ps(a,b); }
inline vec<__m128,float> operator>=(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmpge_ps(a,b); }
inline vec<__m128,float> operator<=(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_cmple_ps(a,b); }
inline vec<__m128,float> operator&(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_and_ps(a,b); }
inline vec<__m128,float> operator|(vec<__m128,float> const &a,vec<__m128,float> const &b) { return _mm_or_ps(a,b); }
#if defined(__SSE3__)
inline float hadd(vec<__m128,float> const &a) {
    __m128 t1 = _mm_hadd_ps(a,a);
    __m128 t2 = _mm_hadd_ps(t1,t1);
    return _mm_cvtss_f32(t2);
    }
#endif
inline int movemask(vec<__m128,float> const &r2) { return _mm_movemask_ps(r2); }
inline vec<__m128,float> sqrt(vec<__m128,float> const &r2) { return _mm_sqrt_ps(r2); }
inline vec<__m128,float> rsqrt(vec<__m128,float> const &r2) {
    vec<__m128,float> r = _mm_rsqrt_ps(r2); /* Approximation */
    return r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
    }
#endif

#if defined(__SSE2__)
/**********************************************************************\
* SSE double precision
\**********************************************************************/

inline vec<__m128d,double>::vec(const double &d) { ymm = _mm_set1_pd(d); }
inline vec<__m128d,double> & vec<__m128d,double>::zero() { ymm = _mm_setzero_pd(); return *this; }
inline vec<__m128d,double> & vec<__m128d,double>::load1(double f) { ymm = _mm_setr_pd(f,0); return *this; }
inline vec<__m128d,double> & vec<__m128d,double>::load(double *pf) { ymm = _mm_loadu_pd(pf); return *this; }
inline const vec<__m128d,double> & vec<__m128d,double>::store(double *pf) const { _mm_storeu_pd(pf,ymm); return *this; }
inline vec<__m128d,double>  operator-(vec<__m128d,double> const &a) {
    return _mm_xor_pd(a,_mm_castsi128_pd(_mm_set1_epi64x(0x8000000000000000)));
    }
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
#if defined(__SSE3__)
inline double hadd(vec<__m128d,double> const &a) {
    return _mm_cvtsd_f64(_mm_hadd_pd(a,a));
    }
#endif
inline int movemask(vec<__m128d,double> const &r2) { return _mm_movemask_pd(r2); }
inline vec<__m128d,double> sqrt(vec<__m128d,double> const &r2) { return _mm_sqrt_pd(r2); }
inline vec<__m128d,double> rsqrt(vec<__m128d,double> const &r2) {
    vec<__m128d,double> r = _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(r2))); /* Approximation */
    r = r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
    return r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
    }
#endif

/**********************************************************************\
* Generic Operators
\**********************************************************************/

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
template<typename v,typename ftype> inline ftype vec<v,ftype>::operator [] (uint32_t idx) const {
    ftype d[width()];
    store(d);
    return d[idx];
    }
#if defined(__AVX__)
typedef vec<__m256,float> fvec;
typedef vec<__m256d,double> dvec;
#elif defined(__SSE__)
typedef vec<__m128,float> fvec;
#if defined(__SSE2__)
typedef vec<__m128d,double> dvec;
#endif
#else
#endif
#endif

#endif
