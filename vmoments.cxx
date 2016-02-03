#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include <math.h>
#include <immintrin.h>
#include "ilc.h"

template<class momFloat>
class vlocReduced {
public:
    momFloat m;
    momFloat x,y,z;
    momFloat xx,yy,xy,xz,yz;
    momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
    momFloat xxxxx,xyyyy,xxxxy,yyyyy,xxxxz,yyyyz,xxxyy,xxyyy,xxxyz,xyyyz,xxyyz;
    };

template<class momFloat>
class vmomReduced {
public:
    momFloat m;
    momFloat xx,yy,xy,xz,yz;
    momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
    };

/******************************************************************************/
class dvec {
    __m256d ymm;
public:
    dvec() {}
    dvec(double d) { ymm = _mm256_set1_pd(d); }
    dvec(__m256d const &d) { ymm = d; }
    operator __m256d() const { return ymm; }
    dvec & zero() { ymm = _mm256_setzero_pd(); return *this; }
    };
static inline dvec operator-(dvec const &a) {
    static const union {
	uint64_t i[4];
	__m256d d;
	} p = {{0x8000000000000000,0x8000000000000000,0x8000000000000000,0x8000000000000000}};
    return _mm256_xor_pd(a,p.d);
    }
static inline double hadd(dvec const &a) {
    __m256d t1 = _mm256_hadd_pd(a,a);
    __m128d t2 = _mm256_extractf128_pd(t1,1);
    __m128d t3 = _mm_add_sd(_mm256_castpd256_pd128(t1),t2);
    return _mm_cvtsd_f64(t3);        
    }
static inline dvec operator*(dvec const &a,dvec const &b) { return _mm256_mul_pd(a,b); }
static inline dvec operator/(dvec const &a,dvec const &b) { return _mm256_div_pd(a,b); }
static inline dvec operator+(dvec const &a,dvec const &b) { return _mm256_add_pd(a,b); }
static inline dvec operator-(dvec const &a,dvec const &b) { return _mm256_sub_pd(a,b); }
static inline dvec & operator+=(dvec &a,dvec const &b) { return a = a + b; }
static inline dvec & operator-=(dvec &a,dvec const &b) { return a = a - b; }
static inline dvec & operator*=(dvec &a,dvec const &b) { return a = a * b; }
static inline dvec & operator/=(dvec &a,dvec const &b) { return a = a / b; }
static inline dvec & zero(dvec &a) { return a = _mm256_setzero_pd(); }
/******************************************************************************/
class fvec {
    __m256 ymm;
public:
    fvec() {}
    fvec(float d) { ymm = _mm256_set1_ps(d); }
    fvec(__m256 const &d) { ymm = d; }
    operator __m256() const { return ymm; }
    fvec & zero() { ymm = _mm256_setzero_ps(); return *this; }
    fvec & load1(float f) {
	ymm = _mm256_setr_ps(f,0,0,0,0,0,0,0);
//1	ymm = _mm256_setzero_ps();
//	ymm = _mm256_castps128_ps256(_mm_load_ss(pf));
//2	ymm = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ss(pf)),_mm_setzero_ps(),1);
	return *this;
	}
    };
static inline fvec operator-(fvec const &a) {
    static const union {
	uint32_t i[8];
	__m256 f;
	} p = {{0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000}};
    return _mm256_xor_ps(a,p.f);
    }
static inline fvec operator*(fvec const &a,fvec const &b) { return _mm256_mul_ps(a,b); }
static inline fvec operator/(fvec const &a,fvec const &b) { return _mm256_div_ps(a,b); }
static inline fvec operator+(fvec const &a,fvec const &b) { return _mm256_add_ps(a,b); }
static inline fvec operator-(fvec const &a,fvec const &b) { return _mm256_sub_ps(a,b); }
static inline fvec & operator+=(fvec &a,fvec const &b) { return a = a + b; }
static inline fvec & operator-=(fvec &a,fvec const &b) { return a = a - b; }
static inline fvec & operator*=(fvec &a,fvec const &b) { return a = a * b; }
static inline fvec & operator/=(fvec &a,fvec const &b) { return a = a / b; }

static inline float hadd(fvec const &a) {
    __m256 t1 = _mm256_hadd_ps(a,a);
    __m256 t2 = _mm256_hadd_ps(t1,t1);
    __m128 t3 = _mm256_extractf128_ps(t2,1);
    return _mm_cvtss_f32(_mm_add_ss(_mm256_castps256_ps128(t2),t3));
    }

static inline fvec rsqrt(fvec const &r2) {
    fvec r = _mm256_rsqrt_ps(r2); /* Approximation */
    return r*(1.5 - 0.5*r*r*r2); /* Newton step correction */
    }

/******************************************************************************/

extern "C"
double momFlocrAddVFmomr5cm(FLOCR *l,float v1,ILC ill,double *tax,double *tay,double *taz) {
    const float onethird = 1.0f/3.0f;
    fvec u2,u3,u4;
    fvec xx,xy,xz,yy,yz,zz;
    fvec xxx,xxy,xyy,yyy,xxz,xyz,yyz;
    fvec R2xx,R2xy,R2xz,R2yy,R2yz,R2x,R2y,R2z,R2,R3xx,R3xy,R3yy,R3xz,R3yz,R3x,R3y,R3z,R3,R4x,R4y,R4z,R4;
    fvec T0,txx,tyy,t1,t1x,t1y,t1z,t1xx,t1yy,t2x,t2y,t2z,t2xx,t2yy,t3x,t3y,t3z,t3xx,t3yy,txxxx,tyyyy;
    fvec t2,vax,vay,vaz;
    fvec llm,llx,lly,llz,llxx,llyy,llxy,llxz,llyz,llxxx,llxyy,llxxy,llyyy,llxxz,llyyz,llxyz;
    fvec llxxxx,llxyyy,llxxxy,llyyyy,llxxxz,llyyyz,llxxyy,llxxyz,llxyyz;
    fvec llxxxxx,llxyyyy,llxxxxy,llyyyyy,llxxxxz,llyyyyz,llxxxyy,llxxyyy,llxxxyz,llxyyyz,llxxyyz;

    llm.load1(l->m);
    llx.load1(l->x);
    lly.load1(l->y);
    llz.load1(l->z);

    llxx.load1(l->xx);
    llyy.load1(l->yy);
    llxy.load1(l->xy);
    llxz.load1(l->xz);
    llyz.load1(l->yz);
    llxxx.load1(l->xxx);
    llxyy.load1(l->xyy);
    llxxy.load1(l->xxy);
    llyyy.load1(l->yyy);
    llxxz.load1(l->xxz);
    llyyz.load1(l->yyz);
    llxyz.load1(l->xyz);
    
    llxxxx.load1(l->xxxx);
    llxyyy.load1(l->xyyy);
    llxxxy.load1(l->xxxy);
    llyyyy.load1(l->yyyy);
    llxxxz.load1(l->xxxz);
    llyyyz.load1(l->yyyz);
    llxxyy.load1(l->xxyy);
    llxxyz.load1(l->xxyz);
    llxyyz.load1(l->xyyz);
    
    llxxxxx.load1(l->xxxxx);
    llxyyyy.load1(l->xyyyy);
    llxxxxy.load1(l->xxxxy);
    llyyyyy.load1(l->yyyyy);
    llxxxxz.load1(l->xxxxz);
    llyyyyz.load1(l->yyyyz);
    llxxxyy.load1(l->xxxyy);
    llxxyyy.load1(l->xxyyy);
    llxxxyz.load1(l->xxxyz);
    llxyyyz.load1(l->xyyyz);
    llxxyyz.load1(l->xxyyz);

    ILCTILE tile;
    ILC_LOOP(ill,tile) {
	ILC_BLK *blk = tile->blk;
	int nBlocks = tile->lstTile.nBlocks;
	int nInLast = tile->lstTile.nInLast;
	int nLeft;
	int j;
	for( j = nInLast; j&SIMD_MASK; j++) {
	    blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
	    blk[nBlocks].m.f[j] = 0.0f;
	    blk[nBlocks].u.f[j] = 0.0f;
	    }

	for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	    int n = ((nLeft ? ILC_PART_PER_BLK : nInLast) + SIMD_MASK) >> SIMD_BITS;
	    for (j=0; j<n; ++j) {
		fvec
		    x = blk->dx.p[j],
		    y = blk->dy.p[j],
		    z = blk->dz.p[j],
		    u = blk->u.p[j],
		    m = blk->m.p[j];
		fvec v=v1;
		fvec dir = rsqrt(x*x + y*y + z*z);

		u *= dir;
		x *= dir;
		y *= dir;
		z *= dir;
		dir = -dir;
		v *= dir;
		u2 = 15.0f*u*u;  /* becomes 15.0f*u2! */
		/*
		** Calculate the funky distance terms.
		*/
		xx = 0.5f*x*x;
		xy = x*y;
		yy = 0.5f*y*y;
		xz = x*z;
		yz = y*z;
		zz = 0.5f*z*z;
		xxx = x*(onethird*xx - zz);
		xxz = z*(xx - onethird*zz);
		yyy = y*(onethird*yy - zz);
		yyz = z*(yy - onethird*zz);
		xx -= zz;
		yy -= zz;
		xxy = y*xx;
		xyy = x*yy;
		xyz = xy*z;

		u3 = u2*u;  /* becomes 5.0f*u3! */

		R2xx = u2*blk->xx.p[j];
		R2xy = u2*blk->xy.p[j];
		R2xz = u2*blk->xz.p[j];
		R2yy = u2*blk->yy.p[j];
		R2yz = u2*blk->yz.p[j];

		u4 = 7.0f*u3*u;  /* becomes 7.0f*5.0f*u4! */

		R2x = x*R2xx + y*R2xy + z*R2xz;
		R2y = x*R2xy + y*R2yy + z*R2yz;
		R2z = x*R2xz + y*R2yz - z*(R2xx + R2yy);

		R3xx = u3*(x*blk->xxx.p[j] + y*blk->xxy.p[j] + z*blk->xxz.p[j]);
		R3xy = u3*(x*blk->xxy.p[j] + y*blk->xyy.p[j] + z*blk->xyz.p[j]);
		R3yy = u3*(x*blk->xyy.p[j] + y*blk->yyy.p[j] + z*blk->yyz.p[j]);
		R3xz = u3*(x*blk->xxz.p[j] + y*blk->xyz.p[j] - z*(fvec(blk->xxx.p[j]) + blk->xyy.p[j]));
		R3yz = u3*(x*blk->xyz.p[j] + y*blk->yyz.p[j] - z*(fvec(blk->xxy.p[j]) + blk->yyy.p[j]));

		R4x = u4*(blk->xxxx.p[j]*xxx + blk->xyyy.p[j]*yyy + blk->xxxy.p[j]*xxy + blk->xxxz.p[j]*xxz + blk->xxyy.p[j]*xyy + blk->xxyz.p[j]*xyz + blk->xyyz.p[j]*yyz);
		R4y = u4*(blk->xyyy.p[j]*xyy + blk->xxxy.p[j]*xxx + blk->yyyy.p[j]*yyy + blk->yyyz.p[j]*yyz + blk->xxyy.p[j]*xxy + blk->xxyz.p[j]*xxz + blk->xyyz.p[j]*xyz);
		R4z = u4*(-blk->xxxx.p[j]*xxz - (fvec(blk->xyyy.p[j]) + blk->xxxy.p[j])*xyz - blk->yyyy.p[j]*yyz + blk->xxxz.p[j]*xxx + blk->yyyz.p[j]*yyy - blk->xxyy.p[j]*(xxz + yyz) + blk->xxyz.p[j]*xxy + blk->xyyz.p[j]*xyy);


		R3x = 0.5f*(x*R3xx + y*R3xy + z*R3xz);
		R3y = 0.5f*(x*R3xy + y*R3yy + z*R3yz);
		R3z = 0.5f*(x*R3xz + y*R3yz - z*(R3xx + R3yy));

		R4 = 0.25f*(x*R4x + y*R4y + z*R4z);

		R2 = 0.5f*(x*R2x + y*R2y + z*R2z);

		R3 = onethird*(x*R3x + y*R3y + z*R3z);

		xx = x*x;
		yy = y*y;

		/*
		** Now we use the 'R's.
		*/
		llm += dir*(m + 0.2f*R2 + R3 + R4);

		dir *= v;
		T0 = -(m + R2 + 7.0f*R3 + 9.0f*R4);

		vax = dir*(T0*x + 0.2f*R2x + R3x + R4x);
		vay = dir*(T0*y + 0.2f*R2y + R3y + R4y);
		vaz = dir*(T0*z + 0.2f*R2z + R3z + R4z);

//		*tax = hadd();
//		*tay = hadd();
//		*taz = hadd();

		llx -= vax;
		lly -= vay;
		llz -= vaz;

		dir *= v;
		T0 = 3.0f*m + 7.0f*(R2 + 9.0f*R3);

		t1 = m + R2 + 7.0f*R3;
		t1x = R2x + 7.0f*R3x;
		t1y = R2y + 7.0f*R3y;
		t1z = R2z + 7.0f*R3z;
		llxx += dir*(T0*xx - t1 - 2.0f*x*t1x + 0.2f*R2xx + R3xx);
		llyy += dir*(T0*yy - t1 - 2.0f*y*t1y + 0.2f*R2yy + R3yy);
		llxy += dir*(T0*xy - y*t1x - x*t1y + 0.2f*R2xy + R3xy);
		llxz += dir*(T0*xz - z*t1x - x*t1z + 0.2f*R2xz + R3xz);
		llyz += dir*(T0*yz - z*t1y - y*t1z + 0.2f*R2yz + R3yz);

		dir *= v;
		T0 = 15.0f*m + 63.0f*R2;
		txx = T0*xx;
		tyy = T0*yy;

		t1 = 3.0f*m + 7.0f*R2;
		t2x = -7.0f*R2x;
		t2y = -7.0f*R2y;
		t2z = -7.0f*R2z;
		t1xx = txx - t1 + 2.0f*x*t2x + R2xx;
		t1yy = tyy - t1 + 2.0f*y*t2y + R2yy;

		llxxx += dir*(x*(txx - 3.0f*(t1 - t2x*x - R2xx)) + 3.0f*R2x);
		llyyy += dir*(y*(tyy - 3.0f*(t1 - t2y*y - R2yy)) + 3.0f*R2y);
		llxxy += dir*(y*t1xx + xx*t2y + R2y + 2.0f*R2xy*x);
		llxxz += dir*(z*t1xx + xx*t2z + R2z + 2.0f*R2xz*x);
		llxyy += dir*(x*t1yy + yy*t2x + R2x + 2.0f*R2xy*y);
		llyyz += dir*(z*t1yy + yy*t2z + R2z + 2.0f*R2yz*y);
		llxyz += dir*(T0*xyz + (yz*t2x + xz*t2y + xy*t2z) + R2xy*z + R2yz*x + R2xz*y);

		dir *= v*m;
		txx = 105.0f*xx;
		tyy = 105.0f*yy;
		t2xx = txx - 90.0f;
		t2yy = tyy - 90.0f;
		llxxxx += dir*(xx*t2xx + 9.0f);
		llyyyy += dir*(yy*t2yy + 9.0f);
		t2xx += 45.0f;
		t2yy += 45.0f;
		llxxxy += dir*xy*t2xx;
		llxxxz += dir*xz*t2xx;
		llxyyy += dir*xy*t2yy;
		llyyyz += dir*yz*t2yy;
		t2xx += 30.0f;
		t2yy += 30.0f;
		llxxyy += dir*(yy*t2xx - xx*15.0f + 3.0f);
		llxxyz += dir*(yz*t2xx);
		llxyyz += dir*(xz*t2yy);

		dir *= v;
		x *= dir;
		y *= dir;
		z *= dir;
		txx = 9.0f*xx - 10.0f;
		tyy = 9.0f*yy - 10.0f;
		xx *= 105.0f;
		yy *= 105.0f;
		xy *= z*105.0f;
		llxxxxx += x*(xx*txx + 225.0f);
		llyyyyy += y*(yy*tyy + 225.0f);
		txx += 4.0f;
		tyy += 4.0f;
		txxxx = xx*txx + 45.0f;
		tyyyy = yy*tyy + 45.0f;
		llxxxxy += y*txxxx;
		llxxxxz += z*txxxx;
		llxyyyy += x*tyyyy;
		llyyyyz += z*tyyyy;
		txx += 3.0f;
		tyy += 3.0f;
		llxxxyz += xy*txx;
		llxyyyz += xy*tyy;
		llxxxyy += x*(yy*txx - xx + 45.0f);
		llxxyyy += y*(xx*tyy - yy + 45.0f);
		tyy += 2.0f;
		llxxyyz += z*(xx*tyy - yy + 15.0f);
		}
	    }
	}

    // 32 terms * 4 instructions = 128 operations
    l->m = hadd(llm);
    l->x = hadd(llx);
    l->y = hadd(lly);
    l->z = hadd(llz);
    l->xx = hadd(llxx);
    l->yy = hadd(llyy);
    l->xy = hadd(llxy);
    l->xz = hadd(llxz);
    l->yz = hadd(llyz);
    l->xxx = hadd(llxxx);
    l->xyy = hadd(llxyy);
    l->xxy = hadd(llxxy);
    l->yyy = hadd(llyyy);
    l->xxz = hadd(llxxz);
    l->yyz = hadd(llyyz);
    l->xyz = hadd(llxyz);
    l->xxxx = hadd(llxxxx);
    l->xyyy = hadd(llxyyy);
    l->xxxy = hadd(llxxxy);
    l->yyyy = hadd(llyyyy);
    l->xxxz = hadd(llxxxz);
    l->yyyz = hadd(llyyyz);
    l->xxyy = hadd(llxxyy);
    l->xxyz = hadd(llxxyz);
    l->xyyz = hadd(llxyyz);
    l->xxxxx = hadd(llxxxxx);
    l->xyyyy = hadd(llxyyyy);
    l->xxxxy = hadd(llxxxxy);
    l->yyyyy = hadd(llyyyyy);
    l->xxxxz = hadd(llxxxxz);
    l->yyyyz = hadd(llyyyyz);
    l->xxxyy = hadd(llxxxyy);
    l->xxyyy = hadd(llxxyyy);
    l->xxxyz = hadd(llxxxyz);
    l->xyyyz = hadd(llxyyyz);
    l->xxyyz = hadd(llxxyyz);

    return(464.0);
}




#if 0
double momLocrAddFmomr5cm(LOCR *l,MOMR *m,double u1,double dx,double dy,double dz,double *tax,double *tay,double *taz) {
    const double onethird = 1.0/3.0;
    dvec u2,u3,u4;
    dvec xx,xy,xz,yy,yz,zz;
    dvec xxx,xxy,xyy,yyy,xxz,xyz,yyz;
    dvec R2xx,R2xy,R2xz,R2yy,R2yz,R2x,R2y,R2z,R2,R3xx,R3xy,R3yy,R3xz,R3yz,R3x,R3y,R3z,R3,R4x,R4y,R4z,R4;
    dvec T0,txx,tyy,t1,t1x,t1y,t1z,t1xx,t1yy,t2x,t2y,t2z,t2xx,t2yy,txxxx,tyyyy;
    dvec v;
    dvec dir;
    dvec x=dx, y=dy, z=dz, u=u1;

    dir = 1.0 / sqrt(dx*dx + dy*dy + dz*dz);

    u *= dir;
    x *= dir;
    y *= dir;
    z *= dir;
    dir = -dir;
    v = dir;
    u2 = 15.0*u*u;  /* becomes 15.0*u2! */
    /*
    ** Calculate the funky distance terms.
    */
    xx = 0.5*x*x;
    xy = x*y;
    yy = 0.5*y*y;
    xz = x*z;
    yz = y*z;
    zz = 0.5*z*z;
    xxx = x*(onethird*xx - zz);
    xxz = z*(xx - onethird*zz);
    yyy = y*(onethird*yy - zz);
    yyz = z*(yy - onethird*zz);
    xx -= zz;
    yy -= zz;
    xxy = y*xx;
    xyy = x*yy;
    xyz = xy*z;

    u3 = u2*u;  /* becomes 5.0*u3! */

    R2xx = u2*m->xx;
    R2xy = u2*m->xy;
    R2xz = u2*m->xz;
    R2yy = u2*m->yy;
    R2yz = u2*m->yz;

    u4 = 7.0*u3*u;  /* becomes 7.0*5.0*u4! */

    R2x = x*R2xx + y*R2xy + z*R2xz;
    R2y = x*R2xy + y*R2yy + z*R2yz;
    R2z = x*R2xz + y*R2yz - z*(R2xx + R2yy);

    R3xx = u3*(x*m->xxx + y*m->xxy + z*m->xxz);
    R3xy = u3*(x*m->xxy + y*m->xyy + z*m->xyz);
    R3yy = u3*(x*m->xyy + y*m->yyy + z*m->yyz);
    R3xz = u3*(x*m->xxz + y*m->xyz - z*(m->xxx + m->xyy));
    R3yz = u3*(x*m->xyz + y*m->yyz - z*(m->xxy + m->yyy));

    R4x = u4*(m->xxxx*xxx + m->xyyy*yyy + m->xxxy*xxy + m->xxxz*xxz + m->xxyy*xyy + m->xxyz*xyz + m->xyyz*yyz);
    R4y = u4*(m->xyyy*xyy + m->xxxy*xxx + m->yyyy*yyy + m->yyyz*yyz + m->xxyy*xxy + m->xxyz*xxz + m->xyyz*xyz);
    R4z = u4*(-m->xxxx*xxz - (m->xyyy + m->xxxy)*xyz - m->yyyy*yyz + m->xxxz*xxx + m->yyyz*yyy - m->xxyy*(xxz + yyz) + m->xxyz*xxy + m->xyyz*xyy);


    R3x = 0.5*(x*R3xx + y*R3xy + z*R3xz);
    R3y = 0.5*(x*R3xy + y*R3yy + z*R3yz);
    R3z = 0.5*(x*R3xz + y*R3yz - z*(R3xx + R3yy));

    R4 = 0.25*(x*R4x + y*R4y + z*R4z);

    R2 = 0.5*(x*R2x + y*R2y + z*R2z);

    R3 = onethird*(x*R3x + y*R3y + z*R3z);

    xx = x*x;
    yy = y*y;

    /*
    ** Now we use the 'R's.
    */
    l->m += dir*(m->m + 0.2*R2 + R3 + R4);

    dir *= v;
    T0 = -(m->m + R2 + 7.0*R3 + 9.0*R4);

    *tax = hadd(dir*(T0*x + 0.2*R2x + R3x + R4x));
    *tay = hadd(dir*(T0*y + 0.2*R2y + R3y + R4y));
    *taz = hadd(dir*(T0*z + 0.2*R2z + R3z + R4z));
    l->x -= *tax;
    l->y -= *tay;
    l->z -= *taz;

    dir *= v;
    T0 = 3.0*m->m + 7.0*(R2 + 9.0*R3);

    t1 = m->m + R2 + 7.0*R3;
    t1x = R2x + 7.0*R3x;
    t1y = R2y + 7.0*R3y;
    t1z = R2z + 7.0*R3z;
    l->xx += dir*(T0*xx - t1 - 2.0*x*t1x + 0.2*R2xx + R3xx);
    l->yy += dir*(T0*yy - t1 - 2.0*y*t1y + 0.2*R2yy + R3yy);
    l->xy += dir*(T0*xy - y*t1x - x*t1y + 0.2*R2xy + R3xy);
    l->xz += dir*(T0*xz - z*t1x - x*t1z + 0.2*R2xz + R3xz);
    l->yz += dir*(T0*yz - z*t1y - y*t1z + 0.2*R2yz + R3yz);

    dir *= v;
    T0 = 15.0*m->m + 63.0*R2;
    txx = T0*xx;
    tyy = T0*yy;

    t1 = 3.0*m->m + 7.0*R2;
    t2x = -7.0*R2x;
    t2y = -7.0*R2y;
    t2z = -7.0*R2z;
    t1xx = txx - t1 + 2.0*x*t2x + R2xx;
    t1yy = tyy - t1 + 2.0*y*t2y + R2yy;

    l->xxx += dir*(x*(txx - 3.0*(t1 - t2x*x - R2xx)) + 3.0*R2x);
    l->yyy += dir*(y*(tyy - 3.0*(t1 - t2y*y - R2yy)) + 3.0*R2y);
    l->xxy += dir*(y*t1xx + xx*t2y + R2y + 2.0*R2xy*x);
    l->xxz += dir*(z*t1xx + xx*t2z + R2z + 2.0*R2xz*x);
    l->xyy += dir*(x*t1yy + yy*t2x + R2x + 2.0*R2xy*y);
    l->yyz += dir*(z*t1yy + yy*t2z + R2z + 2.0*R2yz*y);
    l->xyz += dir*(T0*xyz + (yz*t2x + xz*t2y + xy*t2z) + R2xy*z + R2yz*x + R2xz*y);

    dir *= v*m->m;
    txx = 105.0*xx;
    tyy = 105.0*yy;
    t2xx = txx - 90.0;
    t2yy = tyy - 90.0;
    l->xxxx += dir*(xx*t2xx + 9.0);
    l->yyyy += dir*(yy*t2yy + 9.0);
    t2xx += 45.0;
    t2yy += 45.0;
    l->xxxy += dir*xy*t2xx;
    l->xxxz += dir*xz*t2xx;
    l->xyyy += dir*xy*t2yy;
    l->yyyz += dir*yz*t2yy;
    t2xx += 30.0;
    t2yy += 30.0;
    l->xxyy += dir*(yy*t2xx - xx*15.0 + 3.0);
    l->xxyz += dir*(yz*t2xx);
    l->xyyz += dir*(xz*t2yy);

    dir *= v;
    x *= dir;
    y *= dir;
    z *= dir;
    txx = 9.0*xx - 10.0;
    tyy = 9.0*yy - 10.0;
    xx *= 105.0;
    yy *= 105.0;
    xy *= z*105.0;
    l->xxxxx += x*(xx*txx + 225.0);
    l->yyyyy += y*(yy*tyy + 225.0);
    txx += 4.0;
    tyy += 4.0;
    txxxx = xx*txx + 45.0;
    tyyyy = yy*tyy + 45.0;
    l->xxxxy += y*txxxx;
    l->xxxxz += z*txxxx;
    l->xyyyy += x*tyyyy;
    l->yyyyz += z*tyyyy;
    txx += 3.0;
    tyy += 3.0;
    l->xxxyz += xy*txx;
    l->xyyyz += xy*tyy;
    l->xxxyy += x*(yy*txx - xx + 45.0);
    l->xxyyy += y*(xx*tyy - yy + 45.0);
    tyy += 2.0;
    l->xxyyz += z*(xx*tyy - yy + 15.0);
    return(464.0);
}
#endif


