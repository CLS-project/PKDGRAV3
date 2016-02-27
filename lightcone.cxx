#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include <math.h>
#include <immintrin.h>
#include "pkd.h"

/******************************************************************************/
class dvec {
    __m256d ymm;
public:
    dvec() {}
    dvec(double d) { ymm = _mm256_set1_pd(d); }
    dvec(__m256d const &d) { ymm = d; }
    operator __m256d() const { return ymm; }
    dvec & zero() { ymm = _mm256_setzero_pd(); return *this; }
    dvec & load1(double d) {
	ymm = _mm256_setr_pd(d,0,0,0);
	return *this;
	}
    dvec & load(double *pd) {
	ymm = _mm256_loadu_pd(pd);
	return *this;
	}
    void store(double * pd) const {
        _mm256_storeu_pd(pd, ymm);
	}
    double operator [] (uint32_t idx) const {
	double d[4];
	store(d);
	return d[idx];
	}
    };

static inline dvec operator-(dvec const &a) {
    static const union {
	uint32_t i[8];
	__m256d d;
	} p = {{0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000,0x80000000}};
    return _mm256_xor_pd(a,p.d);
    }
static inline dvec operator*(dvec const &a,dvec const &b) { return _mm256_mul_pd(a,b); }
static inline dvec operator/(dvec const &a,dvec const &b) { return _mm256_div_pd(a,b); }
static inline dvec operator+(dvec const &a,dvec const &b) { return _mm256_add_pd(a,b); }
static inline dvec operator-(dvec const &a,dvec const &b) { return _mm256_sub_pd(a,b); }
static inline dvec operator>(dvec const &a,dvec const &b) { return _mm256_cmp_pd(a,b,_CMP_GT_OQ); }
static inline dvec operator>=(dvec const &a,dvec const &b) { return _mm256_cmp_pd(a,b,_CMP_GE_OQ); }
static inline dvec operator<(dvec const &a,dvec const &b) { return _mm256_cmp_pd(a,b,_CMP_LT_OQ); }
static inline dvec operator<=(dvec const &a,dvec const &b) { return _mm256_cmp_pd(a,b,_CMP_LE_OQ); }
static inline dvec operator&(dvec const &a,dvec const &b) { return _mm256_and_pd(a,b); }
static inline dvec operator|(dvec const &a,dvec const &b) { return _mm256_or_pd(a,b); }
static inline dvec & operator+=(dvec &a,dvec const &b) { return a = a + b; }
static inline dvec & operator-=(dvec &a,dvec const &b) { return a = a - b; }
static inline dvec & operator*=(dvec &a,dvec const &b) { return a = a * b; }
static inline dvec & operator/=(dvec &a,dvec const &b) { return a = a / b; }
static inline dvec & operator&=(dvec &a,dvec const &b) { return a = a & b; }
static inline dvec & operator|=(dvec &a,dvec const &b) { return a = a | b; }
static inline int allzero(dvec const &r2) { return _mm256_movemask_pd(r2); }

static inline dvec sqrt(dvec const &r2) {
    return _mm256_sqrt_pd(r2);
    }
static inline dvec max(dvec const &a,dvec const &b) {
    return _mm256_max_pd(a,b);
    }
static inline dvec min(dvec const &a,dvec const &b) {
    return _mm256_min_pd(a,b);
    }
/******************************************************************************/

extern "C" void addToLightCone(PKD pkd,double *r,PARTICLE *p,int bParticleOutput);

#define NBOX 184

extern "C"
void pkdProcessLightCone(PKD pkd,PARTICLE *p,double dLookbackFac,double dLookbackFacLCP,double dDriftDelta,double dKickDelta) {
    const double dLightSpeed = dLightSpeedSim(pkd->param.dBoxSize);
    const double mrLCP = dLightSpeed*dLookbackFacLCP;

    int nBox = NBOX; /* Check all 184 by default */
    double dxStart;

    dxStart = (dLookbackFac*dLightSpeed - 3.0)/(dKickDelta*dLightSpeed);
    if (dxStart > 1) return;
    else if (dxStart < 0) {
	dxStart = (dLookbackFac*dLightSpeed - 2.0)/(dKickDelta*dLightSpeed);
	if (dxStart >= 0) dxStart = 0;
	else {
	    /*
	    ** Check only 64!
	    */
	    nBox = 64;
	    dxStart = (dLookbackFac*dLightSpeed - 1.0)/(dKickDelta*dLightSpeed);
	    if (dxStart >= 0) dxStart = 0;
	    else {
		/*
		** Check only 8!
		*/
		nBox = 8;
		dxStart = 0;
		}
	    }
	}

    const vel_t *v = pkdVel(pkd,p);
    double r0[3],r1[3];
    int j;
    pkdGetPos1(pkd,p,r0);
    for (j=0;j<3;++j) {
	if (r0[j] < -0.5) r0[j] += 1.0;
	else if (r0[j] >= 0.5) r0[j] -= 1.0;
	}

    struct {
	double dt;
	double fOffset;
	int jPlane;
	} isect[4], temp;

    for (j=0;j<3;++j) {
	isect[j].dt = (0.5 - r0[j])/v[j];
	if (isect[j].dt > 0.0) {
	    /*
	    ** Particle crosses the upper j-coordinate boundary of the unit cell at isect[j].dt.
	    */
	    isect[j].fOffset = -1.0;
	    }
	else {
	    /*
	    ** Particle crosses the lower j-coordinate boundary of the unit cell at isect[j].dt.
	    */
	    isect[j].dt = (-0.5 - r0[j])/v[j];
	    isect[j].fOffset = 1.0;
	    }
	isect[j].jPlane = j;
	}
    isect[3].dt = dDriftDelta;
    isect[3].fOffset = 0.0;
    isect[3].jPlane = 3;
    /*
    ** Sort them!
    */
    if (isect[0].dt>isect[1].dt) { temp = isect[0]; isect[0] = isect[1]; isect[1] = temp; }
    if (isect[2].dt>isect[3].dt) { temp = isect[2]; isect[2] = isect[3]; isect[3] = temp; }
    temp = isect[1]; isect[1] = isect[2]; isect[2] = temp;
    if (isect[0].dt>isect[1].dt) { temp = isect[0]; isect[0] = isect[1]; isect[1] = temp; }
    if (isect[2].dt>isect[3].dt) { temp = isect[2]; isect[2] = isect[3]; isect[3] = temp; }
    if (isect[1].dt>isect[2].dt) { temp = isect[1]; isect[1] = isect[2]; isect[2] = temp; }


    dvec xStart = dxStart;


    nBox /= 4; // SIMD width
    int k;
    for (k=0;k<4;++k) {
	double dtApprox, dt;
	dvec dlbt;

	if (k==0) {
	    /*
	    ** Check lightcone from 0 <= dt < isect[k].dt
	    */
	    dt = isect[k].dt;
	    dtApprox = dt/dDriftDelta*dKickDelta;
	    dlbt = dLookbackFac;
	    }
	else {
	    /*
	    ** Check lightcone from isect[k-1].dt <= dt < isect[k].dt
	    */
	    dt = isect[k].dt - isect[k-1].dt;
	    dtApprox = dt/dDriftDelta*dKickDelta;
	    dlbt = dLookbackFac - dtApprox;
	    }
	double t0 = dlbt*dlbt*dLightSpeed*dLightSpeed;
	double t1 = (dlbt - dtApprox)*(dlbt - dtApprox)*dLightSpeed*dLightSpeed;
	for (j=0;j<3;++j) r1[j] = r0[j] + dt*v[j];
	for(int iOct=0; iOct<nBox; ++iOct) {
	    dvec off0, off1, off2;
	    off0.load(pkd->lcOffset0+iOct*4);
	    off1.load(pkd->lcOffset1+iOct*4);
	    off2.load(pkd->lcOffset2+iOct*4);
	    dvec vrx0 = off0 + r0[0];
	    dvec vry0 = off1 + r0[1];
	    dvec vrz0 = off2 + r0[2];
	    dvec vrx1 = off0 + r1[0];
	    dvec vry1 = off1 + r1[1];
	    dvec vrz1 = off2 + r1[2];
	    dvec mr0 = vrx0*vrx0 + vry0*vry0 + vrz0*vrz0;
	    dvec mr1 = vrx1*vrx1 + vry1*vry1 + vrz1*vrz1;
	    int msk = allzero((t1 <= max(mr0,mr1)) & (t0 >= min(mr0,mr1)));
	    if (msk) {
		mr0 = sqrt(mr0);
		mr1 = sqrt(mr1);
		dvec vx = (dLightSpeed*dlbt - mr0)/(dLightSpeed*dtApprox - mr0 + mr1);
		msk = allzero(vx >= xStart & vx < 1.0);
		if (msk) {
		    dvec vr[3];
		    vr[0] = (1-vx)*vrx0 + vx*vrx1;
		    vr[1] = (1-vx)*vry0 + vx*vry1;
		    vr[2] = (1-vx)*vrz0 + vx*vrz1;
		    for(int j=0; j<4; ++j) {
			if (msk & (1<<j)) {
			    double r[3];
			    r[0] = vr[0][j];
			    r[1] = vr[1][j];
			    r[2] = vr[2][j];
			    /*
			    ** Create a new light cone particle.
			    */
			    double mr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
			    addToLightCone(pkd,r,p,pkd->param.bLightConeParticles && (mr <= mrLCP));
			    }
			}
		    }
		}
	    }
	if (isect[k].jPlane == 3) break;
	/*
	** Now we need to reposition r0 to the new segment.
	*/
	for (j=0;j<3;++j) r0[j] = r1[j];
	r0[isect[k].jPlane] += isect[k].fOffset;
	}
    }
