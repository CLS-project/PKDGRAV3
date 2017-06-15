#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include <math.h>
#include <immintrin.h>
#include "pkd.h"
#include "simd.h"

extern "C" void addToLightCone(PKD pkd,double *r,float fPot,PARTICLE *p,int bParticleOutput);

#define NBOX 184

extern "C"
void pkdProcessLightCone(PKD pkd,PARTICLE *p,float fPot,double dLookbackFac,double dLookbackFacLCP,double dDriftDelta,double dKickDelta) {
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
		if (dLookbackFac >= 0) dxStart = 0;
		else return;   /* Nothing to check */
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
    /*
    ** Figure out at what delta t the particle would cross each of the 
    ** walls of the unit cell. Usually this delta t will be larger than
    ** the actual drift, but if not we need to treat each segment of the
    ** drift that is created by crossing a wall (there can be at most 4 
    ** segments in 3-D).
    */
    for (j=0;j<3;++j) {
	if (v[j] > 0) {
	    /*
	    ** Particle crosses the upper j-coordinate boundary of the unit cell at isect[j].dt.
	    */
	    isect[j].dt = (0.5 - r0[j])/v[j];
	    isect[j].fOffset = -1.0;
	    }
	else if (v[j] < 0) {
	    /*
	    ** Particle crosses the lower j-coordinate boundary of the unit cell at isect[j].dt.
	    */
	    isect[j].dt = (-0.5 - r0[j])/v[j];
	    isect[j].fOffset = 1.0;
	    }
	else {
	    /*
	    ** Particle is not moving in this dimension!
	    */
	    isect[j].dt = 2*dDriftDelta; /* this will be ignored after the sort */
	    isect[j].fOffset = 0.0;
	    }
	assert(isect[j].dt >= 0);
	isect[j].jPlane = j;
	}
    /*
    ** Last option is that the particle does not cross any wall so we make an
    ** entry which contains just the drift of the particle.
    */
    isect[3].dt = dDriftDelta;
    isect[3].fOffset = 0.0;
    isect[3].jPlane = 3;
    /*
    ** Sort them by drift time, we will deal with each segment of the drift in this 
    ** order.
    */
    if (isect[0].dt>isect[1].dt) { temp = isect[0]; isect[0] = isect[1]; isect[1] = temp; }
    if (isect[2].dt>isect[3].dt) { temp = isect[2]; isect[2] = isect[3]; isect[3] = temp; }
    temp = isect[1]; isect[1] = isect[2]; isect[2] = temp;
    if (isect[0].dt>isect[1].dt) { temp = isect[0]; isect[0] = isect[1]; isect[1] = temp; }
    if (isect[2].dt>isect[3].dt) { temp = isect[2]; isect[2] = isect[3]; isect[3] = temp; }
    if (isect[1].dt>isect[2].dt) { temp = isect[1]; isect[1] = isect[2]; isect[2] = temp; }


    dvec xStart = dxStart;

    assert(nBox%dvec::width() == 0);
    nBox /= dvec::width();
    int k;
    for (k=0;k<4;++k) {
	double dtApprox, dt;
	dvec dlbt;
	assert(isect[k].dt <= dDriftDelta);
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
	    /*
	    ** continue the light surface from the end of the previous segment
	    */
	    dlbt = dLookbackFac - isect[k-1].dt/dDriftDelta*dKickDelta;
	    }
	dvec t0 = dlbt*dlbt*dLightSpeed*dLightSpeed;
	dvec t1 = (dlbt - dtApprox)*(dlbt - dtApprox)*dLightSpeed*dLightSpeed;
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
	    dmask msk = (t1 <= max(mr0,mr1)) & (t0 >= min(mr0,mr1));
	    if (!testz(msk)) {
		mr0 = sqrt(mr0);
		mr1 = sqrt(mr1);
		dvec vx = (dLightSpeed*dlbt - mr0)/(dLightSpeed*dtApprox - mr0 + mr1);
		msk = (vx >= xStart & vx < 1.0);
		if (!testz(msk)) {
		    dvec vr[3];
		    vr[0] = (1.0-vx)*vrx0 + vx*vrx1;
		    vr[1] = (1.0-vx)*vry0 + vx*vry1;
		    vr[2] = (1.0-vx)*vrz0 + vx*vrz1;
		    int m = movemask(msk);
		    for(int j=0; j<dvec::width(); ++j) {
			if (m & (1<<j)) {
			    double r[3];
			    r[0] = vr[0][j];
			    r[1] = vr[1][j];
			    r[2] = vr[2][j];
			    /*
			    ** Create a new light cone particle.
			    */
			    double mr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
			    addToLightCone(pkd,r,fPot,p,pkd->param.bLightConeParticles && (mr <= mrLCP));
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
