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

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
#else
    #define PRIu64 "llu"
#endif
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "smoothfcn.h"

#ifdef M43D
/* M43D Creates a 3D kernel by convolution of 3D tophats the way M4(1D) is made in 1D */
#define BALL2(fBall) ((fBall)*(fBall))
#define KERNEL(ak,ar2) { \
        ak = sqrt(ar2); \
        if (ar2 < 1.0) ak = 6.*0.25/350./3. *(1360+ar2*(-2880 \
             +ar2*(3528+ak*(-1890+ak*(-240+ak*(270-6*ar2)))))); \
        else if (ar2 < 4.0) ak = 6.*0.25/350./3. *(7040-1152/ak+ak*(-10080+ak*(2880+ak*(4200 \
             +ak*(-3528+ak*(630+ak*(240+ak*(-90+2*ar2)))))))); \
        else ak = 0.0;\
        }
#define DKERNEL(adk,ar2) { \
        adk = sqrt(ar2); \
        if (ar2 < 1.0) adk = 6.*0.25/350./3. * (-2880*2 \
             +ar2*(3528*4+ adk*(-1890*5 + adk*(-240*6+ adk*(270*7-6*9*ar2))))); \
        else if (ar2 < 4.0) adk = 6.*0.25/350./3. *((1152/ar2-10080)/adk+(2880*2+adk*(4200*3 \
             +adk*(-3528*4+adk*(630*5+adk*(240*6 +adk*(-90*7+2*9*ar2))))))); \
        else adk = 0.0;\
        }

#else
#ifdef HSHRINK
/* HSHRINK M4 Kernel uses an effective h of (pi/6)^(1/3) times h for nSmooth neighbours */
#define dSHRINKFACTOR        0.805995977
#define BALL2(fBall) ((fBall)*(fBall)*(dSHRINKFACTOR*dSHRINKFACTOR))
#define KERNEL(ak,ar2) { \
        ak = 2.0 - sqrt(ar2); \
        if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
        else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
        else ak = 0.0; \
        }
#define DKERNEL(adk,ar2) { \
        adk = sqrt(ar2); \
        if (ar2 < 1.0) { \
            adk = -3 + 2.25*adk; \
            } \
        else if (ar2 < 4.0) { \
            adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
            } \
        else adk = 0.0; \
        }

#else
/* Standard M_4 Kernel */
#define BALL2(fBall) ((fBall)*(fBall))
#define KERNEL(ak,ar2) { \
        ak = 2.0 - sqrt(ar2); \
        if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2); \
        else if (ar2 < 4.0) ak = 0.25*ak*ak*ak; \
        else ak = 0.0;\
        }
#define DKERNEL(adk,ar2) { \
        adk = sqrt(ar2); \
        if (ar2 < 1.0) { \
            adk = -3 + 2.25*adk; \
            } \
        else if (ar2 < 4.0) { \
            adk = -0.75*(2.0-adk)*(2.0-adk)/adk; \
            } \
        else adk = 0.0;\
        }
#endif
#endif

void NullSmooth(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
}

void initBall(void *vpkd, void *p) {
    auto pkd = static_cast<PKD>(vpkd);
    pkdSetBall(pkd,(PARTICLE *)p,0.0);
}

void BallSmooth(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    pkdSetBall(pkd,p,fBall);
}

void initDensity(void *vpkd, void *p) {
    auto pkd = static_cast<PKD>(vpkd);
    pkdSetDensity(pkd,(PARTICLE *)p,0.0);
}

void combDensity(void *vpkd, void *p1,const void *p2) {
    auto pkd = static_cast<PKD>(vpkd);
    pkdSetDensity(pkd,(PARTICLE *)p1,pkdDensity(pkd,(PARTICLE *)p1)+pkdDensity(pkd,(PARTICLE *)p2));
}

void DensityF1(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    double ih2,r2,rs,fDensity,fMass;
    int i;

    ih2 = 1.0/BALL2(fBall);
    fDensity = 0.0;
    for (i=0; i<nSmooth; ++i) {
        fMass = pkdMass(pkd,nnList[i].pPart);
        r2 = nnList[i].fDist2*ih2;
        rs = 1 - r2;
        if (rs < 0) rs = 0.0;
        fDensity += rs*fMass;
    }
    fDensity *= 1.875f*M_1_PI*sqrtf(ih2)*ih2; /* F1 Kernel (15/8) */
    if (smf->pfDensity) *smf->pfDensity = fDensity;
    else pkdSetDensity(pkd,p,fDensity);
}

void DensityM3(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    double ih2,r2,rs,fDensity,fMass;
    int i;
    ih2 = 1.0f/BALL2(fBall);
    fDensity = 0.0;
    for (i=0; i<nSmooth; ++i) {
        fMass = pkdMass(pkd,nnList[i].pPart);
        r2 = nnList[i].fDist2*ih2;
        if (r2 < 1.0) {
            double r = sqrt(r2);
            rs = 1.0f - r;
            rs *= rs*rs; /* rs^3 */
            if (r < 0.5f) {
                double rs2 = 0.5f - r;
                rs2 *= rs2*rs2; /* rs2^3 */
                rs -= 4.0f*rs2;
            }
        }
        else rs = 0.0;
        fDensity += rs*fMass;
    }
    pkdSetDensity(pkd,p,16.0f*M_1_PI*sqrtf(ih2)*ih2*fDensity);
}

void LinkGradientM3(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    double ih2,r2,rs,fMass,fNorm,frho[3], idrho, r2min, dr;
    int i, j;
    ih2 = 1.0/BALL2(fBall);
    fNorm = 16.0f*M_1_PI*ih2*ih2*sqrtf(ih2);
    frho[0] = frho[1] = frho[2] = 0.0;
    for (i=0; i<nSmooth; ++i) {
        fMass = pkdMass(pkd,nnList[i].pPart);
        r2 = nnList[i].fDist2*ih2;
        if (r2 < 1.0) {
            double r = sqrt(r2);
            if (r < 0.5f) {
                rs = -6.0f + 9.0f*r;
            }
            else {
                rs = 1.0f - r;
                rs *= rs; /* rs^2 */
                rs *= -3.0f;
                rs /= r;
            }
        }
        else rs = 0.0;
        rs *= fNorm*fMass;
        rs *= (pkdDensity(pkd,nnList[i].pPart) - pkdDensity(pkd,p))/pkdDensity(pkd,nnList[i].pPart);
        frho[0] -= nnList[i].dx*rs;
        frho[1] -= nnList[i].dy*rs;
        frho[2] -= nnList[i].dz*rs;
    }
    idrho = 1.0/sqrt(frho[0]*frho[0] + frho[1]*frho[1] + frho[2]*frho[2]);
    for (j=0; j<3; ++j) frho[j] *= 0.5*idrho*fBall;
    r2min = HUGE_VALF;
    if (nSmooth==0) pkdSetGroup(pkd, p, -1);
    for (i=0; i<nSmooth; ++i) {
        dr = nnList[i].dx - frho[0];
        r2 = dr*dr;
        dr = nnList[i].dy - frho[1];
        r2 += dr*dr;
        dr = nnList[i].dz - frho[2];
        r2 += dr*dr;
        if (r2 < r2min) {
            r2min = r2;
            smf->hopParticleLink.iPid = nnList[i].iPid;
            smf->hopParticleLink.iIndex = nnList[i].iIndex;
        }
    }
}

void LinkHopChains(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    MDL mdl = pkd->mdl;
    int i, gid1, gid2;
    GHtmpGroupTable *g1, *g2, g;
    gid1 = pkdGetGroup(pkd,p);
    g1 = &pkd->tmpHopGroups[gid1];
    for (i=0; i<nSmooth; ++i) {
        gid2 = pkdGetGroup(pkd,nnList[i].pPart);
        if (nnList[i].iPid==pkd->Self() && gid1==gid2) continue;
        g.iPid = nnList[i].iPid;
        g.iIndex = gid2;
        g2 = static_cast<GHtmpGroupTable *>(mdlAcquire(mdl,CID_GROUP,g.iIndex,g.iPid));

        /* Remote is authoratative. Update myself, but also what I currently link to. */
        if (g1->iPid > g2->iPid || (g1->iPid == g2->iPid && g1->iIndex > g2->iIndex)) {
            smf->bDone = 0;
            g = *g1;
            g1->iPid = g2->iPid;
            g1->iIndex = g2->iIndex;
            mdlRelease(mdl,CID_GROUP,g2);
            g2 = static_cast<GHtmpGroupTable *>(mdlAcquire(mdl,CID_GROUP,g.iIndex,g.iPid));
        }

        /* Update remote (or what we were pointing to) and what it points to if necessary. */
        while (g1->iPid < g2->iPid || (g1->iPid == g2->iPid && g1->iIndex < g2->iIndex) ) {
            smf->bDone = 0;
            g = *g2;
            g2->iPid = g1->iPid;
            g2->iIndex = g1->iIndex;
            mdlRelease(mdl,CID_GROUP,g2);
            g2 = static_cast<GHtmpGroupTable *>(mdlAcquire(mdl,CID_GROUP,g.iIndex,g.iPid));
        }
        mdlRelease(mdl,CID_GROUP,g2);
    }
}

void Density(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    double ih2,r2,rs,fDensity,fMass;
    int i;

    ih2 = 4.0/BALL2(fBall);
    fDensity = 0.0;
    for (i=0; i<nSmooth; ++i) {
        fMass = pkdMass(pkd,nnList[i].pPart);
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        fDensity += rs*fMass;
    }
    pkdSetDensity(pkd,p,M_1_PI*sqrt(ih2)*ih2*fDensity);
}

void DensitySym(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    double fNorm,ih2,r2,rs,fMassQ,fMassP;
    int i;
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(fBall));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        rs *= fNorm;
        q = nnList[i].pPart;
        fMassQ = pkdMass(pkd,q);
        pkdSetDensity(pkd,p,pkdDensity(pkd,p) + rs*fMassQ);
        pkdSetDensity(pkd,q,pkdDensity(pkd,q) + rs*fMassP);
    }
}

void PrintNN(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    int i;

    printf("%" PRIu64 ":",(uint64_t)p->iOrder);
    for (i=0; i<nSmooth; ++i) {
        if (pkdIsActive(pkd,nnList[i].pPart))
            printf("%" PRIu64 " ",(uint64_t)nnList[i].pPart->iOrder);
        else
            printf("\033[7m%" PRIu64 "\033[0m ",(uint64_t)nnList[i].pPart->iOrder);
    }
    printf("\n");
}

void initDenDVDX(void *vpkd, void *p) {
}

void combDenDVDX(void *vpkd, void *p1,const void *p2) {
}

/* Gather only version */
/* JW: What types should dx etc... have -- why is NN using FLOAT ? */
#ifndef OPTIM_REMOVE_UNUSED
void DenDVDX(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double ih2,ih,r2,rs,rs1,fDensity,qMass,fNorm,fNorm1,vFac;
    double dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
    double dvx,dvy,dvz,dx,dy,dz,trace;
    double curlv[3];
    int i;

    assert(pkd->oFieldOffset[oSph]);
    psph = pkdSph(pkd,p);
    ih = 2.0/fBall;  ih2 = ih*ih;
    vFac = (smf->bComove ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot (physical) */
    fNorm = M_1_PI*ih2*ih;
    fNorm1 = fNorm*ih2;
    fDensity = 0.0;
    dvxdx = 0; dvxdy = 0; dvxdz= 0;
    dvydx = 0; dvydy = 0; dvydz= 0;
    dvzdx = 0; dvzdy = 0; dvzdz= 0;

    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        q = nnList[i].pPart;
        qMass = pkdMass(pkd,q);
        qsph = pkdSph(pkd,q);
        KERNEL(rs,r2);
        fDensity += rs*qMass;
        DKERNEL(rs1,r2);
        rs1 *= qMass;
        dx = nnList[i].dx;
        dy = nnList[i].dy;
        dz = nnList[i].dz;
        dvx = (-psph->vPred[0] + qsph->vPred[0])*vFac - dx*smf->H; /* NB: dx = px - qx */
        dvy = (-psph->vPred[1] + qsph->vPred[1])*vFac - dy*smf->H;
        dvz = (-psph->vPred[2] + qsph->vPred[2])*vFac - dz*smf->H;
        dvxdx += dvx*dx*rs1;
        dvxdy += dvx*dy*rs1;
        dvxdz += dvx*dz*rs1;
        dvydx += dvy*dx*rs1;
        dvydy += dvy*dy*rs1;
        dvydz += dvy*dz*rs1;
        dvzdx += dvz*dx*rs1;
        dvzdy += dvz*dy*rs1;
        dvzdz += dvz*dz*rs1;
    }
    fDensity*=fNorm;
    pkdSetDensity(pkd,p,fDensity);
    psph->c = sqrt(smf->gamma*(smf->gamma-1)*psph->uPred);
    fNorm1 /= fDensity;
    trace = dvxdx+dvydy+dvzdz;
    psph->divv =  fNorm1*trace; /* physical */

    psph->BalsaraSwitch=1;
    if (smf->iViscosityLimiter) {
        if (psph->divv!=0.0) {
            curlv[0] = fNorm1*(dvzdy - dvydz);
            curlv[1] = fNorm1*(dvxdz - dvzdx);
            curlv[2] = fNorm1*(dvydx - dvxdy);
            psph->BalsaraSwitch = fabs(psph->divv)/
                                  (fabs(psph->divv)+sqrt(curlv[0]*curlv[0]+
                                          curlv[1]*curlv[1]+
                                          curlv[2]*curlv[2]));
        }
        else {
            psph->BalsaraSwitch = 0;
        }

    }

    if (smf->iDiffusion) {
        double onethirdtrace = (1./3.)*trace;
        /* Build Traceless Strain Tensor (not yet normalized) */
        double sxx = dvxdx - onethirdtrace; /* pure compression/expansion doesn't diffuse */
        double syy = dvydy - onethirdtrace;
        double szz = dvzdz - onethirdtrace;
        double sxy = 0.5*(dvxdy + dvydx); /* pure rotation doesn't diffuse */
        double sxz = 0.5*(dvxdz + dvzdx);
        double syz = 0.5*(dvydz + dvzdy);
        /* diff coeff., nu ~ C L^2 S (add C via dMetalDiffusionConstant, assume L ~ h) */
        if (smf->iDiffusion == 2) psph->diff = 1;
        else psph->diff = fNorm1*fBall*fBall*sqrt(2*(sxx*sxx + syy*syy + szz*szz + 2*(sxy*sxy + sxz*sxz + syz*syz)));
    }
    else psph->diff = 0;

}

/* JW: Do I need to differentiate init of Original Particle and Cached Copy?
   -- YES accel zeroed on cache copy */
void initSphForcesParticle(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    if (pkdIsActive(pkd,p)) {
        SPHFIELDS *psph = pkdSph(pkd,p);
        psph->uDot = 0;
        psph->fMetalsDot = 0;
        assert(0); /* What happens here? */
        // if (!pkd->param.bDoGravity) { /* Normally these are zeroed in Gravity */
        //     p->uNewRung = 0;
        //     pkdAccel(pkd,p)[0] = 0;
        //     pkdAccel(pkd,p)[1] = 0;
        //     pkdAccel(pkd,p)[2] = 0;
        //     }
    }
}

void initSphForces(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    assert(!pkd->bNoParticleOrder);
    if (pkdIsActive(pkd,p)) {
        SPHFIELDS *psph = pkdSph(pkd,p);
        psph->uDot = 0;
        psph->fMetalsDot = 0;
        p->uNewRung = 0;
        pkdAccel(pkd,p)[0] = 0;
        pkdAccel(pkd,p)[1] = 0;
        pkdAccel(pkd,p)[2] = 0; /* JW: Cached copies have zero accel! && rung */
    }
}

void combSphForces(void *vpkd, void *p1,const void *p2) {
    PKD pkd = (PKD) vpkd;
    assert(!pkd->bNoParticleOrder);
    if (pkdIsActive(pkd,p1)) {
        SPHFIELDS *psph1 = pkdSph(pkd,p1);
        const SPHFIELDS *psph2 = pkdSphRO(pkd,p2);
        float *a1 = pkdAccel(pkd,p1);
        const float *a2 = pkdAccelRO(pkd,p2);
        psph1->uDot += psph2->uDot;
        psph1->fMetalsDot += psph2->fMetalsDot;
        a1[0] += a2[0];
        a1[1] += a2[1];
        a1[2] += a2[2];
        if (((PARTICLE *) p2)->uNewRung > ((PARTICLE *) p1)->uNewRung)
            ((PARTICLE *) p1)->uNewRung = ((PARTICLE *) p2)->uNewRung;
    }
}

void SphForces(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {

    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double ih2,r2,rs1,rq,rp;
    double dx,dy,dz,dvx,dvy,dvz,dvdotdr;
    double pPoverRho2,pPoverRho2f,pMass;
    double qPoverRho2,qPoverRho2f;
    double ph,pc,pDensity,visc,hav,absmu,Accp,Accq;
    double fNorm,fNorm1,aFac,vFac,gammainv,dtC,dtMu,dtEst,dt2;
    float *pa,*qa;
    int i,pActive,qActive;
    uint8_t uNewRung;

    assert(!pkd->bNoParticleOrder);

    aFac = (smf->a);        /* comoving acceleration factor */
    vFac = (smf->bComove ? 1./(smf->a*smf->a) : 1.0); /* converts v to xdot */
    gammainv = 1/smf->gamma;
    dtC = (1+0.6*smf->alpha)/(smf->a*smf->dEtaCourant);
    dtMu = (0.6*smf->beta)/(smf->a*smf->dEtaCourant);

    psph = pkdSph(pkd,p);
    pc = psph->c;
    pDensity = pkdDensity(pkd,p);
    pMass = pkdMass(pkd,p);
    pPoverRho2 = gammainv*psph->c*psph->c/pDensity;
    pPoverRho2f = pPoverRho2;
    ph = 0.5*fBall;
    /* JW: Active tests here -- Rung info in pkd */
    pActive = pkdIsActive(pkd,p);
    pa = pkdAccel(pkd,p);

    ih2 = 1/(ph*ph);
    fNorm = 0.5*M_1_PI*ih2/ph;
    fNorm1 = fNorm*ih2; /* converts to physical u */

    for (i=0; i<nSmooth; ++i) {
        q = nnList[i].pPart;
        /* JW: Active tests here -- Rung info in pkd */
        qActive = pkdIsActive(pkd,q);
        if (!pActive && !qActive) continue;
        qsph = pkdSph(pkd,q);
        qa = pkdAccel(pkd,q);

        r2 = nnList[i].fDist2*ih2;
        DKERNEL(rs1,r2);
        rs1 *= fNorm1;
        rp = rs1 * pMass;
        rq = rs1 * pkdMass(pkd,q);

        dx = nnList[i].dx;
        dy = nnList[i].dy;
        dz = nnList[i].dz;
        dvx = psph->vPred[0] - qsph->vPred[0];
        dvy = psph->vPred[1] - qsph->vPred[1];
        dvz = psph->vPred[2] - qsph->vPred[2];
        dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz)
                  + nnList[i].fDist2*smf->H;

        qPoverRho2 = gammainv*qsph->c*qsph->c/pkdDensity(pkd,q);
        qPoverRho2f = qPoverRho2;

#define DIFFUSIONThermal() \
    { double diff = 2*smf->dThermalDiffusionCoeff*(psph->diff+qsph->diff)*(psph->uPred-qsph->uPred) \
        /(pkdDensity(pkd,p)+pkdDensity(pkd,q));     \
      PACTIVE( psph->uDot += diff*rq ); \
      QACTIVE( qsph->uDot -= diff*rp ); }

#define DIFFUSIONMetals() \
    { double diff = 2*smf->dMetalDiffusionCoeff*(psph->diff+qsph->diff)*(psph->fMetals - qsph->fMetals) \
        /(pkdDensity(pkd,p)+pkdDensity(pkd,q));     \
      PACTIVE( psph->fMetalsDot += diff*rq ); \
      QACTIVE( qsph->fMetalsDot -= diff*rp ); }

        /* JW: Star form will need this! */
#define DIFFUSIONMetalsOxygen()
#define DIFFUSIONMetalsIron()
#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch+b->BalsaraSwitch))
#define ALPHA (smf->alpha)
#define BETA (smf->beta)
#define PRES_PDV(a,b) (a)
#define PRES_ACC(a,b) (a+b)

#define SphForcesACTIVECODE() \
        if (dvdotdr>0.0) { \
        PACTIVE( psph->uDot += rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdr; ); \
        QACTIVE( qsph->uDot += rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdr; ); \
        PACTIVE( Accp = (PRES_ACC(pPoverRho2f,qPoverRho2f)); ); \
        QACTIVE( Accq = (PRES_ACC(qPoverRho2f,pPoverRho2f)); ); \
        absmu = 0; \
        } \
        else {  \
        hav=0.5*(ph+0.5*pkdBall(pkd,q));  /* h mean */ \
        absmu = -hav*dvdotdr*smf->a  \
            /(nnList[i].fDist2+0.01*hav*hav); /* mu multiply by a to be consistent with physical c */ \
        visc = SWITCHCOMBINE(psph,qsph)* \
            (ALPHA*(pc + qsph->c) + BETA*2*absmu)  \
            *absmu/(pDensity + pkdDensity(pkd,q));  \
        PACTIVE( psph->uDot += rq*(PRES_PDV(pPoverRho2,qPoverRho2) + 0.5*visc)*dvdotdr; ); \
        QACTIVE( qsph->uDot += rp*(PRES_PDV(qPoverRho2,pPoverRho2) + 0.5*visc)*dvdotdr; ); \
        PACTIVE( Accp = (PRES_ACC(pPoverRho2f,qPoverRho2f) + visc); ); \
        QACTIVE( Accq = (PRES_ACC(qPoverRho2f,pPoverRho2f) + visc); ); \
        } \
        dtEst = ph/(dtC*psph->c+dtMu*absmu);    \
        dt2 = 0.5*pkdBall(pkd,q)/(dtC*qsph->c+dtMu*absmu); \
        if (dt2 < dtEst) dtEst=dt2; \
        uNewRung = pkdDtToRung(dtEst,smf->dDelta,MAX_RUNG); \
            PACTIVE( if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung; );    \
            QACTIVE( if (uNewRung > q->uNewRung ) q->uNewRung = uNewRung; );    \
        PACTIVE( Accp *= rq*aFac; );/* aFac - convert to comoving acceleration */ \
        QACTIVE( Accq *= rp*aFac; ); \
        PACTIVE( pa[0] -= Accp * dx; ); \
        PACTIVE( pa[1] -= Accp * dy; ); \
        PACTIVE( pa[2] -= Accp * dz; ); \
        QACTIVE( qa[0] += Accq * dx; ); \
        QACTIVE( qa[1] += Accq * dy; ); \
        QACTIVE( qa[2] += Accq * dz; ); \
            DIFFUSIONThermal(); \
            DIFFUSIONMetals(); \
            DIFFUSIONMetalsOxygen(); \
            DIFFUSIONMetalsIron();

        if (pActive) {
            if (qActive) {
#define PACTIVE(xxx) xxx
#define QACTIVE(xxx) xxx
                SphForcesACTIVECODE();
            }
            else {
#undef QACTIVE
#define QACTIVE(xxx)
                SphForcesACTIVECODE();
            }
        }
        else if (qActive) {
#undef PACTIVE
#define PACTIVE(xxx)
#undef QACTIVE
#define QACTIVE(xxx) xxx
            SphForcesACTIVECODE();
        }

    } /* end neighbour loop */
}

void initDistDeletedGas(void *vpkd,void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    /*
    * Zero out accumulated quantities.
    */
    if (pkdIsDeleted(pkd,p)) return; /* deleted */

    *((float *) pkdField(p,pkd->oFieldOffset[oMass]))=0;
    pkdSetPos(pkd,p,0,0);
    pkdSetPos(pkd,p,1,0);
    pkdSetPos(pkd,p,2,0);
    pkdVel(pkd,p)[0] = 0;
    pkdVel(pkd,p)[1] = 0;
    pkdVel(pkd,p)[2] = 0;
    pkd_vPred(pkd,p)[0] = 0;
    pkd_vPred(pkd,p)[1] = 0;
    pkd_vPred(pkd,p)[2] = 0;
    pkdAccel(pkd,p)[0] = 0;
    pkdAccel(pkd,p)[1] = 0;
    pkdAccel(pkd,p)[2] = 0;
    pkdSph(pkd,p)->u = 0;
    pkdSph(pkd,p)->uPred = 0;
    pkdSph(pkd,p)->uDot = 0;
    pkdSph(pkd,p)->fMetals = 0;
    pkdSph(pkd,p)->fMetalsDot = 0;
    pkdSph(pkd,p)->fMetalsPred = 0;
}


void combDistDeletedGas(void *vpkd,void *vp1,const void *vp2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = vp1,*p2 = (void *)vp2;
    float *p1mass;
    const float *p2mass;
    double f1,f2,m;
    /*
     * Distribute u, v, and fMetals for particles returning from cache
     * so that everything is conserved nicely.
     */
    if (pkdIsDeleted(pkd,p1)) return; /* deleted */

    p2mass = pkdFieldRO(p2,pkd->oFieldOffset[oMass]);
    if (*p2mass > 0) {
        p1mass = pkdField(p1,pkd->oFieldOffset[oMass]);
        m = (*p1mass+*p2mass);
        f1=(*p1mass)/m;
        f2=(*p2mass)/m;
        *p1mass = m;
        pkdSetPos(pkd,p1,0,f1*pkdPos(pkd,p1,0)+f2*pkdPos(pkd,p2,0));
        pkdSetPos(pkd,p1,1,f1*pkdPos(pkd,p1,1)+f2*pkdPos(pkd,p2,1));
        pkdSetPos(pkd,p1,2,f1*pkdPos(pkd,p1,2)+f2*pkdPos(pkd,p2,2));
        pkdVel(pkd,p1)[0] = f1*pkdVel(pkd,p1)[0]+f2*pkdVel(pkd,p2)[0];
        pkdVel(pkd,p1)[1] = f1*pkdVel(pkd,p1)[1]+f2*pkdVel(pkd,p2)[1];
        pkdVel(pkd,p1)[2] = f1*pkdVel(pkd,p1)[2]+f2*pkdVel(pkd,p2)[2];
        pkd_vPred(pkd,p1)[0] = f1*pkd_vPred(pkd,p1)[0]+f2*pkd_vPred(pkd,p2)[0];
        pkd_vPred(pkd,p1)[1] = f1*pkd_vPred(pkd,p1)[1]+f2*pkd_vPred(pkd,p2)[1];
        pkd_vPred(pkd,p1)[2] = f1*pkd_vPred(pkd,p1)[2]+f2*pkd_vPred(pkd,p2)[2];
        pkdAccel(pkd,p1)[0] = f1*pkdAccel(pkd,p1)[0]+f2*pkdAccel(pkd,p2)[0];
        pkdAccel(pkd,p1)[1] = f1*pkdAccel(pkd,p1)[1]+f2*pkdAccel(pkd,p2)[1];
        pkdAccel(pkd,p1)[2] = f1*pkdAccel(pkd,p1)[2]+f2*pkdAccel(pkd,p2)[2];
        pkdSph(pkd,p1)->u = f1*pkdSph(pkd,p1)->u+f2*pkdSph(pkd,p2)->u;
        pkdSph(pkd,p1)->uPred = f1*pkdSph(pkd,p1)->uPred+f2*pkdSph(pkd,p2)->uPred;
        pkdSph(pkd,p1)->uDot = f1*pkdSph(pkd,p1)->uDot+f2*pkdSph(pkd,p2)->uDot;
        pkdSph(pkd,p1)->fMetals = f1*pkdSph(pkd,p1)->fMetals+f2*pkdSph(pkd,p2)->fMetals;
        pkdSph(pkd,p1)->fMetalsPred = f1*pkdSph(pkd,p1)->fMetalsPred+f2*pkdSph(pkd,p2)->fMetalsPred;
        pkdSph(pkd,p1)->fMetalsDot = f1*pkdSph(pkd,p1)->fMetalsDot+f2*pkdSph(pkd,p2)->fMetalsDot;
    }
}

void DistDeletedGas(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = (PKD) smf->pkd;
    PARTICLE *q;
    double fNorm,ih2,r2,rs,rstot,fp,fq;
    float *pmass, *qmass;
    double m,delta_m;
    int i;

    /* JW: Assert deleted? -- yes once smooth is reliable */
    if (!pkdIsDeleted( pkd, p )) return; /* not deleted */

    pmass = pkdField(p,pkd->oFieldOffset[oMass]);
    if (*pmass <= 0) return;

    ih2 = 4.0/(fBall*fBall);
    rstot = 0;
    for (i=0; i<nSmooth; ++i) {
        q = nnList[i].pPart;
        if (pkdIsDeleted(pkd,q)) continue; /* deleted */
        assert(pkdIsGas(pkd,q));
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        rstot += rs;
    }
    assert(rstot > 0); /* What if all neighbours deleted -- does that happen? */

    fNorm = *pmass/rstot;
    for (i=0; i<nSmooth; ++i) {
        q = nnList[i].pPart;
        if (pkdIsDeleted(pkd,q)) continue; /* deleted */
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);

        delta_m = rs*fNorm;
        qmass = pkdField(q,pkd->oFieldOffset[oMass]);
        m = *qmass + delta_m;
        fp = delta_m/m;
        fq = (*qmass)/m;
        *qmass = m;

        pkdSetPos(pkd,q,0,fq*pkdPos(pkd,q,0)+fp*pkdPos(pkd,p,0));
        pkdSetPos(pkd,q,1,fq*pkdPos(pkd,q,1)+fp*pkdPos(pkd,p,1));
        pkdSetPos(pkd,q,2,fq*pkdPos(pkd,q,2)+fp*pkdPos(pkd,p,2));
        pkdVel(pkd,q)[0] = fq*pkdVel(pkd,q)[0]+fp*pkdVel(pkd,p)[0];
        pkdVel(pkd,q)[1] = fq*pkdVel(pkd,q)[1]+fp*pkdVel(pkd,p)[1];
        pkdVel(pkd,q)[2] = fq*pkdVel(pkd,q)[2]+fp*pkdVel(pkd,p)[2];
        pkd_vPred(pkd,q)[0] = fq*pkd_vPred(pkd,q)[0]+fp*pkd_vPred(pkd,p)[0];
        pkd_vPred(pkd,q)[1] = fq*pkd_vPred(pkd,q)[1]+fp*pkd_vPred(pkd,p)[1];
        pkd_vPred(pkd,q)[2] = fq*pkd_vPred(pkd,q)[2]+fp*pkd_vPred(pkd,p)[2];
        pkdAccel(pkd,q)[0] = fq*pkdAccel(pkd,q)[0]+fp*pkdAccel(pkd,p)[0];
        pkdAccel(pkd,q)[1] = fq*pkdAccel(pkd,q)[1]+fp*pkdAccel(pkd,p)[1];
        pkdAccel(pkd,q)[2] = fq*pkdAccel(pkd,q)[2]+fp*pkdAccel(pkd,p)[2];
        pkdSph(pkd,q)->u = fq*pkdSph(pkd,q)->u+fp*pkdSph(pkd,p)->u;
        pkdSph(pkd,q)->uPred = fq*pkdSph(pkd,q)->uPred+fp*pkdSph(pkd,p)->uPred;
        pkdSph(pkd,q)->uDot = fq*pkdSph(pkd,q)->uDot+fp*pkdSph(pkd,p)->uDot;
        pkdSph(pkd,q)->fMetals = fq*pkdSph(pkd,q)->fMetals+fp*pkdSph(pkd,p)->fMetals;
        pkdSph(pkd,q)->fMetalsPred = fq*pkdSph(pkd,q)->fMetalsPred+fp*pkdSph(pkd,p)->fMetalsPred;
        pkdSph(pkd,q)->fMetalsDot = fq*pkdSph(pkd,q)->fMetalsDot+fp*pkdSph(pkd,p)->fMetalsDot;
    }
    *pmass = 0; /* All distributed */
}

void initDistSNEnergy(void *vpkd,void *vp) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = vp;
    /*
    * Zero out accumulated quantities.
    */
    if (!pkdIsGas(pkd,p)) return; /* not gas */

    *((float *) pkdField(p,pkd->oFieldOffset[oMass]))=0;
    pkdSetPos(pkd,p,0,0);
    pkdSetPos(pkd,p,1,0);
    pkdSetPos(pkd,p,2,0);
    pkdVel(pkd,p)[0] = 0;
    pkdVel(pkd,p)[1] = 0;
    pkdVel(pkd,p)[2] = 0;
    pkd_vPred(pkd,p)[0] = 0;
    pkd_vPred(pkd,p)[1] = 0;
    pkd_vPred(pkd,p)[2] = 0;
    pkdAccel(pkd,p)[0] = 0;
    pkdAccel(pkd,p)[1] = 0;
    pkdAccel(pkd,p)[2] = 0;
    pkdSph(pkd,p)->u = 0;
    pkdSph(pkd,p)->uPred = 0;
    pkdSph(pkd,p)->uDot = 0;
    pkdSph(pkd,p)->fMetals = 0;
    pkdSph(pkd,p)->fMetalsDot = 0;
    pkdSph(pkd,p)->fMetalsPred = 0;
}

void combDistSNEnergy(void *vpkd,void *vp1,const void *vp2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = vp1,*p2=(void *)vp2;
    float *p1mass, *p2mass;
    double f1,f2,m;
    /*
     * Distribute u, v, and fMetals for particles returning from cache
     * so that everything is conserved nicely.
     */
    if (!pkdIsGas(pkd,p1)) return; /* not gas */

    p2mass = pkdField(p2,pkd->oFieldOffset[oMass]);
    if (*p2mass > 0) {
        p1mass = pkdField(p1,pkd->oFieldOffset[oMass]);
        m = (*p1mass+*p2mass);
        f1=(*p1mass)/m;
        f2=(*p2mass)/m;
        *p1mass = m;
        pkdSetPos(pkd,p1,0,f1*pkdPos(pkd,p1,0)+f2*pkdPos(pkd,p2,0));
        pkdSetPos(pkd,p1,1,f1*pkdPos(pkd,p1,1)+f2*pkdPos(pkd,p2,1));
        pkdSetPos(pkd,p1,2,f1*pkdPos(pkd,p1,2)+f2*pkdPos(pkd,p2,2));
        pkdVel(pkd,p1)[0] = f1*pkdVel(pkd,p1)[0]+f2*pkdVel(pkd,p2)[0];
        pkdVel(pkd,p1)[1] = f1*pkdVel(pkd,p1)[1]+f2*pkdVel(pkd,p2)[1];
        pkdVel(pkd,p1)[2] = f1*pkdVel(pkd,p1)[2]+f2*pkdVel(pkd,p2)[2];
        pkd_vPred(pkd,p1)[0] = f1*pkd_vPred(pkd,p1)[0]+f2*pkd_vPred(pkd,p2)[0];
        pkd_vPred(pkd,p1)[1] = f1*pkd_vPred(pkd,p1)[1]+f2*pkd_vPred(pkd,p2)[1];
        pkd_vPred(pkd,p1)[2] = f1*pkd_vPred(pkd,p1)[2]+f2*pkd_vPred(pkd,p2)[2];
        pkdAccel(pkd,p1)[0] = f1*pkdAccel(pkd,p1)[0]+f2*pkdAccel(pkd,p2)[0];
        pkdAccel(pkd,p1)[1] = f1*pkdAccel(pkd,p1)[1]+f2*pkdAccel(pkd,p2)[1];
        pkdAccel(pkd,p1)[2] = f1*pkdAccel(pkd,p1)[2]+f2*pkdAccel(pkd,p2)[2];
        pkdSph(pkd,p1)->u = f1*pkdSph(pkd,p1)->u+f2*pkdSph(pkd,p2)->u;
        pkdSph(pkd,p1)->uPred = f1*pkdSph(pkd,p1)->uPred+f2*pkdSph(pkd,p2)->uPred;
        pkdSph(pkd,p1)->uDot = f1*pkdSph(pkd,p1)->uDot+f2*pkdSph(pkd,p2)->uDot;
        pkdSph(pkd,p1)->fMetals = f1*pkdSph(pkd,p1)->fMetals+f2*pkdSph(pkd,p2)->fMetals;
        pkdSph(pkd,p1)->fMetalsPred = f1*pkdSph(pkd,p1)->fMetalsPred+f2*pkdSph(pkd,p2)->fMetalsPred;
        pkdSph(pkd,p1)->fMetalsDot = f1*pkdSph(pkd,p1)->fMetalsDot+f2*pkdSph(pkd,p2)->fMetalsDot;
    }
}

void DistSNEnergy(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = (PKD) smf->pkd;
    PARTICLE *q,*qmin;
    float qBall;
    double ih2,r2,r2min,fp,fq;
    float *pmass, *qmass;
    double m,im,delta_m,delta_u,delta_Z;
    double dtp,dt,c,dtC,dtNew,Timer;
    int i,uNewRung;

    assert(!pkd->bNoParticleOrder);

    if (!pkdIsStar( pkd, p )) return; /* not a star */

    Timer = *pkd_Timer(pkd,p);
    dtp = smf->dDelta/(1<<p->uRung);/* size of step just completed by star */

    if (smf->dTime-dtp > Timer+smf->SFdtFeedbackDelay) {
        assert(pkdSph(pkd,p)->u == 0); /* Star must have fedback by now */
        return;
    }

    /* OLD: sign of Timer approach    if (Timer > 0) return; */
    assert(pkdSph(pkd,p)->u == 1); /* Star must not have fedback */

    /* Keep star timestep of order Courant time in FB region */
    dtNew=smf->SFdFBFac*fBall;
    uNewRung = pkdDtToRung(dtNew,smf->dDelta,smf->iMaxRung-1);
    if (uNewRung > p->uNewRung) p->uNewRung = uNewRung;

    /* Has it gone off yet ? */
    if (smf->dTime <= Timer + smf->SFdtFeedbackDelay) {
        /* No -- but attempt to warn particles the SN is coming... by lowering timesteps */
        /* Star dt should be suppressed in right way -- see SFdvFB  */
        uNewRung = p->uNewRung;
        for (i=0; i<nSmooth; ++i) {
            q = nnList[i].pPart;
            if (!pkdIsGas(pkd,q)) continue;
            if (pkdIsDeleted(pkd,q)) continue; /* deleted */
            if (uNewRung > q->uNewRung) q->uNewRung = uNewRung;
        }
        return; /* Now exit */
    }

    dtC = (1+0.6*smf->alpha)/(smf->a*smf->dEtaCourant);
    pmass = pkdField(p,pkd->oFieldOffset[oMass]);

    ih2 = 4.0/(fBall*fBall);
    r2min = 1e37f;
    qmin = NULL;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        q = nnList[i].pPart;
        if (!pkdIsGas(pkd,q)) continue;
        if (pkdIsDeleted(pkd,q)) continue; /* deleted */
        if (r2 < r2min) {
            r2min = r2;
            qmin = q;
        }
    }
    assert(qmin!=NULL); /* What if all neighbours invalid -- does that happen? */
    q=qmin;

    {
        qmass = pkdField(q,pkd->oFieldOffset[oMass]);
        delta_m = *pmass*smf->SFdMassLossPerStarMass;
        m = *qmass + delta_m;
        im = 1/m;
        fp = delta_m*im;
        fq = (*qmass)*im;
        delta_u = smf->SFdESNPerStarMass*delta_m*im;
        delta_Z = smf->SFdZMassPerStarMass*delta_m*im;
        pkdStar(pkd,q)->fTimer = smf->dTime+smf->SFdtCoolingShutoff;

        *qmass = m;
        pkdSetPos(pkd,q,0,fq*pkdPos(pkd,q,0)+fp*pkdPos(pkd,p,0));
        pkdSetPos(pkd,q,1,fq*pkdPos(pkd,q,1)+fp*pkdPos(pkd,p,1));
        pkdSetPos(pkd,q,2,fq*pkdPos(pkd,q,2)+fp*pkdPos(pkd,p,2));
        pkdVel(pkd,q)[0] = fq*pkdVel(pkd,q)[0]+fp*pkdVel(pkd,p)[0];
        pkdVel(pkd,q)[1] = fq*pkdVel(pkd,q)[1]+fp*pkdVel(pkd,p)[1];
        pkdVel(pkd,q)[2] = fq*pkdVel(pkd,q)[2]+fp*pkdVel(pkd,p)[2];
        pkd_vPred(pkd,q)[0] = fq*pkd_vPred(pkd,q)[0]+fp*pkd_vPred(pkd,p)[0];
        pkd_vPred(pkd,q)[1] = fq*pkd_vPred(pkd,q)[1]+fp*pkd_vPred(pkd,p)[1];
        pkd_vPred(pkd,q)[2] = fq*pkd_vPred(pkd,q)[2]+fp*pkd_vPred(pkd,p)[2];
        pkdAccel(pkd,q)[0] = fq*pkdAccel(pkd,q)[0]+fp*pkdAccel(pkd,p)[0];
        pkdAccel(pkd,q)[1] = fq*pkdAccel(pkd,q)[1]+fp*pkdAccel(pkd,p)[1];
        pkdAccel(pkd,q)[2] = fq*pkdAccel(pkd,q)[2]+fp*pkdAccel(pkd,p)[2];
        pkdSph(pkd,q)->u = fq*pkdSph(pkd,q)->u+delta_u;
        pkdSph(pkd,q)->uPred = fq*pkdSph(pkd,q)->uPred+delta_u;
        pkdSph(pkd,q)->fMetals = fq*pkdSph(pkd,q)->fMetals+delta_Z;
        pkdSph(pkd,q)->fMetalsPred = fq*pkdSph(pkd,q)->fMetalsPred+delta_Z;

        c = sqrt(smf->gamma*(smf->gamma-1)*pkdSph(pkd,q)->uPred);
        qBall = pkdBall(pkd,q);
        dt = 0.5*qBall/(dtC*c);
        uNewRung = pkdDtToRung(dt,smf->dDelta,MAX_RUNG);

        if (q->uNewRung > uNewRung) uNewRung = q->uNewRung;
    }

    /* The SN is here ... lower timesteps */
    for (i=0; i<nSmooth; ++i) {
        q = nnList[i].pPart;
        if (!pkdIsGas(pkd,q)) continue;
        if (pkdIsDeleted(pkd,q)) continue; /* deleted */
        if (uNewRung > q->uNewRung) q->uNewRung = uNewRung;
    }

    *pmass -= delta_m; /* Lower star mass */
    pkdSph(pkd,p)->u = 0; /* Mark star as fedback */
}

void initMeanVel(void *vpkd, void *pvoid) {
    PKD pkd = (PKD)vpkd;
    PARTICLE *p = pvoid;
    VELSMOOTH *pvel;
    int j;
    assert(pkd);
    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    for (j=0; j<3; ++j) pvel->vmean[j] = 0.0;
}

void combMeanVel(void *vpkd, void *p1void,const void *p2void) {
    PKD pkd = (PKD)vpkd;
    PARTICLE *p1 = p1void;
    const PARTICLE *p2 = p2void;
    VELSMOOTH *p1vel;
    const VELSMOOTH *p2vel;
    int j;

    assert(pkd);
    p1vel = pkdField(p1,pkd->oFieldOffset[oVelSmooth]);
    p2vel = pkdFieldRO(p2,pkd->oFieldOffset[oVelSmooth]);

    for (j=0; j<3; ++j) p1vel->vmean[j] += p2vel->vmean[j];
}

void MeanVel(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    vel_t *qv;
    VELSMOOTH *pvel;
    float v[3];
    double ih2,r2,rs,fMass;
    int i,j;

    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    ih2 = 4.0/BALL2(fBall);
    for (j=0; j<3; ++j) v[j] = 0.0;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        q = nnList[i].pPart;
        fMass = pkdMass(pkd,q);
        qv = pkdVel(pkd,q);
        for (j=0; j<3; ++j) v[j] += rs*fMass/pkdDensity(pkd,q)*qv[j];
    }
    for (j=0; j<3; ++j) pvel->vmean[j] = M_1_PI*sqrt(ih2)*ih2*v[j];
}

void MeanVelSym(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    VELSMOOTH *pvel, *qvel;
    vel_t *pv, *qv;
    double fNorm,ih2,r2,rs,fMassQ,fMassP;
    int i,j;

    pv = pkdVel(pkd,p);
    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(fBall));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        rs *= fNorm;
        q = nnList[i].pPart;
        qv = pkdVel(pkd,q);
        qvel = pkdField(q,pkd->oFieldOffset[oVelSmooth]);
        fMassQ = pkdMass(pkd,q);
        for (j=0; j<3; ++j) {
            pvel->vmean[j] += rs*fMassQ/pkdDensity(pkd,q)*qv[j];
            qvel->vmean[j] += rs*fMassP/pkdDensity(pkd,p)*pv[j];
        }
    }
}



void initDivv(void *vpkd, void *pvoid) {
    PKD pkd = (PKD)vpkd;
    PARTICLE *p = pvoid;
    VELSMOOTH *pvel;
    assert(pkd);
    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    pvel->divv = 0.0;
}

void combDivv(void *vpkd, void *p1void,const void *p2void) {
    PKD pkd = (PKD)vpkd;
    VELSMOOTH *p1vel;
    const VELSMOOTH *p2vel;
    PARTICLE *p1 = p1void;
    const PARTICLE *p2 = p2void;
    assert(pkd);
    p1vel = pkdField(p1,pkd->oFieldOffset[oVelSmooth]);
    p2vel = pkdFieldRO(p2,pkd->oFieldOffset[oVelSmooth]);
    p1vel->divv += p2vel->divv;
}

void Divv(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    vel_t *qv;
    vel_t *pv;
    VELSMOOTH *pvel;
    double fNorm,ih2,r2,rs,fMass,dvdotdr;
    int i;

    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    pv = pkdVel(pkd,p);
    ih2 = 4.0/BALL2(fBall);
    fNorm = M_1_PI*ih2*ih2;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        DKERNEL(rs,r2);
        rs *= fNorm;
        q = nnList[i].pPart;
        fMass = pkdMass(pkd,q);
        qv = pkdVel(pkd,q);
        dvdotdr = (qv[0] - pv[0])*nnList[i].dx + (qv[1] - pv[1])*nnList[i].dy + (qv[2] - pv[2])*nnList[i].dz;
        pvel->divv += rs*fMass/pkdDensity(pkd,q)*dvdotdr;
    }
}

void DivvSym(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    VELSMOOTH *pvel, *qvel;
    vel_t *pv, *qv;
    double fNorm,ih2,r2,rs,fMassQ,fMassP,dvdotdr;
    int i;

    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    pv = pkdVel(pkd,p);
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(fBall));
    fNorm = 0.5*M_1_PI*ih2*ih2;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        DKERNEL(rs,r2);
        rs *= fNorm;
        q = nnList[i].pPart;
        qv = pkdVel(pkd,q);
        qvel = pkdField(q,pkd->oFieldOffset[oVelSmooth]);
        fMassQ = pkdMass(pkd,q);
        dvdotdr = (qv[0] - pv[0])*nnList[i].dx + (qv[1] - pv[1])*nnList[i].dy + (qv[2] - pv[2])*nnList[i].dz;
        pvel->divv += rs*fMassQ/pkdDensity(pkd,q)*dvdotdr;
        qvel->divv += rs*fMassP/pkdDensity(pkd,p)*dvdotdr;
    }
}



void initVelDisp2(void *vpkd, void *pvoid) {
    PKD pkd = (PKD)vpkd;
    PARTICLE *p = pvoid;
    VELSMOOTH *pvel;
    assert(pkd);
    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    pvel->veldisp2 = 0.0;
}

void combVelDisp2(void *vpkd, void *p1void,const void *p2void) {
    PKD pkd = (PKD)vpkd;
    PARTICLE *p1 = p1void;
    const PARTICLE *p2 = p2void;
    VELSMOOTH *p1vel;
    const VELSMOOTH *p2vel;
    assert(pkd);
    p1vel = pkdField(p1,pkd->oFieldOffset[oVelSmooth]);
    p2vel = pkdFieldRO(p2,pkd->oFieldOffset[oVelSmooth]);
    p1vel->veldisp2 += p2vel->veldisp2;
}

void VelDisp2(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    vel_t *qv;
    VELSMOOTH *pvel;
    double fNorm,ih2,r2,rs,fMass,tv,tv2;
    int i;

    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    ih2 = 4.0/BALL2(fBall);
    fNorm = M_1_PI*sqrt(ih2)*ih2;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        rs *= fNorm;
        q = nnList[i].pPart;
        fMass = pkdMass(pkd,q);
        qv = pkdVel(pkd,q);

        tv = qv[0] - pvel->vmean[0] - pvel->divv*nnList[i].dx;
        tv2 = tv*tv;
        tv = qv[1] - pvel->vmean[1] - pvel->divv*nnList[i].dy;
        tv2 += tv*tv;
        tv = qv[2] - pvel->vmean[2] - pvel->divv*nnList[i].dz;
        tv2 += tv*tv;

        pvel->veldisp2 += rs*fMass/pkdDensity(pkd,q)*tv2;
    }
}

void VelDisp2Sym(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    VELSMOOTH *pvel, *qvel;
    vel_t *qv;
    double fNorm,ih2,r2,rs,fMassQ,fMassP,tv,tv2;
    int i;

    pvel = pkdField(p,pkd->oFieldOffset[oVelSmooth]);
    fMassP = pkdMass(pkd,p);
    ih2 = 4.0/(BALL2(fBall));
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0; i<nSmooth; ++i) {
        r2 = nnList[i].fDist2*ih2;
        KERNEL(rs,r2);
        rs *= fNorm;
        q = nnList[i].pPart;
        qv = pkdVel(pkd,q);
        qvel = pkdField(q,pkd->oFieldOffset[oVelSmooth]);
        fMassQ = pkdMass(pkd,q);

        tv = qv[0] - pvel->vmean[0] - pvel->divv*nnList[i].dx;
        tv2 = tv*tv;
        tv = qv[1] - pvel->vmean[1] - pvel->divv*nnList[i].dy;
        tv2 += tv*tv;
        tv = qv[2] - pvel->vmean[2] - pvel->divv*nnList[i].dz;
        tv2 += tv*tv;

        pvel->veldisp2 += rs*fMassQ/pkdDensity(pkd,q)*tv2;
        qvel->veldisp2 += rs*fMassP/pkdDensity(pkd,p)*tv2;
    }
}

void AddRelaxation(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    double  vSigma2,L,fRel,e,pmax,pmin,vMean[3],feps;
    double beta,gamma,rho;
    int i,j;
    PARTICLE *q;
    vel_t *v;
    double *pRelax;

    pRelax = pkdField(p,pkd->oFieldOffset[oRelaxation]);

    beta = 10.0; /* pmax=beta*l, large beta parameter reduces eps dependence  */
    gamma = 0.17; /* gamma scales overall relaxation time */
    feps = 1.0; /* cuttoff for pmin at feps*epsilon or p_90 deg */

    for (j=0; j<3; ++j) {
        vMean[j] = 0.0;
    }
    vSigma2 = 0.0;
    rho = 0.0;
    for (i=0; i<nSmooth; ++i) {
        q = nnList[i].pPart;
        v = pkdVel(pkd,q);
        rho += pkdMass(pkd,q);
        for (j=0; j<3; ++j) {
            vSigma2 += v[j]*v[j];
            vMean[j] += v[j];
        }
    }
    vSigma2 /= nSmooth;
    rho /= 4.18*pow(fBall,3.0);
    for (j=0; j<3; ++j) {
        vMean[j] /= nSmooth;
        vSigma2 -= vMean[j]*vMean[j];
    }
    vSigma2 = 0.33333*vSigma2; /* now its the one dimensional vel.dispersion squared */
    pmin = 2.0*pkdMass(pkd,p)/(6.0*vSigma2);  /* pmin in comoving units */
    e = feps*pkdSoft(pkd,p);
    if (pmin < e) pmin = e; /* pmin is the bigger one of epsilon and p90 */
    pmax = beta*pow(pkdMass(pkd,p)/rho, 0.333333); /* pmax=beta*mean interp. separation */
    if (pmax < pmin || vSigma2 < 0) return;  /* no relaxation if pmin > pmax */
    L = pmax/pmin;
    fRel = 0.5*( log(1+L*L) - L*L/(1+L*L) );
    fRel *= pkdMass(pkd,p)*rho*pow(vSigma2, -1.5)/gamma;
    *pRelax += fRel * smf->dDeltaT;
}

#endif //OPTIM_REMOVE_UNUSED
