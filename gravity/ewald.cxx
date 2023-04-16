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

#include <stdio.h>
#ifdef HAVE_MALLOC_H
    #include <malloc.h>
#endif
#include <new>
#include <math.h>
#include <assert.h>
#include "ewald.h"
#include "pkd.h"
#include "pkd.h"
#include "qeval.h"
#include "moments.h"
#include "gravity/grav.h"
#ifdef USE_SIMD_EWALD
    #include "core/vmath.h"
#endif/*USE_SIMD_EWALD*/
#include "cuda/cudautil.h"

template<class F,class E,class M>
static int evalEwald(
    const E &ew, const M &mom,
    F &ax, F &ay, F &az, F &fPot,
    const F &x, const F &y, const F &z,
    const F &g0, const F &g1, const F &g2, const F &g3,const F &g4, const F &g5) {
    F onethird = 1.0/3.0;

    F xx = 0.5*x*x;
    F xxx = onethird*xx*x;
    F xxy = xx*y;
    F xxz = xx*z;
    F yy = 0.5*y*y;
    F yyy = onethird*yy*y;
    F xyy = yy*x;
    F yyz = yy*z;
    F zz = 0.5*z*z;
    F zzz = onethird*zz*z;
    F xzz = zz*x;
    F yzz = zz*y;
    F xy = x*y;
    F xyz = xy*z;
    F xz = x*z;
    F yz = y*z;
    F Q2mirx = F(mom.xx)*x + F(mom.xy)*y + F(mom.xz)*z;
    F Q2miry = F(mom.xy)*x + F(mom.yy)*y + F(mom.yz)*z;
    F Q2mirz = F(mom.xz)*x + F(mom.yz)*y + F(mom.zz)*z;
    F Q3mirx = F(mom.xxx)*xx + F(mom.xxy)*xy + F(mom.xxz)*xz + F(mom.xyy)*yy + F(mom.xyz)*yz + F(mom.xzz)*zz;
    F Q3miry = F(mom.xxy)*xx + F(mom.xyy)*xy + F(mom.xyz)*xz + F(mom.yyy)*yy + F(mom.yyz)*yz + F(mom.yzz)*zz;
    F Q3mirz = F(mom.xxz)*xx + F(mom.xyz)*xy + F(mom.xzz)*xz + F(mom.yyz)*yy + F(mom.yzz)*yz + F(mom.zzz)*zz;
    F Q4mirx = F(mom.xxxx)*xxx + F(mom.xxxy)*xxy + F(mom.xxxz)*xxz + F(mom.xxyy)*xyy + F(mom.xxyz)*xyz +
               F(mom.xxzz)*xzz + F(mom.xyyy)*yyy + F(mom.xyyz)*yyz + F(mom.xyzz)*yzz + F(mom.xzzz)*zzz;
    F Q4miry = F(mom.xxxy)*xxx + F(mom.xxyy)*xxy + F(mom.xxyz)*xxz + F(mom.xyyy)*xyy + F(mom.xyyz)*xyz +
               F(mom.xyzz)*xzz + F(mom.yyyy)*yyy + F(mom.yyyz)*yyz + F(mom.yyzz)*yzz + F(mom.yzzz)*zzz;
    F Q4mirz = F(mom.xxxz)*xxx + F(mom.xxyz)*xxy + F(mom.xxzz)*xxz + F(mom.xyyz)*xyy + F(mom.xyzz)*xyz +
               F(mom.xzzz)*xzz + F(mom.yyyz)*yyy + F(mom.yyzz)*yyz + F(mom.yzzz)*yzz + F(mom.zzzz)*zzz;
    F Q4x = F(ew.Q4xx)*x + F(ew.Q4xy)*y + F(ew.Q4xz)*z;
    F Q4y = F(ew.Q4xy)*x + F(ew.Q4yy)*y + F(ew.Q4yz)*z;
    F Q4z = F(ew.Q4xz)*x + F(ew.Q4yz)*y + F(ew.Q4zz)*z;
    F Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (F(ew.Q3x)*x + F(ew.Q3y)*y + F(ew.Q3z)*z) + F(ew.Q4);
    F Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
    F Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
    F Qta = g1*F(mom.m) - g2*F(ew.Q2) + g3*Q2mir + g4*Q3mir + g5*Q4mir;
    fPot -= g0*F(mom.m) - g1*F(ew.Q2) + g2*Q2mir + g3*Q3mir + g4*Q4mir;
    ax += g2*(Q2mirx - F(ew.Q3x)) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
    ay += g2*(Q2miry - F(ew.Q3y)) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
    az += g2*(Q2mirz - F(ew.Q3z)) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;

    return COST_FLOP_EWALD;
}

#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
double evalEwaldSIMD( PKD pkd,ewaldSIMD *ews,
                      dvec &ax, dvec &ay, dvec &az, dvec &dPot,
                      dvec x, dvec y, dvec z, dvec r2, const dmask &doerfc ) {
    dvec dir,dir2,a,g0,g1,g2,g3,g4,g5,alphan;
    dvec rerf,rerfc,ex2;
    dvec alpha2x2 = 2.0 * dvec(ews->ewp.alpha2);

    dir = rsqrt(r2);
    dir2 = dir*dir;
    ex2 = exp(-r2*dvec(ews->ewp.alpha2));
    a = ex2 * dvec(ews->ewp.ka) * dir2;

    verf(dvec(ews->ewp.alpha)*r2*dir,dvec(ews->ewp.ialpha)*dir,ex2,rerf,rerfc);

    g0 = dir * mask_mov(-rerf,doerfc,rerfc);
    g1 = g0*dir2 + a;
    alphan = alpha2x2;
    g2 = 3.0*g1*dir2 + alphan*a;
    alphan *= alpha2x2;
    g3 = 5.0*g2*dir2 + alphan*a;
    alphan *= alpha2x2;
    g4 = 7.0*g3*dir2 + alphan*a;
    alphan *= alpha2x2;
    g5 = 9.0*g4*dir2 + alphan*a;

    evalEwald<dvec,struct ewaldSIMD::PEWALDVARS,struct ewaldSIMD::PMOMC>(
            ews->ewp,ews->ewm,ax,ay,az,dPot,x,y,z,g0,g1,g2,g3,g4,g5);

    return COST_FLOP_EWALD * dvec::width();
}
#endif

double pkdParticleEwald(PKD pkd,double *r, float *pa, float *pPot,double *pdFlopSingle, double *pdFlopDouble) {
    struct EwaldVariables &ew = pkd->ew;
    EwaldTable *ewt = &pkd->ewt;
    const MOMC &restrict mom = ew.mom;
    double L,Pot,ax,ay,az,dx,dy,dz,x,y,z,r2;
#ifdef USE_SIMD_EWALD
    dvec dPot,dax,day,daz;
    fvec fPot,fax,fay,faz,fx,fy,fz;
    dvec::array_t px,py,pz,pr2,pInHole;
    int nSIMD = 0;
#endif
    int i,ix,iy,iz;
    int bInHole,bInHolex,bInHolexy;
    double dFlopSingle = 0;
    double dFlopDouble = 0;
    int nLoop = 0;

    L = ew.Lbox;
    dx = r[0] - ew.r[0];
    dy = r[1] - ew.r[1];
    dz = r[2] - ew.r[2];

    ax = ay = az = 0.0;
    Pot = mom.m*ew.k1;
#ifdef USE_SIMD_EWALD
    dPot.zero();
    dax.zero();
    day.zero();
    daz.zero();
#endif
    for (ix=-ew.nEwReps; ix<=ew.nEwReps; ++ix) {
        bInHolex = (abs(ix) <= ew.nReps);
        x = dx + ix*L;
        for (iy=-ew.nEwReps; iy<=ew.nEwReps; ++iy) {
            bInHolexy = (bInHolex && abs(iy) <= ew.nReps);
            y = dy + iy*L;
            for (iz=-ew.nEwReps; iz<=ew.nEwReps; ++iz) {
                bInHole = (bInHolexy && abs(iz) <= ew.nReps);
                z = dz + iz*L;
                r2 = x*x + y*y + z*z;
                if (r2 > ew.fEwCut2 && !bInHole) continue;
                if (r2 < ew.fInner2) {
                    double g0,g1,g2,g3,g4,g5,alphan;
                    /*
                     * For small r, series expand about
                     * the origin to avoid errors caused
                     * by cancellation of large terms.
                     */
                    alphan = ew.ka;
                    r2 *= ew.alpha2;
                    g0 = alphan*((1.0/3.0)*r2 - 1.0);
                    alphan *= 2*ew.alpha2;
                    g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
                    alphan *= 2*ew.alpha2;
                    g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
                    alphan *= 2*ew.alpha2;
                    g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
                    alphan *= 2*ew.alpha2;
                    g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
                    alphan *= 2*ew.alpha2;
                    g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));

                    dFlopDouble += evalEwald<double,struct EwaldVariables,MOMC>(ew,ew.mom,ax,ay,az,Pot,x,y,z,g0,g1,g2,g3,g4,g5);
                }
                else {
#if defined(USE_SIMD_EWALD)
                    px[nSIMD] = x;
                    py[nSIMD] = y;
                    pz[nSIMD] = z;
                    pr2[nSIMD] = r2;
                    pInHole[nSIMD] = bInHole;
//          doerfc.i[nSIMD] = bInHole ? 0 : UINT64_MAX;
                    if (++nSIMD == dvec::width()) {
                        dFlopDouble += evalEwaldSIMD(pkd,&pkd->es,dax,day,daz,dPot,dvec(px),dvec(py),dvec(pz),dvec(pr2),dvec(pInHole) == 0.0);
                        nSIMD = 0;
                    }
#else
                    double dir,dir2,a;
                    double g0,g1,g2,g3,g4,g5,alphan;
                    dir = 1/sqrt(r2);
                    dir2 = dir*dir;
                    a = exp(-r2*ew.alpha2);
                    a *= ew.ka*dir2;


                    if (bInHole) {
                        g0 = -erf(ew.alpha*r2*dir);
                    }
                    else {
                        g0 = erfc(ew.alpha*r2*dir);
                    }
                    g0 *= dir;
                    g1 = g0*dir2 + a;
                    alphan = 2*ew.alpha2;
                    g2 = 3*g1*dir2 + alphan*a;
                    alphan *= 2*ew.alpha2;
                    g3 = 5*g2*dir2 + alphan*a;
                    alphan *= 2*ew.alpha2;
                    g4 = 7*g3*dir2 + alphan*a;
                    alphan *= 2*ew.alpha2;
                    g5 = 9*g4*dir2 + alphan*a;
                    dFlopDouble += evalEwald<double,struct EwaldVariables,MOMC>(ew,ew.mom,ax,ay,az,Pot,x,y,z,g0,g1,g2,g3,g4,g5);
                    //dFlopDouble += evalEwald(ew,&ax,&ay,&az,&Pot,x,y,z,g0,g1,g2,g3,g4,g5);
#endif
                }
                ++nLoop;
            }
        }
    }
#if defined(USE_SIMD_EWALD)
    /* Finish remaining SIMD operations if necessary */
    if (nSIMD) { /* nSIMD can be 0 through 7 */
#define M (~0LL)
#if defined(__AVX512F__)
        static const i64v::array_t keepmask[] = {{0,0,0,0,0,0,0,0},{M,0,0,0,0,0,0,0},{M,M,0,0,0,0,0,0},{M,M,M,0,0,0,0,0},
            {M,M,M,M,0,0,0,0},{M,M,M,M,M,0,0,0},{M,M,M,M,M,M,0,0},{M,M,M,M,M,M,M,0}
        };
#elif defined(__AVX__)
        static const i64v::array_t keepmask[] = {{0,0,0,0},{M,0,0,0},{M,M,0,0},{M,M,M,0}};
#else
        static const i64v::array_t keepmask[] = {{0,0},{M,0}};
#endif
#undef M
        dvec t,tax=0, tay=0, taz=0, tpot=0;
        evalEwaldSIMD(pkd,&pkd->es,tax,tay,taz,tpot,dvec(px),dvec(py),dvec(pz),dvec(pr2),dvec(pInHole) == 0.0);
        dFlopDouble += COST_FLOP_EWALD * nSIMD;
        t = cast_dvec(i64v(keepmask[nSIMD]));
        tax = tax & t;
        tay = tay & t;
        taz = taz & t;
        tpot = tpot & t;
        dax += tax;
        day += tay;
        daz += taz;
        dPot-= tpot;
        nSIMD = 0;
    }
#endif

#ifdef USE_SIMD_EWALD
    /* h-loop is done in float precision */
    fax = cvt_fvec(dax);
    fay = cvt_fvec(day);
    faz = cvt_fvec(daz);
    fPot = cvt_fvec(dPot);

    fx = dx;
    fy = dy;
    fz = dz;

    nLoop = (ew.nEwhLoop+fvec::mask()) / fvec::width();
    i = 0;
    do {
        fvec hdotx,s,c,t;
        hdotx = fvec(ewt->hx.p[i])*fx + fvec(ewt->hy.p[i])*fy + fvec(ewt->hz.p[i])*fz;
        fvec svec,cvec;

        sincosf(fvec(hdotx),svec,cvec);
        s = svec; c = cvec;

        fPot += fvec(ewt->hSfac.p[i])*s + fvec(ewt->hCfac.p[i])*c;
        s *= fvec(ewt->hCfac.p[i]);
        c *= fvec(ewt->hSfac.p[i]);
        t = s - c;
        fax += fvec(ewt->hx.p[i])*t;
        fay += fvec(ewt->hy.p[i])*t;
        faz += fvec(ewt->hz.p[i])*t;
    } while (++i < nLoop);
    dFlopSingle += ew.nEwhLoop*COST_FLOP_HLOOP;

    ax += hadd(fax);
    ay += hadd(fay);
    az += hadd(faz);
    Pot += hadd(fPot);
#else
    /*
    ** Scoring for the h-loop (+,*)
    **  Without trig = (10,14)
    **      Trig est.    = 2*(6,11)  same as 1/sqrt scoring.
    **      Total        = (22,36)
    **                   = 58
    */
    for (i=0; i<ew.nEwhLoop; ++i) {
        double hdotx,s,c,t;
        hdotx = ewt->hx.f[i]*dx + ewt->hy.f[i]*dy + ewt->hz.f[i]*dz;
        c = cos(hdotx);
        s = sin(hdotx);
        Pot += ewt->hCfac.f[i]*c + ewt->hSfac.f[i]*s;
        t = ewt->hCfac.f[i]*s - ewt->hSfac.f[i]*c;
        ax += ewt->hx.f[i]*t;
        ay += ewt->hy.f[i]*t;
        az += ewt->hz.f[i]*t;
    }
    dFlopDouble += ew.nEwhLoop*COST_FLOP_HLOOP;
#endif
    pa[0] += ax;
    pa[1] += ay;
    pa[2] += az;
    *pPot += Pot;

    *pdFlopSingle += dFlopSingle;
    *pdFlopDouble += dFlopDouble;
    return dFlopDouble + dFlopSingle;
}

void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut) {
    struct EwaldVariables *const ew = &pkd->ew;
    EwaldTable *const ewt = &pkd->ewt;
    const MOMC *restrict mom = &ew->mom;
    int i,hReps,hx,hy,hz,h2;
    double k4,L;
    double gam[6],mfacc,mfacs;
    double ax,ay,az;
    const int iOrder = 4;

    L = pkd->fPeriod[0];
    ew->Lbox = L;
    /*
    ** Create SIMD versions of the moments.
    */
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    pkd->es.ewm.m = dvec(mom->m);
    pkd->es.ewm.xx = dvec(mom->xx);
    pkd->es.ewm.yy = dvec(mom->yy);
    pkd->es.ewm.xy = dvec(mom->xy);
    pkd->es.ewm.xz = dvec(mom->xz);
    pkd->es.ewm.yz = dvec(mom->yz);
    pkd->es.ewm.xxx = dvec(mom->xxx);
    pkd->es.ewm.xyy = dvec(mom->xyy);
    pkd->es.ewm.xxy = dvec(mom->xxy);
    pkd->es.ewm.yyy = dvec(mom->yyy);
    pkd->es.ewm.xxz = dvec(mom->xxz);
    pkd->es.ewm.yyz = dvec(mom->yyz);
    pkd->es.ewm.xyz = dvec(mom->xyz);
    pkd->es.ewm.xxxx = dvec(mom->xxxx);
    pkd->es.ewm.xyyy = dvec(mom->xyyy);
    pkd->es.ewm.xxxy = dvec(mom->xxxy);
    pkd->es.ewm.yyyy = dvec(mom->yyyy);
    pkd->es.ewm.xxxz = dvec(mom->xxxz);
    pkd->es.ewm.yyyz = dvec(mom->yyyz);
    pkd->es.ewm.xxyy = dvec(mom->xxyy);
    pkd->es.ewm.xxyz = dvec(mom->xxyz);
    pkd->es.ewm.xyyz = dvec(mom->xyyz);
    pkd->es.ewm.zz = dvec(mom->zz);
    pkd->es.ewm.xzz = dvec(mom->xzz);
    pkd->es.ewm.yzz = dvec(mom->yzz);
    pkd->es.ewm.zzz = dvec(mom->zzz);
    pkd->es.ewm.xxzz = dvec(mom->xxzz);
    pkd->es.ewm.xyzz = dvec(mom->xyzz);
    pkd->es.ewm.xzzz = dvec(mom->xzzz);
    pkd->es.ewm.yyzz = dvec(mom->yyzz);
    pkd->es.ewm.yzzz = dvec(mom->yzzz);
    pkd->es.ewm.zzzz = dvec(mom->zzzz);
#endif

    /*
    ** Set up traces of the complete multipole moments.
    */
    ew->Q4xx = 0.5*(mom->xxxx + mom->xxyy + mom->xxzz);
    ew->Q4xy = 0.5*(mom->xxxy + mom->xyyy + mom->xyzz);
    ew->Q4xz = 0.5*(mom->xxxz + mom->xyyz + mom->xzzz);
    ew->Q4yy = 0.5*(mom->xxyy + mom->yyyy + mom->yyzz);
    ew->Q4yz = 0.5*(mom->xxyz + mom->yyyz + mom->yzzz);
    ew->Q4zz = 0.5*(mom->xxzz + mom->yyzz + mom->zzzz);
    ew->Q4 = 0.25*(ew->Q4xx + ew->Q4yy + ew->Q4zz);
    ew->Q3x = 0.5*(mom->xxx + mom->xyy + mom->xzz);
    ew->Q3y = 0.5*(mom->xxy + mom->yyy + mom->yzz);
    ew->Q3z = 0.5*(mom->xxz + mom->yyz + mom->zzz);
    ew->Q2 = 0.5*(mom->xx + mom->yy + mom->zz);
    ew->nReps = nReps;
    ew->nEwReps = ceil(fEwCut);
    ew->nEwReps = ew->nEwReps > nReps ? ew->nEwReps : nReps;
    ew->fEwCut2 = fEwCut*fEwCut*L*L;
    ew->fInner2 = 1.2e-3*L*L;
    ew->alpha = 2.0/L;
    ew->ialpha = 0.5 * L;
    ew->alpha2 = ew->alpha*ew->alpha;
    ew->k1 = M_PI/(ew->alpha2*L*L*L);
    ew->ka = 2.0*ew->alpha/sqrt(M_PI);
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    pkd->es.ewp.Q4xx = dvec(ew->Q4xx);
    pkd->es.ewp.Q4xy = dvec(ew->Q4xy);
    pkd->es.ewp.Q4xz = dvec(ew->Q4xz);
    pkd->es.ewp.Q4yy = dvec(ew->Q4yy);
    pkd->es.ewp.Q4yz = dvec(ew->Q4yz);
    pkd->es.ewp.Q4zz = dvec(ew->Q4zz);
    pkd->es.ewp.Q4 = dvec(ew->Q4);
    pkd->es.ewp.Q3x = dvec(ew->Q3x);
    pkd->es.ewp.Q3y = dvec(ew->Q3y);
    pkd->es.ewp.Q3z = dvec(ew->Q3z);
    pkd->es.ewp.Q2 = dvec(ew->Q2);
    pkd->es.ewp.fEwCut2 = dvec(ew->fEwCut2);
    pkd->es.ewp.fInner2 = dvec(ew->fInner2);
    pkd->es.ewp.alpha = dvec(ew->alpha);
    pkd->es.ewp.ialpha = dvec(ew->ialpha);
    pkd->es.ewp.alpha2 = dvec(ew->alpha2);
    pkd->es.ewp.k1 = dvec(ew->k1);
    pkd->es.ewp.ka = dvec(ew->ka);
#endif


    /*
    ** Now setup stuff for the h-loop.
    */
    hReps = ceil(fhCut);
    k4 = M_PI*M_PI/(ew->alpha*ew->alpha*L*L);

    i = (int)pow(1+2*hReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE__)
    i = (i + fvec::mask()) & ~fvec::mask();
#endif
    if ( i>ew->nMaxEwhLoop ) {
        ew->nMaxEwhLoop = i;
        ewt->hx.f = new (std::align_val_t(sizeof(fvec))) ewaldFloatType[ew->nMaxEwhLoop];
        assert(ewt->hx.f != NULL);
        ewt->hy.f = new (std::align_val_t(sizeof(fvec))) ewaldFloatType[ew->nMaxEwhLoop];
        assert(ewt->hy.f != NULL);
        ewt->hz.f = new (std::align_val_t(sizeof(fvec))) ewaldFloatType[ew->nMaxEwhLoop];
        assert(ewt->hz.f != NULL);
        ewt->hCfac.f = new (std::align_val_t(sizeof(fvec))) ewaldFloatType[ew->nMaxEwhLoop];
        assert(ewt->hCfac.f != NULL);
        ewt->hSfac.f = new (std::align_val_t(sizeof(fvec))) ewaldFloatType[ew->nMaxEwhLoop];
        assert(ewt->hSfac.f != NULL);
    }
    ew->nEwhLoop = i;
    i = (int)pow(1+2*ew->nEwReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    i = (i + fvec::mask()) & ~fvec::mask();
#endif
    i = 0;
    for (hx=-hReps; hx<=hReps; ++hx) {
        for (hy=-hReps; hy<=hReps; ++hy) {
            for (hz=-hReps; hz<=hReps; ++hz) {
                h2 = hx*hx + hy*hy + hz*hz;
                if (h2 == 0) continue;
                if (h2 > fhCut*fhCut) continue;
                assert (i < ew->nMaxEwhLoop);
                gam[0] = exp(-k4*h2)/(M_PI*h2*L);
                gam[1] = 2*M_PI/L*gam[0];
                gam[2] = -2*M_PI/L*gam[1];
                gam[3] = 2*M_PI/L*gam[2];
                gam[4] = -2*M_PI/L*gam[3];
                gam[5] = 2*M_PI/L*gam[4];
                gam[1] = 0.0;
                gam[3] = 0.0;
                gam[5] = 0.0;
                ax = 0.0;
                ay = 0.0;
                az = 0.0;
                mfacc = 0.0;
                QEVAL(iOrder,ew->mom,gam,hx,hy,hz,ax,ay,az,mfacc);
                gam[0] = exp(-k4*h2)/(M_PI*h2*L);
                gam[1] = 2*M_PI/L*gam[0];
                gam[2] = -2*M_PI/L*gam[1];
                gam[3] = 2*M_PI/L*gam[2];
                gam[4] = -2*M_PI/L*gam[3];
                gam[5] = 2*M_PI/L*gam[4];
                gam[0] = 0.0;
                gam[2] = 0.0;
                gam[4] = 0.0;
                ax = 0.0;
                ay = 0.0;
                az = 0.0;
                mfacs = 0.0;
                QEVAL(iOrder,ew->mom,gam,hx,hy,hz,ax,ay,az,mfacs);
                ewt->hx.f[i] = 2*M_PI/L*hx;
                ewt->hy.f[i] = 2*M_PI/L*hy;
                ewt->hz.f[i] = 2*M_PI/L*hz;
                ewt->hCfac.f[i] = mfacc;
                ewt->hSfac.f[i] = mfacs;
                ++i;
            }
        }
    }
    ew->nEwhLoop = i;
    while (i<ew->nMaxEwhLoop) {
        ewt->hx.f[i] = 0;
        ewt->hy.f[i] = 0;
        ewt->hz.f[i] = 0;
        ewt->hCfac.f[i] = 0;
        ewt->hSfac.f[i] = 0;
        ++i;
    }
#ifdef USE_CL
    clEwaldInit(pkd->mdl->clCtx,ew,ewt);
    mdlThreadBarrier(pkd->mdl);
#endif
#ifdef USE_CUDA
    auto cuda = reinterpret_cast<CudaClient *>(pkd->cudaClient);
    // Only one thread needs to transfer the tables to the GPU
    if (pkd->mdl->Core()==0) cuda->setupEwald(ew,ewt);
    pkd->mdl->ThreadBarrier();
#endif
}
