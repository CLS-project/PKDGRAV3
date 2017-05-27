#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include "ewald.h"
#include "pkd.h"
#include "pkd.h"
#include "qeval.h"
#include "moments.h"
#include "grav.h"
#include "vmath.h"

static int evalEwald(struct EwaldVariables *ew,double *ax, double *ay, double *az, double *fPot,
    double x, double y, double z,
    double g0, double g1, double g2, double g3,double g4, double g5) {
    const MOMC * restrict mom = &ew->mom;
    double onethird = 1.0/3.0;
    double xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
    double Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
    double Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;

    xx = 0.5*x*x;
    xxx = onethird*xx*x;
    xxy = xx*y;
    xxz = xx*z;
    yy = 0.5*y*y;
    yyy = onethird*yy*y;
    xyy = yy*x;
    yyz = yy*z;
    zz = 0.5*z*z;
    zzz = onethird*zz*z;
    xzz = zz*x;
    yzz = zz*y;
    xy = x*y;
    xyz = xy*z;
    xz = x*z;
    yz = y*z;
    Q2mirx = mom->xx*x + mom->xy*y + mom->xz*z;
    Q2miry = mom->xy*x + mom->yy*y + mom->yz*z;
    Q2mirz = mom->xz*x + mom->yz*y + mom->zz*z;
    Q3mirx = mom->xxx*xx + mom->xxy*xy + mom->xxz*xz + mom->xyy*yy + mom->xyz*yz + mom->xzz*zz;
    Q3miry = mom->xxy*xx + mom->xyy*xy + mom->xyz*xz + mom->yyy*yy + mom->yyz*yz + mom->yzz*zz;
    Q3mirz = mom->xxz*xx + mom->xyz*xy + mom->xzz*xz + mom->yyz*yy + mom->yzz*yz + mom->zzz*zz;
    Q4mirx = mom->xxxx*xxx + mom->xxxy*xxy + mom->xxxz*xxz + mom->xxyy*xyy + mom->xxyz*xyz +
	mom->xxzz*xzz + mom->xyyy*yyy + mom->xyyz*yyz + mom->xyzz*yzz + mom->xzzz*zzz;
    Q4miry = mom->xxxy*xxx + mom->xxyy*xxy + mom->xxyz*xxz + mom->xyyy*xyy + mom->xyyz*xyz +
	mom->xyzz*xzz + mom->yyyy*yyy + mom->yyyz*yyz + mom->yyzz*yzz + mom->yzzz*zzz;
    Q4mirz = mom->xxxz*xxx + mom->xxyz*xxy + mom->xxzz*xxz + mom->xyyz*xyy + mom->xyzz*xyz +
	mom->xzzz*xzz + mom->yyyz*yyy + mom->yyzz*yyz + mom->yzzz*yzz + mom->zzzz*zzz;
    Q4x = ew->Q4xx*x + ew->Q4xy*y + ew->Q4xz*z;
    Q4y = ew->Q4xy*x + ew->Q4yy*y + ew->Q4yz*z;
    Q4z = ew->Q4xz*x + ew->Q4yz*y + ew->Q4zz*z;
    Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ew->Q3x*x + ew->Q3y*y + ew->Q3z*z) + ew->Q4;
    Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
    Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
    Qta = g1*mom->m - g2*ew->Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
    *fPot -= g0*mom->m - g1*ew->Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
    *ax += g2*(Q2mirx - ew->Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
    *ay += g2*(Q2miry - ew->Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
    *az += g2*(Q2mirz - ew->Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;

    return COST_FLOP_EWALD;
    }

#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
static const struct CONSTS {
    vdouble onequarter,onethird,half,one,two,three,five,seven,nine;
    } consts = {
        {SIMD_DCONST(0.25)},
        {SIMD_DCONST(1.0/3.0)},
        {SIMD_DCONST(0.5)},
        {SIMD_DCONST(1.0)},
        {SIMD_DCONST(2.0)},
        {SIMD_DCONST(3.0)},
        {SIMD_DCONST(5.0)},
        {SIMD_DCONST(7.0)},
        {SIMD_DCONST(9.0)},
    };

#define vrsqrt SIMD_DRSQRT_EXACT

double evalEwaldSIMD( PKD pkd,ewaldSIMD *ews, dvec &ax, dvec &ay, dvec &az, dvec &dPot, v_df Ix, v_df Iy, v_df Iz, v_df Ir2, dmask doerfc ) {
    dvec dir,dir2,a,g0,g1,g2,g3,g4,g5,alphan;
    dvec xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
    dvec Qta,Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
    dvec Q3mirx,Q3miry,Q3mirz,Q3mir,Q2mirx,Q2miry,Q2mirz,Q2mir;
    dvec rerf,rerfc,ex2,t,tx,ty,tz,tpot;
    dvec alpha2x2 = 2.0 * dvec(ews->ewp.alpha2.p);
    static const double onethird = 1.0 / 3.0;
    struct EwaldVariables * const ew = &pkd->ew;
    const MOMC * restrict mom = &ew->mom;
    dvec x=Ix, y=Iy, z=Iz, r2=Ir2;

    dir = rsqrt(r2);
    dir2 = dir*dir;
    ex2 = exp(-r2*ews->ewp.alpha2.p);
    a = ex2 * ews->ewp.ka.p * dir2;

    verf(ews->ewp.alpha.p*r2*dir,ews->ewp.ialpha.p*dir,ex2,rerf,rerfc);

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

    xx = 0.5*x*x;
    xxx = onethird*xx*x;
    xxy = xx*y;
    xxz = xx*z;
    yy = 0.5*y*y;
    yyy = onethird*yy*y;
    xyy = yy*x;
    yyz = yy*z;
    zz = 0.5*z*z;
    zzz = onethird*zz*z;
    xzz = zz*x;
    yzz = zz*y;
    xy = x*y;
    xyz = xy*z;
    xz = x*z;
    yz = y*z;

#if 0
    Q2mirx = mom->xx*x + mom->xy*y + mom->xz*z;
    Q2miry = mom->xy*x + mom->yy*y + mom->yz*z;
    Q2mirz = mom->xz*x + mom->yz*y + mom->zz*z;
    Q3mirx = mom->xxx*xx + mom->xxy*xy + mom->xxz*xz + mom->xyy*yy + mom->xyz*yz + mom->xzz*zz;
    Q3miry = mom->xxy*xx + mom->xyy*xy + mom->xyz*xz + mom->yyy*yy + mom->yyz*yz + mom->yzz*zz;
    Q3mirz = mom->xxz*xx + mom->xyz*xy + mom->xzz*xz + mom->yyz*yy + mom->yzz*yz + mom->zzz*zz;
    Q4mirx = mom->xxxx*xxx + mom->xxxy*xxy + mom->xxxz*xxz + mom->xxyy*xyy + mom->xxyz*xyz +
	mom->xxzz*xzz + mom->xyyy*yyy + mom->xyyz*yyz + mom->xyzz*yzz + mom->xzzz*zzz;
    Q4miry = mom->xxxy*xxx + mom->xxyy*xxy + mom->xxyz*xxz + mom->xyyy*xyy + mom->xyyz*xyz +
	mom->xyzz*xzz + mom->yyyy*yyy + mom->yyyz*yyz + mom->yyzz*yzz + mom->yzzz*zzz;
    Q4mirz = mom->xxxz*xxx + mom->xxyz*xxy + mom->xxzz*xxz + mom->xyyz*xyy + mom->xyzz*xyz +
	mom->xzzz*xzz + mom->yyyz*yyy + mom->yyzz*yyz + mom->yzzz*yzz + mom->zzzz*zzz;
    Q4x = ew->Q4xx*x + ew->Q4xy*y + ew->Q4xz*z;
    Q4y = ew->Q4xy*x + ew->Q4yy*y + ew->Q4yz*z;
    Q4z = ew->Q4xz*x + ew->Q4yz*y + ew->Q4zz*z;
    Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ew->Q3x*x + ew->Q3y*y + ew->Q3z*z) + ew->Q4;
    Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
    Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
    Qta = g1*mom->m - g2*ew->Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;

    dPot = g0*mom->m - g1*ew->Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
    ax = g2*(Q2mirx - ew->Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
    ay = g2*(Q2miry - ew->Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
    az = g2*(Q2mirz - ew->Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
#elif 1
    Q2mirx = ews->ewm.xx.p*x + ews->ewm.xy.p*y + ews->ewm.xz.p*z;
    Q2miry = ews->ewm.xy.p*x + ews->ewm.yy.p*y + ews->ewm.yz.p*z;
    Q2mirz = ews->ewm.xz.p*x + ews->ewm.yz.p*y + ews->ewm.zz.p*z;
    Q3mirx = ews->ewm.xxx.p*xx + ews->ewm.xxy.p*xy + ews->ewm.xxz.p*xz + ews->ewm.xyy.p*yy + ews->ewm.xyz.p*yz + ews->ewm.xzz.p*zz;
    Q3miry = ews->ewm.xxy.p*xx + ews->ewm.xyy.p*xy + ews->ewm.xyz.p*xz + ews->ewm.yyy.p*yy + ews->ewm.yyz.p*yz + ews->ewm.yzz.p*zz;
    Q3mirz = ews->ewm.xxz.p*xx + ews->ewm.xyz.p*xy + ews->ewm.xzz.p*xz + ews->ewm.yyz.p*yy + ews->ewm.yzz.p*yz + ews->ewm.zzz.p*zz;
    Q4mirx = ews->ewm.xxxx.p*xxx + ews->ewm.xxxy.p*xxy + ews->ewm.xxxz.p*xxz + ews->ewm.xxyy.p*xyy + ews->ewm.xxyz.p*xyz +
	ews->ewm.xxzz.p*xzz + ews->ewm.xyyy.p*yyy + ews->ewm.xyyz.p*yyz + ews->ewm.xyzz.p*yzz + ews->ewm.xzzz.p*zzz;
    Q4miry = ews->ewm.xxxy.p*xxx + ews->ewm.xxyy.p*xxy + ews->ewm.xxyz.p*xxz + ews->ewm.xyyy.p*xyy + ews->ewm.xyyz.p*xyz +
	ews->ewm.xyzz.p*xzz + ews->ewm.yyyy.p*yyy + ews->ewm.yyyz.p*yyz + ews->ewm.yyzz.p*yzz + ews->ewm.yzzz.p*zzz;
    Q4mirz = ews->ewm.xxxz.p*xxx + ews->ewm.xxyz.p*xxy + ews->ewm.xxzz.p*xxz + ews->ewm.xyyz.p*xyy + ews->ewm.xyzz.p*xyz +
	ews->ewm.xzzz.p*xzz + ews->ewm.yyyz.p*yyy + ews->ewm.yyzz.p*yyz + ews->ewm.yzzz.p*yzz + ews->ewm.zzzz.p*zzz;
    Q4x = ews->ewp.Q4xx.p*x + ews->ewp.Q4xy.p*y + ews->ewp.Q4xz.p*z;
    Q4y = ews->ewp.Q4xy.p*x + ews->ewp.Q4yy.p*y + ews->ewp.Q4yz.p*z;
    Q4z = ews->ewp.Q4xz.p*x + ews->ewp.Q4yz.p*y + ews->ewp.Q4zz.p*z;
    Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (ews->ewp.Q3x.p*x + ews->ewp.Q3y.p*y + ews->ewp.Q3z.p*z) + ews->ewp.Q4.p;
    Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
    Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
    Qta = g1*ews->ewm.m.p - g2*ews->ewp.Q2.p + g3*Q2mir + g4*Q3mir + g5*Q4mir;

    dPot = g0*ews->ewm.m.p - g1*ews->ewp.Q2.p + g2*Q2mir + g3*Q3mir + g4*Q4mir;
    ax = g2*(Q2mirx - ews->ewp.Q3x.p) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
    ay = g2*(Q2miry - ews->ewp.Q3y.p) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
    az = g2*(Q2mirz - ews->ewp.Q3z.p) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
#else
    Q2mirx = SIMD_DMUL(ews->ewm.xx.p,x);
    Q2mirx = SIMD_DMADD(ews->ewm.xy.p,y,Q2mirx);
    Q2mirx = SIMD_DMADD(ews->ewm.xz.p,z,Q2mirx);
    Q2miry = SIMD_DMUL(ews->ewm.xy.p,x);
    Q2miry = SIMD_DMADD(ews->ewm.yy.p,y,Q2miry);
    Q2miry = SIMD_DMADD(ews->ewm.yz.p,z,Q2miry);
    Q2mirz = SIMD_DMUL(ews->ewm.xz.p,x);
    Q2mirz = SIMD_DMADD(ews->ewm.yz.p,y,Q2mirz);
    Q2mirz = SIMD_DMADD(ews->ewm.zz.p,z,Q2mirz);
    Q3mirx = SIMD_DMUL(ews->ewm.xxx.p,xx);
    Q3mirx = SIMD_DMADD(ews->ewm.xxy.p,xy,Q3mirx);
    Q3mirx = SIMD_DMADD(ews->ewm.xxz.p,xz,Q3mirx);
    Q3mirx = SIMD_DMADD(ews->ewm.xyy.p,yy,Q3mirx);
    Q3mirx = SIMD_DMADD(ews->ewm.xyz.p,yz,Q3mirx);
    Q3mirx = SIMD_DMADD(ews->ewm.xzz.p,zz,Q3mirx);
    Q3miry = SIMD_DMUL(ews->ewm.xxy.p,xx);
    Q3miry = SIMD_DMADD(ews->ewm.xyy.p,xy,Q3miry);
    Q3miry = SIMD_DMADD(ews->ewm.xyz.p,xz,Q3miry);
    Q3miry = SIMD_DMADD(ews->ewm.yyy.p,yy,Q3miry);
    Q3miry = SIMD_DMADD(ews->ewm.yyz.p,yz,Q3miry);
    Q3miry = SIMD_DMADD(ews->ewm.yzz.p,zz,Q3miry);
    Q3mirz = SIMD_DMUL(ews->ewm.xxz.p,xx);
    Q3mirz = SIMD_DMADD(ews->ewm.xyz.p,xy,Q3mirz);
    Q3mirz = SIMD_DMADD(ews->ewm.xzz.p,xz,Q3mirz);
    Q3mirz = SIMD_DMADD(ews->ewm.yyz.p,yy,Q3mirz);
    Q3mirz = SIMD_DMADD(ews->ewm.yzz.p,yz,Q3mirz);
    Q3mirz = SIMD_DMADD(ews->ewm.zzz.p,zz,Q3mirz);
    Q4mirx = SIMD_DMUL(ews->ewm.xxxx.p,xxx);
    Q4mirx = SIMD_DMADD(ews->ewm.xxxy.p,xxy,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xxxz.p,xxz,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xxyy.p,xyy,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xxyz.p,xyz,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xxzz.p,xzz,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xyyy.p,yyy,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xyyz.p,yyz,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xyzz.p,yzz,Q4mirx);
    Q4mirx = SIMD_DMADD(ews->ewm.xzzz.p,zzz,Q4mirx);
    Q4miry = SIMD_DMUL(ews->ewm.xxxy.p,xxx);
    Q4miry = SIMD_DMADD(ews->ewm.xxyy.p,xxy,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.xxyz.p,xxz,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.xyyy.p,xyy,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.xyyz.p,xyz,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.xyzz.p,xzz,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.yyyy.p,yyy,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.yyyz.p,yyz,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.yyzz.p,yzz,Q4miry);
    Q4miry = SIMD_DMADD(ews->ewm.yzzz.p,zzz,Q4miry);
    Q4mirz = SIMD_DMUL(ews->ewm.xxxz.p,xxx);
    Q4mirz = SIMD_DMADD(ews->ewm.xxyz.p,xxy,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.xxzz.p,xxz,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.xyyz.p,xyy,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.xyzz.p,xyz,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.xzzz.p,xzz,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.yyyz.p,yyy,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.yyzz.p,yyz,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.yzzz.p,yzz,Q4mirz);
    Q4mirz = SIMD_DMADD(ews->ewm.zzzz.p,zzz,Q4mirz);
    Q4x = SIMD_DMUL(ews->ewp.Q4xx.p,x);
    Q4x = SIMD_DMADD(ews->ewp.Q4xy.p,y,Q4x);
    Q4x = SIMD_DMADD(ews->ewp.Q4xz.p,z,Q4x);
    Q4y = SIMD_DMUL(ews->ewp.Q4xy.p,x);
    Q4y = SIMD_DMADD(ews->ewp.Q4yy.p,y,Q4y);
    Q4y = SIMD_DMADD(ews->ewp.Q4yz.p,z,Q4y);
    Q4z = SIMD_DMUL(ews->ewp.Q4xz.p,x);
    Q4z = SIMD_DMADD(ews->ewp.Q4yz.p,y,Q4z);
    Q4z = SIMD_DMADD(ews->ewp.Q4zz.p,z,Q4z);
    Q2mir = SIMD_DMUL(Q2mirx,x);
    Q2mir = SIMD_DMADD(Q2miry,y,Q2mir);
    Q2mir = SIMD_DMADD(Q2mirz,z,Q2mir);
    Q2mir = SIMD_DMUL(Q2mir,consts.half.p);
    Q2mir = SIMD_DNMADD(ews->ewp.Q3x.p,x,Q2mir);
    Q2mir = SIMD_DNMADD(ews->ewp.Q3y.p,y,Q2mir);
    Q2mir = SIMD_DNMADD(ews->ewp.Q3z.p,z,Q2mir);
    Q2mir = SIMD_DADD(Q2mir,ews->ewp.Q4.p);
    Q3mir = SIMD_DMUL(Q3mirx,x);
    Q3mir = SIMD_DMADD(Q3miry,y,Q3mir);
    Q3mir = SIMD_DMADD(Q3mirz,z,Q3mir);
    Q3mir = SIMD_DMUL(Q3mir,consts.onethird.p);
    t =             SIMD_DMUL(Q4x,x);
    t = SIMD_DMADD(Q4y,y,t);
    t = SIMD_DMADD(Q4z,z,t);
    Q3mir = SIMD_DNMADD(t,consts.half.p,Q3mir);
    Q4mir = SIMD_DMUL(Q4mirx,x);
    Q4mir = SIMD_DMADD(Q4miry,y,Q4mir);
    Q4mir = SIMD_DMADD(Q4mirz,z,Q4mir);
    Q4mir = SIMD_DMUL(Q4mir,consts.onequarter.p);
    Qta = SIMD_DMUL(g1,ews->ewm.m.p);
    Qta = SIMD_DNMADD(g2,ews->ewp.Q2.p,Qta);
    Qta = SIMD_DMADD(g3,Q2mir,Qta);
    Qta = SIMD_DMADD(g4,Q3mir,Qta);
    Qta = SIMD_DMADD(g5,Q4mir,Qta);

    tpot = SIMD_DMUL(g0,ews->ewm.m.p);
    tpot = SIMD_DNMADD(g1,ews->ewp.Q2.p,tpot);
    tpot = SIMD_DMADD(g2,Q2mir,tpot);
    tpot = SIMD_DMADD(g3,Q3mir,tpot);
    tpot = SIMD_DMADD(g4,Q4mir,tpot);
    tx = SIMD_DMUL(g2,SIMD_DSUB(Q2mirx,ews->ewp.Q3x.p));
    tx = SIMD_DMADD(g3,SIMD_DSUB(Q3mirx,Q4x),tx);
    tx = SIMD_DMADD(g4,Q4mirx,tx);
    tx = SIMD_DNMADD(x,Qta,tx);
    ty = SIMD_DMUL(g2,SIMD_DSUB(Q2miry,ews->ewp.Q3y.p));
    ty = SIMD_DMADD(g3,SIMD_DSUB(Q3miry,Q4y),ty);
    ty = SIMD_DMADD(g4,Q4miry,ty);
    ty = SIMD_DNMADD(y,Qta,ty);
    tz = SIMD_DMUL(g2,SIMD_DSUB(Q2mirz,ews->ewp.Q3z.p));
    tz = SIMD_DMADD(g3,SIMD_DSUB(Q3mirz,Q4z),tz);
    tz = SIMD_DMADD(g4,Q4mirz,tz);
    tz = SIMD_DNMADD(z,Qta,tz);

    dPot = tpot;
    ax = tx;
    ay = ty;
    az = tz;

#endif


    return COST_FLOP_EWALD * SIMD_DWIDTH;
    }
#endif

extern "C"
double pkdParticleEwald(PKD pkd,double *r, float *pa, float *pPot,double *pdFlopSingle, double *pdFlopDouble) {
    struct EwaldVariables *ew = &pkd->ew;
    EwaldTable *ewt = &pkd->ewt;
    const MOMC * restrict mom = &ew->mom;
    double L,Pot,ax,ay,az,dx,dy,dz,x,y,z,r2;
#ifdef USE_SIMD_EWALD
    dvec dPot,dax,day,daz;
    fvec fPot,fax,fay,faz,fx,fy,fz;
    vdouble px, py, pz, pr2,pInHole;
    int nSIMD = 0;
#endif
    int i,ix,iy,iz;
    int bInHole,bInHolex,bInHolexy;
    double dFlopSingle = 0;
    double dFlopDouble = 0;
    int nLoop = 0;

    L = ew->Lbox;
    dx = r[0] - ew->r[0];
    dy = r[1] - ew->r[1];
    dz = r[2] - ew->r[2];

    ax = ay = az = 0.0;
    Pot = mom->m*ew->k1;
#ifdef USE_SIMD_EWALD
    dPot.zero();
    dax.zero();
    day.zero();
    daz.zero();
#endif
    for (ix=-ew->nEwReps;ix<=ew->nEwReps;++ix) {
	bInHolex = (abs(ix) <= ew->nReps);
	x = dx + ix*L;
	for (iy=-ew->nEwReps;iy<=ew->nEwReps;++iy) {
	    bInHolexy = (bInHolex && abs(iy) <= ew->nReps);
	    y = dy + iy*L;
	    for (iz=-ew->nEwReps;iz<=ew->nEwReps;++iz) {
		bInHole = (bInHolexy && abs(iz) <= ew->nReps);
		z = dz + iz*L;
		r2 = x*x + y*y + z*z;
		if (r2 > ew->fEwCut2 && !bInHole) continue;
		if (r2 < ew->fInner2) {
		    double g0,g1,g2,g3,g4,g5,alphan;
		    /*
		     * For small r, series expand about
		     * the origin to avoid errors caused
		     * by cancellation of large terms.
		     */
		    alphan = ew->ka;
		    r2 *= ew->alpha2;
		    g0 = alphan*((1.0/3.0)*r2 - 1.0);
		    alphan *= 2*ew->alpha2;
		    g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
		    alphan *= 2*ew->alpha2;
		    g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
		    alphan *= 2*ew->alpha2;
		    g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
		    alphan *= 2*ew->alpha2;
		    g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
		    alphan *= 2*ew->alpha2;
		    g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));

		    dFlopDouble += evalEwald(ew,&ax,&ay,&az,&Pot,x,y,z,g0,g1,g2,g3,g4,g5);
		    }
		else {
#if defined(USE_SIMD_EWALD)
		    px.d[nSIMD] = x;
		    py.d[nSIMD] = y;
		    pz.d[nSIMD] = z;
		    pr2.d[nSIMD] = r2;
		    pInHole.d[nSIMD] = bInHole;
//		    doerfc.i[nSIMD] = bInHole ? 0 : UINT64_MAX;
		    if (++nSIMD == SIMD_DWIDTH) {
		        dvec tax, tay, taz, tpot;
			dFlopDouble += evalEwaldSIMD(pkd,&pkd->es,tax,tay,taz,tpot,px.p,py.p,pz.p,pr2.p,dvec(pInHole.p) == 0.0);
 			dax += tax;
			day += tay;
			daz += taz;
			dPot-= tpot;
			nSIMD = 0;
			}
#else
		    double dir,dir2,a;
		    double g0,g1,g2,g3,g4,g5,alphan;
		    dir = 1/sqrt(r2);
		    dir2 = dir*dir;
		    a = exp(-r2*ew->alpha2);
		    a *= ew->ka*dir2;


		    if (bInHole) {
			g0 = -erf(ew->alpha*r2*dir);
			}
		    else {
			g0 = erfc(ew->alpha*r2*dir);
			}
		    g0 *= dir;
		    g1 = g0*dir2 + a;
		    alphan = 2*ew->alpha2;
		    g2 = 3*g1*dir2 + alphan*a;
		    alphan *= 2*ew->alpha2;
		    g3 = 5*g2*dir2 + alphan*a;
		    alphan *= 2*ew->alpha2;
		    g4 = 7*g3*dir2 + alphan*a;
		    alphan *= 2*ew->alpha2;
		    g5 = 9*g4*dir2 + alphan*a;
		    dFlopDouble += evalEwald(ew,&ax,&ay,&az,&Pot,x,y,z,g0,g1,g2,g3,g4,g5);
#endif
		    }
		++nLoop;
		}
	    }
	}
#if defined(USE_SIMD_EWALD)
    /* Finish remaining SIMD operations if necessary */
    if (nSIMD) { /* nSIMD can be 0 through 7 */
#define M 0xffffffffffffffff
#if defined(__AVX512F__)
	static const vint64 keepmask[] = {{0,0,0,0,0,0,0,0},{M,0,0,0,0,0,0,0},{M,M,0,0,0,0,0,0},{M,M,M,0,0,0,0,0},
	    {M,M,M,M,0,0,0,0},{M,M,M,M,M,0,0,0},{M,M,M,M,M,M,0,0},{M,M,M,M,M,M,M,0}};
#elif defined(__AVX__)
	static const vint64 keepmask[] = {{0,0,0,0},{M,0,0,0},{M,M,0,0},{M,M,M,0}};
#else
	static const vint64 keepmask[] = {{0,0},{M,0}};
#endif
#undef M
	dvec t, tax, tay, taz, tpot;
	evalEwaldSIMD(pkd,&pkd->es,tax,tay,taz,tpot,px.p,py.p,pz.p,pr2.p,dvec(pInHole.p) == 0.0);
	dFlopDouble += COST_FLOP_EWALD * nSIMD;
	t = keepmask[nSIMD].pd;
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
    fax = SIMD_D2F(dax);
    fay = SIMD_D2F(day);
    faz = SIMD_D2F(daz);
    fPot = SIMD_D2F(dPot);

    fx = SIMD_SPLAT(dx);
    fy = SIMD_SPLAT(dy);
    fz = SIMD_SPLAT(dz);

    nLoop = (ew->nEwhLoop+SIMD_MASK) >> SIMD_BITS;
    i = 0;
    do {
	v_sf hdotx,s,c,t;
	hdotx = SIMD_MADD(ewt->hz.p[i],fz,SIMD_MADD(ewt->hy.p[i],fy,SIMD_MUL(ewt->hx.p[i],fx)));
	fvec svec,cvec;

	sincosf(fvec(hdotx),svec,cvec);
	s = svec; c = cvec;

	fPot = SIMD_ADD(fPot,SIMD_ADD(SIMD_MUL(ewt->hSfac.p[i],s),SIMD_MUL(ewt->hCfac.p[i],c)));
	s = SIMD_MUL(ewt->hCfac.p[i],s);
	c = SIMD_MUL(ewt->hSfac.p[i],c);
	t = SIMD_SUB(s,c);
	fax = SIMD_ADD(fax,SIMD_MUL(ewt->hx.p[i],t));
	fay = SIMD_ADD(fay,SIMD_MUL(ewt->hy.p[i],t));
	faz = SIMD_ADD(faz,SIMD_MUL(ewt->hz.p[i],t));
	} while(++i < nLoop);
    dFlopSingle += ew->nEwhLoop*COST_FLOP_HLOOP;

    ax += hadd(fax);
    ay += hadd(fay);
    az += hadd(faz);
    Pot += hadd(fPot);
#else
    /*
    ** Scoring for the h-loop (+,*)
    ** 	Without trig = (10,14)
    **	    Trig est.	 = 2*(6,11)  same as 1/sqrt scoring.
    **		Total        = (22,36)
    **					 = 58
    */
    for (i=0;i<ew->nEwhLoop;++i) {
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
    dFlopDouble += ew->nEwhLoop*COST_FLOP_HLOOP;
#endif
    pa[0] += ax;
    pa[1] += ay;
    pa[2] += az;
    *pPot += Pot;

    *pdFlopSingle += dFlopSingle;
    *pdFlopDouble += dFlopDouble;
    return dFlopDouble + dFlopSingle;
    }

extern "C"
void pkdEwaldInit(PKD pkd,int nReps,double fEwCut,double fhCut) {
    struct EwaldVariables * const ew = &pkd->ew;
    EwaldTable * const ewt = &pkd->ewt;
    const MOMC * restrict mom = &ew->mom;
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
    pkd->es.ewm.m.p = SIMD_DSPLAT(mom->m);
    pkd->es.ewm.xx.p = SIMD_DSPLAT(mom->xx);
    pkd->es.ewm.yy.p = SIMD_DSPLAT(mom->yy);
    pkd->es.ewm.xy.p = SIMD_DSPLAT(mom->xy);
    pkd->es.ewm.xz.p = SIMD_DSPLAT(mom->xz);
    pkd->es.ewm.yz.p = SIMD_DSPLAT(mom->yz);
    pkd->es.ewm.xxx.p = SIMD_DSPLAT(mom->xxx);
    pkd->es.ewm.xyy.p = SIMD_DSPLAT(mom->xyy);
    pkd->es.ewm.xxy.p = SIMD_DSPLAT(mom->xxy);
    pkd->es.ewm.yyy.p = SIMD_DSPLAT(mom->yyy);
    pkd->es.ewm.xxz.p = SIMD_DSPLAT(mom->xxz);
    pkd->es.ewm.yyz.p = SIMD_DSPLAT(mom->yyz);
    pkd->es.ewm.xyz.p = SIMD_DSPLAT(mom->xyz);
    pkd->es.ewm.xxxx.p = SIMD_DSPLAT(mom->xxxx);
    pkd->es.ewm.xyyy.p = SIMD_DSPLAT(mom->xyyy);
    pkd->es.ewm.xxxy.p = SIMD_DSPLAT(mom->xxxy);
    pkd->es.ewm.yyyy.p = SIMD_DSPLAT(mom->yyyy);
    pkd->es.ewm.xxxz.p = SIMD_DSPLAT(mom->xxxz);
    pkd->es.ewm.yyyz.p = SIMD_DSPLAT(mom->yyyz);
    pkd->es.ewm.xxyy.p = SIMD_DSPLAT(mom->xxyy);
    pkd->es.ewm.xxyz.p = SIMD_DSPLAT(mom->xxyz);
    pkd->es.ewm.xyyz.p = SIMD_DSPLAT(mom->xyyz);
    pkd->es.ewm.zz.p = SIMD_DSPLAT(mom->zz);
    pkd->es.ewm.xzz.p = SIMD_DSPLAT(mom->xzz);
    pkd->es.ewm.yzz.p = SIMD_DSPLAT(mom->yzz);
    pkd->es.ewm.zzz.p = SIMD_DSPLAT(mom->zzz);
    pkd->es.ewm.xxzz.p = SIMD_DSPLAT(mom->xxzz);
    pkd->es.ewm.xyzz.p = SIMD_DSPLAT(mom->xyzz);
    pkd->es.ewm.xzzz.p = SIMD_DSPLAT(mom->xzzz);
    pkd->es.ewm.yyzz.p = SIMD_DSPLAT(mom->yyzz);
    pkd->es.ewm.yzzz.p = SIMD_DSPLAT(mom->yzzz);
    pkd->es.ewm.zzzz.p = SIMD_DSPLAT(mom->zzzz);
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
    ew->nEwReps = d2i(ceil(fEwCut));
    ew->nEwReps = ew->nEwReps > nReps ? ew->nEwReps : nReps;
    ew->fEwCut2 = fEwCut*fEwCut*L*L;
    ew->fInner2 = 1.2e-3*L*L;
    ew->alpha = 2.0/L;
    ew->ialpha = 0.5 * L;
    ew->alpha2 = ew->alpha*ew->alpha;
    ew->k1 = M_PI/(ew->alpha2*L*L*L);
    ew->ka = 2.0*ew->alpha/sqrt(M_PI);
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    pkd->es.ewp.Q4xx.p = SIMD_DSPLAT(ew->Q4xx);
    pkd->es.ewp.Q4xy.p = SIMD_DSPLAT(ew->Q4xy);
    pkd->es.ewp.Q4xz.p = SIMD_DSPLAT(ew->Q4xz);
    pkd->es.ewp.Q4yy.p = SIMD_DSPLAT(ew->Q4yy);
    pkd->es.ewp.Q4yz.p = SIMD_DSPLAT(ew->Q4yz);
    pkd->es.ewp.Q4zz.p = SIMD_DSPLAT(ew->Q4zz);
    pkd->es.ewp.Q4.p = SIMD_DSPLAT(ew->Q4);
    pkd->es.ewp.Q3x.p = SIMD_DSPLAT(ew->Q3x);
    pkd->es.ewp.Q3y.p = SIMD_DSPLAT(ew->Q3y);
    pkd->es.ewp.Q3z.p = SIMD_DSPLAT(ew->Q3z);
    pkd->es.ewp.Q2.p = SIMD_DSPLAT(ew->Q2);
    pkd->es.ewp.fEwCut2.p = SIMD_DSPLAT(ew->fEwCut2);
    pkd->es.ewp.fInner2.p = SIMD_DSPLAT(ew->fInner2);
    pkd->es.ewp.alpha.p = SIMD_DSPLAT(ew->alpha);
    pkd->es.ewp.ialpha.p = SIMD_DSPLAT(ew->ialpha);
    pkd->es.ewp.alpha2.p = SIMD_DSPLAT(ew->alpha2);
    pkd->es.ewp.k1.p = SIMD_DSPLAT(ew->k1);
    pkd->es.ewp.ka.p = SIMD_DSPLAT(ew->ka);
#endif


    /*
    ** Now setup stuff for the h-loop.
    */
    hReps = d2i(ceil(fhCut));
    k4 = M_PI*M_PI/(ew->alpha*ew->alpha*L*L);

    i = (int)pow(1+2*hReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE__)
    i = (i + SIMD_MASK) & ~SIMD_MASK;
#endif
    if ( i>ew->nMaxEwhLoop ) {
	ew->nMaxEwhLoop = i;
	ewt->hx.f = (float *)SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hx.f));
	assert(ewt->hx.f != NULL);
	ewt->hy.f = (float *)SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hy.f));
	assert(ewt->hy.f != NULL);
	ewt->hz.f = (float *)SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hz.f));
	assert(ewt->hz.f != NULL);
	ewt->hCfac.f = (float *)SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hCfac.f));
	assert(ewt->hCfac.f != NULL);
	ewt->hSfac.f = (float *)SIMD_malloc(ew->nMaxEwhLoop*sizeof(ewt->hSfac.f));
	assert(ewt->hSfac.f != NULL);
	}
    ew->nEwhLoop = i;
    i = (int)pow(1+2*ew->nEwReps,3);
#if defined(USE_SIMD_EWALD) && defined(__SSE2__)
    i = (i + SIMD_MASK) & ~SIMD_MASK;
#endif
    i = 0;
    for (hx=-hReps;hx<=hReps;++hx) {
	for (hy=-hReps;hy<=hReps;++hy) {
	    for (hz=-hReps;hz<=hReps;++hz) {
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
    while(i<ew->nMaxEwhLoop) {
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
    cudaEwaldInit(pkd->mdl->cudaCtx,ew,ewt);
    mdlThreadBarrier(pkd->mdl);
#endif


    }
