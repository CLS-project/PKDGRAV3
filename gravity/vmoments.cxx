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
#include <stdint.h>
#include <math.h>
#include "ilc.h"

#include "core/simd.h"

template<typename F>
struct ResultLOCR {
    F m;
    F x,y,z;
    F xx,yy,xy,xz,yz;
    F xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    F xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
    F xxxxx,xyyyy,xxxxy,yyyyy,xxxxz,yyyyz,xxxyy,xxyyy,xxxyz,xyyyz,xxyyz;
    F vdirLsum,vnormLsum;
    void zero() {memset(this,0,sizeof(ResultLOCR));}
    ResultLOCR operator+=(const ResultLOCR &rhs) {
        m += rhs.m;
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        xx += rhs.xx;
        yy += rhs.yy;
        xy += rhs.xy;
        xz += rhs.xz;
        yz += rhs.yz;
        xxx += rhs.xxx;
        xyy += rhs.xyy;
        xxy += rhs.xxy;
        yyy += rhs.yyy;
        xxz += rhs.xxz;
        yyz += rhs.yyz;
        xyz += rhs.xyz;
        xxxx += rhs.xxxx;
        xyyy += rhs.xyyy;
        xxxy += rhs.xxxy;
        yyyy += rhs.yyyy;
        xxxz += rhs.xxxz;
        yyyz += rhs.yyyz;
        xxyy += rhs.xxyy;
        xxyz += rhs.xxyz;
        xyyz += rhs.xyyz;
        xxxxx += rhs.xxxxx;
        xyyyy += rhs.xyyyy;
        xxxxy += rhs.xxxxy;
        yyyyy += rhs.yyyyy;
        xxxxz += rhs.xxxxz;
        yyyyz += rhs.yyyyz;
        xxxyy += rhs.xxxyy;
        xxyyy += rhs.xxyyy;
        xxxyz += rhs.xxxyz;
        xyyyz += rhs.xyyyz;
        xxyyz += rhs.xxyyz;
        return *this;
    }
};

template<typename T>
ResultLOCR<T> FlocrSetVmomr5cm(T v,T u,T x,T y,T z,T ax,T ay,T az,T m,T imaga,
                               T Mxx,T Myy,T Mxy,T Mxz,T Myz,
                               T Mxxx,T Mxyy,T Mxxy,T Myyy,T Mxxz,T Myyz,T Mxyz,
                               T Mxxxx,T Mxyyy,T Mxxxy,T Myyyy,T Mxxxz,T Myyyz,T Mxxyy,T Mxxyz,T Mxyyz) {
    ResultLOCR<fvec> result;
    const T onethird = 1.0 / 3.0;

    auto dir = rsqrt(x*x + y*y + z*z);
    auto sdir = dir;
    u *= dir;
    x *= dir;
    y *= dir;
    z *= dir;
    dir = -dir;
    v *= dir;
    auto u2 = 15.0f*u*u;  // becomes 15.0f*u2!

    auto xx = 0.5f*x*x;
    auto xy = x*y;
    auto yy = 0.5f*y*y;
    auto xz = x*z;
    auto yz = y*z;
    auto zz = 0.5f*z*z;
    auto xxx = x*(onethird*xx - zz);
    auto xxz = z*(xx - onethird*zz);
    auto yyy = y*(onethird*yy - zz);
    auto yyz = z*(yy - onethird*zz);
    xx -= zz;
    yy -= zz;
    auto xxy = y*xx;
    auto xyy = x*yy;
    auto xyz = xy*z;

    auto u3 = u2*u;  // becomes 5.0f*u3!

    auto R2xx = u2*Mxx;
    auto R2xy = u2*Mxy;
    auto R2xz = u2*Mxz;
    auto R2yy = u2*Myy;
    auto R2yz = u2*Myz;

    auto u4 = 7.0f*u3*u;  // becomes 7.0f*5.0f*u4!

    auto R2x = x*R2xx + y*R2xy + z*R2xz;
    auto R2y = x*R2xy + y*R2yy + z*R2yz;
    auto R2z = x*R2xz + y*R2yz - z*(R2xx + R2yy);

    auto R3xx = u3*(x*Mxxx + y*Mxxy + z*Mxxz);
    auto R3xy = u3*(x*Mxxy + y*Mxyy + z*Mxyz);
    auto R3yy = u3*(x*Mxyy + y*Myyy + z*Myyz);
    auto R3xz = u3*(x*Mxxz + y*Mxyz - z*(Mxxx + Mxyy));
    auto R3yz = u3*(x*Mxyz + y*Myyz - z*(Mxxy + Myyy));

    auto R4x = u4*(Mxxxx*xxx + Mxyyy*yyy + Mxxxy*xxy + Mxxxz*xxz + Mxxyy*xyy + Mxxyz*xyz + Mxyyz*yyz);
    auto R4y = u4*(Mxyyy*xyy + Mxxxy*xxx + Myyyy*yyy + Myyyz*yyz + Mxxyy*xxy + Mxxyz*xxz + Mxyyz*xyz);
    auto R4z = u4*(-Mxxxx*xxz - (Mxyyy + Mxxxy)*xyz - Myyyy*yyz + Mxxxz*xxx + Myyyz*yyy - Mxxyy*(xxz + yyz) + Mxxyz*xxy + Mxyyz*xyy);

    auto R3x = 0.5f*(x*R3xx + y*R3xy + z*R3xz);
    auto R3y = 0.5f*(x*R3xy + y*R3yy + z*R3yz);
    auto R3z = 0.5f*(x*R3xz + y*R3yz - z*(R3xx + R3yy));

    auto R4 = 0.25f*(x*R4x + y*R4y + z*R4z);

    auto R2 = 0.5f*(x*R2x + y*R2y + z*R2z);

    auto R3 = onethird*(x*R3x + y*R3y + z*R3z);

    xx = x*x;
    yy = y*y;

    //Now we use the 'R's.
    result.m = dir*(m + 0.2f*R2 + R3 + R4);

    dir *= v;
    auto T0 = -(m + R2 + 7.0f*R3 + 9.0f*R4);

    auto vax = dir*(T0*x + 0.2f*R2x + R3x + R4x);
    auto vay = dir*(T0*y + 0.2f*R2y + R3y + R4y);
    auto vaz = dir*(T0*z + 0.2f*R2z + R3z + R4z);

    auto adotai = ax*vax + ay*vay + az*vaz;
    adotai = maskz_mov(adotai > 0.0f,adotai);
    adotai *= imaga;
    auto vd2 = adotai * adotai;
    result.vdirLsum = sdir * vd2;
    result.vnormLsum = vd2;

    result.x = -vax;
    result.y = -vay;
    result.z = -vaz;

    dir *= v;
    T0 = 3.0f*m + 7.0f*(R2 + 9.0f*R3);

    auto t1 = m + R2 + 7.0f*R3;
    auto t1x = R2x + 7.0f*R3x;
    auto t1y = R2y + 7.0f*R3y;
    auto t1z = R2z + 7.0f*R3z;
    result.xx = dir*(T0*xx - t1 - 2.0f*x*t1x + 0.2f*R2xx + R3xx);
    result.yy = dir*(T0*yy - t1 - 2.0f*y*t1y + 0.2f*R2yy + R3yy);
    result.xy = dir*(T0*xy - y*t1x - x*t1y + 0.2f*R2xy + R3xy);
    result.xz = dir*(T0*xz - z*t1x - x*t1z + 0.2f*R2xz + R3xz);
    result.yz = dir*(T0*yz - z*t1y - y*t1z + 0.2f*R2yz + R3yz);

    dir *= v;
    T0 = 15.0f*m + 63.0f*R2;
    auto txx = T0*xx;
    auto tyy = T0*yy;

    t1 = 3.0f*m + 7.0f*R2;
    auto t2x = -7.0f*R2x;
    auto t2y = -7.0f*R2y;
    auto t2z = -7.0f*R2z;
    auto t1xx = txx - t1 + 2.0f*x*t2x + R2xx;
    auto t1yy = tyy - t1 + 2.0f*y*t2y + R2yy;

    result.xxx = dir*(x*(txx - 3.0f*(t1 - t2x*x - R2xx)) + 3.0f*R2x);
    result.yyy = dir*(y*(tyy - 3.0f*(t1 - t2y*y - R2yy)) + 3.0f*R2y);
    result.xxy = dir*(y*t1xx + xx*t2y + R2y + 2.0f*R2xy*x);
    result.xxz = dir*(z*t1xx + xx*t2z + R2z + 2.0f*R2xz*x);
    result.xyy = dir*(x*t1yy + yy*t2x + R2x + 2.0f*R2xy*y);
    result.yyz = dir*(z*t1yy + yy*t2z + R2z + 2.0f*R2yz*y);
    result.xyz = dir*(T0*xyz + (yz*t2x + xz*t2y + xy*t2z) + R2xy*z + R2yz*x + R2xz*y);

    dir *= v*m;
    txx = 105.0f*xx;
    tyy = 105.0f*yy;
    auto t2xx = txx - 90.0f;
    auto t2yy = tyy - 90.0f;
    result.xxxx = dir*(xx*t2xx + 9.0f);
    result.yyyy = dir*(yy*t2yy + 9.0f);
    t2xx += 45.0f;
    t2yy += 45.0f;
    result.xxxy = dir*xy*t2xx;
    result.xxxz = dir*xz*t2xx;
    result.xyyy = dir*xy*t2yy;
    result.yyyz = dir*yz*t2yy;
    t2xx += 30.0f;
    t2yy += 30.0f;
    result.xxyy = dir*(yy*t2xx - xx*15.0f + 3.0f);
    result.xxyz = dir*(yz*t2xx);
    result.xyyz = dir*(xz*t2yy);

    dir *= v;
    x *= dir;
    y *= dir;
    z *= dir;
    txx = 9.0f*xx - 10.0f;
    tyy = 9.0f*yy - 10.0f;
    xx *= 105.0f;
    yy *= 105.0f;
    xy *= z*105.0f;
    result.xxxxx = x*(xx*txx + 225.0f);
    result.yyyyy = y*(yy*tyy + 225.0f);
    txx += 4.0f;
    tyy += 4.0f;
    auto txxxx = xx*txx + 45.0f;
    auto tyyyy = yy*tyy + 45.0f;
    result.xxxxy = y*txxxx;
    result.xxxxz = z*txxxx;
    result.xyyyy = x*tyyyy;
    result.yyyyz = z*tyyyy;
    txx += 3.0f;
    tyy += 3.0f;
    result.xxxyz = xy*txx;
    result.xyyyz = xy*tyy;
    result.xxxyy = x*(yy*txx - xx + 45.0f);
    result.xxyyy = y*(xx*tyy - yy + 45.0f);
    tyy += 2.0f;
    result.xxyyz = z*(xx*tyy - yy + 15.0f);

    return result;
}

template<typename T>
auto EvalBlock(int n,ilcBlock &blk,T v,T ax,T ay,T az,T imaga) {
    // Sentinal values
    while (n&fvec::mask()) {
        blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e18f;
        blk.m.s[n] = 0.0f;
        blk.u.s[n] = 0.0f;
        ++n;
    }
    n /= fvec::width(); // Now number of blocks
    ResultLOCR<T> result;
    result.zero();
    for (auto i=0; i<n; ++i) {
        result += FlocrSetVmomr5cm(v,blk.u.v[i],blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],ax,ay,az,blk.m.v[i],imaga,
                                   blk.xx.v[i],blk.yy.v[i],blk.xy.v[i],blk.xz.v[i],blk.yz.v[i],
                                   blk.xxx.v[i],blk.xyy.v[i],blk.xxy.v[i],blk.yyy.v[i],blk.xxz.v[i],blk.yyz.v[i],blk.xyz.v[i],
                                   blk.xxxx.v[i],blk.xyyy.v[i],blk.xxxy.v[i],blk.yyyy.v[i],blk.xxxz.v[i],blk.yyyz.v[i],blk.xxyy.v[i],blk.xxyz.v[i],blk.xyyz.v[i]);
    }
    return result;
}

double momFlocrSetVFmomr5cm(FLOCR *l,float v1,ilcList &ill,const float *a,float *pfdirLsum, float *pfnormLsum) {
    float dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0.0) {
        dimaga = 1.0f / sqrtf(dimaga);
    }
    fvec imaga = dimaga;

    ResultLOCR<fvec> result;
    result.zero();

    fvec v(v1), ax(a[0]), ay(a[1]), az(a[2]);
    for (auto tile : ill) {
        auto nBlocks = tile.count() / tile.width;
        for (auto iBlock=0; iBlock<nBlocks; ++iBlock) {
            result += EvalBlock(tile.width,tile[iBlock],v,ax,ay,az,imaga);
        }
        result += EvalBlock(tile.count() - nBlocks*tile.width,tile[nBlocks],v,ax,ay,az,imaga);
    }

    // 32 terms * 4 instructions = 128 operations
    l->m = hadd(result.m);
    l->x = hadd(result.x);
    l->y = hadd(result.y);
    l->z = hadd(result.z);
    l->xx = hadd(result.xx);
    l->yy = hadd(result.yy);
    l->xy = hadd(result.xy);
    l->xz = hadd(result.xz);
    l->yz = hadd(result.yz);
    l->xxx = hadd(result.xxx);
    l->xyy = hadd(result.xyy);
    l->xxy = hadd(result.xxy);
    l->yyy = hadd(result.yyy);
    l->xxz = hadd(result.xxz);
    l->yyz = hadd(result.yyz);
    l->xyz = hadd(result.xyz);
    l->xxxx = hadd(result.xxxx);
    l->xyyy = hadd(result.xyyy);
    l->xxxy = hadd(result.xxxy);
    l->yyyy = hadd(result.yyyy);
    l->xxxz = hadd(result.xxxz);
    l->yyyz = hadd(result.yyyz);
    l->xxyy = hadd(result.xxyy);
    l->xxyz = hadd(result.xxyz);
    l->xyyz = hadd(result.xyyz);
    l->xxxxx = hadd(result.xxxxx);
    l->xyyyy = hadd(result.xyyyy);
    l->xxxxy = hadd(result.xxxxy);
    l->yyyyy = hadd(result.yyyyy);
    l->xxxxz = hadd(result.xxxxz);
    l->yyyyz = hadd(result.yyyyz);
    l->xxxyy = hadd(result.xxxyy);
    l->xxyyy = hadd(result.xxyyy);
    l->xxxyz = hadd(result.xxxyz);
    l->xyyyz = hadd(result.xyyyz);
    l->xxyyz = hadd(result.xxyyz);

    *pfdirLsum += hadd(result.vdirLsum);
    *pfnormLsum += hadd(result.vnormLsum);

    // OpCount(*,+,-) = (269,160,43,1) = 472.0 + 7.0
    // mov: 147.0
    return (479.0 * ill.count());
}