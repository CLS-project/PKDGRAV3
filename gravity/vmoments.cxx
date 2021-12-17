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
#include <immintrin.h>
#include "ilc.h"

#include "core/simd.h"

double momFlocrSetVFmomr5cm(FLOCR *l,float v1,ILC ill,const float *a,float *pfdirLsum, float *pfnormLsum) {
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
    fvec vdirLsum, vnormLsum;

    float dimaga = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (dimaga > 0.0) {
        dimaga = 1.0f / sqrtf(dimaga);
    }
    fvec vimaga = dimaga;

    llm.zero();
    llx.zero();
    lly.zero();
    llz.zero();
    llxx.zero();
    llyy.zero();
    llxy.zero();
    llxz.zero();
    llyz.zero();
    llxxx.zero();
    llxyy.zero();
    llxxy.zero();
    llyyy.zero();
    llxxz.zero();
    llyyz.zero();
    llxyz.zero();
    llxxxx.zero();
    llxyyy.zero();
    llxxxy.zero();
    llyyyy.zero();
    llxxxz.zero();
    llyyyz.zero();
    llxxyy.zero();
    llxxyz.zero();
    llxyyz.zero();
    llxxxxx.zero();
    llxyyyy.zero();
    llxxxxy.zero();
    llyyyyy.zero();
    llxxxxz.zero();
    llyyyyz.zero();
    llxxxyy.zero();
    llxxyyy.zero();
    llxxxyz.zero();
    llxyyyz.zero();
    llxxyyz.zero();

    vdirLsum.zero();
    vnormLsum.zero();

    int nLoops = 0;
    ILCTILE tile;
    ILC_LOOP(ill,tile) {
        ILC_BLK *blk = tile->blk;
        int nBlocks = tile->lstTile.nBlocks;
        int nInLast = tile->lstTile.nInLast;
        int nLeft;
        int j;
        for ( j = nInLast; j&fvec::mask(); j++) {
            blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
            blk[nBlocks].m.f[j] = 0.0f;
            blk[nBlocks].u.f[j] = 0.0f;
        }

        for ( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
            int n = ((nLeft ? ILC_PART_PER_BLK : nInLast) + fvec::mask()) / fvec::width();
            for (j=0; j<n; ++j) {
                fvec
                x = fvec(blk->dx.p[j]),
                y = fvec(blk->dy.p[j]),
                z = fvec(blk->dz.p[j]),
                u = fvec(blk->u.p[j]),
                m = fvec(blk->m.p[j]);
                fvec v=v1;
                fvec dir = rsqrt(x*x + y*y + z*z);
                fvec sdir = dir;

                ++nLoops;

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

                R2xx = u2*fvec(blk->xx.p[j]);
                R2xy = u2*fvec(blk->xy.p[j]);
                R2xz = u2*fvec(blk->xz.p[j]);
                R2yy = u2*fvec(blk->yy.p[j]);
                R2yz = u2*fvec(blk->yz.p[j]);

                u4 = 7.0f*u3*u;  /* becomes 7.0f*5.0f*u4! */

                R2x = x*R2xx + y*R2xy + z*R2xz;
                R2y = x*R2xy + y*R2yy + z*R2yz;
                R2z = x*R2xz + y*R2yz - z*(R2xx + R2yy);

                R3xx = u3*(x*fvec(blk->xxx.p[j]) + y*fvec(blk->xxy.p[j]) + z*fvec(blk->xxz.p[j]));
                R3xy = u3*(x*fvec(blk->xxy.p[j]) + y*fvec(blk->xyy.p[j]) + z*fvec(blk->xyz.p[j]));
                R3yy = u3*(x*fvec(blk->xyy.p[j]) + y*fvec(blk->yyy.p[j]) + z*fvec(blk->yyz.p[j]));
                R3xz = u3*(x*fvec(blk->xxz.p[j]) + y*fvec(blk->xyz.p[j]) - z*(fvec(blk->xxx.p[j]) + fvec(blk->xyy.p[j])));
                R3yz = u3*(x*fvec(blk->xyz.p[j]) + y*fvec(blk->yyz.p[j]) - z*(fvec(blk->xxy.p[j]) + fvec(blk->yyy.p[j])));

                R4x = u4*(fvec(blk->xxxx.p[j])*xxx + fvec(blk->xyyy.p[j])*yyy + fvec(blk->xxxy.p[j])*xxy + fvec(blk->xxxz.p[j])*xxz + fvec(blk->xxyy.p[j])*xyy + fvec(blk->xxyz.p[j])*xyz + fvec(blk->xyyz.p[j])*yyz);
                R4y = u4*(fvec(blk->xyyy.p[j])*xyy + fvec(blk->xxxy.p[j])*xxx + fvec(blk->yyyy.p[j])*yyy + fvec(blk->yyyz.p[j])*yyz + fvec(blk->xxyy.p[j])*xxy + fvec(blk->xxyz.p[j])*xxz + fvec(blk->xyyz.p[j])*xyz);
                R4z = u4*(-fvec(blk->xxxx.p[j])*xxz - (fvec(blk->xyyy.p[j]) + fvec(blk->xxxy.p[j]))*xyz - fvec(blk->yyyy.p[j])*yyz + fvec(blk->xxxz.p[j])*xxx + fvec(blk->yyyz.p[j])*yyy - fvec(blk->xxyy.p[j])*(xxz + yyz) + fvec(blk->xxyz.p[j])*xxy + fvec(blk->xyyz.p[j])*xyy);


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

                fvec adotai = a[0]*vax + a[1]*vay + a[2]*vaz;
                adotai = maskz_mov(adotai > 0.0f,adotai);
//      adotai &= adotai > 0.0f;
                adotai *= vimaga;
                fvec vd2 = adotai * adotai;
                vdirLsum += sdir * vd2;
                vnormLsum += vd2;

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

    *pfdirLsum += hadd(vdirLsum);
    *pfnormLsum += hadd(vnormLsum);

    // OpCount(*,+,-) = (269,160,43,1) = 472.0 + 7.0
    // mov: 147.0
    return (479.0 * 8 * nLoops);
}
