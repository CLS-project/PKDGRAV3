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
#ifndef xxUSE_SIMD_PC
#define MPICH_SKIP_MPICXX
#include "core/simd.h"
#include "pkd.h"
#include "pc.h"

template<typename BLOCK> struct ilist::EvalBlock<ResultPC<fvec>,BLOCK> {
    typedef ResultPC<fvec> result_type;
    const fvec fx,fy,fz,pSmooth2,Pax,Pay,Paz,imaga;

    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy,fvec fz,fvec pSmooth2,fvec Pax,fvec Pay,fvec Paz,fvec imaga)
        : fx(fx),fy(fy),fz(fz),pSmooth2(pSmooth2),Pax(Pax),Pay(Pay),Paz(Paz),imaga(imaga) {}

    result_type operator()(int n,BLOCK &blk) {
        // Sentinal values
        while (n&fvec::mask()) {
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e18f;
            blk.m.s[n] = 0.0f;
            blk.u.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        ResultPC<fvec> result;
        result.zero();
        for (auto i=0; i<n; ++i) {
            result += EvalPC<fvec,fmask,true>(fx, fy, fz, pSmooth2,
                                              blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],blk.m.v[i],blk.u.v[i],
                                              blk.xxxx.v[i], blk.xxxy.v[i], blk.xxxz.v[i], blk.xxyz.v[i], blk.xxyy.v[i],
                                              blk.yyyz.v[i], blk.xyyz.v[i], blk.xyyy.v[i], blk.yyyy.v[i],
                                              blk.xxx.v[i], blk.xyy.v[i], blk.xxy.v[i], blk.yyy.v[i], blk.xxz.v[i], blk.yyz.v[i], blk.xyz.v[i],
                                              blk.xx.v[i], blk.xy.v[i], blk.xz.v[i], blk.yy.v[i], blk.yz.v[i],
#ifdef USE_DIAPOLE
                                              blk.x.v[i], blk.y.v[i], blk.z.v[i],
#endif
                                              Pax, Pay, Paz,imaga);
        }
        return result;
    }
};

void pkdGravEvalPC(const PINFOIN &Part, ilcTile &tile,  PINFOOUT &Out ) {
    const float *a = Part.a;
    float a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    fvec imaga = a2 > 0.0f ? 1.0f / sqrtf(a2) : 0.0f;
    ilist::EvalBlock<ResultPC<fvec>,ilcBlock> eval(
        Part.r[0],Part.r[1],Part.r[2],Part.fSmooth2,a[0],a[1],a[2],imaga);;
    auto result = EvalTile(tile,eval);
    Out.a[0] += hadd(result.ax);
    Out.a[1] += hadd(result.ay);
    Out.a[2] += hadd(result.az);
    Out.fPot += hadd(result.pot);
    Out.dirsum += hadd(result.ir);
    Out.normsum += hadd(result.norm);
}
#endif/*USE_SIMD_PC*/
