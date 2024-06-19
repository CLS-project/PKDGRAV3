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
#ifndef xxxUSE_SIMD_PP
#define MPICH_SKIP_MPICXX
#include "core/simd.h"
#include "pkd.h"
#include "pp.h"
#include <algorithm>

template<typename BLOCK> struct ilist::EvalBlock<ResultPP<fvec>,BLOCK> {
    typedef ResultPP<fvec> result_type;
    const fvec fx,fy,fz,pSmooth2,Pax,Pay,Paz,imaga;

    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy,fvec fz,fvec pSmooth2,fvec Pax,fvec Pay,fvec Paz,fvec imaga)
        : fx(fx),fy(fy),fz(fz),pSmooth2(pSmooth2),Pax(Pax),Pay(Pay),Paz(Paz),imaga(imaga) {}

    result_type operator()(int n,BLOCK &blk) {
        // Sentinal values
        while (n&fvec::mask()) {
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e18f;
            blk.m.s[n] = 0.0f;
            blk.fourh2.s[n] = 1e-18f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i=0; i<n; ++i) {
            result += EvalPP<fvec,fmask>(fx,fy,fz,pSmooth2,blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],blk.fourh2.v[i],blk.m.v[i],Pax,Pay,Paz,imaga);
        }
        return result;
    }
};

void pkdGravEvalPP(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out ) {
    float a2 = blitz::dot(Part.a,Part.a);
    fvec imaga = a2 > 0.0f ? 1.0f / sqrtf(a2) : 0.0f;

    ilist::EvalBlock<ResultPP<fvec>,ilpBlock> eval(
        Part.r[0],Part.r[1],Part.r[2],Part.fSmooth2,Part.a[0],Part.a[1],Part.a[2],imaga);;
    auto result = EvalTile(tile,eval);
    Out.a[0] += hadd(result.ax);
    Out.a[1] += hadd(result.ay);
    Out.a[2] += hadd(result.az);
    Out.fPot += hadd(result.pot);
    Out.dirsum += hadd(result.ir);
    Out.normsum += hadd(result.norm);
}

template<typename BLOCK> struct ilist::EvalBlock<ResultDensity<fvec>,BLOCK> {
    typedef ResultDensity<fvec> result_type;
    const fvec fx,fy,fz,fBall,iMat;
    const SPHOptions *const SPHoptions;

    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy,fvec fz,fvec fBall,fvec iMat,SPHOptions *SPHoptions)
        : fx(fx),fy(fy),fz(fz),fBall(fBall),iMat(iMat),SPHoptions(SPHoptions) {}

    result_type operator()(int n,BLOCK &blk) {
        // Sentinal values
        while (n&fvec::mask()) {
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e14f;
            blk.m.s[n] = 0.0f;
            blk.iMat.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i=0; i<n; ++i) {
            result += EvalDensity<fvec,fmask>(fx,fy,fz,fBall,iMat,1.0f,
                                              blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],blk.m.v[i],blk.iMat.v[i],1.0f,
                                              SPHoptions->kernelType,SPHoptions->doInterfaceCorrection);
        }
        return result;
    }
};

void pkdDensityEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions ) {
    ilist::EvalBlock<ResultDensity<fvec>,BlockPP<ILC_PART_PER_BLK>> eval(
                Part.r[0],Part.r[1],Part.r[2],Part.fBall,Part.iMat,SPHoptions);
    auto result = EvalTile(tile,eval);
    Out.rho += hadd(result.rho);
    Out.drhodfball += hadd(result.drhodfball);
    Out.nden += hadd(result.nden);
    Out.dndendfball += hadd(result.dndendfball);
    Out.nSmooth += hadd(result.nSmooth);
    Out.imbalanceX += hadd(result.imbalanceX);
    Out.imbalanceY += hadd(result.imbalanceY);
    Out.imbalanceZ += hadd(result.imbalanceZ);
}

template<typename BLOCK> struct ilist::EvalBlock<ResultDensityCorrection<fvec>,BLOCK> {
    typedef ResultDensityCorrection<fvec> result_type;
    const fvec fx,fy,fz,fBall;
    const SPHOptions *const SPHoptions;

    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy,fvec fz,fvec fBall,SPHOptions *SPHoptions)
        : fx(fx),fy(fy),fz(fz),fBall(fBall),SPHoptions(SPHoptions) {}

    result_type operator()(int n,BLOCK &blk) {
        // Sentinal values
        while (n&fvec::mask()) {
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e14f;
            blk.T.s[n] = 0.0f;
            blk.P.s[n] = 0.0f;
            blk.expImb2.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i=0; i<n; ++i) {
            result += EvalDensityCorrection<fvec,fmask>(fx,fy,fz,fBall,1.0f,
                      blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],blk.T.v[i],blk.P.v[i],blk.expImb2.v[i],1.0f,
                      SPHoptions->kernelType);
        }
        return result;
    }
};

void pkdDensityCorrectionEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions ) {
    ilist::EvalBlock<ResultDensityCorrection<fvec>,BlockPP<ILC_PART_PER_BLK>> eval(
                Part.r[0],Part.r[1],Part.r[2],Part.fBall,SPHoptions);
    auto result = EvalTile(tile,eval);
    Out.corrT += hadd(result.corrT);
    Out.corrP += hadd(result.corrP);
    Out.corr += hadd(result.corr);
}

template<typename BLOCK, bool doShearStrengthModel> struct ilist::EvalBlock<ResultSPHForces<fvec,doShearStrengthModel>,BLOCK> {
    typedef ResultSPHForces<fvec,doShearStrengthModel> result_type;
    const fvec fx,fy,fz,fBall,Omega,vx,vy,vz,rho,P,c;
    const fvec Sxx, Syy, Sxy, Sxz, Syz;
    const SPHOptions *const SPHoptions;
    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy,fvec fz,fvec fBall,fvec Omega,fvec vx,fvec vy,fvec vz,
              fvec rho,fvec P,fvec c,
              fvec Sxx, fvec Syy, fvec Sxy, fvec Sxz, fvec Syz, SPHOptions *SPHoptions)
        : fx(fx),fy(fy),fz(fz),fBall(fBall),Omega(Omega),vx(vx),vy(vy),vz(vz),
          rho(rho),P(P),c(c),
          Sxx(Sxx),Syy(Syy),Sxy(Sxy),Sxz(Sxz),Syz(Syz), SPHoptions(SPHoptions) {}

    result_type operator()(int n,BLOCK &blk) {
        // Sentinal values
        while (n&fvec::mask()) {
            // This is a little trick to speed up the calculation. By setting
            // unused entries in the list to have a zero mass, the resulting
            // forces are zero. Setting the distance to a large value avoids
            // softening the non-existent forces which is slightly faster.
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e14f;
            blk.fBall.s[n] = blk.Omega.s[n] = blk.rho.s[n] = 1.0f;
            blk.vx.s[n] = blk.vy.s[n] = blk.vz.s[n] = 0.0f;
            blk.m.s[n] = 0.0f;
            blk.P.s[n] = 0.0f;
            blk.c.s[n] = 0.0f;
            blk.uRung.s[n] = 0.0f;
            blk.Sxx.s[n] = 0.0f;
            blk.Syy.s[n] = 0.0f;
            blk.Sxy.s[n] = 0.0f;
            blk.Sxz.s[n] = 0.0f;
            blk.Syz.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i=0; i<n; ++i) {
            result += EvalSPHForces<fvec,fmask,doShearStrengthModel>(
                          fx,fy,fz,fBall,Omega,
                          vx,vy,vz,rho,P,c,1.0f,
                          Sxx, Syy, Sxy, Sxz, Syz,
                          blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],blk.m.v[i],blk.fBall.v[i],blk.Omega.v[i],
                          blk.vx.v[i],blk.vy.v[i],blk.vz.v[i],blk.rho.v[i],blk.P.v[i],blk.c.v[i],blk.uRung.v[i],1.0f,
                          blk.Sxx.v[i], blk.Syy.v[i], blk.Sxy.v[i], blk.Sxz.v[i], blk.Syz.v[i],
                          SPHoptions->kernelType,SPHoptions->epsilon,SPHoptions->alpha,SPHoptions->beta,
                          SPHoptions->EtaCourant,SPHoptions->a,SPHoptions->H,SPHoptions->useIsentropic);
        }
        return result;
    }
};

template<bool doShearStrengthModel>
void combineSPHForcesResult(PINFOOUT &Out, ResultSPHForces<fvec,doShearStrengthModel> &result, SPHOptions *SPHoptions) {
    Out.uDot += hadd(result.uDot);
    Out.a[0] += hadd(result.ax);
    Out.a[1] += hadd(result.ay);
    Out.a[2] += hadd(result.az);
    Out.divv += hadd(result.divv);
    // This should be a horizontal minimum for an fvec, resulting in a float containing the smallest float in the fvec
    for (int k = 0; k < result.dtEst.width(); k++) {
        Out.dtEst = std::min(Out.dtEst,result.dtEst[k]);
    }
    for (int k = 0; k < result.maxRung.width(); k++) {
        Out.maxRung = std::max(Out.maxRung,result.maxRung[k]);
    }
    assert(Out.dtEst > 0);
    if (doShearStrengthModel) {
        Out.dvxdx += hadd(result.dvxdx);
        Out.dvxdy += hadd(result.dvxdy);
        Out.dvxdz += hadd(result.dvxdz);
        Out.dvydx += hadd(result.dvydx);
        Out.dvydy += hadd(result.dvydy);
        Out.dvydz += hadd(result.dvydz);
        Out.dvzdx += hadd(result.dvzdx);
        Out.dvzdy += hadd(result.dvzdy);
        Out.dvzdz += hadd(result.dvzdz);
        Out.Cinvxx += hadd(result.Cinvxx);
        Out.Cinvxy += hadd(result.Cinvxy);
        Out.Cinvxz += hadd(result.Cinvxz);
        Out.Cinvyx += hadd(result.Cinvyx);
        Out.Cinvyy += hadd(result.Cinvyy);
        Out.Cinvyz += hadd(result.Cinvyz);
        Out.Cinvzx += hadd(result.Cinvzx);
        Out.Cinvzy += hadd(result.Cinvzy);
        Out.Cinvzz += hadd(result.Cinvzz);
    }
};

void pkdSPHForcesEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions ) {
    if (SPHoptions->doShearStrengthModel) {
        ilist::EvalBlock<ResultSPHForces<fvec,true>,BlockPP<ILC_PART_PER_BLK>> eval(
                    Part.r[0],Part.r[1],Part.r[2],Part.fBall,Part.Omega,
                    Part.v[0],Part.v[1],Part.v[2],Part.rho,Part.P,Part.cs,
                    Part.Sxx, Part.Syy, Part.Sxy, Part.Sxz, Part.Syz,
                    SPHoptions);
        auto result = EvalTile(tile,eval);
        combineSPHForcesResult<true>(Out, result, SPHoptions);
    }
    else {
        ilist::EvalBlock<ResultSPHForces<fvec,false>,BlockPP<ILC_PART_PER_BLK>> eval(
                    Part.r[0],Part.r[1],Part.r[2],Part.fBall,Part.Omega,
                    Part.v[0],Part.v[1],Part.v[2],Part.rho,Part.P,Part.cs,
                    Part.Sxx, Part.Syy, Part.Sxy, Part.Sxz, Part.Syz,
                    SPHoptions);
        auto result = EvalTile(tile,eval);
        combineSPHForcesResult<false>(Out, result, SPHoptions);
    }
}

#endif/*USE_SIMD_PP*/
