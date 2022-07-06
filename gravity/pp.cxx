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
    const float *a = Part.a;
    float a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    fvec imaga = a2 > 0.0f ? 1.0f / sqrtf(a2) : 0.0f;

    ilist::EvalBlock<ResultPP<fvec>,ilpBlock> eval(
        Part.r[0],Part.r[1],Part.r[2],Part.fSmooth2,a[0],a[1],a[2],imaga);;
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
    const fvec fx,fy,fz,fBall;
    const uint64_t kernelType;

    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy,fvec fz,fvec fBall,uint64_t kernelType)
        : fx(fx),fy(fy),fz(fz),fBall(fBall),kernelType(kernelType) {}

    result_type operator()(int n,BLOCK &blk) {
        // Sentinal values
        while (n&fvec::mask()) {
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e18f;
            blk.m.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i=0; i<n; ++i) {
            result += EvalDensity<fvec,fmask>(fx,fy,fz,blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],blk.m.v[i],fBall,kernelType);
        }
        return result;
    }
};

void pkdDensityEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions ) {
    ilist::EvalBlock<ResultDensity<fvec>,BlockPP<ILC_PART_PER_BLK>> eval(
                Part.r[0],Part.r[1],Part.r[2],Part.fBall,SPHoptions->kernelType);
    auto result = EvalTile(tile,eval);
    Out.rho += hadd(result.arho);
    Out.drhodfball += hadd(result.adrhodfball);
    Out.nden += hadd(result.anden);
    Out.dndendfball += hadd(result.adndendfball);
    Out.nSmooth += hadd(result.anSmooth);
}

template<typename BLOCK> struct ilist::EvalBlock<ResultSPHForces<fvec>,BLOCK> {
    typedef ResultSPHForces<fvec> result_type;
    const fvec fx,fy,fz,fBall,Omega,vx,vy,vz,rho,P,c;
    const i32v species;
    const SPHOptions *const SPHoptions;
    EvalBlock() = default;
    EvalBlock(fvec fx, fvec fy,fvec fz,fvec fBall,fvec Omega,fvec vx,fvec vy,fvec vz,
              fvec rho,fvec P,fvec c,i32v species,SPHOptions *SPHoptions)
        : fx(fx),fy(fy),fz(fz),fBall(fBall),Omega(Omega),vx(vx),vy(vy),vz(vz),
          rho(rho),P(P),c(c),species(species), SPHoptions(SPHoptions) {}

    result_type operator()(int n,BLOCK &blk) {
        // Sentinal values
        while (n&fvec::mask()) {
            // This is a little trick to speed up the calculation. By setting
            // unused entries in the list to have a zero mass, the resulting
            // forces are zero. Setting the distance to a large value avoids
            // softening the non-existent forces which is slightly faster.
            blk.dx.s[n] = blk.dy.s[n] = blk.dz.s[n] = 1e18f;
            blk.fBall.s[n] = blk.Omega.s[n] = blk.rho.s[n] = 1e18f;
            blk.vx.s[n] = blk.vy.s[n] = blk.vz.s[n] = 0.0f;
            blk.m.s[n] = 0.0f;
            blk.P.s[n] = 0.0f;
            blk.c.s[n] = 0.0f;
            blk.uRung.s[n] = 0.0f;
            ++n;
        }
        n /= fvec::width(); // Now number of blocks
        result_type result;
        result.zero();
        for (auto i=0; i<n; ++i) {
            result += EvalSPHForces<fvec,fmask,i32v>(
                          fx,fy,fz,fBall,Omega,vx,vy,vz,rho,P,c,species,
                          blk.dx.v[i],blk.dy.v[i],blk.dz.v[i],blk.m.v[i],blk.fBall.v[i],blk.Omega.v[i],
                          blk.vx.v[i],blk.vy.v[i],blk.vz.v[i],blk.rho.v[i],blk.P.v[i],blk.c.v[i],blk.species.v[i],blk.uRung.v[i],
                          SPHoptions->kernelType,SPHoptions->epsilon,SPHoptions->alpha,SPHoptions->beta,
                          SPHoptions->EtaCourant,SPHoptions->a,SPHoptions->H,SPHoptions->useIsentropic);
        }
        return result;
    }
};

void pkdSPHForcesEval(const PINFOIN &Part, ilpTile &tile,  PINFOOUT &Out, SPHOptions *SPHoptions ) {
    ilist::EvalBlock<ResultSPHForces<fvec>,BlockPP<ILC_PART_PER_BLK>> eval(
                Part.r[0],Part.r[1],Part.r[2],Part.fBall,Part.Omega,
                Part.v[0],Part.v[1],Part.v[2],Part.rho,Part.P,Part.cs,Part.species,
                SPHoptions);
    auto result = EvalTile(tile,eval);
    Out.uDot += hadd(result.uDot);
    Out.a[0] += hadd(result.ax);
    Out.a[1] += hadd(result.ay);
    Out.a[2] += hadd(result.az);
    Out.divv += hadd(result.divv);
    // This should be a horizontal minimum for an fvec, resulting in a float containing the smallest float in the fvec
    for (int k = 0; k < result.dtEst.width(); k++) {
        Out.dtEst = fmin(Out.dtEst,result.dtEst[k]);
    }
    for (int k = 0; k < result.maxRung.width(); k++) {
        Out.maxRung = fmax(Out.maxRung,result.maxRung[k]);
    }
    assert(Out.dtEst > 0);
}

#endif/*USE_SIMD_PP*/
