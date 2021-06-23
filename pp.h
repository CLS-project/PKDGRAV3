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

#ifndef CUDA_DEVICE
#define CUDA_DEVICE
#endif
template<class F,class M,bool bGravStep>
CUDA_DEVICE void EvalPP(
	const F &Pdx, const F &Pdy, const F &Pdz, const F &Psmooth2,     // Particle
	const F &Idx, const F &Idy, const F &Idz, const F &fourh2, const F &Im, // Interaction(s)
	F &ax, F &ay, F &az, F &pot,         // results
	const F &Pax, const F &Pay, const F &Paz,const F &imaga, F &ir, F &norm) {
    static const float minSoftening = 1e-18f;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;
    M vcmp = d2 < fourh2;
    F td2 = max(minSoftening,max(d2,fourh2));
    F dir = rsqrt(td2);
    F dir2 = dir * dir;

    td2 = dir2 * d2; /* for SOFTENED */
    dir2 = dir2 * dir;

    /* dir and dir2 are valid now for both softened and unsoftened particles */
    /* Now we apply the fix to softened particles only */
    if (!testz(vcmp)) {
	td2 = maskz_mov(vcmp,F(1.0f - td2));
	dir *= 1.0f + td2*(0.5f + td2*(3.0f/8.0f + td2*(45.0f/32.0f)));
	dir2 *= 1.0f + td2*(1.5f + td2*(135.0f/16.0f));
	}
    dir2 *= -Im;
    ax = dx * dir2;
    ay = dy * dir2;
    az = dz * dir2;
    pot = -Im*dir;

    /* Calculations for determining the timestep. */
    if (bGravStep) {
	F adotai = Pax*ax + Pay*ay + Paz*az;
	adotai = maskz_mov(adotai>0.0f & d2>=Psmooth2,adotai) * imaga;
	norm = adotai * adotai;
	ir = dir * norm;
	}
    }

template<class F,class M,bool bGravStep>
CUDA_DEVICE void EvalDensity(
	const F &Pdx, const F &Pdy, const F &Pdz,     // Particle
	const F &Idx, const F &Idy, const F &Idz, const F &Im, const F & fBall, // Interaction(s)
	F &arho, F &adrhodfball, F &anden, F& adndendfball, F &anSmooth,         // results
    SPHOptions SPHoptions) {
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;

    F r, w, dwdr, dWdfball;
    F ifBall, C;
    F t1, t2, t3;
    M mask1;

    int kernelType = SPHoptions.kernelType;

    ifBall = 1.0f / fBall;
    r = sqrt(d2) * ifBall;

    M r_lt_one = r < 1.0f;

    if (!testz(r_lt_one)) {
        // There is some work to do
        SPHKERNEL_INIT(r, ifBall, C, t1, mask1, kernelType);
        SPHKERNEL(r, w, t1, t2, t3, r_lt_one, mask1, kernelType);
        DSPHKERNEL_DR(r, dwdr, t1, t2, t3, r_lt_one, mask1, kernelType);
        DSPHKERNEL_DFBALL(r, ifBall, w, dwdr, C, dWdfball, kernelType);

        // return the density
        anden = C * w;
        arho = Im * anden;

        // return the density derivative
        adndendfball = dWdfball;
        adrhodfball = Im * adndendfball;

        // return the number of particles used
        anSmooth = maskz_mov(r_lt_one,1.0f);
    } else {
        // No work to do
        arho = 0.0f;
        adrhodfball = 0.0f;
        anden = 0.0f;
        adndendfball = 0.0f;
        anSmooth = 0.0f;
    }
    }

template<class F,class M,class Ivec,bool bGravStep>
CUDA_DEVICE void EvalSPHForces(
	const F &Pdx, const F &Pdy, const F &Pdz, const F &PfBall, const F &POmega,     // Particle
    const F &Pvx, const F &Pvy, const F &Pvz, const F &Prho, const F &PP, const F &Pc, const Ivec & Pspecies,
	const F &Idx, const F &Idy, const F &Idz, const F & Im, const F &IfBall, const F &IOmega,     // Interactions
    const F &Ivx, const F &Ivy, const F &Ivz, const F &Irho, const F &IP, const F &Ic, const Ivec & Ispecies,
    SPHOptions SPHoptions) {
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;

    // for (int i=0; i<8;i++) {
        // printf("element: %d, Pdx = %f, Pdy = %f, Pdz = %f, PfBall = %f, POmega = %f, Pvx = %f, Pvy = %f, Pvz = %f, Prho = %f, PP = %f, Pc = %f, Pspecies = %d, Idx = %f, Idy = %f, Idz = %f, Im = %f, IfBall = %f, IOmega = %f, Ivx = %f, Ivy = %f, Ivz = %f, Irho = %f, IP = %f, Ic = %f, Ispecies = %d\n",
        // i,Pdx[i],Pdy[i],Pdz[i],PfBall[i],POmega[i],Pvx[i],Pvy[i],Pvz[i],Prho[i],PP[i],Pc[i],Pspecies[i],Idx[i],Idy[i],Idz[i],Im[i],IfBall[i],IOmega[i],Ivx[i],Ivy[i],Ivz[i],Irho[i],IP[i],Ic[i],Ispecies[i]);
    // }

    F r;
    F ifBall;

    ifBall = 1.0f / PfBall;
    r = sqrt(d2) * ifBall;

    M r_lt_one = r < 1.0f;

    if (!testz(r_lt_one)) {
        // There is some work to do
    } else {
        // No work to do
    }
    }