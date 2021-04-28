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
	F &arho, F &adrhodfball         // results
    ) {
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;

    F p, w, dwdp;
    F ifBall;
    F t1, t2, t3;

    F one = 1.0f;
    F two = 2.0f;
    F three = 3.0f;
    F six = 6.0f;
    F eight = 8.0f;

    ifBall = one / fBall;
    p = sqrt(d2) * ifBall;

    M pltone = p < one;

    if (!testz(pltone)) {
        // There is some work to do
        M pltonehalf = p < (one / two);
        t1 = p - one;

        // Evaluate the kernel
        t2 = one + six * p * p * t1;
        t3 = - two * t1 * t1 * t1;
        w = maskz_mov(pltone,t3);
        w = mask_mov(w,pltonehalf,t2);

        // Evaluate the kernel derivative
        t2 = six * p * (three * p - two);
        t3 = - six * t1 * t1;
        dwdp = maskz_mov(pltone,t3);
        dwdp = mask_mov(dwdp,pltonehalf,t2);

        // Normalization factor
        F C = eight * M_1_PI * ifBall * ifBall * ifBall;

        // return the density
        arho = Im * C * w;

        // return the density derivative
        adrhodfball = - Im * C * ifBall * (three * w + dwdp * p);
    } else {
        // No work to do
        arho = 0.0f;
        adrhodfball = 0.0f;
    }
    }
