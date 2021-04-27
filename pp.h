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
	F &arho, F &adrhodh         // results
    ) {
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;

    F ar2, ak, adk;
    F ifBall2;

    F onefourth = 1.0f/4.0f;
    F threefourths = 3.0f/4.0f;
    F one = 1.0f;
    F two = 2.0f;
    F four = 4.0f;
    F three = 3.0f;
    F ninefourths = 9.0f/4.0f;

    F t1 = 0.0f;
    F t2 = 0.0f;
    F t3 = 0.0f;

    ifBall2 = four / (fBall * fBall);
    ar2 = d2 * ifBall2;

    M ar2stone = ar2 < one;
    M ar2stfour = ar2 < four;

    if (!testz(ar2stfour)){
    // There is some work to do

    // Evaluate the kernel

    //ak = 2.0 - sqrt(ar2)
    t1 = two - sqrt(ar2);

    //if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2)
    t2 = (one - threefourths * t1 * ar2);

    //else if (ar2 < 4.0) ak = 0.25*ak*ak*ak
    t3 = onefourth * t1 * t1 * t1;

    // else ak = 0;

    ak = maskz_mov(ar2stfour,t3);
    ak = mask_mov(ak,ar2stone,t2);

    // Evaluate the kernel derivative

    // adk = sqrt(ar2)
    t1 = sqrt(ar2);

    // if (ar2 < 1.0) adk = -3 + 2.25*adk;
    t2 = - three + ninefourths * t1;

    // else if (ar2 < 4.0) adk = -0.75*(2.0-adk)*(2.0-adk)/adk;
    t3 = - threefourths * (two - t1) * (two - t1) / t1;

    // else adk = 0;

    adk = maskz_mov(ar2stfour,t3);
    adk = mask_mov(adk,ar2stone,t2);

    // return the density scaled with the volume
    arho = M_1_PI*sqrt(ifBall2)*ifBall2 * Im * ak;

    // return the derivative of the density wrt fBall
    // the four is actually an 8 but it also should be adk/2
    adrhodh = M_1_PI*sqrt(ifBall2)*ifBall2 * Im * adk * (- four) * d2 / (fBall * fBall * fBall);
    } else {
    // No work to do
    arho = 0.0f;
    adrhodh = 0.0f;
    }
    }
