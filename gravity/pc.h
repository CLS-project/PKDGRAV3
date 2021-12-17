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
CUDA_DEVICE void EvalPC(
    const F &Pdx, const F &Pdy, const F &Pdz, const F &Psmooth2, // Particle
    const F &Idx, const F &Idy, const F &Idz, const F &Im, const F &Iu, // Interaction(s)
    const F &Ixxxx,const F &Ixxxy,const F &Ixxxz,const F &Ixxyz,const F &Ixxyy,const F &Iyyyz,const F &Ixyyz,const F &Ixyyy,const F &Iyyyy,
    const F &Ixxx,const F &Ixyy,const F &Ixxy,const F &Iyyy,const F &Ixxz,const F &Iyyz,const F &Ixyz,
    const F &Ixx,const F &Ixy,const F &Ixz,const F &Iyy,const F &Iyz,
#ifdef USE_DIAPOLE
    const F &Ix, const F &Iy, const F &Iz,
#endif
    F &ax, F &ay, F &az, F &pot,     // Results
    const F &Pax, const F &Pay, const F &Paz,const F &imaga,F &ir, F &norm) {
    const F onethird = 1.0f/3.0f;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;
    F dir = rsqrt(d2);
    F u = Iu*dir;
    F g1 = dir*u;
    F g2 = 3.0f*g1*u;
    F g3 = 5.0f*g2*u;
    F g4 = 7.0f*g3*u;
    /*
    ** Calculate the funky distance terms.
    */
    F x = dx*dir;
    F y = dy*dir;
    F z = dz*dir;
    F xx = 0.5f*x*x;
    F xy = x*y;
    F xz = x*z;
    F yy = 0.5f*y*y;
    F yz = y*z;
    F zz = 0.5f*z*z;
    F xxx = x*(onethird*xx - zz);
    F xxz = z*(xx - onethird*zz);
    F yyy = y*(onethird*yy - zz);
    F yyz = z*(yy - onethird*zz);
    xx -= zz;
    yy -= zz;
    F xxy = y*xx;
    F xyy = x*yy;
    F xyz = xy*z;

    /*
    ** Now calculate the interaction up to Hexadecapole order.
    */
    F tx = g4*(Ixxxx*xxx + Ixyyy*yyy + Ixxxy*xxy + Ixxxz*xxz + Ixxyy*xyy + Ixxyz*xyz + Ixyyz*yyz);
    F ty = g4*(Ixyyy*xyy + Ixxxy*xxx + Iyyyy*yyy + Iyyyz*yyz + Ixxyy*xxy + Ixxyz*xxz + Ixyyz*xyz);
    F tz = g4*(-Ixxxx*xxz - (Ixyyy + Ixxxy)*xyz - Iyyyy*yyz + Ixxxz*xxx + Iyyyz*yyy - Ixxyy*(xxz + yyz) + Ixxyz*xxy + Ixyyz*xyy);
    g4 = 0.25f*(tx*x + ty*y + tz*z);
    xxx = g3*(Ixxx*xx + Ixyy*yy + Ixxy*xy + Ixxz*xz + Ixyz*yz);
    xxy = g3*(Ixyy*xy + Ixxy*xx + Iyyy*yy + Iyyz*yz + Ixyz*xz);
    xxz = g3*(-(Ixxx + Ixyy)*xz - (Ixxy + Iyyy)*yz + Ixxz*xx + Iyyz*yy + Ixyz*xy);
    g3 = onethird*(xxx*x + xxy*y + xxz*z);
    xx = g2*(Ixx*x + Ixy*y + Ixz*z);
    xy = g2*(Iyy*y + Ixy*x + Iyz*z);
    xz = g2*(-(Ixx + Iyy)*z + Ixz*x + Iyz*y);
    g2 = 0.5f*(xx*x + xy*y + xz*z);
    F g0 = dir * Im;
    pot = -(g0 + g2 + g3 + g4);
    g0 += 5.0f*g2 + 7.0f*g3 + 9.0f*g4;
#ifdef USE_DIAPOLE
    yy = g1*Ix;
    yz = g1*Iy;
    zz = g1*Iz;
    g1 = yy*x + yz*y + zz*z;
    pot -= g1;
    g0 += 3.0f*g1;
#else
    yy = 0.0f;
    yz = 0.0f;
    zz = 0.0f;
#endif
    ax = dir*(yy + xx + xxx + tx - x*g0);
    ay = dir*(yz + xy + xxy + ty - y*g0);
    az = dir*(zz + xz + xxz + tz - z*g0);

    /* Calculations for determining the timestep. */
    if (bGravStep) {
        F adotai = Pax*ax + Pay*ay + Paz*az;
        adotai = maskz_mov(adotai>0.0f & d2>Psmooth2,adotai) * imaga;
        norm = adotai * adotai;
        ir = dir * norm;
    }
}
