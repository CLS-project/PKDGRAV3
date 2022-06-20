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

#ifdef __CUDACC__
    #define PC_CUDA_BOTH __host__ __device__
#else
    #define PC_CUDA_BOTH
#endif
template<class F=float>
struct ResultPC {
    F ax, ay, az, pot;
    F ir, norm;
    PC_CUDA_BOTH void zero() {ax=ay=az=pot=ir=norm=0;}
    PC_CUDA_BOTH ResultPC<F> operator+=(const ResultPC<F> rhs) {
        ax += rhs.ax;
        ay += rhs.ay;
        az += rhs.az;
        pot += rhs.pot;
        ir += rhs.ir;
        norm += rhs.norm;
        return *this;
    }
};
template<class F,class M,bool bGravStep>
PC_CUDA_BOTH ResultPC<F> EvalPC(
    F Pdx, F Pdy, F Pdz, F Psmooth2, // Particle
    F Idx, F Idy, F Idz, F Im, F Iu, // Interaction(s)
    F Ixxxx,F Ixxxy,F Ixxxz,F Ixxyz,F Ixxyy,F Iyyyz,F Ixyyz,F Ixyyy,F Iyyyy,
    F Ixxx,F Ixyy,F Ixxy,F Iyyy,F Ixxz,F Iyyz,F Ixyz,F Ixx,F Ixy,F Ixz,F Iyy,F Iyz,
#ifdef USE_DIAPOLE
    F Ix, F Iy, F Iz,
#endif
    F Pax, F Pay, F Paz,F imaga) {
    ResultPC<F> result;
    const F onethird = 1.0f/3.0f;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;
    F dir = rsqrt(d2);
    F x = dx*dir;
    F y = dy*dir;
    F z = dz*dir;
    {
        // Order 0
        F g = -dir * Im;
        result.pot= g;
        result.ax = x*g;
        result.ay = y*g;
        result.az = z*g;
    }
    F u = Iu*dir;
    F g1 = dir*u;
    F g2 = 3.0f*g1*u;
    F g3 = 5.0f*g2*u;
    F g4 = 7.0f*g3*u;

    F xx = 0.5f*x*x;
    F yy = 0.5f*y*y;
    F zz = 0.5f*z*z;

    // Now calculate the interactions up to Hexadecapole order.
    {
        // Order 4
        F xxx = x*(onethird*xx - zz);
        F tx = Ixxxx*xxx;
        F ty = Ixxxy*xxx;
        F tz = Ixxxz*xxx;
        F xxz = z*(xx - onethird*zz);
        tx += Ixxxz*xxz;
        ty += Ixxyz*xxz;
        tz -= Ixxxx*xxz + Ixxyy*xxz;
        F yyz = z*(yy - onethird*zz);
        tx += Ixyyz*yyz;
        ty += Iyyyz*yyz;
        tz -= Ixxyy*yyz + Iyyyy*yyz;
        F yyy = y*(onethird*yy - zz);
        tx += Ixyyy*yyy;
        ty += Iyyyy*yyy;
        tz += Iyyyz*yyy;

        xx -= zz;
        yy -= zz;

        F xyy = x*yy;
        tx += Ixxyy*xyy;
        ty += Ixyyy*xyy;
        tz += Ixyyz*xyy;
        F xxy = y*xx;
        tx += Ixxxy*xxy;
        ty += Ixxyy*xxy;
        tz += Ixxyz*xxy;
        F xyz = x*y*z;
        tx += Ixxyz*xyz;
        ty += Ixyyz*xyz;
        tz -= (Ixyyy + Ixxxy)*xyz;

        result.ax += (tx *= g4);
        result.ay += (ty *= g4);
        result.az += (tz *= g4);
        F g = 0.25f*(tx*x + ty*y + tz*z);
        result.pot -= g;
        result.ax -= x*9.0f*g;
        result.ay -= y*9.0f*g;
        result.az -= z*9.0f*g;
    }
    {
        // Order 3
        F tx = Ixxx*xx + Ixyy*yy;
        F ty = Ixxy*xx + Iyyy*yy;
        F tz = Ixxz*xx + Iyyz*yy;
        F xy = x*y;
        tx += Ixxy*xy;
        ty += Ixyy*xy;
        tz += Ixyz*xy;
        F xz = x*z;
        tx += Ixxz*xz;
        ty += Ixyz*xz;
        tz -= (Ixxx + Ixyy)*xz;
        F yz = y*z;
        tx += Ixyz*yz;
        ty += Iyyz*yz;
        tz -= (Ixxy + Iyyy)*yz;
        result.ax += (tx *= g3); // xxx
        result.ay += (ty *= g3); // xxy
        result.az += (tz *= g3); // xxz

        F g = onethird*(tx*x + ty*y + tz*z);
        result.pot -= g;
        result.ax -= x*7.0f*g;
        result.ay -= y*7.0f*g;
        result.az -= z*7.0f*g;
    }
    {
        // Order 2
        F t,g;
        result.ax += (t = g2*(Ixx*x + Ixy*y + Ixz*z)); // xx
        g = t*x;
        result.ay += (t = g2*(Iyy*y + Ixy*x + Iyz*z)); // xy
        g += t*y;
        result.az += (t = g2*(-(Ixx + Iyy)*z + Ixz*x + Iyz*y)); // xz
        g += t*z;
        g *= 0.5f;
        result.pot -= g;
        result.ax -= x*5.0f*g;
        result.ay -= y*5.0f*g;
        result.az -= z*5.0f*g;
    }

#ifdef USE_DIAPOLE
    {
        F yy = g1*Ix;
        result.ax += yy;
        F yz = g1*Iy;
        result.ay += yz;
        F zz = g1*Iz;
        result.az += zz;
        F g = yy*x + yz*y + zz*z;
        result.pot -= g;
        result.ax -= x*3.0f*g;
        result.ay -= y*3.0f*g;
        result.az -= z*3.0f*g;
    }
#endif
    result.ax *= dir;
    result.ay *= dir;
    result.az *= dir;

    /* Calculations for determining the timestep. */
    if (bGravStep) {
        F adotai = Pax*result.ax + Pay*result.ay + Paz*result.az;
        adotai = maskz_mov((adotai>0.0f) & (d2>Psmooth2),adotai) * imaga;
        result.norm = adotai * adotai;
        result.ir = dir * result.norm;
    }
    return result;
}
#undef PC_CUDA_BOTH
