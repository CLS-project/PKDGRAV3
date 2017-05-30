#ifndef CUDA_DEVICE
#define CUDA_DEVICE
#endif
template<class F,class M,bool bGravStep>
CUDA_DEVICE void EvalPC(
	F Pdx, F Pdy, F Pdz, F Psmooth2, // Particle
	F Idx, F Idy, F Idz, F Im, F Iu, // Interaction(s)
	F Ixxxx,F Ixxxy,F Ixxxz,F Ixxyz,F Ixxyy,F Iyyyz,F Ixyyz,F Ixyyy,F Iyyyy,
	F Ixxx,F Ixyy,F Ixxy,F Iyyy,F Ixxz,F Iyyz,F Ixyz,
	F Ixx,F Ixy,F Ixz,F Iyy,F Iyz,
	F &ax, F &ay, F &az, F &pot,     // Results
	F Pax, F Pay, F Paz,F imaga,F &ir, F &norm) {
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
