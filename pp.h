#ifndef CUDA_DEVICE
#define CUDA_DEVICE
#endif
template<class F,class M,bool bGravStep>
CUDA_DEVICE void EvalPP(
	F Pdx, F Pdy, F Pdz, F Psmooth2,     // Particle
	F Idx, F Idy, F Idz, F fourh2, F Im, // Interaction(s)
	F &ax, F &ay, F &az, F &pot,         // results
	F Pax, F Pay, F Paz,F imaga, F &ir, F &norm) {
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
