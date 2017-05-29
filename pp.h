template<class F,class M,bool bGravStep>
void EvalPP(
	F dx, F dy, F dz,F fourh2, F m,F smooth2,
	F &ax, F &ay, F &az, F &pot,
	F iax, F iay, F iaz,F imaga,
	F &ir, F &norm) {
    static const float minSoftening = 1e-18f;
    F d2 = dx*dx + dy*dy + dz*dz;
    M vcmp = d2 < fourh2;
    F td2 = max(minSoftening,max(d2,fourh2));
    F pir = rsqrt(td2);
    F pir2 = pir * pir;

    td2 = pir2 * d2; /* for SOFTENED */
    pir2 = pir2 * pir;

    /* pir and pir2 are valid now for both softened and unsoftened particles */
    /* Now we apply the fix to softened particles only */
    if (!testz(vcmp)) {
	td2 = maskz_mov(vcmp,F(1.0f - td2));
	pir *= 1.0f + td2*(0.5f + td2*(3.0f/8.0f + td2*(45.0f/32.0f)));
	pir2 *= 1.0f + td2*(1.5f + td2*(135.0f/16.0f));
	}
    pir2 *= -m;
    ax = dx * pir2;
    ay = dy * pir2;
    az = dz * pir2;
    pot = -m*pir;
    if (bGravStep) {
	/* Time stepping criteria stuff */
	F padotai = iax*ax + iay*ay + iaz*az;
	vcmp = (padotai>0.0f) & (d2>=smooth2);
	maskz_mov(vcmp,padotai);
	padotai= padotai * imaga;

	ir = pir * padotai;
	norm = pir * td2;
	}
}
