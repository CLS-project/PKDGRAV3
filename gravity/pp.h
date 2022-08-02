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
    #define PP_CUDA_BOTH __host__ __device__
#else
    #define PP_CUDA_BOTH
#endif
template<class F=float>
struct ResultPP {
    F ax, ay, az, pot;
    F ir, norm;
    PP_CUDA_BOTH void zero() { ax=ay=az=pot=ir=norm=0; }
    PP_CUDA_BOTH ResultPP<F> operator+=(const ResultPP<F> rhs) {
        ax += rhs.ax;
        ay += rhs.ay;
        az += rhs.az;
        pot += rhs.pot;
        ir += rhs.ir;
        norm += rhs.norm;
        return *this;
    }
};
template<class F,class M>
PP_CUDA_BOTH ResultPP<F> EvalPP(
    F Pdx, F Pdy, F Pdz, F Psmooth2,       // Particle
    F Idx, F Idy, F Idz, F fourh2, F Im) { // Interaction(s)
    ResultPP<F> result;
    constexpr float minSoftening = 1e-18f;
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
    result.ax = dx * dir2;
    result.ay = dy * dir2;
    result.az = dz * dir2;
    result.pot = -Im*dir;
    result.ir = dir;
    result.norm = d2;
    return result;
}

// Calculate additional terms for GravStep
template<class F,class M>
PP_CUDA_BOTH ResultPP<F> EvalPP(
    F Pdx, F Pdy, F Pdz, F Psmooth2,     // Particle
    F Idx, F Idy, F Idz, F fourh2, F Im, // Interaction(s)
    F Pax, F Pay, F Paz, F imaga) {
    ResultPP<F> result = EvalPP<F,M>(Pdx,Pdy,Pdz,Psmooth2,Idx,Idy,Idz,fourh2,Im);
    F adotai = Pax*result.ax + Pay*result.ay + Paz*result.az;
    adotai = maskz_mov((adotai>0.0f) & (result.norm>=Psmooth2),adotai) * imaga;
    result.norm = adotai * adotai;
    result.ir *= result.norm;
    return result;
}

template<class F=float>
struct ResultDensity {
    F arho, adrhodfball, anden, adndendfball, anSmooth;
    void zero() { arho=adrhodfball=anden=adndendfball=anSmooth=0; }
    ResultDensity<F> operator+=(const ResultDensity<F> rhs) {
        arho += rhs.arho;
        adrhodfball += rhs.adrhodfball;
        anden += rhs.anden;
        adndendfball += rhs.adndendfball;
        anSmooth += rhs.anSmooth;
        return *this;
    }
};
template<class F,class M>
PP_CUDA_BOTH ResultDensity<F> EvalDensity(
    F Pdx, F Pdy, F Pdz,     // Particle
    F Idx, F Idy, F Idz, F Im, F fBall,  // Interaction(s)
    int kernelType) {
    ResultDensity<F> result;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;

    F r, w, dwdr, dWdfball;
    F ifBall, C;
    F t1, t2, t3;
    M mask1;

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
        result.anden = C * w;
        result.arho = Im * result.anden;

        // return the density derivative
        result.adndendfball = dWdfball;
        result.adrhodfball = Im * result.adndendfball;

        // return the number of particles used
        result.anSmooth = maskz_mov(r_lt_one,1.0f);
    }
    else {
        result.zero(); // No work to do
    }
    return result;
}

template<class F=float>
struct ResultSPHForces {
    F uDot, ax, ay, az, divv, dtEst, maxRung;
    void zero() { uDot=ax=ay=az=divv=maxRung=0; dtEst=HUGE_VALF; }
    ResultSPHForces<F> operator+=(const ResultSPHForces<F> rhs) {
        uDot += rhs.uDot;
        ax += rhs.ax;
        ay += rhs.ay;
        az += rhs.az;
        divv += rhs.divv;
        dtEst = min(dtEst,rhs.dtEst); // CAREFUL! We use "min" here
        maxRung = max(maxRung,rhs.maxRung);
        return *this;
    }
};
template<class F,class M,class Ivec>
PP_CUDA_BOTH ResultSPHForces<F> EvalSPHForces(
    F Pdx, F Pdy, F Pdz, F PfBall, F POmega,     // Particle
    F Pvx, F Pvy, F Pvz, F Prho, F PP, F Pc, Ivec Pspecies,
    F Idx, F Idy, F Idz, F Im, F IfBall, F IOmega,      // Interactions
    F Ivx, F Ivy, F Ivz, F Irho, F IP, F Ic, Ivec Ispecies, F uRung,
    int kernelType, float epsilon, float alpha, float beta,
    float EtaCourant,float a,float H,bool useIsentropic) {
    ResultSPHForces<F> result;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;
    F dvx = Pvx - Ivx;
    F dvy = Pvy - Ivy;
    F dvz = Pvz - Ivz;

    F d, Pr, Ir;
    F t1, t2, t3;
    F PifBall, IifBall, PC, IC, Pdwdr, Idwdr, PdWdr, IdWdr;
    F PdWdx, PdWdy, PdWdz, IdWdx, IdWdy, IdWdz, dWdx, dWdy, dWdz;
    F cij, rhoij, fBallij, dvdotdx, muij, Piij;
    F PPoverRho2, IPoverRho2, dtC, dtMu;
    F vFac, aFac;
    M Pr_lt_one, Ir_lt_one, mask1, dvdotdx_st_zero;


    PifBall = 1.0f / PfBall;
    IifBall = 1.0f / IfBall;
    d = sqrt(d2);
    Pr = d * PifBall;
    Ir = d * IifBall;

    Pr_lt_one = Pr < 1.0f;
    Ir_lt_one = Ir < 1.0f;

    if (!testz(Pr_lt_one) || !testz(Ir_lt_one)) {
        // There is some work to do

        // First convert velocities
        aFac = a;
        vFac = 1.0f / (a * a);
        dvx = vFac * dvx;
        dvy = vFac * dvy;
        dvz = vFac * dvz;

        // Kernel derivatives
        SPHKERNEL_INIT(Pr, PifBall, PC, t1, mask1, kernelType);
        DSPHKERNEL_DR(Pr, Pdwdr, t1, t2, t3, Pr_lt_one, mask1, kernelType);
        PdWdr = PC * Pdwdr;
        SPHKERNEL_INIT(Ir, IifBall, IC, t1, mask1, kernelType);
        DSPHKERNEL_DR(Ir, Idwdr, t1, t2, t3, Ir_lt_one, mask1, kernelType);
        IdWdr = IC * Idwdr;

        // Kernel gradients, separate at the moment, as i am not sure if we need them separately
        // at some point. If we don't, can save some operations, by combining earlier.
        t1 = PdWdr * PifBall / d;
        mask1 = Pr > 0.0f;
        t1 = maskz_mov(mask1,t1);
        PdWdx = t1 * dx;
        PdWdy = t1 * dy;
        PdWdz = t1 * dz;
        t1 = IdWdr * IifBall / d;
        mask1 = Ir > 0.0f;
        t1 = maskz_mov(mask1,t1);
        IdWdx = t1 * dx;
        IdWdy = t1 * dy;
        IdWdz = t1 * dz;
        dWdx = 0.5f * (PdWdx + IdWdx);
        dWdy = 0.5f * (PdWdy + IdWdy);
        dWdz = 0.5f * (PdWdz + IdWdz);

        // PoverRho2
        PPoverRho2 = PP / (POmega * Prho * Prho);
        IPoverRho2 = IP / (IOmega * Irho * Irho);

        // Artificial viscosity
        cij = 0.5f * (Pc + Ic);
        rhoij = 0.5f * (Prho + Irho);
        fBallij = 0.25f * (PfBall + IfBall); // here we need h, not fBall
        dvdotdx = dvx * dx + dvy * dy + dvz * dz + H * d2;
        dvdotdx_st_zero = dvdotdx < 0.0f;
        muij = fBallij * dvdotdx * aFac / (d2 + epsilon * fBallij * fBallij);
        muij = maskz_mov(dvdotdx_st_zero,muij);
        Piij = (-alpha * cij * muij + beta * muij * muij) / rhoij;

        // du/dt
        if (useIsentropic) {
            result.uDot = 0.5f * Piij * Im * (dvx * dWdx + dvy * dWdy + dvz * dWdz);
        }
        else {
            result.uDot = (PPoverRho2 + 0.5f * Piij) * Im * (dvx * dWdx + dvy * dWdy + dvz * dWdz);
        }

        // acceleration
        result.ax = - Im * (PPoverRho2 * PdWdx + IPoverRho2 * IdWdx + Piij * dWdx) * aFac;
        result.ay = - Im * (PPoverRho2 * PdWdy + IPoverRho2 * IdWdy + Piij * dWdy) * aFac;
        result.az = - Im * (PPoverRho2 * PdWdz + IPoverRho2 * IdWdz + Piij * dWdz) * aFac;

        // divv
        result.divv = Im / Irho * (dvx * dWdx + dvy * dWdy + dvz * dWdz);

        // timestep
        dtC = (1.0f + 0.6f * alpha) / (aFac * EtaCourant);
        dtMu = 0.6f * beta / (aFac * EtaCourant);
        result.dtEst = 0.5f * PfBall / (dtC * Pc - dtMu * muij);
        mask1 = Pr_lt_one | Ir_lt_one;
        result.dtEst = mask_mov(HUGE_VALF,mask1,result.dtEst);
        result.maxRung = maskz_mov(mask1,uRung);

        // for (int index = 0; index<8;index++) {
        // if (ax[index] != ax[index]) {
        // printf("ax = %.15e ax = %.15e ax = %.15e\n",ax[index],ay[index],az[index]);
        // printf("Pdx = %.15e Pdy = %.15e Pdz = %.15e\n",Pdx[index],Pdy[index],Pdz[index]);
        // printf("PfBall = %.15e POmega = %.15e Prho = %.15e\n",PfBall[index],POmega[index],Prho[index]);
        // printf("Pvx = %.15e Pvy = %.15e Pvy = %.15e\n",Pvx[index],Pvy[index],Pvz[index]);
        // printf("PP = %.15e Pc = %.15e Pspecies = %d\n",PP[index],Pc[index],Pspecies[index]);
        // printf("Idx = %.15e Idy = %.15e Idz = %.15e\n",Idx[index],Idy[index],Idz[index]);
        // printf("Im = %.15e IfBall = %.15e IOmega = %.15e\n Irho = %.15e\n",Im[index],IfBall[index],IOmega[index],Irho[index]);
        // printf("Ivx = %.15e Ivy = %.15e Ivz = %.15e\n",Ivx[index],Ivy[index],Ivz[index]);
        // printf("IP = %.15e Ic = %.15e Ispecies = %d\n",IP[index],Ic[index],Ispecies[index]);
        // printf("Im = %.15e PPoverRho2 = %.15e IPoverRho2 = %.15e Piij = %.15e dWdz = %.15e aFac = %.15e\n",Im[index],PPoverRho2[index], IPoverRho2[index], Piij[index], dWdz[index],aFac[index]);
        // assert(0);
        // }
        // }
    }
    else {
        result.zero();
    }
    return result;
}
#undef PP_CUDA_BOTH
