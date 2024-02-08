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
    typedef F type;
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
    F rho, drhodfball, nden, dndendfball, nSmooth, imbalanceX, imbalanceY, imbalanceZ;
    PP_CUDA_BOTH void zero() { rho=drhodfball=nden=dndendfball=nSmooth=imbalanceX=imbalanceY=imbalanceZ=0; }
    PP_CUDA_BOTH ResultDensity<F> operator+=(const ResultDensity<F> rhs) {
        rho += rhs.rho;
        drhodfball += rhs.drhodfball;
        nden += rhs.nden;
        dndendfball += rhs.dndendfball;
        nSmooth += rhs.nSmooth;
        imbalanceX += rhs.imbalanceX;
        imbalanceY += rhs.imbalanceY;
        imbalanceZ += rhs.imbalanceZ;
        return *this;
    }
};
template<class F,class M>
PP_CUDA_BOTH ResultDensity<F> EvalDensity(
    F Pdx, F Pdy, F Pdz, F fBall, F PiMat,     // Particle
    F Idx, F Idy, F Idz, F Im, F IiMat,  // Interaction(s)
    int kernelType, bool doInterfaceCorrection) {
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
        result.nden = C * w;
        result.rho = Im * result.nden;

        // return the density derivative
        result.dndendfball = dWdfball;
        result.drhodfball = Im * result.dndendfball;

        // return the number of particles used
        result.nSmooth = maskz_mov(r_lt_one,F(1.0f));

        // Calculate the imbalance values
        if (doInterfaceCorrection) {
            M mask2 = PiMat == IiMat;
            F plus_one = 1.0f;
            F minus_one = -1.0f;
            F kappa = mask_mov(minus_one,mask2,plus_one);
            result.imbalanceX = kappa * dx * result.rho;
            result.imbalanceY = kappa * dy * result.rho;
            result.imbalanceZ = kappa * dz * result.rho;
        }
        else {
            result.imbalanceX = 0.0f;
            result.imbalanceY = 0.0f;
            result.imbalanceZ = 0.0f;
        }
    }
    else {
        result.zero(); // No work to do
    }
    return result;
}

template<class F=float>
struct ResultDensityCorrection {
    F corrT, corrP, corr;
    PP_CUDA_BOTH void zero() { corrT=corrP=corr=0; }
    PP_CUDA_BOTH ResultDensityCorrection<F> operator+=(const ResultDensityCorrection<F> rhs) {
        corrT += rhs.corrT;
        corrP += rhs.corrP;
        corr += rhs.corr;
        return *this;
    }
};
template<class F,class M>
PP_CUDA_BOTH ResultDensityCorrection<F> EvalDensityCorrection(
    F Pdx, F Pdy, F Pdz, F fBall,     // Particle
    F Idx, F Idy, F Idz, F T, F P, F expImb2,  // Interaction(s)
    int kernelType) {
    ResultDensityCorrection<F> result;
    F dx = Idx + Pdx;
    F dy = Idy + Pdy;
    F dz = Idz + Pdz;
    F d2 = dx*dx + dy*dy + dz*dz;

    F r, w;
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

        result.corr = expImb2 * C * w;
        result.corrT = T * result.corr;
        result.corrP = P * result.corr;
    }
    else {
        result.zero(); // No work to do
    }
    return result;
}

template<class F,bool doShearStrengthModel>
struct ResultSPHForces {
    F uDot, ax, ay, az, divv, dtEst, maxRung;
    F dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
    F Cinvxx, Cinvxy, Cinvxz, Cinvyx, Cinvyy, Cinvyz, Cinvzx, Cinvzy, Cinvzz;
    PP_CUDA_BOTH void zero() {
        uDot=ax=ay=az=divv=maxRung=0.0f;
        dtEst=1e14f;
        if (doShearStrengthModel) {
            dvxdx=dvxdy=dvxdz=dvydx=dvydy=dvydz=dvzdx=dvzdy=dvzdz=0.0f;
            Cinvxx=Cinvxy=Cinvxz=Cinvyx=Cinvyy=Cinvyz=Cinvzx=Cinvzy=Cinvzz=0.0f;
        }
    }
    PP_CUDA_BOTH ResultSPHForces<F,doShearStrengthModel> operator+=(const ResultSPHForces<F,doShearStrengthModel> rhs) {
        uDot += rhs.uDot;
        ax += rhs.ax;
        ay += rhs.ay;
        az += rhs.az;
        divv += rhs.divv;
        dtEst = min(dtEst,rhs.dtEst); // CAREFUL! We use "min" here
        maxRung = max(maxRung,rhs.maxRung);
        if (doShearStrengthModel) {
            dvxdx += rhs.dvxdx;
            dvxdy += rhs.dvxdy;
            dvxdz += rhs.dvxdz;
            dvydx += rhs.dvydx;
            dvydy += rhs.dvydy;
            dvydz += rhs.dvydz;
            dvzdx += rhs.dvzdx;
            dvzdy += rhs.dvzdy;
            dvzdz += rhs.dvzdz;
            Cinvxx += rhs.Cinvxx;
            Cinvxy += rhs.Cinvxy;
            Cinvxz += rhs.Cinvxz;
            Cinvyx += rhs.Cinvyx;
            Cinvyy += rhs.Cinvyy;
            Cinvyz += rhs.Cinvyz;
            Cinvzx += rhs.Cinvzx;
            Cinvzy += rhs.Cinvzy;
            Cinvzz += rhs.Cinvzz;
        }
        return *this;
    }
};
template<class F,class M, bool doShearStrengthModel>
PP_CUDA_BOTH ResultSPHForces<F,doShearStrengthModel> EvalSPHForces(
    F Pdx, F Pdy, F Pdz, F PfBall, F POmega,     // Particle
    F Pvx, F Pvy, F Pvz, F Prho, F PP, F Pc,
    F PSxx, F PSyy, F PSxy, F PSxz, F PSyz,
    F Idx, F Idy, F Idz, F Im, F IfBall, F IOmega,      // Interactions
    F Ivx, F Ivy, F Ivz, F Irho, F IP, F Ic, F uRung,
    F ISxx, F ISyy, F ISxy, F ISxz, F ISyz,
    int kernelType, float epsilon, float alpha, float beta,
    float EtaCourant,float a,float H,bool useIsentropic) {
    ResultSPHForces<F,doShearStrengthModel> result;
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
    F cij, rhoij, hij, dvdotdx, muij, Piij;
    F POneOverRho2, IOneOverRho2, minusImOverRho;

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
        t1 = PdWdr * PifBall / (d * POmega);
        mask1 = Pr > 0.0f;
        t1 = maskz_mov(mask1,t1);
        PdWdx = t1 * dx;
        PdWdy = t1 * dy;
        PdWdz = t1 * dz;
        t1 = IdWdr * IifBall / (d * IOmega);
        mask1 = Ir > 0.0f;
        t1 = maskz_mov(mask1,t1);
        IdWdx = t1 * dx;
        IdWdy = t1 * dy;
        IdWdz = t1 * dz;
        // The symmetrized kernel does for now not use the Omega correction, so we multiply it out again
        dWdx = 0.5f * (PdWdx * POmega + IdWdx * IOmega);
        dWdy = 0.5f * (PdWdy * POmega + IdWdy * IOmega);
        dWdz = 0.5f * (PdWdz * POmega + IdWdz * IOmega);

        // OneOverRho2
        POneOverRho2 = 1.0f / (Prho * Prho);
        IOneOverRho2 = 1.0f / (Irho * Irho);

        // Artificial viscosity
        cij = 0.5f * (Pc + Ic);
        rhoij = 0.5f * (Prho + Irho);
        hij = 0.25f * (PfBall + IfBall); // here we need h, not fBall
        dvdotdx = dvx * dx + dvy * dy + dvz * dz + H * d2;
        dvdotdx_st_zero = dvdotdx < 0.0f;
        muij = hij * dvdotdx * aFac / (d2 + epsilon * hij * hij);
        muij = maskz_mov(dvdotdx_st_zero,muij);
        Piij = (-alpha * cij * muij + beta * muij * muij) / rhoij;

        // du/dt
        if (useIsentropic) {
            result.uDot = 0.5f * Piij * Im * (dvx * dWdx + dvy * dWdy + dvz * dWdz);
        }
        else {
            result.uDot = Im * (PP * POneOverRho2 * (dvx * PdWdx + dvy * PdWdy + dvz * PdWdz) + 0.5f * Piij * (dvx * dWdx + dvy * dWdy + dvz * dWdz));
        }

        // acceleration
        result.ax = - Im * (PP * POneOverRho2 * PdWdx + IP * IOneOverRho2 * IdWdx + Piij * dWdx) * aFac;
        result.ay = - Im * (PP * POneOverRho2 * PdWdy + IP * IOneOverRho2 * IdWdy + Piij * dWdy) * aFac;
        result.az = - Im * (PP * POneOverRho2 * PdWdz + IP * IOneOverRho2 * IdWdz + Piij * dWdz) * aFac;

        // No transformation back into cosmology units, as strength makes no sense in cosmology. This saves multiplications.
        if (doShearStrengthModel) {
            result.ax += Im * (POneOverRho2 * (PSxx * PdWdx + PSxy * PdWdy + PSxz * PdWdz) + IOneOverRho2 * (ISxx * IdWdx + ISxy * IdWdy + ISxz * IdWdz));
            result.ay += Im * (POneOverRho2 * (PSxy * PdWdx + PSyy * PdWdy + PSyz * PdWdz) + IOneOverRho2 * (ISxy * IdWdx + ISyy * IdWdy + ISyz * IdWdz));
            result.az += Im * (POneOverRho2 * (PSxz * PdWdx + PSyz * PdWdy - (PSxx + PSyy) * PdWdz) + IOneOverRho2 * (ISxz * IdWdx + ISyz * IdWdy - (ISxx + ISyy) * IdWdz));

            result.uDot -= Im * POneOverRho2 * (dvx * (PSxx * PdWdx + PSxy * PdWdy + PSxz * PdWdz) + dvy * (PSxy * PdWdx + PSyy * PdWdy + PSyz * PdWdz) + dvz * (PSxz * PdWdx + PSyz * PdWdy - (PSxx + PSyy) * PdWdz));

            minusImOverRho = - Im / Irho;

            result.dvxdx = minusImOverRho * dvx * PdWdx;
            result.dvxdy = minusImOverRho * dvx * PdWdy;
            result.dvxdz = minusImOverRho * dvx * PdWdz;
            result.dvydx = minusImOverRho * dvy * PdWdx;
            result.dvydy = minusImOverRho * dvy * PdWdy;
            result.dvydz = minusImOverRho * dvy * PdWdz;
            result.dvzdx = minusImOverRho * dvz * PdWdx;
            result.dvzdy = minusImOverRho * dvz * PdWdy;
            result.dvzdz = minusImOverRho * dvz * PdWdz;

            result.Cinvxx = minusImOverRho * dx * PdWdx;
            result.Cinvxy = minusImOverRho * dx * PdWdy;
            result.Cinvxz = minusImOverRho * dx * PdWdz;
            result.Cinvyx = minusImOverRho * dy * PdWdx;
            result.Cinvyy = minusImOverRho * dy * PdWdy;
            result.Cinvyz = minusImOverRho * dy * PdWdz;
            result.Cinvzx = minusImOverRho * dz * PdWdx;
            result.Cinvzy = minusImOverRho * dz * PdWdy;
            result.Cinvzz = minusImOverRho * dz * PdWdz;
        }

        // divv
        result.divv = Im / Irho * (dvx * dWdx + dvy * dWdy + dvz * dWdz);

        // timestep
        result.dtEst = aFac * EtaCourant * 0.5f * PfBall / ((1.0f + 0.6f * alpha) * Pc - 0.6f * beta * muij);
        mask1 = Pr_lt_one | Ir_lt_one;
        result.dtEst = mask_mov(1e14f,mask1,result.dtEst);
        result.maxRung = maskz_mov(mask1,uRung);

        /*for (int index = 0; index < result.dtEst.width(); index++) {
            if (!isfinite(result.ax[index]) || !isfinite(result.ay[index]) || !isfinite(result.az[index]) || !isfinite(result.uDot[index]) || (!isfinite(result.dtEst[index]) && result.dtEst[index] != 1e14f) || result.dtEst[index] <= 0.0f) {
                printf("An error occured in EvalSPHForces:\n\
                ax = %.5e, ay = %.5e, az = %.5e, uDot = %.5e, dtEst = %.5e\n\
                Pdx = %.5e, Pdy = %.5e, Pdz = %.5e, Pvx = %.5e, Pvy = %.5e, Pvz = %.5e\n\
                PfBall = %.5e, POmega = %.5e, Prho = %.5e, PP = %.5e, Pc = %.5e\n\
                Idx = %.5e, Idy = %.5e, Idz = %.5e, Ivx = %.5e, Ivy = %.5e, Ivz = %.5e\n\
                IfBall = %.5e, IOmega = %.5e, Irho = %.5e, IP = %.5e, Ic = %.5e, Im = %.5e\n\
                PPoverRho2 = %.5e, IPoverRho2 = %.5e, Piij = %.5e, muij = %.5e\n\
                PdWdx = %.5e, PdWdy = %.5e, PdWdz = %.5e\n\
                IdWdx = %.5e, IdWdy = %.5e, IdWdz = %.5e\n\
                dWdx = %.5e, dWdy = %.5e, dWdz = %.5e\n\
                ",result.ax[index],result.ay[index],result.az[index],result.uDot[index],result.dtEst[index],
                       Pdx[index],Pdy[index],Pdz[index],Pvx[index],Pvy[index],Pvz[index],
                       PfBall[index],POmega[index],Prho[index],PP[index],Pc[index],
                       Idx[index],Idy[index],Idz[index],Ivx[index],Ivy[index],Ivz[index],
                       IfBall[index],IOmega[index],Irho[index],IP[index],Ic[index],Im[index],
                       PPoverRho2[index],IPoverRho2[index],Piij[index],muij[index],
                       PdWdx[index],PdWdy[index],PdWdz[index],
                       IdWdx[index],IdWdy[index],IdWdz[index],
                       dWdx[index],dWdy[index],dWdz[index]);
            }
        }*/
    }
    else {
        result.zero();
    }
    return result;
}
#undef PP_CUDA_BOTH
