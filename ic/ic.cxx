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

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "pst.h"
#include "ic.h"
#include "whitenoise.hpp"
#include "core/gridinfo.hpp"
#ifdef BLACKHOLES
    #include "blackhole/evolve.h"
#endif
using namespace gridinfo;
using namespace blitz;
using namespace mdl;
static const std::complex<float> I(0,1);

typedef blitz::Array<basicParticle,3> basicParticleArray;
static basicParticleArray getOutputArray(PKD pkd,GridInfo &G,real_array_t &R) {
    void *pData = mdlSetArray(pkd->mdl,0,0,pkd->particles);
    if (blitz::product(R.shape())) {
        basicParticleArray fullOutput(
            reinterpret_cast<basicParticle *>(pData),
            blitz::shape(G.n1r(),G.n2(),G.nz()), blitz::neverDeleteData,
            RegularArray(G.sz()));
        basicParticleArray output = fullOutput(
                                        blitz::Range::all(),
                                        blitz::Range(R.base(1),R.base(1)+R.extent(1)-1),
                                        blitz::Range(G.sz(),G.ez()-1));
        output.reindexSelf(dimension_t(0,R.base(1),G.sz()));
        return output;
    }
    else return basicParticleArray();
}



// A blitz++ friendly wrap function. returns "ik" given array index
// Range: (-iNyquist,iNyquist] where iNyquist = m/2
static float fwrap(float v,float m) {
    return v - (v > m*0.5 ? m : 0);
}
BZ_DEFINE_BINARY_FUNC(Fn_fwrap,fwrap)
BZ_DECLARE_ARRAY_ET_BINARY(fwrap,     Fn_fwrap)
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(fwrap,     Fn_fwrap, float)

class Power {
    typedef struct {
        class Power *power;
        double r;
    } varianceParameters;
    static double variance_integrand(double ak, void *params);
public:
    virtual ~Power() = default;
    virtual double getAmplitude(double k) = 0;
    double variance(double dRadius,double k0,double k1);
};

double Power::variance_integrand(double ak, void *params) {
    varianceParameters *vprm = reinterpret_cast<varianceParameters *>(params);
    double x, w;
    /* Window function for spherical tophat of given radius (e.g., 8 Mpc/h) */
    x = ak * vprm->r;
    w = 3.0*(sin(x)-x*cos(x))/(x*x*x);
    return vprm->power->getAmplitude(ak)*ak*ak*w*w*4.0*M_PI;
}

double Power::variance(double dRadius,double k0,double k1) {
    varianceParameters vprm;
    gsl_function F;
    double result, error;
    gsl_integration_workspace *W = gsl_integration_workspace_alloc (1000);
    vprm.power = this;
    vprm.r = dRadius; /* 8 Mpc/h for example */
    F.function = &variance_integrand;
    F.params = &vprm;
    gsl_integration_qag(&F, exp(k0), exp(k1),
                        0.0, 1e-6, 1000, GSL_INTEG_GAUSS61, W, &result, &error);
    gsl_integration_workspace_free(W);
    return result;
}

class PowerTransfer : public Power {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    double *tk, *tf;
    double spectral;
    double normalization;
    int nTf;
    double rise0, rise1;
public:
    virtual double getAmplitude(double k);
    PowerTransfer(CSM csm, double a,int nTf, double *tk, double *tf);
    virtual ~PowerTransfer();
    double variance(double dRadius);
};

double PowerTransfer::getAmplitude(double k) {
    double lk = log(k);
    double T;

    if (lk > tk[nTf-1]) /* Extrapolate beyond kmax */
        T =  tf[nTf-1] + (lk - tk[nTf-1]) * rise1;
    else if (lk < tk[0]) /* Extrapolate beyond kmin */
        T =  tf[0] + (lk - tk[0]) * rise0;
    else
        T = gsl_spline_eval(spline,lk,acc);
    T = exp(T);
    return pow(k,spectral) * normalization * T * T;
}

double PowerTransfer::variance(double dRadius) {
    return Power::variance(dRadius,tk[0], tk[nTf-1]);
}

PowerTransfer::~PowerTransfer() {
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
}

PowerTransfer::PowerTransfer(CSM csm, double a,int nTf, double *tk, double *tf) {
    if (csm->val.classData.bClass)
        return;
    double D1_0, D2_0, D1_a, D2_a;
    double f1_0, f2_0, f1_a, f2_a;
    csmComoveGrowth(csm, 1.0, &D1_0, &D2_0, &f1_0, &f2_0);
    csmComoveGrowth(csm, a, &D1_a, &D2_a, &f1_a, &f2_a);
    normalization = 1.0;
    spectral = csm->val.dSpectral;
    this->nTf = nTf;
    this->tk = tk;
    this->tf = tf;
    rise0 = (tf[0] - tf[1]) / (tk[0] - tk[1]);
    rise1 = (tf[nTf-1] - tf[nTf-2]) / (tk[nTf-1] - tk[nTf-2]);
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc (gsl_interp_cspline, nTf);
    gsl_spline_init(spline, tk, tf, nTf);
    float dSigma8 = csm->val.dSigma8;
    if (dSigma8 > 0) {
        dSigma8 *= D1_a/D1_0;
        normalization *= dSigma8*dSigma8 / variance(8.0);
    }
    else if (csm->val.dNormalization > 0) {
        normalization = csm->val.dNormalization * D1_a/D1_0;
        dSigma8 = sqrt(variance(8.0));
    }
}

int pkdGenerateIC(PKD pkd, MDLFFT fft, int iSeed, int bFixed, float fPhase, int nGrid, int b2LPT, double dBoxSize,
                  double a, int nTf, double *tk, double *tf, double *noiseMean, double *noiseCSQ) {
    double pi = 4.*atan(1.);
    CSM csm = pkd->csm;
    int bClass = csm->val.classData.bClass;
    mdlClass *mdl = reinterpret_cast<mdlClass *>(pkd->mdl);
    NoiseGenerator ng(iSeed, bFixed, fPhase);
    /* Prepare growth factors and transfer function */
    double D1_0, D2_0, f1_0, f2_0;
    double D1_a, D2_a, f1_a, f2_a;
    csmComoveGrowth(csm, 1., &D1_0, &D2_0, &f1_0, &f2_0);
    csmComoveGrowth(csm,  a, &D1_a, &D2_a, &f1_a, &f2_a);
    double dOmega = csm->val.dOmega0/(a*a*a*pow(csmExp2Hub(csm, a)/csm->val.dHubble0, 2.0));
    double f1_approx =    pow(dOmega, 5./ 9.);
    double f2_approx = 2.*pow(dOmega, 6./11.);
    PowerTransfer transfer(csm, a, nTf, tk, tf);
    if (bClass) {
        double D1_internal, D2_internal, f1_internal, f2_internal;
        csm->val.classData.bClassGrowth = 0;
        csmComoveGrowth(csm, a, &D1_internal, &D2_internal, &f1_internal, &f2_internal);
        csm->val.classData.bClassGrowth = 1;
        if (mdl->Self() == 0) {
            printf("f1 = %.12g (CLASS), %.12g (internal but CLASS background), %.12g (approx)\n", f1_a, f1_internal, f1_approx);
            if (b2LPT)
                printf("f2 = %.12g (CLASS), %.12g (internal but CLASS background), %.12g (approx)\n", f2_a, f2_internal, f2_approx);
        }
    }
    else {
        if (mdl->Self() == 0) {
            printf("f1 = %.12g (internal), %.12g (approx)\n", f1_a, f1_approx);
            if (b2LPT)
                printf("f2 = %.12g (internal), %.12g (approx)\n", f2_a, f2_approx);
        }
    }
    fflush(stdout);
    /* Set up grids.
    ** Note that "data" points to the same block for all threads.
    ** Particles will overlap K[0] through K[5] eventually.
    */
    GridInfo G(pkd->mdl, fft);
    complex_array_t K[10];
    real_array_t    R[10];
    auto data = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->particles));
    for (auto i=0; i<10; ++i) {
        G.setupArray(data, K[i]);
        G.setupArray(data, R[i]);
        data += fft->rgrid->nLocal;
    }
    /* Create a local view of our part of the output array */
    basicParticleArray output = getOutputArray(pkd, G, R[0]);
    int nLocal = blitz::product(R[7].shape());
    /* Do 2LPT before 1LPT */
    float fft_normalization = dBoxSize/pow(nGrid, 3.);
    float k_fundamental = 2*pi/dBoxSize;
    float transfer_factor = pow(k_fundamental, 1.5)/dBoxSize;
    float vel_factor1 = a*a*csmExp2Hub(csm, a)*f1_a;
    float vel_factor2 = a*a*csmExp2Hub(csm, a)*f2_a;
    float kx, ky, kz, k2, x, amp;
    if (b2LPT) {
        /* Generate primordial white noise for 2LPT */
        if (mdl->Self()==0) {
            printf("Generating primordial noise\n"); fflush(stdout);
        }
        ng.FillNoise(K[6], nGrid, noiseMean, noiseCSQ);
        /* 2LPT particle displacements */
        if (mdl->Self() == 0) {
            printf("Imprinting 2LPT density spectrum\n"); fflush(stdout);
        }
        complex<float> psi_x, psi_y, psi_z;
        for (auto index = K[6].begin(); index != K[6].end(); index++) {
            auto pos = index.position();
            kx = fwrap(pos[0], fft->rgrid->n1)*k_fundamental;
            ky = fwrap(pos[1], fft->rgrid->n2)*k_fundamental;
            kz = fwrap(pos[2], fft->rgrid->n3)*k_fundamental;
            k2 = kx*kx + ky*ky + kz*kz;
            if (k2 == 0.) {
                for (auto i = 0; i < 6; i++)
                    K[i](pos) = 0.;
                continue;
            }
            if (bClass)
                amp = -csmDelta_m(csm, a, sqrt(k2));  // delta < 0 in CLASS convention
            else
                amp = transfer_factor*sqrt(transfer.getAmplitude(sqrt(k2)));
            /* 1LPT grid values */
            psi_x = (*index)*amp*kx/k2*(-I);
            psi_y = (*index)*amp*ky/k2*(-I);
            psi_z = (*index)*amp*kz/k2*(-I);
            /* 2LPT grids */
            K[0](pos) = psi_x*kx*I;  // \psi_x,x
            K[1](pos) = psi_y*ky*I;  // \psi_y,y
            K[2](pos) = psi_z*kz*I;  // \psi_z,z
            K[3](pos) = psi_x*ky*I;  // \psi_x,y = \psi_y,x
            K[4](pos) = psi_y*kz*I;  // \psi_y,z = \psi_z,y
            K[5](pos) = psi_z*kx*I;  // \psi_z,x = \psi_x,z
        }
        char terms[6][3] = {"xx", "yy", "zz", "xy", "yz", "zx"};
        for (auto i = 0; i < 6; i++) {
            if (mdl->Self() == 0) {
                printf("Generating 2LPT %s term\n", terms[i]); fflush(stdout);
            }
            mdl->IFFT(fft, (FFTW3(complex) *)K[i].dataFirst());
        }
        /* Calculate the source term */
        if (mdl->Self()==0) {
            printf("Generating 2LPT source term\n"); fflush(stdout);
        }
        amp = D2_a/(D1_a*D1_a)*fft_normalization;
        for (auto index = R[6].begin(); index != R[6].end(); index++) {
            auto pos = index.position();
            *index = amp*(
                         + R[0](pos)*R[1](pos)  // \psi_x,x * \psi_y,y
                         + R[1](pos)*R[2](pos)  // \psi_y,y * \psi_z,z
                         + R[2](pos)*R[0](pos)  // \psi_z,z * \psi_x,x
                         - R[3](pos)*R[3](pos)  // \psi_x,y * \psi_y,x = (\psi_x,y)^2
                         - R[4](pos)*R[4](pos)  // \psi_y,z * \psi_z,y = (\psi_y,z)^2
                         - R[5](pos)*R[5](pos)  // \psi_z,x * \psi_x,z = (\psi_z,x)^2
                     );
        }
        mdl->FFT(fft, R[6].dataFirst());
        for (auto index = K[6].begin(); index != K[6].end(); index++) {
            auto pos = index.position();
            kx = fwrap(pos[0], fft->rgrid->n1)*k_fundamental;
            ky = fwrap(pos[1], fft->rgrid->n2)*k_fundamental;
            kz = fwrap(pos[2], fft->rgrid->n3)*k_fundamental;
            k2 = kx*kx + ky*ky + kz*kz;
            if (k2 == 0.) {
                for (auto i = 0; i < 3; i++)
                    K[7 + i](pos) = 0.;
                continue;
            }
            K[7](pos) = (*index)*kx/k2*(-I);
            K[8](pos) = (*index)*ky/k2*(-I);
            K[9](pos) = (*index)*kz/k2*(-I);
        }
        for (auto i = 0; i < 3; i++) {
            if (mdl->Self() == 0) {
                printf("Generating 2LPT %c displacements\n", 'x' + i); fflush(stdout);
            }
            mdl->IFFT(fft, (FFTW3(complex) *)K[7 + i].dataFirst());
        }
        /* Set particle positions and velocities to the 2LPT values */
        if (mdl->Self() == 0 ) {
            printf("Displacing particle positions and boosting velocities (2LPT)\n"); fflush(stdout);
        }
        for (auto index = output.begin(); index != output.end(); index++) {
            auto pos = index.position();
            for (auto i = 0; i < 3; i++) {
                x = R[7 + i](pos);
                index->dr[i] = x;
                index-> v[i] = x*vel_factor2;
            }
        }
    }  // done with 2LPT
    /* Generate primordial white noise for 1LPT */
    if (mdl->Self() == 0) {
        printf("Generating primordial noise\n"); fflush(stdout);
    }
    ng.FillNoise(K[6], nGrid, noiseMean, noiseCSQ);
    /* 1LPT particle displacements and velocities */
    for (auto quantity = 0; quantity < 1 + bClass; quantity++) {  // 0: positions, 1: velocities
        if (mdl->Self() == 0) {
            printf("Imprinting %s spectrum\n", quantity == 0 ? "density" : "velocity");
            fflush(stdout);
        }
        for (auto index = K[6].begin(); index != K[6].end(); index++) {
            auto pos = index.position();
            kx = fwrap(pos[0], fft->rgrid->n1)*k_fundamental;
            ky = fwrap(pos[1], fft->rgrid->n2)*k_fundamental;
            kz = fwrap(pos[2], fft->rgrid->n3)*k_fundamental;
            k2 = kx*kx + ky*ky + kz*kz;
            if (k2 == 0.) {
                for (auto i = 0; i < 3; i++)
                    K[7 + i](pos) = 0.;
                continue;
            }
            if (bClass) {
                if (quantity == 0)
                    amp = -csmDelta_m(csm, a, sqrt(k2));  // delta < 0 in CLASS convention
                else
                    amp =  csmTheta_m(csm, a, sqrt(k2));
            }
            else
                /* The same transfer function is used for both
                ** delta and theta when using back-scaling.
                */
                amp = transfer_factor*sqrt(transfer.getAmplitude(sqrt(k2)));
            K[7](pos) = (*index)*amp*kx/k2*(-I);
            K[8](pos) = (*index)*amp*ky/k2*(-I);
            K[9](pos) = (*index)*amp*kz/k2*(-I);
        }
        for (auto i = 0; i < 3; i++) {
            if (mdl->Self() == 0) {
                printf("Generating %c %s\n", 'x' + i, quantity == 0 ? "displacements" : "velocities");
                fflush(stdout);
            }
            mdl->IFFT(fft, (FFTW3(complex) *)K[7 + i].dataFirst());
        }
        /* Displace particle positions or/and boost velocities */
        if (mdl->Self() == 0) {
            if (bClass) {
                if (quantity == 0)
                    printf("Displacing particle positions\n");
                else
                    printf("Boosting particle velocities\n");
            }
            else
                printf("Displacing particle positions and boosting velocities\n");
            fflush(stdout);
        }
        for (auto index = output.begin(); index != output.end(); index++) {
            auto pos = index.position();
            for (auto i = 0; i < 3; i++) {
                if (quantity == 0) {
                    if (b2LPT) index->dr[i] += R[7 + i](pos);
                    else       index->dr[i]  = R[7 + i](pos);
                }
                else if (b2LPT) index->v[i] += R[7 + i](pos);
                else            index->v[i]  = R[7 + i](pos);
                if (!bClass) {
                    if (b2LPT)  index->v[i] += vel_factor1*R[7 + i](pos);
                    else        index->v[i]  = vel_factor1*R[7 + i](pos);
                }
            }
        }
    }
    return nLocal;
}

int pltMoveIC(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inMoveIC *in = reinterpret_cast<struct inMoveIC *>(vin);
    int i;
    mdlClass *mdl = reinterpret_cast<mdlClass *>(pst->mdl);

    mdlassert(mdl,nIn == sizeof(struct inMoveIC));
    assert(pstOnNode(pst)); /* We pass around pointers! */
    if (pstNotCore(pst)) {
        struct inMoveIC icUp;

        icUp.pBase = in->pBase;
        icUp.nMove = pst->nUpper * in->nMove / pst->nLeaves;
        in->nMove -= icUp.nMove;
        icUp.iStart = in->iStart + in->nMove;
        icUp.fMass = in->fMass;
        icUp.fSoft = in->fSoft;
        icUp.nGrid = in->nGrid;
        icUp.bICgas = in->bICgas;
        icUp.dInitialT = in->dInitialT;
        icUp.dInitialH = in->dInitialH;
#ifdef HAVE_HELIUM
        icUp.dInitialHe = in->dInitialHe;
#endif
#ifdef HAVE_CARBON
        icUp.dInitialC = in->dInitialC;
#endif
#ifdef HAVE_NITROGEN
        icUp.dInitialN = in->dInitialN;
#endif
#ifdef HAVE_OXYGEN
        icUp.dInitialO = in->dInitialO;
#endif
#ifdef HAVE_NEON
        icUp.dInitialNe = in->dInitialNe;
#endif
#ifdef HAVE_MAGNESIUM
        icUp.dInitialMg = in->dInitialMg;
#endif
#ifdef HAVE_SILICON
        icUp.dInitialSi = in->dInitialSi;
#endif
#ifdef HAVE_IRON
        icUp.dInitialFe = in->dInitialFe;
#endif
#ifdef HAVE_METALLICITY
        icUp.dInitialMetallicity = in->dInitialMetallicity;
#endif
        icUp.dExpansion = in->dExpansion;
        icUp.dOmegaRate = in->dOmegaRate;
        icUp.dTuFac = in->dTuFac;

        int rID = mdl->ReqService(pst->idUpper,PLT_MOVEIC,&icUp,nIn);
        mdl->GetReply(rID);
        pltMoveIC(pst->pstLower,in,nIn,NULL,0);
    }
    else {
        PKD pkd = plcl->pkd;
        pkd->particles.clearClasses();
        double inGrid = 1.0 / in->nGrid;
        float fGasMass, fDarkMass, fGasSoft, fDarkSoft;
        if (in->bICgas) {
            fGasMass = 2.0*in->fMass*in->dOmegaRate;
            fDarkMass = 2.0*in->fMass*(1.0 - in->dOmegaRate);
            fGasSoft = in->fSoft * sqrt(in->dOmegaRate);
            fDarkSoft = in->fSoft;
        }
        else {
            fDarkMass = in->fMass;
            fDarkSoft = in->fSoft;
        }
        for (i=in->nMove-1; i>=0; --i) {
            auto p = pkd->particles[i];
            auto &Vel = p.velocity();
            // If we have no particle order convert directly to Integerized positions.
            // We do this to save space as an "Integer" particle is small.
            if (pkd->bIntegerPosition && pkd->bNoParticleOrder) {
                integerParticle *b = ((integerParticle *)in->pBase) + in->iStart + i;
                integerParticle temp;
                memcpy(&temp,b,sizeof(temp));
                Vel[2] = temp.v[2];
                Vel[1] = temp.v[1];
                Vel[0] = temp.v[0];
                auto &r = p.raw_position<int32_t>();
                r[2] = temp.r[2];
                r[1] = temp.r[1];
                r[0] = temp.r[0];
            }
            else {
                expandParticle *b = ((expandParticle *)in->pBase) + in->iStart + i;
                expandParticle temp;
                memcpy(&temp,b,sizeof(temp));
                Vel[2] = temp.v[2];
                Vel[1] = temp.v[1];
                Vel[0] = temp.v[0];
                blitz::TinyVector<double,3> r(  temp.dr[0] + (temp.ix+0.5) * inGrid - 0.5,
                                                temp.dr[1] + (temp.iy+0.5) * inGrid - 0.5,
                                                temp.dr[2] + (temp.iz+0.5) * inGrid - 0.5);
                p.set_position(r);
                if (pkd->particles.present(PKD_FIELD::oParticleID)) {
                    auto &ID = p.ParticleID();
                    ID = temp.ix + in->nGrid*(temp.iy + 1ul*in->nGrid*temp.iz);
                }
                if (!pkd->bNoParticleOrder)
                    p.set_order(temp.ix + in->nGrid*(temp.iy + 1ul*in->nGrid*temp.iz));
            }
            pkd->particles.setClass(fDarkMass,fDarkSoft,0,FIO_SPECIES_DARK,&p);
            p.set_marked(true);
            p.set_rung(0);
            if (pkd->bNoParticleOrder) p.set_group(0);
            else p.set_new_rung(0);
            if (pkd->particles.present(PKD_FIELD::oPotential)) p.potential() = 0;
            if (in->bICgas) {
                auto pgas = pkd->particles[i+in->nMove];
                pgas = p;
                pkd->particles.setClass(fGasMass, fGasSoft, 0, FIO_SPECIES_SPH, &pgas);

                if (pkd->particles.present(PKD_FIELD::oParticleID)) {
                    auto &ID = pgas.ParticleID();
                    ID += in->nGrid*in->nGrid*in->nGrid;
                }
                if (!pkd->bNoParticleOrder)
                    pgas.set_order(pgas.order() + in->nGrid*in->nGrid*in->nGrid);
                pgas.set_position(pgas.position() + inGrid*0.5);

                auto &sph = pgas.sph();
                sph.pNeighborList = NULL;

                auto &VelGas = pgas.velocity();
                // Change the scale factor dependency
                double a_m1 = 1./in->dExpansion;
                VelGas *= a_m1;

                /* Fill the SPHFIELDS with some initial values */
                double u = in->dInitialT * in->dTuFac;
                assert(pkd->particles.present(PKD_FIELD::oSph));
                auto &Sph = pgas.sph();
                Sph.afElemMass[ELEMENT_H]  = in->dInitialH  * fGasMass;
#ifdef HAVE_HELIUM
                Sph.afElemMass[ELEMENT_He] = in->dInitialHe * fGasMass;
#endif
#ifdef HAVE_CARBON
                Sph.afElemMass[ELEMENT_C]  = in->dInitialC  * fGasMass;
#endif
#ifdef HAVE_NITROGEN
                Sph.afElemMass[ELEMENT_N]  = in->dInitialN  * fGasMass;
#endif
#ifdef HAVE_OXYGEN
                Sph.afElemMass[ELEMENT_O]  = in->dInitialO  * fGasMass;
#endif
#ifdef HAVE_NEON
                Sph.afElemMass[ELEMENT_Ne] = in->dInitialNe * fGasMass;
#endif
#ifdef HAVE_MAGNESIUM
                Sph.afElemMass[ELEMENT_Mg] = in->dInitialMg * fGasMass;
#endif
#ifdef HAVE_SILICON
                Sph.afElemMass[ELEMENT_Si] = in->dInitialSi * fGasMass;
#endif
#ifdef HAVE_IRON
                Sph.afElemMass[ELEMENT_Fe] = in->dInitialFe * fGasMass;
#endif
#ifdef HAVE_METALLICITY
                Sph.fMetalMass = in->dInitialMetallicity * fGasMass;
#endif
                Sph.vPred = VelGas;
                Sph.Frho = 0.0;
                Sph.Fmom = 0.0;
                Sph.Fene = 0.0;
                Sph.E = u + 0.5*blitz::dot(Sph.vPred,Sph.vPred);
                Sph.E *= fGasMass;
                Sph.Uint = u*fGasMass;
                assert(Sph.E>0);
                Sph.mom = fGasMass*VelGas;
                Sph.lastMom = 0.;
                Sph.lastE = Sph.E;
#ifdef ENTROPY_SWITCH
                Sph.S = 0.0;
                Sph.lastS = 0.0;
                Sph.maxEkin = 0.0;
#endif
                Sph.lastUint = Sph.Uint;
                Sph.lastHubble = 0.0;
                Sph.lastMass = fGasMass;
                Sph.lastAcc = 0.;
#ifndef USE_MFM
                Sph.lastDrDotFrho = 0.;
                Sph.drDotFrho = 0.;
#endif
                //Sph.fLastBall = 0.0;
                Sph.lastUpdateTime = -1.;
                // Sph.nLastNeighs = 100;
#ifdef STAR_FORMATION
                Sph.SFR = 0.;
#endif
#ifdef COOLING
                Sph.lastCooling = 0.;
                Sph.cooling_dudt = 0.;
#endif
#if defined(FEEDBACK) || defined(BLACKHOLES)
                Sph.fAccFBEnergy = 0.;
#endif
                Sph.uWake = 0;
                Sph.omega = 0.0;
#ifdef BLACKHOLES
                Sph.BHAccretor.iIndex = NOT_ACCRETED;
                Sph.BHAccretor.iPid   = NOT_ACCRETED;
#endif
            }
        }
        if (in->bICgas) pkd->SetLocal( pkd->nActive = in->nMove * 2);
        else pkd->SetLocal(pkd->nActive = in->nMove);
    }
    return 0;
}

int pstMoveIC(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    struct inGenerateIC *in = reinterpret_cast<struct inGenerateIC *>(vin);
    mdlClass *mdl = reinterpret_cast<mdlClass *>(pst->mdl);

    if (pstOffNode(pst)) {
        int rID = mdl->ReqService(pst->idUpper,PST_MOVEIC,vin,nIn);
        pstMoveIC(pst->pstLower,vin,nIn,NULL,0);
        mdl->GetReply(rID);
    }
    else {
        MDLFFT fft = pkd->fft;
        int myProc = mdl->Proc();
        int iProc;
        uint64_t iUnderBeg=0, iOverBeg=0;
        uint64_t iUnderEnd, iOverEnd;
        uint64_t iBeg, iEnd;

        int *scount = reinterpret_cast<int *>(malloc(sizeof(int)*mdl->Procs())); assert(scount!=NULL);
        int *sdisps = reinterpret_cast<int *>(malloc(sizeof(int)*mdl->Procs())); assert(sdisps!=NULL);
        int *rcount = reinterpret_cast<int *>(malloc(sizeof(int)*mdl->Procs())); assert(rcount!=NULL);
        int *rdisps = reinterpret_cast<int *>(malloc(sizeof(int)*mdl->Procs())); assert(rdisps!=NULL);

        assert(pstAmNode(pst));
        assert(fft != NULL);

        uint64_t nPerNode = (uint64_t)mdl->Cores() * pkd->FreeStore();
        uint64_t nLocal = (int64_t)fft->rgrid->rn[myProc] * in->nGrid*in->nGrid;

        /* Calculate how many slots are free (under) and how many need to be sent (over) before my rank */
        iUnderBeg = iOverBeg = 0;
        for (iProc=0; iProc<myProc; ++iProc) {
            uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
            if (nOnNode>nPerNode) iOverBeg += nOnNode - nPerNode;
            else iUnderBeg += nPerNode - nOnNode;
        }
        size_t nSize = sizeof(expandParticle);
        if (pkd->bIntegerPosition && pkd->bNoParticleOrder) nSize = sizeof(integerParticle);
        char *pBase = (char *)pkd->particles.Element(0);
        char *pRecv = pBase + nSize*nLocal;
        //char *eBase;
        if (nLocal > nPerNode) {      /* Too much here: send extra particles to other nodes */
            //eBase = pBase + nSize*nPerNode;
            iUnderBeg = 0;
            iOverEnd = iOverBeg + nLocal - nPerNode;
            for (iProc=0; iProc<mdl->Procs(); ++iProc) {
                rcount[iProc] = rdisps[iProc] = 0; // We cannot receive anything
                uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
                if (nOnNode<nPerNode) {
                    iUnderEnd = iUnderBeg + nPerNode - nOnNode;
                    /* The transfer condition */
                    if (iUnderEnd>iOverBeg && iUnderBeg<iOverEnd) {
                        iBeg = iOverBeg>iUnderBeg ? iOverBeg : iUnderBeg;
                        iEnd = iOverEnd<iUnderEnd ? iOverEnd : iUnderEnd;
                        scount[iProc] = (iEnd-iBeg);
                        sdisps[iProc] = (iBeg-iOverBeg);
                        nLocal -= iEnd-iBeg;
                    }
                    else scount[iProc] = sdisps[iProc] = 0;
                    iUnderBeg = iUnderEnd;
                }
                else scount[iProc] = sdisps[iProc] = 0;
            }
            assert(nLocal == nPerNode);
        }
        else if (nLocal < nPerNode) { /* We have room: *maybe* receive particles from other nodes */
            //eBase = pBase + nSize*nLocal;
            iOverBeg = 0;
            iUnderEnd = iUnderBeg + nPerNode - nLocal;
            for (iProc=0; iProc<mdl->Procs(); ++iProc) {
                scount[iProc] = sdisps[iProc] = 0; // We have nothing to send
                uint64_t nOnNode = fft->rgrid->rn[iProc] * in->nGrid*in->nGrid;
                if (nOnNode>nPerNode) {
                    iOverEnd = iOverBeg + nOnNode - nPerNode;
                    if (iOverEnd>iUnderBeg && iOverBeg<iUnderEnd) {
                        iBeg = iOverBeg>iUnderBeg ? iOverBeg : iUnderBeg;
                        iEnd = iOverEnd<iUnderEnd ? iOverEnd : iUnderEnd;
                        rcount[iProc] = (iEnd-iBeg);
                        rdisps[iProc] = (iBeg-iUnderBeg);
                        nLocal += iEnd-iBeg;
                    }
                    else rcount[iProc] = rdisps[iProc] = 0;
                    iOverBeg = iOverEnd;
                }
                else rcount[iProc] = rdisps[iProc] = 0;
            }
            assert(nLocal <= nPerNode);
        }
        else {
            for (iProc=0; iProc<mdl->Procs(); ++iProc) {
                rcount[iProc] = rdisps[iProc] = 0; // We cannot receive anything
                scount[iProc] = sdisps[iProc] = 0; // We have nothing to send
            }
        }
        mdl->Alltoallv(nSize,
                       pBase + nSize*nPerNode, scount, sdisps,
                       pRecv,            rcount, rdisps);
        free(scount);
        free(sdisps);
        free(rcount);
        free(rdisps);
        mdlFFTNodeFinish(pst->mdl,fft);
        pkd->fft = NULL;

        /* We need to relocate the particles */
        struct inMoveIC move;
        move.pBase = (overlayedParticle *)pkd->particles.Element(0);
        move.iStart = 0;
        move.nMove = nLocal;
        move.fMass = in->dBoxMass;
        move.fSoft = 1.0 / (50.0*in->nGrid);
        move.nGrid = in->nGrid;
        move.bICgas = in->bICgas;
        move.dInitialT = in->dInitialT;
        move.dInitialH = in->dInitialH;
#ifdef HAVE_HELIUM
        move.dInitialHe = in->dInitialHe;
#endif
#ifdef HAVE_CARBON
        move.dInitialC = in->dInitialC;
#endif
#ifdef HAVE_NITROGEN
        move.dInitialN = in->dInitialN;
#endif
#ifdef HAVE_OXYGEN
        move.dInitialO = in->dInitialO;
#endif
#ifdef HAVE_NEON
        move.dInitialNe = in->dInitialNe;
#endif
#ifdef HAVE_MAGNESIUM
        move.dInitialMg = in->dInitialMg;
#endif
#ifdef HAVE_SILICON
        move.dInitialSi = in->dInitialSi;
#endif
#ifdef HAVE_IRON
        move.dInitialFe = in->dInitialFe;
#endif
#ifdef HAVE_METALLICITY
        move.dInitialMetallicity = in->dInitialMetallicity;
#endif
        move.dExpansion = in->dExpansion;
        move.dOmegaRate = in->dOmegaRate;
        move.dTuFac = in->dTuFac;
        pltMoveIC(pst,&move,sizeof(move),NULL,0);
    }
    return 0;
}

/* NOTE: only called when on-node -- pointers are passed around. */
int pltGenerateIC(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inGenerateICthread *tin = reinterpret_cast<struct inGenerateICthread *>(vin);
    struct inGenerateIC *in = tin->ic;
    struct outGenerateIC *out = reinterpret_cast<struct outGenerateIC *>(vout), outUp;
    mdlClass *mdl = reinterpret_cast<mdlClass *>(pst->mdl);
    mdlassert(mdl,nIn == sizeof(struct inGenerateICthread));
    mdlassert(mdl,vout != NULL);
    assert(pstOnNode(pst)); /* We pass around pointers! */

    if (pstNotCore(pst)) {
        int rID = mdl->ReqService(pst->idUpper,PLT_GENERATEIC,vin,nIn);
        pltGenerateIC(pst->pstLower,vin,nIn,vout,nOut);
        mdl->GetReply(rID,&outUp);
        out->N += outUp.N;
        out->noiseMean += outUp.noiseMean;
        out->noiseCSQ += outUp.noiseCSQ;
    }
    else {
        out->N = pkdGenerateIC(plcl->pkd, tin->fft, in->iSeed, in->bFixed, in->fPhase,
                               in->nGrid, in->b2LPT, in->dBoxSize, in->dExpansion, in->nTf,
                               in->k, in->tf, &out->noiseMean, &out->noiseCSQ);
        out->dExpansion = in->dExpansion;
    }

    return sizeof(struct outGenerateIC);
}

int pstGenerateIC(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    PKD pkd = plcl->pkd;
    mdlClass *mdl = reinterpret_cast<mdlClass *>(pst->mdl);
    struct inGenerateIC *in = reinterpret_cast<struct inGenerateIC *>(vin);
    struct outGenerateIC *out = reinterpret_cast<struct outGenerateIC *>(vout), outUp;
    int64_t i;

    mdlassert(pst->mdl,nIn == sizeof(struct inGenerateIC));
    mdlassert(pst->mdl,vout != NULL);

    if (pstAmNode(pst)) {
        struct inGenerateICthread tin;
        MDLFFT fft = mdlFFTNodeInitialize(pst->mdl,in->nGrid,in->nGrid,in->nGrid,0,0);
        tin.ic = reinterpret_cast<struct inGenerateIC *>(vin);
        tin.fft = fft;
        pltGenerateIC(pst,&tin,sizeof(tin),vout,nOut);

        int myProc = mdlProc(pst->mdl);
        uint64_t nLocal = (int64_t)fft->rgrid->rn[myProc] * in->nGrid*in->nGrid;

        /* Expand the particles by adding an iOrder */
        assert(sizeof(expandParticle) >= sizeof(basicParticle));
        overlayedParticle   *pbBase = (overlayedParticle *)pkd->particles.Element(0);
        int iz = fft->rgrid->rs[myProc] + fft->rgrid->rn[myProc];
        int iy=0, ix=0;
        float inGrid = 1.0 / in->nGrid;
        for (i=nLocal-1; i>=0; --i) {
            basicParticle  *b = &pbBase->b + i;
            basicParticle temp;
            memcpy(&temp,b,sizeof(temp));
            if (ix>0) --ix;
            else {
                ix = in->nGrid-1;
                if (iy>0) --iy;
                else {
                    iy = in->nGrid-1;
                    --iz;
                    assert(iz>=0);
                }
            }
            // If we have no particle order convert directly to Integerized positions.
            // We do this to save space as an "Integer" particle is small.
            if (pkd->bIntegerPosition && pkd->bNoParticleOrder) {
                integerParticle *p = &pbBase->i + i;
                p->v[2] = temp.v[2];
                p->v[1] = temp.v[1];
                p->v[0] = temp.v[0];
                p->r[2] = pkdDblToIntPos(pkd,temp.dr[2] + (iz+0.5) * inGrid - 0.5);
                p->r[1] = pkdDblToIntPos(pkd,temp.dr[1] + (iy+0.5) * inGrid - 0.5);
                p->r[0] = pkdDblToIntPos(pkd,temp.dr[0] + (ix+0.5) * inGrid - 0.5);
            }
            else {
                expandParticle *p = &pbBase->e + i;
                p->v[2] = temp.v[2];
                p->v[1] = temp.v[1];
                p->v[0] = temp.v[0];
                p->dr[2] = temp.dr[2];
                p->dr[1] = temp.dr[1];
                p->dr[0] = temp.dr[0];
                p->ix = ix;
                p->iy = iy;
                p->iz = iz;
            }
        }
        assert(ix==0 && iy==0 && iz==fft->rgrid->rs[myProc]);
        /* Now we need to move excess particles between nodes so nStore is obeyed. */
        pkd->fft = fft; /* This is freed in pstMoveIC() */
    }
    else if (pstNotCore(pst)) {
        int rID = mdl->ReqService(pst->idUpper,PST_GENERATEIC,in,nIn);
        pstGenerateIC(pst->pstLower,in,nIn,vout,nOut);
        mdl->GetReply(rID,&outUp);
        out->N += outUp.N;
        out->noiseMean += outUp.noiseMean;
        out->noiseCSQ += outUp.noiseCSQ;
    }
    return sizeof(struct outGenerateIC);
}
