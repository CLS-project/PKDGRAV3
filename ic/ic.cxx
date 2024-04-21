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
    if (acc != NULL)
        gsl_interp_accel_free(acc);
    if (spline != NULL)
        gsl_spline_free(spline);
}

PowerTransfer::PowerTransfer(CSM csm, double a, int nTf, double *tk, double *tf) {
    acc = NULL;
    spline = NULL;
    if (csm->val.classData.bClass)
        return;
    double D1_0, D2_0, D3a_0, D3b_0, D3c_0, f1_0, f2_0, f3a_0, f3b_0, f3c_0;
    double D1_a, D2_a, D3a_a, D3b_a, D3c_a, f1_a, f2_a, f3a_a, f3b_a, f3c_a;
    csmComoveGrowth(csm, 1.0,
                    &D1_0, &D2_0, &D3a_0, &D3b_0, &D3c_0,
                    &f1_0, &f2_0, &f3a_0, &f3b_0, &f3c_0);
    csmComoveGrowth(csm, a,
                    &D1_a, &D2_a, &D3a_a, &D3b_a, &D3c_a,
                    &f1_a, &f2_a, &f3a_a, &f3b_a, &f3c_a);
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







// !!!
typedef struct {
    int term1;
    int term2;
    int term3a;
    int term3b;
    int term3c;
} enabledTerms;
enabledTerms checkEnabledTerms(PKD pkd) {
    int enable_term1 = 1;
    int enable_term2 = 1;
    int enable_term3a = 1;
    int enable_term3b = 1;
    int enable_term3c = 1;
    char *term_str;
    term_str = getenv("term1");
    if (term_str != NULL && strcmp(term_str, "False") == 0)
        enable_term1 = 0;
    term_str = getenv("term2");
    if (term_str != NULL && strcmp(term_str, "False") == 0)
        enable_term2 = 0;
    term_str = getenv("term3a");
    if (term_str != NULL && strcmp(term_str, "False") == 0)
        enable_term3a = 0;
    term_str = getenv("term3b");
    if (term_str != NULL && strcmp(term_str, "False") == 0)
        enable_term3b = 0;
    term_str = getenv("term3c");
    if (term_str != NULL && strcmp(term_str, "False") == 0)
        enable_term3c = 0;
    if (pkd->Self() == 0) {
        printf("\nEnabled terms: ");
        if (enable_term1) printf("1 ");
        if (enable_term2) printf("2 ");
        if (enable_term3a) printf("3a ");
        if (enable_term3b) printf("3b ");
        if (enable_term3c) printf("3c ");
        printf("\n\n");
    }
    enabledTerms enabled;
    enabled.term1 = enable_term1;
    enabled.term2 = enable_term2;
    enabled.term3a = enable_term3a;
    enabled.term3b = enable_term3b;
    enabled.term3c = enable_term3c;
    return enabled;
}





/* Function for out-of-place differentiation (once of twice)
** of a Fourier grid, followed by an inverse Fourier transform.
*/
void diffIFFT(PKD pkd, MDLFFT fft, complex_array_t &gridIn, complex_array_t &gridOut,
              float kFundamental, int i, int j = -1) {
    uint32_t nvec[] = {fft->rgrid->n1, fft->rgrid->n2, fft->rgrid->n3};
    if (j == -1) {
        /* Differentiate along dimension i */
        for (auto index = gridOut.begin(); index != gridOut.end(); index++) {
            auto pos = index.position();
            float ki = fwrap(pos[i], nvec[i])*kFundamental;
            *index = gridIn(pos)*ki*I;  // diff. along dimension i: ki*sqrt(-1)
        }
    }
    else {
        /* Differentiate along dimensions i and j */
        for (auto index = gridOut.begin(); index != gridOut.end(); index++) {
            auto pos = index.position();
            float ki = fwrap(pos[i], nvec[i])*kFundamental;
            float kj = fwrap(pos[j], nvec[j])*kFundamental;
            *index = -gridIn(pos)*ki*kj;  // diff. along dimension i and j: -ki*kj
        }
    }
    /* Transform to real space */
    pkd->mdl->IFFT(fft, (FFTW3(complex) *)gridOut.dataFirst());
}

/* Function which applies the Fourier transform of a grid
** in real space, followed by an in-place inverse Laplacian
** (possily with an additional factor).
*/
void FFTLaplacianInverse(PKD pkd, MDLFFT fft, real_array_t &grid_R, complex_array_t &grid_K,
                         float kFundamental, float factor = 1) {
    /* Transform to Fourier space */
    pkd->mdl->FFT(fft, grid_R.dataFirst());
    /* Apply inverse Laplacian */
    for (auto index = grid_K.begin(); index != grid_K.end(); index++) {
        auto pos = index.position();
        auto kx = fwrap(pos[0], fft->rgrid->n1);
        auto ky = fwrap(pos[1], fft->rgrid->n2);
        auto kz = fwrap(pos[2], fft->rgrid->n3);
        auto k2 = kx*kx + ky*ky + kz*kz;
        if (k2 == 0) {
            *index = 0.;
            continue;
        }
        float k2Phys = k2*kFundamental*kFundamental;
        *index *= factor/(-k2Phys);  // inverse Laplacian: 1/(-k^2)
    }
}

/* Helper functionality for building up LPT potentials */
typedef struct {
    double D1;
    double f1;
    double D2;
    double f2;
    double D3a;
    double f3a;
    double D3b;
    double f3b;
    double D3c;
    double f3c;
} growthFactors;
typedef struct {
    int gridIndex = -1;
    int diffs[2] = {-1, -1};
} differentiatedPotential;
enum class TMP_STATE {
    INACTIVE,  // contents of temporary grid not relevant
    ACTIVE,    // contents of temporary grid may be reused
};
typedef struct {
    int gridIndex = -1;
    TMP_STATE state = TMP_STATE::INACTIVE;
    differentiatedPotential diffPot;
} tmpNote;
enum class LPT_TERM_OP {
    EQ,      //  =
    PLS_EQ,  // +=
    MIN_EQ,  // -=
};
tmpNote *getTmpNote(std::queue<tmpNote *> &notes) {
    auto note = notes.front();
    /* Rotate */
    notes.pop();
    notes.push(note);
    return note;
}
int compareDiffPots(differentiatedPotential diffPot0, differentiatedPotential diffPot1) {
    return diffPot0.gridIndex == diffPot1.gridIndex
           && diffPot0.diffs[0] == diffPot1.diffs[0]
           && diffPot0.diffs[1] == diffPot1.diffs[1];
}
int reusableTmp(tmpNote note, differentiatedPotential diffPot) {
    return note.state == TMP_STATE::ACTIVE && compareDiffPots(note.diffPot, diffPot);
}
void handleLPTTerm(
    PKD pkd, MDLFFT fft, real_array_t *R, complex_array_t *K, std::queue<tmpNote *> &notes, float kFundamental,
    int gridIndexOut, LPT_TERM_OP op, float factor, std::vector<differentiatedPotential> diffPots
) {
    assert(diffPots.size() >= 2);
    /* Sort differentiations to come (maximizing reusability) */
    for (auto &diffPot : diffPots) {
        if (diffPot.diffs[0] > diffPot.diffs[1]) {
            auto diff = diffPot.diffs[0];
            diffPot.diffs[0] = diffPot.diffs[1];
            diffPot.diffs[1] = diff;
        }
    }
    /* The notes keep track of the current state of the
    ** temporary grids. If the first needed diff. potential is present
    ** amongst the temporary grids, move that to the front
    ** so that it will be reused as the base grid.
    */
    tmpNote *note;
    auto diffPot = diffPots[0];
    int reuseFirst = 0;
    for (auto _ = 0; _ < notes.size(); _++) {
        note = notes.front();
        if (reusableTmp(*note, diffPot)) {
            /* First needed diff. potential present amongst
            ** temporary grids. It has been moved to the front.
            */
            if (pkd->Self() == 0) printf("!!! Move reusable tmp to front\n");
            reuseFirst = 1;
            break;
        }
        note = getTmpNote(notes);  // just for the rotation
    }
    if (!reuseFirst) {
        /* The first needed diff. potential was not found amongst the
        ** temporary grids. Instead search for the second (unique)
        ** diff. potential. If found, place in back of the queue,
        ** so that it is not used as the base grid.
        */
        for (auto n = 1; n < diffPots.size(); n++) {
            if (diffPot.gridIndex != diffPots[n].gridIndex
                    || diffPot.diffs[0] != diffPots[n].diffs[0]
                    || diffPot.diffs[1] != diffPots[n].diffs[1]
               ) {
                /* Found unique second diff. potential */
                diffPot = diffPots[n];
                for (auto _ = 0; _ < notes.size(); _++) {
                    note = getTmpNote(notes);
                    if (reusableTmp(*note, diffPot)) {
                        /* Second needed diff. potential present amongst
                        ** temporary grids. It has been moved to the back.
                        */
                        if (pkd->Self() == 0) printf("!!! Move reusable tmp to back\n");
                        break;
                    }
                }
                break;
            }
        }
    }
    /* Obtain the first potential, differentiated as needed.
    ** This will be stored in the base grid.
    */
    diffPot = diffPots[0];
    auto noteBase = getTmpNote(notes);
    auto gridBase_R = R[noteBase->gridIndex];
    auto gridBase_K = K[noteBase->gridIndex];
    if (!reusableTmp(*noteBase, diffPot)) {
        if (pkd->Self() == 0) printf("!!! Constructing base grid\n");
        diffIFFT(pkd, fft, K[diffPot.gridIndex], gridBase_K, kFundamental,
                 diffPot.diffs[0], diffPot.diffs[1]);
    }
    else {
        /* Reuse base grid as is */
        if (pkd->Self() == 0) printf("!!! Reusing base grid\n");
    }
    noteBase->state = TMP_STATE::INACTIVE;  // flag as being non-reusable
    /* Obtain grid for temporary storage */
    auto noteTmp = getTmpNote(notes);
    auto gridTmp_R = R[noteTmp->gridIndex];
    auto gridTmp_K = K[noteTmp->gridIndex];
    /* Multiply remaining grids into base grid */
    for (auto n = 1; n < diffPots.size(); n++) {
        diffPot = diffPots[n];
        /* Obtain the next grid */
        if (n == 1 && compareDiffPots(diffPot, diffPots[0])) {
            /* Squared sub-expression at the beginning.
            ** Copy base grid into temporary grid. We could avoid this
            ** copying, but doing it allows for later reusage.
            */
            if (pkd->Self() == 0) printf("!!! Square: Copying base into tmp\n");
            for (auto index = gridBase_K.begin(); index != gridBase_K.end(); index++) {
                auto pos = index.position();
                gridTmp_K(pos) = *index;
            }
        }
        else if (!reusableTmp(*noteTmp, diffPot)) {
            /* Construct grid */
            if (pkd->Self() == 0) printf("!!! Constructing tmp grid\n");
            diffIFFT(pkd, fft, K[diffPot.gridIndex], gridTmp_K, kFundamental,
                     diffPot.diffs[0], diffPot.diffs[1]);
        }
        else {
            /* Reuse temporary grid as is */
            if (pkd->Self() == 0) printf("!!! Reusing tmp\n");
        }
        /* Note down the contents of the temporary grid */
        noteTmp->diffPot = diffPot;
        noteTmp->state = TMP_STATE::ACTIVE;
        /* Multiply into base grid */
        for (auto index = gridBase_R.begin(); index != gridBase_R.end(); index++) {
            auto pos = index.position();
            *index *= gridTmp_R(pos);
        }
    }
    /* Move constructed term from base grid
    ** to output grid and apply factor.
    */
    auto gridOut_R = R[gridIndexOut];
    if (factor == 1.) {
        switch (op) {
        case LPT_TERM_OP::EQ:
            for (auto index = gridOut_R.begin(); index != gridOut_R.end(); index++) {
                auto pos = index.position();
                *index = gridBase_R(pos);
            }
            break;
        case LPT_TERM_OP::PLS_EQ:
            for (auto index = gridOut_R.begin(); index != gridOut_R.end(); index++) {
                auto pos = index.position();
                *index += gridBase_R(pos);
            }
            break;
        case LPT_TERM_OP::MIN_EQ:
            for (auto index = gridOut_R.begin(); index != gridOut_R.end(); index++) {
                auto pos = index.position();
                *index -= gridBase_R(pos);
            }
            break;
        }
    }
    else {
        switch (op) {
        case LPT_TERM_OP::EQ:
            for (auto index = gridOut_R.begin(); index != gridOut_R.end(); index++) {
                auto pos = index.position();
                *index = factor*gridBase_R(pos);
            }
            break;
        case LPT_TERM_OP::PLS_EQ:
            for (auto index = gridOut_R.begin(); index != gridOut_R.end(); index++) {
                auto pos = index.position();
                *index += factor*gridBase_R(pos);
            }
            break;
        case LPT_TERM_OP::MIN_EQ:
            for (auto index = gridOut_R.begin(); index != gridOut_R.end(); index++) {
                auto pos = index.position();
                *index -= factor*gridBase_R(pos);
            }
            break;
        }
    }
}

/* Function for carrying out 1LPT (Zeldovich). Besides setting particle
** positions and velocities, the \Psi1 potential (in Fourier space)
** will be available after this function returns.
*/
void carryout1LPT(
    PKD pkd, MDLFFT fft, basicParticleArray output, real_array_t *R, complex_array_t *K,
    int gridIndexPhi1, std::queue<tmpNote *> &notes,
    growthFactors growth, int iSeed, int bFixed, float fPhase, int nGrid, double dBoxSize,
    double a, int nTf, double *tk, double *tf, double *noiseMean, double *noiseCSQ,
    int printIndent = 0
) {

    // !!!
    auto enabled = checkEnabledTerms(pkd);

    float kFundamental = 2.*M_PI/dBoxSize;
    CSM csm = pkd->csm;
    int bClass = csm->val.classData.bClass;
    int gridIndexTmp0 = getTmpNote(notes)->gridIndex;
    int gridIndexTmp1 = getTmpNote(notes)->gridIndex;
    int onlyOneTmpGrid = (gridIndexTmp0 == gridIndexTmp1);
    auto Phi1_K = K[gridIndexPhi1];
    auto tmp0_R = R[gridIndexTmp0];
    auto tmp0_K = K[gridIndexTmp0];
    auto tmp1_K = K[gridIndexTmp1];
    float transferFactor = pow(kFundamental, 1.5)/dBoxSize;
    float velocityFactor = a*a*csmExp2Hub(csm, a)*growth.f1;
    /* Prepare transfer function (only used when not using CLASS) */
    PowerTransfer transfer(csm, a, nTf, tk, tf);
    /* Set particle velocities and positions (in that order).
    ** By doing the positions last, the \Phi1 grid ends up
    ** containing the actual \Phi1 values.
    */
    int noiseGenerated = 0;
    for (int variable = bClass; variable >= 0; variable--) {  // first \theta, then \delta
        /* Generate primordial white noise */
        auto &noise_K = tmp1_K;  // primordial noise
        if (!noiseGenerated || onlyOneTmpGrid) {
            if (pkd->Self() == 0) printf("%*sGenerating primordial noise\n", printIndent, "");
            NoiseGenerator ng(iSeed, bFixed, fPhase);
            ng.FillNoise(noise_K, nGrid, noiseMean, noiseCSQ);
            noiseGenerated = 1;
        }
        /* With variable == 1 (\theta):
        **   Create potential d\Phi1/dtau in Fourier space, defined by
        **     \nabla^2 d\Phi1/dtau = \theta
        **   The boost field is then
        **     u1 = \nabla d\Phi1/dtau
        ** With variable == 0 (\delta):
        **   Create potential \Phi1 in Fourier space, defined by
        **     \nabla^2 \Phi1 = -\delta
        **   The displacement field is then
        **     \Psi1 = \nabla \Phi1
        **   For back-scaling (!bClass):
        **     Here we furhter use
        **       \theta ~= - a*H*f*\delta
        **       => d\Phi1/dtau ~= a H f \Phi1
        **       => u1 ~= a H f \nabla \Phi1
        **     The velocities used in PKDGRAV has an extra factor a,
        **     so really we need to use (a^2 H f) as the prefactor.
        **     Note that this extra a is already in the \theta stored
        **     as CLASS data.
        */
        if (pkd->Self() == 0) {
            if (!bClass) {
                printf("%*sBuilding potential\n", printIndent, "");
            }
            else {
                printf("%*sBuilding potential (%s)\n", printIndent, "",
                       variable ? "velocities" : "positions");
            }
        }
        for (auto index = noise_K.begin(); index != noise_K.end(); index++) {
            auto pos = index.position();
            auto kx = fwrap(pos[0], fft->rgrid->n1);
            auto ky = fwrap(pos[1], fft->rgrid->n2);
            auto kz = fwrap(pos[2], fft->rgrid->n3);
            auto k2 = kx*kx + ky*ky + kz*kz;
            if (k2 == 0) {
                Phi1_K(pos) = 0.;
                continue;
            }
            float k2Phys = k2*kFundamental*kFundamental;
            float amp;
            if (!bClass) {
                amp = transferFactor*sqrt(transfer.getAmplitude(sqrt(k2Phys)));
            }
            else if (variable == 0) {
                amp = -csmDelta_m(csm, a, sqrt(k2Phys));
            }
            else {    // variable == 1
                amp = csmTheta_m(csm, a, sqrt(k2Phys));
            }
            Phi1_K(pos) = (*index)*amp/(-k2Phys);  // inverse Laplacian: 1/(-k^2)
        }
        /* Construct displacement field from potential,
        ** displace positions and boost velocities.
        */
        if (pkd->Self() == 0) {
            if (!bClass) {
                printf("%*sDisplacing positions and boosting velocities\n", printIndent, "");
            }
            else if (variable == 0) {
                printf("%*sDisplacing positions\n", printIndent, "");
            }
            else {    // variable == 1
                printf("%*sBoosting velocities\n", printIndent, "");
            }
        }
        auto &Psi1_R = tmp0_R;  // displacement field
        auto &Psi1_K = tmp0_K;

        // !!!
        if (!enabled.term1) continue;

        for (auto i = 0; i < 3; i++) {
            diffIFFT(pkd, fft, Phi1_K, Psi1_K, kFundamental, i);
            for (auto index = output.begin(); index != output.end(); index++) {
                auto pos = index.position();
                if (!bClass) {
                    index->dr[i] = Psi1_R(pos);
                    index-> v[i] = velocityFactor*Psi1_R(pos);
                }
                else if (variable == 0) {
                    index->dr[i] = Psi1_R(pos);
                }
                else {    // variable == 1
                    index-> v[i] = Psi1_R(pos);
                }
            }
        }
    }
}

/* Function for carrying out 2LPT. Besides updating particle positions
** and velocities, the \Psi2 potential (in Fourier space) will be
** available after this function returns. A populated \Psi1 (in Fourier
** space) is expected as input.
*/
void carryout2LPT(
    PKD pkd, MDLFFT fft, basicParticleArray output, real_array_t *R, complex_array_t *K,
    int gridIndexPhi1, int gridIndexPhi2, std::queue<tmpNote *> &notes,
    growthFactors growth, int nGrid, double dBoxSize, double a,
    int printIndent = 0
) {
    float kFundamental = 2.*M_PI/dBoxSize;
    float fftFactor = dBoxSize/pow(nGrid, 3.);
    auto Phi2_R = R[gridIndexPhi2];
    auto Phi2_K = K[gridIndexPhi2];
    float velocityFactor = a*a*csmExp2Hub(pkd->csm, a)*growth.f2;
    /* The 2LPT potential looks like
    **   \nabla^2 \Phi2 = D2/D1^2 (
    **       - \Phi1,00 \Phi1,11
    **       - \Phi1,11 \Phi1,22
    **       - \Phi1,22 \Phi1,00
    **       + \Phi1,01 \Phi1,01
    **       + \Phi1,12 \Phi1,12
    **       + \Phi1,20 \Phi1,20
    **   )
    ** where all growth factors are taken to be positive.
    */
    if (pkd->Self() == 0) printf("%*sBuilding potential\n", printIndent, "");
    float potentialFactor = growth.D2/pow(growth.D1, 2.);
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi2, LPT_TERM_OP::EQ,    -1., {{gridIndexPhi1, 0, 0}, {gridIndexPhi1, 1, 1}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi2, LPT_TERM_OP::MIN_EQ, 1., {{gridIndexPhi1, 1, 1}, {gridIndexPhi1, 2, 2}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi2, LPT_TERM_OP::MIN_EQ, 1., {{gridIndexPhi1, 2, 2}, {gridIndexPhi1, 0, 0}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi2, LPT_TERM_OP::PLS_EQ, 1., {{gridIndexPhi1, 0, 1}, {gridIndexPhi1, 0, 1}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi2, LPT_TERM_OP::PLS_EQ, 1., {{gridIndexPhi1, 1, 2}, {gridIndexPhi1, 1, 2}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental, gridIndexPhi2, LPT_TERM_OP::PLS_EQ, 1.,
    {{gridIndexPhi1, 2, 0}, {gridIndexPhi1, 2, 0}});
    FFTLaplacianInverse(pkd, fft, Phi2_R, Phi2_K, kFundamental, fftFactor*potentialFactor);
    /* Construct displacement field from potential,
    ** displace positions and boost velocities.
    */
    if (pkd->Self() == 0) printf("%*sDisplacing positions and boosting velocities\n", printIndent, "");
    auto note = getTmpNote(notes);
    note->state = TMP_STATE::INACTIVE;  // flag as being non-reusable
    auto Psi2_R = R[note->gridIndex];
    auto Psi2_K = K[note->gridIndex];

    // !!!
    auto enabled = checkEnabledTerms(pkd);
    if (!enabled.term2) return;

    for (auto i = 0; i < 3; i++) {
        diffIFFT(pkd, fft, Phi2_K, Psi2_K, kFundamental, i);
        for (auto index = output.begin(); index != output.end(); index++) {
            auto pos = index.position();
            index->dr[i] += Psi2_R(pos);
            index-> v[i] += velocityFactor*Psi2_R(pos);
        }
    }
}

/* Function for carrying out 3LPT ('a' term only). Populated \Psi1
** and \Psi2 (both in Fourier space) are expected as input.
*/
void carryout3aLPT(
    PKD pkd, MDLFFT fft, basicParticleArray output, real_array_t *R, complex_array_t *K,
    int gridIndexPhi1, int gridIndexPhi2, int gridIndexPhi3a, std::queue<tmpNote *> &notes,
    growthFactors growth, int nGrid, double dBoxSize, double a,
    int printIndent = 0
) {
    float kFundamental = 2.*M_PI/dBoxSize;
    float fftFactor = dBoxSize/pow(nGrid, 3.);
    auto Phi3a_R = R[gridIndexPhi3a];
    auto Phi3a_K = K[gridIndexPhi3a];
    float velocityFactor = a*a*csmExp2Hub(pkd->csm, a)*growth.f3a;
    /* The 3LPT (a) potential looks like
    **   \nabla^2 \Phi3a = D3a/D1^3 (
    **       +   \Phi1,20 \Phi1,20 \Phi1,11
    **       -   \Phi1,11 \Phi1,22 \Phi1,00
    **       +   \Phi1,00 \Phi1,12 \Phi1,12
    **       - 2 \Phi1,12 \Phi1,20 \Phi1,01
    **       +   \Phi1,01 \Phi1,01 \Phi1,22
    **   )
    ** where all growth factors are taken to be positive.
    */
    if (pkd->Self() == 0) printf("%*sBuilding potential (a)\n", printIndent, "");
    float potentialFactor = growth.D3a/pow(growth.D1, 3.);
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi3a, LPT_TERM_OP::EQ,     1., {{gridIndexPhi1, 2, 0}, {gridIndexPhi1, 2, 0}, {gridIndexPhi1, 1, 1}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi3a, LPT_TERM_OP::MIN_EQ, 1., {{gridIndexPhi1, 1, 1}, {gridIndexPhi1, 2, 2}, {gridIndexPhi1, 0, 0}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi3a, LPT_TERM_OP::PLS_EQ, 1., {{gridIndexPhi1, 0, 0}, {gridIndexPhi1, 1, 2}, {gridIndexPhi1, 1, 2}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi3a, LPT_TERM_OP::MIN_EQ, 2., {{gridIndexPhi1, 1, 2}, {gridIndexPhi1, 2, 0}, {gridIndexPhi1, 0, 1}});
    handleLPTTerm(pkd, fft, R, K, notes, kFundamental,
    gridIndexPhi3a, LPT_TERM_OP::PLS_EQ, 1., {{gridIndexPhi1, 0, 1}, {gridIndexPhi1, 0, 1}, {gridIndexPhi1, 2, 2}});
    FFTLaplacianInverse(pkd, fft, Phi3a_R, Phi3a_K, kFundamental, fftFactor*potentialFactor);
    /* Construct displacement field from potential,
    ** displace positions and boost velocities.
    */
    if (pkd->Self() == 0) printf("%*sDisplacing positions and boosting velocities\n", printIndent, "");
    auto note = getTmpNote(notes);
    note->state = TMP_STATE::INACTIVE;  // flag as being non-reusable
    auto Psi3a_R = R[note->gridIndex];
    auto Psi3a_K = K[note->gridIndex];

    // !!!
    auto enabled = checkEnabledTerms(pkd);
    if (!enabled.term3a) return;

    for (auto i = 0; i < 3; i++) {
        diffIFFT(pkd, fft, Phi3a_K, Psi3a_K, kFundamental, i);
        for (auto index = output.begin(); index != output.end(); index++) {
            auto pos = index.position();
            index->dr[i] += Psi3a_R(pos);
            index-> v[i] += velocityFactor*Psi3a_R(pos);
        }
    }
}

/* Function responsible for generating particle initial conditions
** through Lagrangian perturbation theory.
*/
int pkdGenerateIC_new(  // !!! remove _new
    PKD pkd, MDLFFT fft, int iSeed, int bFixed, float fPhase, int nGrid, int iLPT, double dBoxSize,
    double a, int nTf, double *tk, double *tf, double *noiseMean, double *noiseCSQ,
    int printIndent = 4
) {
    int bExtraTmpGridFor1LPT = 1;  // allow one additional grid when doing 1LPT only?
    if ((iLPT < 1) || (iLPT > 3)) {
        fprintf(stderr, "pkdGenerateIC() implements 1LPT, 2LPT, 3LPT, but iLPT = %d\n", iLPT);
        abort();
    }
    /* Prepare growth factors and rates */
    growthFactors growth;
    csmComoveGrowth(pkd->csm, a,
                    &growth.D1, &growth.D2, &growth.D3a, &growth.D3b, &growth.D3c,
                    &growth.f1, &growth.f2, &growth.f3a, &growth.f3b, &growth.f3c);
    /* Set up grids.
    ** Note that "data" points to the same block for all threads.
    ** Particles will overlap K[0] through K[5] eventually.
    */
    int gridIndexPhi1 = -1;   // 1LPT potential (persistent)
    int gridIndexPhi2 = -1;   // 2LPT potential (persistent)
    int gridIndexPhi3 = -1;   // 3LPT potentials
    int gridIndexTmp0 = -1;   // temporary grid
    int gridIndexTmp1 = -1;   // temporary grid
    int nGrids = 6; // particle positions and velocities
    if (iLPT >= 1) {
        /* For 1LPT we need 1 potential grid and one temporary grid */
        gridIndexPhi1 = nGrids++;
        gridIndexTmp0 = nGrids++;
        /* We can save time by having one additional temporary grid */
        gridIndexTmp1 = gridIndexTmp0;
        if (iLPT == 1 && bExtraTmpGridFor1LPT) {
            gridIndexTmp1 = nGrids++;
        }
    }
    if (iLPT >= 2) {
        /* For 2LPT we need 1 additional potential grid
        ** and one additional temporary grid.
        */
        gridIndexPhi2 = nGrids++;
        gridIndexTmp1 = nGrids++;
    }
    if (iLPT >= 3) {
        /* For 3LPT we need 1 additional potential grid */
        gridIndexPhi3 = nGrids++;
    }
    real_array_t    R[nGrids];
    complex_array_t K[nGrids];
    GridInfo G(pkd->mdl, fft);
    auto data = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl, 0, 0, pkd->particles));
    for (auto i = 0; i < nGrids; i++) {
        G.setupArray(data, K[i]);
        G.setupArray(data, R[i]);
        data += fft->rgrid->nLocal;
    }
    /* Create a local view of our part of the output array */
    basicParticleArray output = getOutputArray(pkd, G, R[0]);
    int nLocal = blitz::product(R[0].shape());
    /* The potentials are constructed from combinations of
    ** differentiated lower-order potentials. During potential
    ** constructions we may sometimes reuse old results stored in the
    ** temporary grids. The notes are used to keep a record of what the
    ** temporary grids currently holds (if anything).
    */
    std::queue<tmpNote *> notes;
    tmpNote note0 {gridIndexTmp0};
    tmpNote note1 {gridIndexTmp1};
    notes.push(&note0);
    notes.push(&note1);
    /* Carry out the LPT IC generation, one order at a time */
    if (iLPT >= 1) {
        /* 1LPT (Zeldovich) */
        if (pkd->Self() == 0) printf("%*sCarrying out 1LPT\n", printIndent, "");
        carryout1LPT(
            pkd, fft, output, R, K, gridIndexPhi1, notes,
            growth, iSeed, bFixed, fPhase, nGrid, dBoxSize, a, nTf, tk, tf, noiseMean, noiseCSQ,
            printIndent + 4
        );
    }
    if (iLPT >= 2) {
        /* 2LPT */
        if (pkd->Self() == 0) printf("%*sCarrying out 2LPT\n", printIndent, "");
        carryout2LPT(
            pkd, fft, output, R, K, gridIndexPhi1, gridIndexPhi2, notes,
            growth, nGrid, dBoxSize, a,
            printIndent + 4
        );
    }
    if (iLPT >= 3) {
        if (pkd->Self() == 0) printf("%*sCarrying out 3LPT\n", printIndent, "");
        /* 3LPT ('a' term) */
        carryout3aLPT(
            pkd, fft, output, R, K, gridIndexPhi1, gridIndexPhi2, gridIndexPhi3, notes,
            growth, nGrid, dBoxSize, a,
            printIndent + 4
        );
    }
    /* Done with LPT */
    return nLocal;
}




// !!! ORI BELOW
int pkdGenerateIC(PKD pkd, MDLFFT fft, int iSeed, int bFixed, float fPhase, int nGrid, int iLPT, double dBoxSize,
                  double a, int nTf, double *tk, double *tf, double *noiseMean, double *noiseCSQ) {


    // !!!
    char *newlptscheme_str;
    newlptscheme_str = getenv("newlptscheme");
    if (newlptscheme_str != NULL && strcmp(newlptscheme_str, "True") == 0)
        return pkdGenerateIC_new(pkd, fft, iSeed, bFixed, fPhase, nGrid, iLPT, dBoxSize,
                                 a, nTf, tk, tf, noiseMean, noiseCSQ);


    // !!!
    if ((iLPT < 1) || (iLPT > 3)) {
        fprintf(stderr, "Old pkdGenerateIC not suitable for iLPT = %d\n", iLPT);
        abort();
    }

    // !!!
    auto enabled = checkEnabledTerms(pkd);


    double pi = 4.*atan(1.);
    CSM csm = pkd->csm;
    int bClass = csm->val.classData.bClass;
    mdlClass *mdl = reinterpret_cast<mdlClass *>(pkd->mdl);
    NoiseGenerator ng(iSeed, bFixed, fPhase);
    /* Prepare growth factors and transfer function */
    double D1_0, D2_0, D3a_0, D3b_0, D3c_0, f1_0, f2_0, f3a_0, f3b_0, f3c_0;
    double D1_a, D2_a, D3a_a, D3b_a, D3c_a, f1_a, f2_a, f3a_a, f3b_a, f3c_a;
    csmComoveGrowth(csm, 1., &D1_0, &D2_0, &D3a_0, &D3b_0, &D3c_0, &f1_0, &f2_0, &f3a_0, &f3b_0, &f3c_0);
    csmComoveGrowth(csm,  a, &D1_a, &D2_a, &D3a_a, &D3b_a, &D3c_a, &f1_a, &f2_a, &f3a_a, &f3b_a, &f3c_a);
    double dOmega = csm->val.dOmega0/(a*a*a*pow(csmExp2Hub(csm, a)/csm->val.dHubble0, 2.0));
    double f1_approx =    pow(dOmega, 5./ 9.);
    double f2_approx = 2.*pow(dOmega, 6./11.);
    PowerTransfer transfer(csm, a, nTf, tk, tf);
    if (bClass) {
        double D1_internal, D2_internal, D3a_internal, D3b_internal, D3c_internal;
        double f1_internal, f2_internal, f3a_internal, f3b_internal, f3c_internal;
        csm->val.classData.bClassGrowth = 0;
        csmComoveGrowth(csm, a,
                        &D1_internal, &D2_internal, &D3a_internal, &D3b_internal, &D3c_internal,
                        &f1_internal, &f2_internal, &f3a_internal, &f3b_internal, &f3c_internal);
        csm->val.classData.bClassGrowth = 1;
        if (mdl->Self() == 0) {
            printf("f1 = %.12g (CLASS), %.12g (internal but CLASS background), %.12g (approx)\n", f1_a, f1_internal, f1_approx);
            if (iLPT)
                printf("f2 = %.12g (CLASS), %.12g (internal but CLASS background), %.12g (approx)\n", f2_a, f2_internal, f2_approx);
        }
    }
    else {
        if (mdl->Self() == 0) {
            printf("f1 = %.12g (internal), %.12g (approx)\n", f1_a, f1_approx);
            if (iLPT)
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
    if (iLPT) {
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
        amp = -D2_a/(D1_a*D1_a)*fft_normalization; // !!! negated due to old convention
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

                // !!!
                if (!enabled.term2) index->dr[i] = 0.0;

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

            // !!!
            if (!enabled.term1) break;

            auto pos = index.position();
            for (auto i = 0; i < 3; i++) {
                if (quantity == 0) {
                    if (iLPT) index->dr[i] += R[7 + i](pos);
                    else       index->dr[i]  = R[7 + i](pos);
                }
                else if (iLPT) index->v[i] += R[7 + i](pos);
                else            index->v[i]  = R[7 + i](pos);
                if (!bClass) {
                    if (iLPT)  index->v[i] += vel_factor1*R[7 + i](pos);
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
        icUp.dBaryonFraction = in->dBaryonFraction;
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
        int nBucket;
        if (in->bICgas) {
            fGasMass = in->fMass*in->dBaryonFraction;
            fDarkMass = in->fMass*(1.0 - in->dBaryonFraction);
            fGasSoft = in->fSoft * pow(in->dBaryonFraction, 1./3.);
            fDarkSoft = in->fSoft * pow(1.0 - in->dBaryonFraction, 1./3.);
            nBucket = in->nBucket;
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
                blitz::TinyVector<double,3> r(temp.dr[0] + (temp.ix+0.5) * inGrid - 0.5,
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

                // Use N-GenIC trick to keep the centre of mass of the original
                // particle at its original position. This is needed for displaced
                // dark matter and gas particles in case the tree is built down to
                // one particle per leaf. Here we displace the particles by a length
                // proportional to half the grid cell size and weighted by
                // the particle mass to preserve the centre of mass and conserve
                // linear and angular momentum and kinetic energy.
                p.set_position(p.position() - inGrid*0.5*in->dBaryonFraction);
                pgas.set_position(pgas.position() + inGrid*0.5*(1.0 - in->dBaryonFraction));

                auto &VelGas = pgas.velocity();
                // Change the scale factor dependency
                double a_m1 = 1./in->dExpansion;
                VelGas *= a_m1;

                /* Fill the meshless::FIELDS with some initial values */
                double u = in->dInitialT * in->dTuFac;
                assert(pkd->particles.present(PKD_FIELD::oSph));
                auto &Sph = pgas.sph();
                Sph.ElemMass[ELEMENT_H]  = in->dInitialH  * fGasMass;
#ifdef HAVE_HELIUM
                Sph.ElemMass[ELEMENT_He] = in->dInitialHe * fGasMass;
#endif
#ifdef HAVE_CARBON
                Sph.ElemMass[ELEMENT_C]  = in->dInitialC  * fGasMass;
#endif
#ifdef HAVE_NITROGEN
                Sph.ElemMass[ELEMENT_N]  = in->dInitialN  * fGasMass;
#endif
#ifdef HAVE_OXYGEN
                Sph.ElemMass[ELEMENT_O]  = in->dInitialO  * fGasMass;
#endif
#ifdef HAVE_NEON
                Sph.ElemMass[ELEMENT_Ne] = in->dInitialNe * fGasMass;
#endif
#ifdef HAVE_MAGNESIUM
                Sph.ElemMass[ELEMENT_Mg] = in->dInitialMg * fGasMass;
#endif
#ifdef HAVE_SILICON
                Sph.ElemMass[ELEMENT_Si] = in->dInitialSi * fGasMass;
#endif
#ifdef HAVE_IRON
                Sph.ElemMass[ELEMENT_Fe] = in->dInitialFe * fGasMass;
#endif
#ifdef HAVE_METALLICITY
                Sph.fMetalMass = in->dInitialMetallicity * fGasMass;
#endif
                Sph.Frho = 0.0;
                Sph.Fmom = 0.0;
                Sph.Fene = 0.0;
                Sph.E = u + 0.5*blitz::dot(VelGas,VelGas);
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
        move.nBucket = in->nBucket;
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
        move.dBaryonFraction = in->dBaryonFraction;
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
                               in->nGrid, in->iLPT, in->dBoxSize, in->dExpansion, in->nTf,
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
