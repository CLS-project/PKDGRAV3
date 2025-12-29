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
#include "pkd_config.h"
#include "pmforces.h"
#include "pkd.h"

#include <vector>
#include "core/aweights.hpp"
#include "ic/whitenoise.hpp"
using namespace gridinfo;
using namespace blitz;

auto ServiceGridCreateFFT::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == 0);
    pst->lcl()->pkd->GridCreateFFT(in->nGrid);
    return 0;
}

void pkdContext::GridCreateFFT(int nGrid) {
    assert(fft == nullptr);
    fft = mdl->FFTInitialize(nGrid, nGrid, nGrid, 0, 0);
}

auto ServiceGridDeleteFFT::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    assert(vin.size() == 0);
    assert(vout.size() == 0);
    pst->lcl()->pkd->GridDeleteFFT();
    return 0;
}

void pkdContext::GridDeleteFFT() {
    assert(fft != nullptr);
    mdl->FFTFinish(fft);
    fft = nullptr;
}

auto ServiceLinearKick::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == 0);
    pst->lcl()->pkd->LinearKick(in->dtOpen, in->dtClose);
    return 0;
}

// A blitz++ friendly wrap function. returns "ik" given array index
// Range: (-iNyquist,iNyquist] where iNyquist = m/2
static float fwrap(float v, float m) {
    return v - (v > m * 0.5 ? m : 0);
}
BZ_DEFINE_BINARY_FUNC(Fn_fwrap, fwrap)
BZ_DECLARE_ARRAY_ET_BINARY(fwrap, Fn_fwrap)
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(fwrap, Fn_fwrap, float)

/********************************************************************************\
*
* GRID Construction
*
\********************************************************************************/

void pkdContext::GenerateLinearGrid(double a, double a_next, double Lbox, int iSeed, int bFixed, float fPhase, int bRho) {
    /* If bRho == 0, we generate the \delta field,
    ** otherwise we generate the \delta\rho field.
    ** If a == a_next, we generate the field at this a.
    ** Otherwize, we generate the weighted average of the
    ** field over the interval [a, a_next]. Note that
    ** averaging is only implemented for bRho == 1.
    */
    if (!bRho && a != a_next) {
        fprintf(stderr, "WARNING: In GenerateLinearGrid(): Averaging of \\delta fields not implemented\n");
        abort();
    }
    /* For the sake of performance we precompute the (possibly averaged)
    ** field at the |k|'s at which the perturbations are tabulated.
    */
    size_t size = csm->val.classData.perturbations.size_k;
    double* logk = (double*)malloc(sizeof(double) * size);
    double* field = (double*)malloc(sizeof(double) * size);
    size_t i;
    double k;
    gsl_interp_accel* acc;
    gsl_spline* spline;
    for (i = 0; i < size; i++) {
        k = csm->val.classData.perturbations.k[i];
        logk[i] = log(k);
        if (bRho) /* Use the \delta\rho field, possibly averaged */
            field[i] = csmDeltaRho_lin(csm, a, a_next, k);
        else /* Use the \delta field */ field[i] = csmDelta_lin(csm, a, k);
        /* The cubic splining is much better without the analytic zeta */
        field[i] /= csmZeta(csm, k);
    }
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(spline, logk, field, size);
    /* Generate grid */
    double ix, iy, iz, i2;
    double iLbox = 2 * M_PI / Lbox;

    GridInfo G(mdl, fft);
    complex_array_t K;
    G.setupArray(mdl->SetArray(Ephemeral<FFTW3(real)>()), K);
    NoiseGenerator ng(iSeed, bFixed, fPhase);
    ng.FillNoise(K, fft->rgrid->n3);
    for (auto index = K.begin(); index != K.end(); ++index) {
        auto pos = index.position();
        iz = fwrap(pos[2], fft->rgrid->n3); // Range: (-iNyquist,iNyquist]
        iy = fwrap(pos[1], fft->rgrid->n2);
        ix = fwrap(pos[0], fft->rgrid->n1);
        i2 = ix * ix + iy * iy + iz * iz;
        if (i2 > 0) {
            k = sqrt((double)i2) * iLbox;
            *index *= csmZeta(csm, k) * gsl_spline_eval(spline, log(k), acc);
        }
        else *index = 0.0;
    }
    /* Cleanup */
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
    free(logk);
    free(field);
}

static double deconvolveLinWindow(int i, int nGrid) {
    double win = M_PI * i / nGrid;
    if (win > 0.1) win = win / sin(win);
    else win = 1.0 / (1.0 - win * win / 6.0 * (1.0 - win * win / 20.0 * (1.0 - win * win / 76.0)));
#if defined(USE_NGP_LIN)
    return win;
#elif defined(USE_CIC_LIN)
    return win * win;
#elif defined(USE_TSC_LIN)
    return win * win * win;
#else
    return win * win * win * win;
#endif
}

/*
 * Green Function for the Laplacian operator in Fourier space */
static double green(int i, int jj, int kk, int nGrid) {
    double g = pow2(sin(M_PI * i / (1.0 * nGrid)));
    g += pow2(sin(M_PI * jj / (1.0 * nGrid)));
    g += pow2(sin(M_PI * kk / (1.0 * nGrid)));
    g *= 4 * nGrid * nGrid;
    if (g == 0.0) return 0.0;
    else return -1.0 / g;
}

auto ServiceSetLinearGrid::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() == 0);
    pst->lcl()->pkd->SetLinearGrid(in->a0, in->a, in->a1, in->dBSize, in->nGrid, in->iSeed, in->bFixed, in->fPhase);
    return 0;
}

void pkdContext::SetLinearGrid(double a0, double a, double a1, double dBSize, int nGrid, int iSeed, int bFixed, float fPhase) {
    /* Grid coordinates in real space :      [0, nGrid].[0, nGrid].[0, nGrid] */
    mdlGridCoord rfirst, rlast;
    /* Grid coordinates in Fourier space : [O, Nyquist].[0, nGrid].[0, nGrid] */
    mdlGridCoord kfirst, klast, kindex;

    /* Scale factors and normalization */
    const double dNormalization = a * a * a * dBSize;

    mdlGridCoordFirstLast(mdl, fft->rgrid, &rfirst, &rlast, 1);
    mdlGridCoordFirstLast(mdl, fft->kgrid, &kfirst, &klast, 0);

    /* Imprint the density grid of the linear species */
    int bRho = 1; /* Generate the \delta\rho field */
    GenerateLinearGrid(a0, a1, dBSize, iSeed, bFixed, fPhase, bRho);

    /*
     * Define the grid arrays : only 3 grids are stored
     * in memory, the other ones are defined in order to
     * have an explicit naming in the code
     */
    auto cDelta_lin_field = mdl->SetArray(Ephemeral<FFTW3(complex)>(), klast.i);
    /* Remember, the grid is now transposed to x,z,y (from x,y,z) */
    auto cForceY = mdl->SetArray(cDelta_lin_field + fft->kgrid->nLocal, klast.i);
    auto cForceZ = mdl->SetArray(cForceY + fft->kgrid->nLocal, klast.i);

    int idx, i, j, jj, k, kk;
    const int iNyquist = nGrid / 2;
    double rePotential, imPotential;
    double dDifferentiate, dPoissonSolve;
    double win_j, win_k;
    /* Here starts the Poisson solver */
    i = j = k = -1;
    for (kindex = kfirst; !mdlGridCoordCompare(&kindex, &klast); mdlGridCoordIncrement(&kindex)) {
        idx = kindex.i;
        if (j != kindex.z) {
            j = kindex.z;
            jj = j > iNyquist ? j - nGrid : j;
            win_j = deconvolveLinWindow(jj, nGrid);
        }
        if (k != kindex.y) {
            k = kindex.y;
            kk = k > iNyquist ? k - nGrid : k;
            win_k = deconvolveLinWindow(kk, nGrid);
        }
        i = kindex.x;
        double win = deconvolveLinWindow(i, nGrid) * win_j * win_k;
        /* Green Function for a discrete Laplacian operator */
        dPoissonSolve = 4 * M_PI * green(i, jj, kk, nGrid) * dNormalization * win;
        /* Solve Poisson equation */

        rePotential = cDelta_lin_field[idx][0] * dPoissonSolve;
        imPotential = cDelta_lin_field[idx][1] * dPoissonSolve;
        /* Differentiaite in Y direction */
        dDifferentiate = nGrid * sin(2 * M_PI * jj / (1.0 * nGrid));
        cForceY[idx][0] = dDifferentiate * imPotential;
        cForceY[idx][1] = -dDifferentiate * rePotential;

        /* Differentiate in Z direction */
        dDifferentiate = nGrid * sin(2 * M_PI * kk / (1.0 * nGrid));
        cForceZ[idx][0] = dDifferentiate * imPotential;
        cForceZ[idx][1] = -dDifferentiate * rePotential;

        /*
         * Differentiate in X direction (over-write the
         * delta_lin field)
         */
        dDifferentiate = nGrid * sin(2 * M_PI * i / (1.0 * nGrid));
        cDelta_lin_field[idx][0] = dDifferentiate * imPotential;
        cDelta_lin_field[idx][1] = -dDifferentiate * rePotential;
    }
    mdl->IFFT(fft, cForceY);
    mdl->IFFT(fft, cForceZ);
    mdl->IFFT(fft, cDelta_lin_field);
}

using float3_t = blitz::TinyVector<float, 3>;

struct tree_node : public KDN {
    bool is_cell() { return iLower != 0; }
    bool is_bucket() { return iLower == 0; }
};

template <int Order, typename F>
static float interpolate(base_type::force_array_t& forces, const F r[3]) {
    AssignmentWeights<Order, F> Hx(r[0]), Hy(r[1]), Hz(r[2]);
    float force = 0;
    for (int i = 0; i <= Order; ++i) {
        for (int j = 0; j <= Order; ++j) {
            for (int k = 0; k <= Order; ++k) {
                force += forces(Hx.i + i, Hy.i + j, Hz.i + k) * Hx.H[i] * Hy.H[j] * Hz.H[k];
            }
        }
    }
    return force;
}

template <typename F>
static float force_interpolate(base_type::force_array_t& forces, const F r[3], int iAssignment = 3) {
    float f;
    switch (iAssignment) {
    case 0:
        f = interpolate<0, F>(forces, r);
        break;
    case 1:
        f = interpolate<1, F>(forces, r);
        break;
    case 2:
        f = interpolate<2, F>(forces, r);
        break;
    case 3:
        f = interpolate<3, F>(forces, r);
        break;
    default:
        f = 0;
        assert(iAssignment >= 0 && iAssignment <= 3);
        abort();
    }
    return f;
}

// Fetch all forces into this local grid. Could be optimized by bulk fetches.

void pkdContext::fetch_forces(int cid, int nGrid, base_type::force_array_t& forces, const base_type::shape_vector_t& lower) {
    auto wrap = [&nGrid](int i) {
        if (i >= nGrid) i -= nGrid;
        else if (i < 0) i += nGrid;
        return i;
    };
    for (auto i = forces.begin(); i != forces.end(); ++i) {
        base_type::shape_vector_t loc = i.position() + lower;
        loc[0] = wrap(loc[0]);
        loc[1] = wrap(loc[1]);
        loc[2] = wrap(loc[2]);
        auto id = mdlFFTrId(mdl, fft, loc[0], loc[1], loc[2]);
        auto idx = mdlFFTrIdx(mdl, fft, loc[0], loc[1], loc[2]);
        auto p = mdl->Fetch<base_type::force_t>(cid, idx, id);
        *i = *p;
    }
}

void pkdContext::LinearKick(vel_t dtOpen, vel_t dtClose, int iAssignment) {
    const std::size_t maxSize = 100000; // We would like this to remain in L2 cache
    std::vector<float> dataX, dataY, dataZ;
    dataX.reserve(maxSize);
    dataY.reserve(maxSize);
    dataZ.reserve(maxSize);
    base_type::shape_vector_t index;
    base_type::position_vector_t ifPeriod = 1.0 / fPeriod;
    int nGrid = fft->rgrid->n1;
    assert(iAssignment >= 0 && iAssignment <= 3);

    int iLocal = mdl->Core() ? 0 : fft->rgrid->nLocal;
    auto forceX = mdl->SetArray(Ephemeral<FFTW3(real)>(), iLocal);
    auto forceY = mdl->SetArray(forceX + fft->rgrid->nLocal, iLocal);
    auto forceZ = mdl->SetArray(forceY + fft->rgrid->nLocal, iLocal);
    mdl->CacheInitialize(CID_GridLinFx, std::make_shared<mdl::arrayCACHEhelper>(forceX, iLocal));
    mdl->CacheInitialize(CID_GridLinFy, std::make_shared<mdl::arrayCACHEhelper>(forceY, iLocal));
    mdl->CacheInitialize(CID_GridLinFz, std::make_shared<mdl::arrayCACHEhelper>(forceZ, iLocal));
    auto pad = (iAssignment + 1) / 2;
    std::vector<std::uint32_t> stack;
    stack.push_back(ROOT);
    while (!stack.empty()) {
        auto kdn = tree[stack.back()];
        stack.pop_back(); // Go to the next node in the tree
        auto bnd = kdn->bound();
        base_type::shape_vector_t ilower = base_type::shape_vector_t(floor((bnd.lower() * ifPeriod + 0.5) * nGrid)) - pad;
        base_type::shape_vector_t iupper = base_type::shape_vector_t(floor((bnd.upper() * ifPeriod + 0.5) * nGrid)) + pad;
        base_type::shape_vector_t ishape = iupper - ilower + 1;
        float3_t flower = ilower;
        std::size_t size = blitz::product(ishape);

        if (size > maxSize) {       // This cell is too large, so we split it and move on
            assert(kdn->is_cell()); // At the moment we cannot handle enormous buckets
            stack.push_back(kdn->rchild());
            stack.push_back(kdn->lchild());
        }
        else {                  // Assign the mass for this range of particles
            dataX.resize(size); // Hold the right number of masses
            dataY.resize(size); // Hold the right number of masses
            dataZ.resize(size); // Hold the right number of masses
            base_type::force_array_t forcesX(dataX.data(), ishape, blitz::neverDeleteData, blitz::ColumnMajorArray<3>());
            base_type::force_array_t forcesY(dataY.data(), ishape, blitz::neverDeleteData, blitz::ColumnMajorArray<3>());
            base_type::force_array_t forcesZ(dataZ.data(), ishape, blitz::neverDeleteData, blitz::ColumnMajorArray<3>());
            fetch_forces(CID_GridLinFx, nGrid, forcesX, ilower);
            fetch_forces(CID_GridLinFy, nGrid, forcesY, ilower);
            fetch_forces(CID_GridLinFz, nGrid, forcesZ, ilower);
            for (auto& p : *kdn) { // All particles in this tree cell
                auto& v = p.velocity();
                float3_t r = p.position();
                r = (r * ifPeriod + 0.5) * nGrid - flower; // Scale and shift to fit in subcube
                v[0] += (dtOpen + dtClose) * force_interpolate(forcesX, r.data(), iAssignment);
                v[1] += (dtOpen + dtClose) * force_interpolate(forcesY, r.data(), iAssignment);
                v[2] += (dtOpen + dtClose) * force_interpolate(forcesZ, r.data(), iAssignment);
            }
        }
    }
    mdl->FinishCache(CID_GridLinFx);
    mdl->FinishCache(CID_GridLinFy);
    mdl->FinishCache(CID_GridLinFz);
}

void pkdContext::MeasureLinearPk(int nGrid, double dA, double dBoxSize, int nBins, int iSeed, int bFixed, float fPhase, double* fK, double* fPower,
                                 uint64_t* nPower) {
    mdlGridCoord first, last, index;
    double ak;
    int i, j, k, idx, ks;
    int iNyquist = nGrid / 2;

    mdlGridCoordFirstLast(mdl, fft->rgrid, &first, &last, 1);
    /* Generate the grid of the linear species again,
    ** this should not be done this way for performances.
    ** - Well, for the gravity computation what is generated
    **   on the grid is \delta\rho averaged over one time step.
    **   Here we want the \delta power spectrum at a specific
    **   time, and so we need to generate the \delta field at this
    **   time. We thus have to re-generate the grid.
    */
    int bRho = 0;
    double a_next = dA; /* Generate the \delta field at dA (no averaging) */
    GenerateLinearGrid(dA, a_next, dBoxSize, iSeed, bFixed, fPhase, bRho);

    /* Remember, the grid is now transposed to x,z,y (from x,y,z) */
    mdlGridCoordFirstLast(mdl, fft->kgrid, &first, &last, 0);
    auto fftDataK = mdl->SetArray(Ephemeral<FFTW3(complex)>(), last.i);
    for (i = 0; i < nBins; i++) {
        fK[i] = 0.0;
        fPower[i] = 0.0;
        nPower[i] = 0;
    }
    /*
    ** Calculate which slabs to process.  Because of zero padding,
    ** there may be nothing to process here, in which case ey will
    ** be less than sy.  This is fine.
    */
#ifdef LINEAR_PK
    double scale = nBins * 1.0 / iNyquist;
#else
    double scale = nBins * 1.0 / log(iNyquist + 1);
#endif
    double dBox2 = dBoxSize * dBoxSize;
    int jj, kk;
    i = j = k = -1;
    for (index = first; !mdlGridCoordCompare(&index, &last); mdlGridCoordIncrement(&index)) {
        if (j != index.z) {
            j = index.z;
            jj = j > iNyquist ? nGrid - j : j;
        }
        if (k != index.y) {
            k = index.y;
            kk = k > iNyquist ? nGrid - k : k;
        }
        i = index.x;
        ak = sqrt(i * i + jj * jj + kk * kk);
        ks = ak;
        if (ks >= 1 && ks <= iNyquist) {
#ifdef LINEAR_PK
            ks = floor((ks - 1.0) * scale);
#else
            ks = floor(log(ks) * scale);
#endif
            assert(ks >= 0 && ks < nBins);
            idx = index.i;
            double delta2 = dBox2 * (pow2(fftDataK[idx][0]) + pow2(fftDataK[idx][1]));
            fK[ks] += ak;
            fPower[ks] += delta2;
            nPower[ks] += 1;
        }
    }
}

auto ServiceMeasureLinearPk::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->lcl()->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    assert(vin.size() == sizeof(input));
    assert(vout.size() >= sizeof(output) * in->nBins);
    pkd->MeasureLinearPk(in->nGrid, in->dA, in->dBoxSize, in->nBins, in->iSeed, in->bFixed, in->fPhase, out->fK, out->fPower, out->nPower);
    std::fill(out->fK + in->nBins, out->fK + max_k_bins, 0.0);
    std::fill(out->fPower + in->nBins, out->fPower + max_k_bins, 0.0);
    std::fill(out->nPower + in->nBins, out->nPower + max_k_bins, 0);
    return sizeof(output);
}
auto ServiceMeasureLinearPk::Combine(std::span<std::byte> vout, std::span<const std::byte> vin, PST) -> int {
    auto out = reinterpret_cast<output*>(vout.data());
    auto in = reinterpret_cast<const output*>(vin.data());
    for (std::size_t i = 0; i < max_k_bins; i++) {
        out->fK[i] += in->fK[i];
        out->fPower[i] += in->fPower[i];
        out->nPower[i] += in->nPower[i];
    }
    return sizeof(output);
}