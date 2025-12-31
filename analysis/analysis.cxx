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
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include "pkd.h"
#include "analysis.h"
#include "ic/ic.h"

// #define SHAPES
#ifdef SHAPES
    #define USE_PCS     /* USE_PCS USE_TSC USE_CIC USE_NGP */
    #define USE_PCS_LIN /* USE_PCS_LIN USE_TSC_LIN USE_CIC_LIN USE_NGP_LIN */
static void transpose(double mat[3][3], double trans_mat[3][3]) {
    int i, j;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            trans_mat[i][j] = mat[j][i];
        }
    }
}

    #define ROTATE(a, i, j, k, l)                                                                                                                              \
        g = a[i][j];                                                                                                                                           \
        h = a[k][l];                                                                                                                                           \
        a[i][j] = g - s * (h + g * tau);                                                                                                                       \
        a[k][l] = h + s * (g - h * tau);

static void jacobi(double a[4][4], int n, double d[4], double v[4][4], int* nrot) {
    int j, iq, ip, i;
    double tresh, theta, tau, t, sm, s, h, g, c;
    double b[4], z[4];

    for (ip = 1; ip <= n; ip++) {
        for (iq = 1; iq <= n; iq++) v[ip][iq] = 0.0;
        v[ip][ip] = 1.0;
    }
    for (ip = 1; ip <= n; ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    *nrot = 0;
    for (i = 1; i <= 100; i++) {
        sm = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++) sm += fabs(a[ip][iq]);
        }
        if (sm <= 1.0e-12) {
            return;
        }
        if (i < 4) tresh = 0.2 * sm / (n * n);
        else tresh = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) && fabs(d[iq]) + g == fabs(d[iq])) a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if (fabs(h) + g == fabs(h)) t = (a[ip][iq]) / h;
                    else {
                        theta = 0.5 * h / (a[ip][iq]);
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0) t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 1; j <= ip - 1; j++) {
                        ROTATE(a, j, ip, j, iq)
                    }
                    for (j = ip + 1; j <= iq - 1; j++) {
                        ROTATE(a, ip, j, j, iq)
                    }
                    for (j = iq + 1; j <= n; j++) {
                        ROTATE(a, ip, j, iq, j)
                    }
                    for (j = 1; j <= n; j++) {
                        ROTATE(v, j, ip, j, iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip = 1; ip <= n; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    printf("<error in jacobi>\n");
}
#endif

/*
** Combiner cache functions
*/

class HelperProfileBins1 : public mdl::arrayCACHEhelper {
  protected:
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto pBin1 = static_cast<PROFILEBIN*>(dst);
        const auto pBin2 = static_cast<const PROFILEBIN*>(src);
        pBin1->dMassInBin += pBin2->dMassInBin;
        pBin1->nParticles += pBin2->nParticles;
        pBin1->L += pBin2->L;
        pBin1->vel_radial += pBin2->vel_radial;
        pBin1->vel_radial_sigma += pBin2->vel_radial_sigma;
    }
    virtual void init(void* dst) override final {
        auto pBin = static_cast<PROFILEBIN*>(dst);
        pBin->dMassInBin = 0.0;
        pBin->L = 0.0;
        pBin->vel_radial = 0.0;
        pBin->vel_radial_sigma = 0.0;
        pBin->nParticles = 0;
    }

  public:
    HelperProfileBins1(PROFILEBIN* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

class HelperProfileBins2 : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto pBin = static_cast<PROFILEBIN*>(dst);
        pBin->vel_tang_sigma = 0.0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto pBin1 = static_cast<PROFILEBIN*>(dst);
        const auto pBin2 = static_cast<const PROFILEBIN*>(src);
        pBin1->vel_tang_sigma += pBin2->vel_tang_sigma;
    }

  public:
    HelperProfileBins2(PROFILEBIN* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};

#ifdef SHAPES
class HelperShapesBins1 : public mdl::arrayCACHEhelper {
  protected:
    virtual void init(void* dst) override final {
        auto pBin = static_cast<SHAPESBIN*>(dst);
        pBin->dMassEnclosed = 0.0;
        pBin->com = 0.0;
        pBin->dInertia = 0.0;
    }
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto pBin1 = static_cast<SHAPESBIN*>(dst);
        const auto pBin2 = static_cast<const SHAPESBIN*>(src);
        pBin1->com += pBin2->com;
        pBin1->dInertia += pBin2->dInertia;
    }

  public:
    HelperShapesBins1(SHAPESBIN* pArray, uint32_t nArray) : arrayCACHEhelper(pArray, nArray, true) {}
};
#endif

/*
** This function will calculate the distance between a particle and a
** reference point.  If a periodic boundary is in effect then the smallest
** possible distance is returned.
*/
double pkdGetDistance2(PKD pkd, particleStore::Particle& p, blitz::TinyVector<double, 3> dCenter, int bPeriodic) {
    blitz::TinyVector<double, 3> dx = p.position() - dCenter;
    // If a periodic wrap results in a smaller distance, then use that.
    if (bPeriodic) {
        for (auto j = 0; j < 3; j++) {
            auto dx2 = dx[j] + (dx[j] < 0.0 ? pkd->fPeriod[j] : -pkd->fPeriod[j]);
            if (dx2 * dx2 < dx[j] * dx[j]) dx[j] = dx2;
        }
    }
    return blitz::dot(dx, dx);
}

struct distance {
    float d2;
    uint32_t i;
    bool operator<(const distance& other) const { return d2 < other.d2; }
    bool operator<(float rhs) const { return d2 < rhs; }
    friend bool operator<(float lhs, const distance& rhs) { return lhs < rhs.d2; }
};

/*
** Use the Ephemeral structure to calculate the distance to each particle
** Sort by distance when finished.
*/
void pkdCalcDistance(PKD pkd, blitz::TinyVector<double, 3> dCenter, int bPeriodic) {
    auto ephemeral = pkd->Ephemeral<distance>();
    int i;

    /*
    ** Initialize the temporary particles.
    */
    for (i = 0; i < pkd->Local(); ++i) {
        auto p = pkd->particles[i];
        ephemeral[i].d2 = pkdGetDistance2(pkd, p, dCenter, bPeriodic);
        ephemeral[i].i = i;
    }
    std::sort(ephemeral, ephemeral + pkd->Local());
}

/*
** Return the mass weighted center of mass and velocity
*/
void pkdCalcCOM(PKD pkd, const blitz::TinyVector<double, 3>& dCenter, double dRadius, int bPeriodic, blitz::TinyVector<double, 3>& com,
                blitz::TinyVector<double, 3>& vcm, blitz::TinyVector<double, 3>& L, double& M, uint64_t& N) {
    using blitz::cross;
    com = 0.0;
    vcm = 0.0;
    L = 0.0;
    M = 0.0;
    N = 0;
    auto dRadius2 = dRadius * dRadius;
    for (auto& p : pkd->particles) {
        double m = p.mass();
        const blitz::TinyVector<double, 3> v = p.velocity();
        const auto r = p.position();
        auto d2 = pkdGetDistance2(pkd, p, dCenter, bPeriodic);
        if (dRadius < 0 || d2 < dRadius2) {
            M += m;
            com += m * r;
            vcm += m * v;
            auto T = cross(r, v);
            L += m * T;
            N++;
        }
    }
}

/*
** Return the mass of all particles
*/
void pkdCalcMtot(PKD pkd, double* M, uint64_t* N) {
    *N = pkd->Local();
    *M = std::accumulate(pkd->particles.begin(), pkd->particles.end(), 0.0, [](double a, auto& p) { return a + p.mass(); });
}

void pkdResetCOM(PKD pkd, blitz::TinyVector<double, 3> r_com, blitz::TinyVector<double, 3> v_com) {
    for (auto& p : pkd->particles) {
        auto r = p.position();
        p.set_position(p.position() - r_com);
        auto& v = p.velocity();
        v -= v_com;
    }
}

// Count the number of elements that are interior to r2 (between r2i and r2o).
int ServiceCountDistance::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->lcl()->pkd;
    assert(vin.size() == sizeof(input));
    assert(vout.size() == sizeof(output));
    decltype(distance::d2) r2i = in->dRadius2Inner, r2o = in->dRadius2Outer; // Avoid potential type conversions below
    auto ephemeral = pkd->Ephemeral<distance>();
    auto upper = std::lower_bound(ephemeral, ephemeral + pkd->Local(), r2o); // First element >= r2o
    auto lower = std::lower_bound(ephemeral, upper, r2i);                    // First element >= r2i
    *out = std::distance(lower, upper);
    return sizeof(output);
}

#ifdef SHAPES

double ell_distance2(const double* r, SHAPESBIN* pShape, double ba, double ca) {
    double dx[3], dx_rot[3];
    int i;

    for (i = 0; i < 3; i++) dx[i] = r[i] - pShape->ell_center[i];
    matrix_vector_mult(dx_rot, pShape->ell_matrix, dx);
    return dx_rot[0] * dx_rot[0] + dx_rot[1] * dx_rot[1] / ba / ba + dx_rot[2] * dx_rot[2] / ca / ca;
}

/*
** NB: We assume that shells do NOT cross.  In other words, if a particle is part
**     of bin n, then it is always part of bin n+1, n+2, ... etc.  In principle
**     this is not necessarily true, but it would be extremely odd for the centers
**     to be off by that much.  This assumption greatly improves the performance,
**     taking the algorithm from Nbin * Npart down to Npart.
*/
static void CalculateInertia(PKD pkd, int nBins, const double* dRadii, SHAPESBIN* shapesBins) {
    auto ephemeral = pkd->Ephemeral<distance>();
    local_t n = pkd->Local();
    SHAPESBIN* pShape;
    int i, j, k;
    int iBin;
    double ell_matrix[3][3], ell_matrix_inv[3][3];

    pkd->mdl->CacheInitialize(CID_SHAPES, std::make_shared<HelperShapesBins1>(shapesBins, pkd->Self() == 0 ? nBins : 0));

    /*
    ** The most efficient way to handle this is to do the calculations for all bins
    */
    for (i = iBin = 0; i < n; i++) {
        auto p = pkd->particles[ephemeral[i].i];
        auto r = p.position();
        double m = p.mass();

        // double r = dRadii[iBin];

        /* Find the bin: Assume that the last particle was close to the correct bin */
        while (ephemeral[i].d2 < dRadii[iBin] * dRadii[iBin] && iBin < nBins) ++iBin;
        while (ephemeral[i].d2 > dRadii[iBin] * dRadii[iBin]) --iBin;

        pShape = pkd->mdl->Acquire<SHAPESBIN>(CID_SHAPES, iBin, 0);
        pShape->dMassEnclosed += m;
        for (j = 0; j < 3; j++) {
            pShape->com[j] += m * r[j];
            for (k = 0; k <= j; k++) {
                pShape->dInertia[j][k] += m * r[j] * r[k];
            }
        }
        pkd->mdl->Release(CID_SHAPES, pShape);
    }
    pkd->mdl->FinishCache(CID_SHAPES);

    if (pkd->Self() == 0) {
        double inertia_cm[4][4];
        double evectors[4][4];
        double evalues[4];
        double VecProd[4], ScalProd;
        int nrot;
        int ia, ib, ic;
        for (i = iBin = 0; iBin < nBins; iBin++) {
            pShape = &shapesBins[iBin];

            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    if (k <= j)
                        inertia_cm[j + 1][k + 1] = ((pShape->dInertia[j][k] - pShape->com[j] * pShape->com[k] / pShape->dMassEnclosed) / pShape->dMassEnclosed);
                    else
                        inertia_cm[j + 1][k + 1] = ((pShape->dInertia[k][j] - pShape->com[k] * pShape->com[j] / pShape->dMassEnclosed) / pShape->dMassEnclosed);
                }
                evalues[j + 1] = inertia_cm[j + 1][j + 1];
            }
            jacobi(inertia_cm, 3, evalues, evectors, &nrot);
            if (evalues[1] >= evalues[2] && evalues[1] >= evalues[3]) {
                ia = 1;
                if (evalues[2] >= evalues[3]) {
                    ib = 2;
                    ic = 3;
                }
                else {
                    ib = 3;
                    ic = 2;
                }
            }
            else if (evalues[2] > evalues[1] && evalues[2] >= evalues[3]) {
                ia = 2;
                if (evalues[1] >= evalues[3]) {
                    ib = 1;
                    ic = 3;
                }
                else {
                    ib = 3;
                    ic = 1;
                }
            }
            else {
                ia = 3;
                if (evalues[1] >= evalues[2]) {
                    ib = 1;
                    ic = 2;
                }
                else {
                    ib = 2;
                    ic = 1;
                }
            }

            /* Check if Eigenvectors are righthanded in 3D :
               ev[ib] x ev[ic] = ev[ia] */

            VecProd[1] = evectors[2][ib] * evectors[3][ic] - evectors[3][ib] * evectors[2][ic];
            VecProd[2] = -evectors[1][ib] * evectors[3][ic] + evectors[3][ib] * evectors[1][ic];
            VecProd[3] = evectors[1][ib] * evectors[2][ic] - evectors[2][ib] * evectors[1][ic];
            ScalProd = evectors[1][ia] * VecProd[1] + evectors[2][ia] * VecProd[2] + evectors[3][ia] * VecProd[3];
            if (ScalProd < 0.0) {
                for (i = 0; i < 3; i++) evectors[i + 1][ia] = -evectors[i + 1][ia];
            }

            // double ba = sqrt((double)(evalues[ib]/evalues[ia])) ;
            // double ca = sqrt((double)(evalues[ic]/evalues[ia])) ;

            /* euler angles for a zyz rotation */
            // double theta = 180. / M_PI * acos((double) evectors[3][ic]);
            double phi = 180. / M_PI * acos((double)evectors[1][ic] / sqrt(evectors[1][ic] * evectors[1][ic] + evectors[2][ic] * evectors[2][ic]));
            double psi = 180. / M_PI *
                         acos((double)(-evectors[2][ic] * evectors[1][ib] + evectors[1][ic] * evectors[2][ib]) /
                              sqrt(evectors[1][ic] * evectors[1][ic] + evectors[2][ic] * evectors[2][ic]));

            /* inverse acos is only defined between 0 and pi therefore we must
               deal with pi to 2*pi */
            if (evectors[2][ic] < 0.0) phi = 360. - phi; /* phi always positive */
            if (evectors[3][ib] < 0.0) psi = 360. - psi; /* psi always positive */

            /*
            for(i = 0; i<3; i++){
            ell_center[i] = pShape->com[i] / pShape->dMassEnclosed;
              }
                      */
            for (i = 0; i < 3; i++) {
                ell_matrix_inv[i][0] = evectors[i + 1][ia];
                ell_matrix_inv[i][1] = evectors[i + 1][ib];
                ell_matrix_inv[i][2] = evectors[i + 1][ic];
            }
            transpose(ell_matrix_inv, ell_matrix);
        }
    }
}

/*
** Calculate shapes for existing profile bins
*/
void pkdShapes(PKD pkd, int nBins, const double* dCenter, const double* dRadii) {
    SHAPESBIN* shapesBins;
    int i, j, k;

    if (pkd->Self() == 0) {
        shapesBins = pkd->mdl->Malloc<SHAPESBIN>(nBins * sizeof(SHAPESBIN));
        assert(shapesBins != NULL);
        /* Start with the given center for every bin */
        for (i = 0; i < nBins; i++) {
            for (j = 0; j < 3; j++) {
                shapesBins[i].ell_center[j] = dCenter[j];
                for (k = 0; k <= j; k++) {
                    shapesBins[i].ell_matrix[j][k] = (j == k) ? 1.0 : 0.0;
                }
            }
        }
    }
    else shapesBins = NULL;

    pkd->mdl->CacheInitialize(CID_BIN, std::make_shared<mdl::arrayCACHEhelper>(pkd->profileBins, pkd->Self() ? 0 : nBins));

    CalculateInertia(pkd, nBins, dRadii, shapesBins);
    pkd->mdl->FinishCache(CID_BIN);

    if (pkd->Self() == 0) {
        pkd->mdl->Free(shapesBins);
    }
}
#endif

/*
** Density Profile
*/
void pkdProfile(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, const blitz::TinyVector<double, 3>& dCenter, const double* dRadii, int nBins,
                const blitz::TinyVector<double, 3>& com, const blitz::TinyVector<double, 3>& vcm, const blitz::TinyVector<double, 3>& L) {
    using blitz::cross;
    using blitz::dot;
    using blitz::TinyVector;
    auto ephemeral = pkd->Ephemeral<distance>();
    local_t n = pkd->Local();
    double r0, r, r2;
    int i, iBin;
    PROFILEBIN* pBin;

    if (pkd->Self() == 0) {
        if (pkd->profileBins != NULL) pkd->mdl->Free(pkd->profileBins);
        pkd->profileBins = pkd->mdl->Malloc<PROFILEBIN>(nBins * sizeof(PROFILEBIN));
        assert(pkd->profileBins != NULL);
        r0 = 0.0;
        for (iBin = 0; iBin < nBins; iBin++) {
            pBin = &pkd->profileBins[iBin];
            r = dRadii[iBin];
            r2 = r * r;
            assert(r > r0);
            pBin->nParticles = 0;
            pBin->dMassInBin = 0.0;
            pBin->dRadius = r;
            pBin->dVolume = (4.0 / 3.0) * M_PI * (r * r2 - r0 * r0 * r0);
            pBin->L[0] = pBin->L[1] = pBin->L[2] = 0.0;
            pBin->vel_radial = pBin->vel_radial_sigma = 0.0;
            pBin->vel_tang_sigma = 0.0;
            r0 = r;
        }
    }
    pkd->mdl->CacheInitialize(CID_BIN, std::make_shared<HelperProfileBins1>(pkd->profileBins, pkd->Self() == 0 ? nBins : 0));

    /*
    ** Now we add all of the particles to the appropriate bin.  NOTE than both
    ** the bins, and the particles are already sorted by distance from the center.
    */
    r0 = 0.0;
    i = 0;
    for (iBin = 0; iBin < nBins; iBin++) {
        pBin = pkd->mdl->Acquire<PROFILEBIN>(CID_BIN, iBin, 0);
        r = pBin->dRadius;
        r = dRadii[iBin];
        r2 = r * r;
        assert(r > r0);

        while (ephemeral[i].d2 <= r2 && i < n) {
            auto P = pkd->particles[ephemeral[i].i];
            auto m = P.mass();
            auto v = P.velocity();
            const auto r = P.position();
            TinyVector<double, 3> delta_x, delta_v, ang_mom;
            double dx2, vel;

            pBin->dMassInBin += m;
            pBin->nParticles++;
            delta_x = r - com;
            delta_v = v - vcm;
            ang_mom = cross(delta_x, delta_v);
            pBin->L += m * ang_mom;
            dx2 = dot(delta_x, delta_x);
            if (dx2 != 0.0) vel = dot(delta_x, delta_v) / sqrt(dx2);
            else vel = sqrt(dot(delta_v, delta_v));
            pBin->vel_radial += m * vel;
            pBin->vel_radial_sigma += m * vel * vel;
            i++;
        }
        r0 = r;
        pkd->mdl->Release(CID_BIN, pBin);
    }
    pkd->mdl->FinishCache(CID_BIN);

    /* We need angular momentum to calculate tangental velocity sigma, so we reopen the cache */
    pkd->mdl->CacheInitialize(CID_BIN, std::make_shared<HelperProfileBins2>(pkd->profileBins, pkd->Self() == 0 ? nBins : 0));
    r0 = 0.0;
    i = 0;
    for (iBin = 0; iBin < nBins; iBin++) {
        pBin = pkd->mdl->Acquire<PROFILEBIN>(CID_BIN, iBin, 0);
        r = dRadii[iBin];
        r2 = r * r;
        assert(r > r0);
        while (ephemeral[i].d2 <= r2 && i < n) {
            auto P = pkd->particles[ephemeral[i].i];
            double m = P.mass();
            auto v = P.velocity();
            const auto r = P.position();
            double dx2;
            TinyVector<double, 3> delta_x, delta_v, vel_tang, vel_shell, vel_tang_pec;

            delta_x = r - com;
            delta_v = v - vcm;
            dx2 = dot(delta_x, delta_x);

            if (dx2 != 0.0) {
                vel_tang = delta_v - dot(delta_v, delta_x) / dx2 * delta_x;
                vel_shell = cross(pBin->L, delta_x);
                vel_tang_pec = vel_tang - 1.0 / dx2 * vel_shell;
                pBin->vel_tang_sigma += m * dot(vel_tang_pec, vel_tang_pec);
            }
            i++;
        }
        r0 = r;
        pkd->mdl->Release(CID_BIN, pBin);
    }
    pkd->mdl->FinishCache(CID_BIN);
}
