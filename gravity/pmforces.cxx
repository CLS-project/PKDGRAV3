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

#include <vector>
#include "pst.h"
#include "core/aweights.hpp"
#include "ic/whitenoise.hpp"
using namespace gridinfo;
using namespace blitz;

// A blitz++ friendly wrap function. returns "ik" given array index
// Range: (-iNyquist,iNyquist] where iNyquist = m/2
static float fwrap(float v,float m) {
    return v - (v > m*0.5 ? m : 0);
}
BZ_DEFINE_BINARY_FUNC(Fn_fwrap,fwrap)
BZ_DECLARE_ARRAY_ET_BINARY(fwrap,     Fn_fwrap)
BZ_DECLARE_ARRAY_ET_BINARY_SCALAR(fwrap,     Fn_fwrap, float)

/********************************************************************************\
*
* GRID Construction
*
\********************************************************************************/

void pkdGenerateLinGrid(PKD pkd, MDLFFT fft, double a, double a_next, double Lbox, int iSeed,
                        int bFixed, float fPhase, int bRho) {
    /* If bRho == 0, we generate the \delta field,
    ** otherwise we generate the \delta\rho field.
    ** If a == a_next, we generate the field at this a.
    ** Otherwize, we generate the weighted average of the
    ** field over the interval [a, a_next]. Note that
    ** averaging is only implemented for bRho == 1.
    */
    if (!bRho && a != a_next) {
        fprintf(stderr,
                "WARNING: In pkdGenerateLinGrid(): Averaging of \\delta fields not implemented\n");
        abort();
    }
    /* For the sake of performance we precompute the (possibly averaged)
    ** field at the |k|'s at which the perturbations are tabulated.
    */
    size_t size = pkd->csm->val.classData.perturbations.size_k;
    double *logk = (double *)malloc(sizeof(double)*size);
    double *field = (double *)malloc(sizeof(double)*size);
    size_t i;
    double k;
    gsl_interp_accel *acc;
    gsl_spline *spline;
    for (i = 0; i < size; i++) {
        k = pkd->csm->val.classData.perturbations.k[i];
        logk[i] = log(k);
        if (bRho)
            /* Use the \delta\rho field, possibly averaged */
            field[i] = csmDeltaRho_lin(pkd->csm, a, a_next, k);
        else
            /* Use the \delta field */
            field[i] = csmDelta_lin(pkd->csm, a, k);
        /* The cubic splining is much better without the analytic zeta */
        field[i] /= csmZeta(pkd->csm, k);
    }
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(spline, logk, field, size);
    /* Generate grid */
    double ix, iy, iz, i2;
    double iLbox = 2*M_PI/Lbox;

    GridInfo G(pkd->mdl,fft);
    complex_array_t K;
    G.setupArray((FFTW3(real) *)mdlSetArray(pkd->mdl,0,0,pkd->pLite),K);
    NoiseGenerator ng(iSeed,bFixed,fPhase);
    ng.FillNoise(K,fft->rgrid->n3);
    for ( auto index=K.begin(); index!=K.end(); ++index ) {
        auto pos = index.position();
        iz = fwrap(pos[2],fft->rgrid->n3); // Range: (-iNyquist,iNyquist]
        iy = fwrap(pos[1],fft->rgrid->n2);
        ix = fwrap(pos[0],fft->rgrid->n1);
        i2 = ix*ix + iy*iy + iz*iz;
        if (i2>0) {
            k = sqrt((double)i2)*iLbox;
            *index *= csmZeta(pkd->csm, k)*gsl_spline_eval(spline, log(k), acc);
        }
        else
            *index = 0.0;
    }
    /* Cleanup */
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
    free(logk);
    free(field);
}

static double deconvolveLinWindow(int i,int nGrid) {
    double win = M_PI * i / nGrid;
    if (win>0.1) win = win / sin(win);
    else win=1.0 / (1.0-win*win/6.0*(1.0-win*win/20.0*(1.0-win*win/76.0)));
#if defined(USE_NGP_LIN)
    return win;
#elif defined(USE_CIC_LIN)
    return win*win;
#elif defined(USE_TSC_LIN)
    return win*win*win;
#else
    return win*win*win*win;
#endif
}

/*
 * Green Function for the Laplacian operator in Fourier space */
static double green(int i, int jj, int kk, int nGrid) {
    double g = pow2(sin(M_PI*i/(1.0*nGrid)));
    g += pow2(sin(M_PI*jj/(1.0*nGrid)));
    g += pow2(sin(M_PI*kk/(1.0*nGrid)));
    g *= 4*nGrid*nGrid;
    if (g ==0.0)
        return 0.0;
    else
        return -1.0/g;
}

extern "C"
void pkdSetLinGrid(PKD pkd, double a0, double a, double a1, double dBSize, int nGrid, int iSeed,
                   int bFixed, float fPhase) {
    MDLFFT fft = pkd->fft;
    /* Grid coordinates in real space :      [0, nGrid].[0, nGrid].[0, nGrid] */
    mdlGridCoord rfirst, rlast;
    /* Grid coordinates in Fourier space : [O, Nyquist].[0, nGrid].[0, nGrid] */
    mdlGridCoord kfirst, klast, kindex;
    /*
     * Define the grid arrays : only 3 grids are stored
     * in memory, the other ones are defined in order to
     * have an explicit naming in the code
     */
    FFTW3(real) *rForceX, *rForceY, *rForceZ;
    FFTW3(complex) *cDelta_lin_field, *cForceY, *cForceZ;

    /* Scale factors and normalization */
    const double dNormalization = a*a*a * dBSize;

    mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&rfirst,&rlast,1);
    mdlGridCoordFirstLast(pkd->mdl,fft->kgrid,&kfirst,&klast,0);

    /* Imprint the density grid of the linear species */
    int bRho = 1;  /* Generate the \delta\rho field */
    pkdGenerateLinGrid(pkd, fft, a0, a1, dBSize, iSeed, bFixed, fPhase, bRho);
    cDelta_lin_field = (FFTW3(complex) *)mdlSetArray(pkd->mdl, klast.i, sizeof(FFTW3(complex)), pkd->pLite);

    /* Remember, the grid is now transposed to x,z,y (from x,y,z) */
    cForceY = (FFTW3(complex) *)mdlSetArray(pkd->mdl,klast.i,sizeof(FFTW3(complex)),cDelta_lin_field + fft->kgrid->nLocal);
    cForceZ = (FFTW3(complex) *)mdlSetArray(pkd->mdl,klast.i, sizeof(FFTW3(complex)),cForceY + fft->kgrid->nLocal);

    int idx, i, j, jj, k, kk;
    const int iNyquist = nGrid / 2 ;
    double rePotential, imPotential;
    double dDifferentiate, dPoissonSolve;
    double win_j, win_k;
    /* Here starts the Poisson solver */
    i = j = k = -1;
    for ( kindex=kfirst; !mdlGridCoordCompare(&kindex,&klast); mdlGridCoordIncrement(&kindex) ) {
        idx = kindex.i;
        if ( j != kindex.z ) {
            j = kindex.z;
            jj = j>iNyquist ? j - nGrid : j;
            win_j = deconvolveLinWindow(jj,nGrid);
        }
        if ( k != kindex.y ) {
            k = kindex.y;
            kk = k>iNyquist ? k - nGrid : k;
            win_k = deconvolveLinWindow(kk,nGrid);
        }
        i = kindex.x;
        double win = deconvolveLinWindow(i,nGrid)*win_j*win_k;
        /* Green Function for a discrete Laplacian operator */
        dPoissonSolve=4*M_PI*green(i,jj,kk,nGrid)*dNormalization*win;
        /* Solve Poisson equation */

        rePotential = cDelta_lin_field[idx][0] * dPoissonSolve;
        imPotential = cDelta_lin_field[idx][1] * dPoissonSolve;
        /* Differentiaite in Y direction */
        dDifferentiate = nGrid*sin(2*M_PI*jj/(1.0*nGrid));
        cForceY[idx][0] =  dDifferentiate * imPotential;
        cForceY[idx][1] = -dDifferentiate * rePotential;

        /* Differentiate in Z direction */
        dDifferentiate = nGrid*sin(2*M_PI*kk/(1.0*nGrid));
        cForceZ[idx][0] =  dDifferentiate * imPotential;
        cForceZ[idx][1] = -dDifferentiate * rePotential;

        /*
         * Differentiate in X direction (over-write the
         * delta_lin field)
         */
        dDifferentiate = nGrid*sin(2*M_PI*i/(1.0*nGrid));
        cDelta_lin_field[idx][0] =  dDifferentiate * imPotential;
        cDelta_lin_field[idx][1] = -dDifferentiate * rePotential;
    }
    mdlIFFT(pkd->mdl, fft, cForceY);
    rForceY = (FFTW3(real) *)mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),cForceY);

    mdlIFFT(pkd->mdl, fft, cForceZ);
    rForceZ = (FFTW3(real) *)mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)),cForceZ);

    mdlIFFT(pkd->mdl, fft, cDelta_lin_field);
    rForceX = (FFTW3(real) *)mdlSetArray(pkd->mdl,rlast.i,sizeof(FFTW3(real)), cDelta_lin_field);
}

extern "C"
int pstSetLinGrid(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inSetLinGrid *in = reinterpret_cast<struct inSetLinGrid *>(vin);
    assert (nIn==sizeof(struct inSetLinGrid) );
    if (pstNotCore(pst)) {
        int rID = mdlReqService(pst->mdl, pst->idUpper, PST_SETLINGRID, vin, nIn);
        pstSetLinGrid(pst->pstLower, vin, nIn, NULL, 0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdSetLinGrid(plcl->pkd, in->a0, in->a, in->a1,
                      in->dBSize, in->nGrid,
                      in ->iSeed, in->bFixed, in->fPhase);
    }
    return 0;
}

typedef blitz::Array<float,3> force_array_t;
typedef blitz::TinyVector<int,3> shape_t;
typedef blitz::TinyVector<double,3> position_t;
typedef blitz::TinyVector<float,3> float3_t;

struct tree_node : public KDN {
    bool is_cell()   { return iLower!=0; }
    bool is_bucket() { return iLower==0; }
};

template<int Order,typename F>
static float interpolate(force_array_t &forces, const F r[3]) {
    AssignmentWeights<Order,F> Hx(r[0]),Hy(r[1]),Hz(r[2]);
    float force = 0;
    for (int i=0; i<Order; ++i) {
        for (int j=0; j<Order; ++j) {
            for (int k=0; k<Order; ++k) {
                force += forces(Hx.i+i,Hy.i+j,Hz.i+k) * Hx.H[i]*Hy.H[j]*Hz.H[k];
            }
        }
    }
    return force;
}

template<typename F>
static float force_interpolate(force_array_t &forces, const F r[3],int iAssignment=4) {
    float f;
    switch (iAssignment) {
    case 1: f=interpolate<1,F>(forces,r); break;
    case 2: f=interpolate<2,F>(forces,r); break;
    case 3: f=interpolate<3,F>(forces,r); break;
    case 4: f=interpolate<4,F>(forces,r); break;
    default: f=0; assert(iAssignment>=1 && iAssignment<=4); abort();
    }
    return f;
}

// Fetch all forces into this local grid. Could be optimized by bulk fetches.
static void fetch_forces(PKD pkd,int cid,int nGrid,force_array_t &forces, const shape_t &lower) {
    auto wrap = [&nGrid](int i) { if (i>=nGrid) i-=nGrid; else if (i<0) i+=nGrid; return i; };
    for (auto i=forces.begin(); i!=forces.end(); ++i) {
        shape_t loc = i.position() + lower;
        loc[0] = wrap(loc[0]); loc[1] = wrap(loc[1]); loc[2] = wrap(loc[2]);
        auto id = mdlFFTrId(pkd->mdl,pkd->fft,loc[0],loc[1],loc[2]);
        auto idx = mdlFFTrIdx(pkd->mdl,pkd->fft,loc[0],loc[1],loc[2]);
        auto p = reinterpret_cast<float *>(mdlFetch(pkd->mdl,cid,idx,id));
        *i = *p;
    }
}

void pkdLinearKick(PKD pkd,vel_t dtOpen,vel_t dtClose, int iAssignment=4) {
    const std::size_t maxSize = 100000; // We would like this to remain in L2 cache
    std::vector<float> dataX, dataY, dataZ;
    dataX.reserve(maxSize);
    dataY.reserve(maxSize);
    dataZ.reserve(maxSize);
    shape_t index;
    position_t fPeriod(pkd->fPeriod), ifPeriod = 1.0 / fPeriod;
    int nGrid = pkd->fft->rgrid->n1;
    assert(iAssignment>=1 && iAssignment<=4);

    int iLocal = mdlCore(pkd->mdl) ? 0 : pkd->fft->rgrid->nLocal;
    FFTW3(real)* forceX = reinterpret_cast<FFTW3(real) *>(mdlSetArray(pkd->mdl,iLocal,sizeof(FFTW3(real)),pkd->pLite));
    FFTW3(real)* forceY = reinterpret_cast<FFTW3(real) *>(mdlSetArray(pkd->mdl,iLocal,sizeof(FFTW3(real)),forceX + pkd->fft->rgrid->nLocal));
    FFTW3(real)* forceZ = reinterpret_cast<FFTW3(real) *>(mdlSetArray(pkd->mdl,iLocal,sizeof(FFTW3(real)),forceY + pkd->fft->rgrid->nLocal));
    mdlROcache(pkd->mdl,CID_GridLinFx,NULL, forceX, sizeof(FFTW3(real)),iLocal);
    mdlROcache(pkd->mdl,CID_GridLinFy,NULL, forceY, sizeof(FFTW3(real)),iLocal );
    mdlROcache(pkd->mdl,CID_GridLinFz,NULL, forceZ, sizeof(FFTW3(real)),iLocal );

    std::vector<std::uint32_t> stack;
    stack.push_back(ROOT);
    while ( !stack.empty()) {
        tree_node *kdn = reinterpret_cast<tree_node *>(pkdTreeNode(pkd,stack.back()));
        stack.pop_back(); // Go to the next node in the tree
        Bound bnd = pkdNodeGetBnd(pkd, kdn);
        shape_t ilower = shape_t(floor((bnd.lower() * ifPeriod + 0.5) * nGrid)) - iAssignment/2;
        shape_t iupper = shape_t(floor((bnd.upper() * ifPeriod + 0.5) * nGrid)) + iAssignment/2;
        shape_t ishape = iupper - ilower + 1;
        float3_t flower = ilower;
        std::size_t size = blitz::product(ishape);

        if (size > maxSize) { // This cell is too large, so we split it and move on
            assert(kdn->is_cell()); // At the moment we cannot handle enormous buckets
            stack.push_back(kdn->iLower+1);
            stack.push_back(kdn->iLower);
        }
        else { // Assign the mass for this range of particles
            dataX.resize(size); // Hold the right number of masses
            dataY.resize(size); // Hold the right number of masses
            dataZ.resize(size); // Hold the right number of masses
            force_array_t forcesX(dataX.data(),ishape,blitz::neverDeleteData,blitz::ColumnMajorArray<3>());
            force_array_t forcesY(dataY.data(),ishape,blitz::neverDeleteData,blitz::ColumnMajorArray<3>());
            force_array_t forcesZ(dataZ.data(),ishape,blitz::neverDeleteData,blitz::ColumnMajorArray<3>());
            fetch_forces(pkd,CID_GridLinFx,nGrid,forcesX,ilower);
            fetch_forces(pkd,CID_GridLinFy,nGrid,forcesY,ilower);
            fetch_forces(pkd,CID_GridLinFz,nGrid,forcesZ,ilower);
            for ( int i=kdn->pLower; i<=kdn->pUpper; ++i) { // All particles in this tree cell
                auto p = pkdParticle(pkd,i);
                auto v = pkdVel(pkd,p);
                position_t dr; pkdGetPos1(pkd,p,dr.data()); // Centered on 0 with period fPeriod
                float3_t r(dr);
                r = (r * ifPeriod + 0.5) * nGrid - flower; // Scale and shift to fit in subcube
                v[0] += (dtOpen + dtClose) * force_interpolate(forcesX, r.data(), iAssignment);
                v[1] += (dtOpen + dtClose) * force_interpolate(forcesY, r.data(), iAssignment);
                v[2] += (dtOpen + dtClose) * force_interpolate(forcesZ, r.data(), iAssignment);
            }
        }
    }
    mdlFinishCache(pkd->mdl,CID_GridLinFx);
    mdlFinishCache(pkd->mdl,CID_GridLinFy);
    mdlFinishCache(pkd->mdl,CID_GridLinFz);
}

extern "C"
int pstLinearKick(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inLinearKick *in = reinterpret_cast<struct inLinearKick *>(vin);
    assert( nIn==sizeof(struct inLinearKick) );

    if (pstNotCore(pst)) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LINEARKICK,vin,nIn);
        pstLinearKick(pst->pstLower,vin,nIn,NULL,0);
        mdlGetReply(pst->mdl,rID,NULL,NULL);
    }
    else {
        pkdLinearKick(plcl->pkd,in->dtOpen,in->dtClose);
    }
    return 0;
}

void pkdMeasureLinPk(PKD pkd, int nGrid, double dA, double dBoxSize,
                     int nBins,  int iSeed, int bFixed, float fPhase,
                     double *fK, double *fPower, uint64_t *nPower) {
    MDLFFT fft = pkd->fft;
    mdlGridCoord first, last, index;
    FFTW3(complex) *fftDataK;
    double ak;
    int i,j,k, idx, ks;
    int iNyquist;

    /* Sort the particles into optimal "cell" order */
    /* Use tree order: QSORT(pkdParticleSize(pkd),pkdParticle(pkd,0),pkd->nLocal,qsort_lt); */

    iNyquist = nGrid / 2;

    mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&first,&last,1);
    /* Generate the grid of the linear species again,
    ** this should not be done this way for performances.
    ** - Well, for the gravity computation what is generated
    **   on the grid is \delta\rho averaged over one time step.
    **   Here we want the \delta power spectrum at a specific
    **   time, and so we need to generate the \delta field at this
    **   time. We thus have to re-generate the grid.
    */
    int bRho = 0; double a_next = dA; /* Generate the \delta field at dA (no averaging) */
    pkdGenerateLinGrid(pkd, fft, dA, a_next, dBoxSize, iSeed, bFixed, fPhase, bRho);

    /* Remember, the grid is now transposed to x,z,y (from x,y,z) */
    mdlGridCoordFirstLast(pkd->mdl,fft->kgrid,&first,&last,0);
    fftDataK = (FFTW3(complex) *)mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(complex)),pkd->pLite);

    for ( i=0; i<nBins; i++ ) {
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
    double scale = nBins * 1.0 / log(iNyquist+1);
#endif
    double dBox2 = dBoxSize * dBoxSize;
    int jj, kk;
    i = j = k = -1;
    for ( index=first; !mdlGridCoordCompare(&index,&last); mdlGridCoordIncrement(&index) ) {
        if ( j != index.z ) {
            j = index.z;
            jj = j>iNyquist ? nGrid - j : j;
        }
        if ( k != index.y ) {
            k = index.y;
            kk = k>iNyquist ? nGrid - k : k;
        }
        i = index.x;
        ak = sqrt(i*i + jj*jj + kk*kk);
        ks = ak;
        if ( ks >= 1 && ks <= iNyquist ) {
#ifdef LINEAR_PK
            ks = floor((ks-1.0) * scale);
#else
            ks = floor(log(ks) * scale);
#endif
            assert(ks>=0 && ks <nBins);
            idx = index.i;
            double delta2 = dBox2*(pow2(fftDataK[idx][0]) + pow2(fftDataK[idx][1]));
            fK[ks] += ak;
            fPower[ks] += delta2;
            nPower[ks] += 1;
        }
    }
}

extern "C"
int pstMeasureLinPk(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inMeasureLinPk *in = reinterpret_cast<struct inMeasureLinPk *>(vin);
    struct outMeasureLinPk *out = reinterpret_cast<struct outMeasureLinPk *>(vout);
    struct outMeasureLinPk *outUpper;
    int i;

    assert( nIn==sizeof(struct inMeasureLinPk) );
    if (pstNotCore(pst)) {
        int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEASURELINPK,vin,nIn);
        pstMeasureLinPk(pst->pstLower,vin,nIn,vout,nOut);
        outUpper = reinterpret_cast<struct outMeasureLinPk *>(malloc(sizeof(struct outMeasureLinPk)));
        assert(outUpper != NULL);
        mdlGetReply(pst->mdl,rID,outUpper,&nOut);
        assert(nOut==sizeof(struct outMeasureLinPk));

        for (i=0; i<in->nBins; i++) {
            out->fK[i] += outUpper->fK[i];
            out->fPower[i] += outUpper->fPower[i];
            out->nPower[i] += outUpper->nPower[i];
        }
        free(outUpper);
    }
    else {
        pkdMeasureLinPk(plcl->pkd, in->nGrid, in->dA, in->dBoxSize,
                        in->nBins, in->iSeed, in->bFixed, in->fPhase,
                        out->fK, out->fPower, out->nPower);
    }
    return sizeof(struct outMeasureLinPk);
}
