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

#include "pst.h"
#include "whitenoise.hpp"
using namespace gridinfo;
using namespace blitz;

static const std::complex<float> I(0,1);

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
    int bFixed, float fPhase, int bRho){
    /* If bRho == 0, we generate the \delta field,
    ** otherwise we generate the \delta\rho field.
    ** If a == a_next, we generate the field at this a.
    ** Otherwize, we generate the weighted average of the
    ** field over the interval [a, a_next]. Note that
    ** averaging is only implemented for bRho == 1.
    */
    if (!bRho && a != a_next){
        fprintf(stderr,
            "WARNING: In pkdGenerateLinGrid(): Averaging of \\delta fields not implemented\n");
        abort();
    }
    /* For the sake of performance we precompute the (possibly averaged)
    ** field at the |k|'s at which the perturbations are tabulated.
    */
    size_t size = pkd->param.csm->val.classData.perturbations.size_k;
    double *logk = (double*)malloc(sizeof(double)*size);
    double *field = (double*)malloc(sizeof(double)*size);
    size_t i;
    double k;
    gsl_interp_accel *acc;
    gsl_spline *spline;
    for (i = 0; i < size; i++){
        k = pkd->param.csm->val.classData.perturbations.k[i];
        logk[i] = log(k);
        if (bRho)
            /* Use the \delta\rho field, possibly averaged */
            field[i] = csmDeltaRho_lin(pkd->param.csm, a, a_next, k);
        else
            /* Use the \delta field */
            field[i] = csmDelta_lin(pkd->param.csm, a, k);
        /* The cubic splining is much better without the analytic zeta */
        field[i] /= csmZeta(pkd->param.csm, k);
    }
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(spline, logk, field, size);
    /* Generate grid */
    double noiseMean, noiseCSQ;
    int idx;
    double ix, iy, iz, i2;
    double iLbox = 2*M_PI/Lbox;
    int iNuquist = fft->rgrid->n3 / 2;

    GridInfo G(pkd->mdl,fft);
    complex_array_t K;
    G.setupArray((FFTW3(real) *)mdlSetArray(pkd->mdl,0,0,pkd->pLite),K);
    pkdGenerateNoise(pkd, iSeed, bFixed, fPhase, fft, K, &noiseMean, &noiseCSQ);
    for( auto index=K.begin(); index!=K.end(); ++index ) {
	auto pos = index.position();
	iz = fwrap(pos[2],fft->rgrid->n3); // Range: (-iNyquist,iNyquist]
	iy = fwrap(pos[1],fft->rgrid->n2);
	ix = fwrap(pos[0],fft->rgrid->n1);
        i2 = ix*ix + iy*iy + iz*iz;
        if (i2>0){
            k = sqrt((double)i2)*iLbox;
            *index *= csmZeta(pkd->param.csm, k)*gsl_spline_eval(spline, log(k), acc);
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
    if(win>0.1) win = win / sin(win);
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
static double green(int i, int jj, int kk, int nGrid){
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
void pkdSetLinGrid(PKD pkd, double dTime, double dTime_next, double dBSize, int nGrid, int iSeed,
    int bFixed, float fPhase) {
        MDLFFT fft = pkd->Linfft;
        /* Grid coordinates in real space :      [0, nGrid].[0, nGrid].[0, nGrid] */
        mdlGridCoord rfirst, rlast, rindex;
        /* Grid coordinates in Fourier space : [O, Nyquist].[0, nGrid].[0, nGrid] */
        mdlGridCoord kfirst, klast, kindex;
        /* 
         * Define the grid arrays : only 3 grids are stored
         * in memory, the other ones are defined in order to
         * have an explicit naming in the code
         */
        FFTW3(real) *rForceX, *rForceY, *rForceZ;
        FFTW3(complex) *cDelta_lin_field, *cForceY, *cForceZ;
#ifdef USE_ITT
        __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
        __itt_string_handle* shMyTask = __itt_string_handle_create("AssignMass_DGrid");
        __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
        mdlThreadBarrier(pkd->mdl);
        __itt_resume();
#endif
        /* Scale factors and normalization */
        const double a      = csmTime2Exp(pkd->param.csm, dTime);
        const double a_next = csmTime2Exp(pkd->param.csm, dTime_next);
        const double dNormalization = a*a*a * dBSize;

        mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&rfirst,&rlast,1);
        mdlGridCoordFirstLast(pkd->mdl,fft->kgrid,&kfirst,&klast,0);

        /* Imprint the density grid of the linear species */
        int bRho = 1;  /* Generate the \delta\rho field */
        pkdGenerateLinGrid(pkd, fft, a, a_next, dBSize, iSeed, bFixed, fPhase, bRho);
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
        for( kindex=kfirst; !mdlGridCoordCompare(&kindex,&klast); mdlGridCoordIncrement(&kindex) ) {
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


#ifdef USE_ITT
        mdlThreadBarrier(pkd->mdl);
        __itt_pause();
        __itt_task_end(domain);
        mdlThreadBarrier(pkd->mdl);
#endif
}

extern "C"
void pstSetLinGrid(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
        LCL *plcl = pst->plcl;
        struct inSetLinGrid *in = reinterpret_cast<struct inSetLinGrid *>(vin);
        assert (nIn==sizeof(struct inSetLinGrid) );
        if (pstNotCore(pst)) {
            int rID = mdlReqService(pst->mdl, pst->idUpper, PST_SETLINGRID, vin, nIn);
            pstSetLinGrid(pst->pstLower, vin, nIn, vout, pnOut);
            mdlGetReply(pst->mdl,rID, vout,pnOut);
        }
        else {
            plcl->pkd->Linfft = mdlFFTInitialize(pst->mdl, in->nGrid, in->nGrid, in->nGrid, 0,0);
            pkdSetLinGrid(plcl->pkd, in->dTime, in->dTime_next,
                in->dBSize, in->nGrid, 
                in ->iSeed, in->bFixed, in->fPhase);
        }
    }

/********************************************************************************\
*
* Force Interpolation
*
\********************************************************************************/

static int wrap(int i,int w) {
    if (i>=w) i -= w;
    else if (i<0) i += w;
    return i;
    }

static inline float pow2(float x) {
    return x*x;
    }

static void cic_weights(int ii[3][2],float H[3][2],const double r[3], int nGrid) {
    int d;
    float rr, h;
    for(d=0; d<3; ++d) {
	rr = r[d] * (float)(nGrid);                 /* coordinates in subcube units [0,NGRID] */
	ii[d][0]  = (int)(rr);                      /* index of nearest grid point [0,NGRID] */
	if (ii[d][0]==nGrid) ii[d][0] = nGrid-1;    /* If very close to 1.0, it could round up, so correct */
	h = rr - (float)ii[d][0];             /* distance to nearest grid point */
	ii[d][1]=wrap(ii[d][0]+1,nGrid);            /* keep track of periodic boundaries */
	H[d][0] = 1.0 - h;              /* calculate CIC weights */
	H[d][1] = h;
	}
    }

static void tsc_weights(int ii[3][3],float H[3][3],const double r[3], int nGrid) {
    int d;
    float rr, h;
    for(d=0; d<3; ++d) {
	rr = r[d] * (float)(nGrid);                 /* coordinates in subcube units [0,NGRID] */
	ii[d][1]  = (int)(rr);                      /* index of nearest grid point [0,NGRID] */
	if (ii[d][1]==nGrid) ii[d][1] = nGrid-1;    /* If very close to 1.0, it could round up, so correct */
	h = (rr-0.5) - (float)ii[d][1];             /* distance to nearest grid point */
	ii[d][2]=wrap(ii[d][1]+1,nGrid);            /* keep track of periodic boundaries */
	ii[d][0]=wrap(ii[d][1]-1,nGrid);
	H[d][0] = 0.5 * pow2(0.5 - h);              /* calculate TSC weights */
	H[d][1] = 0.75 - pow2(h);
	H[d][2] = 0.5 * pow2(0.5 + h);
	}
    }

static inline float pow3(float x) {
    return x*x*x;
    }

static void pcs_weights(int ii[3][4],float H[3][4],const double r[3], int nGrid) {
    int d,i;
    float rr, h;
    for(d=0; d<3; ++d) {
	rr = r[d] * (float)(nGrid);                 /* coordinates in subcube units [0,NGRID] */
	int g = (int)(rr);                          /* index of nearest grid point [0,NGRID] */
	if (g==nGrid) g = nGrid-1;                  /* If very close to 1.0, it could round up, so correct */
	h = (rr-0.5) - (float)g;                    /* distance to nearest grid point */
	int b = h > 0.0 ? -1 : -2;                  /* the kernel is 4x4x4, so choose the correct start cell */
	for(i=0; i<4; ++i) {
	    float s = fabs(i + b - h );
	    ii[d][i] = wrap(i + b + g,nGrid);       /* keep track of periodic boundaries */
	    if ( s < 1.0f ) H[d][i] = 1.0f/6.0f * ( 4.0f - 6.0f*s*s + 3.0f*s*s*s);
	    else if ( s < 2.0f ) H[d][i] = 1.0f/6.0f * pow3(2.0f - s);
	    else H[d][i] = 0.0f;
	    }
	}
    }

static void force_accumulate(PKD pkd, MDLFFT fft, int cid, int x, int y, int z, float *force, float w){
        /* If the weight is non zero */
        if (w > 0){
                /* 
                 * Find the value on the force grid (which
                 * corresponds to cid) at position (x, y, z)
                 */ 
                FFTW3(real)* p;
                int id, idx;
                id = mdlFFTrId(pkd->mdl, fft, x, y, z);
                idx = mdlFFTrIdx(pkd->mdl, fft, x, y, z);
                p = (FFTW3(real)*)mdlFetch(pkd->mdl, cid, idx, id);
                /* Accumulate the value in *force */
                *force += *p * w;
        }
}

#if defined(USE_NGP_LIN)
static void ngp_addForce(PKD pkd, MDLFFT fft,int cid, int nGrid,
                double x, double y, double z, float* force){
        int ix, iy, iz;
        ix = (int)(x * nGrid);
        iy = (int)(y * nGrid);
        iz = (int)(z * nGrid);
        if (ix==nGrid) ix = nGrid-1;
        if (iy==nGrid) iy = nGrid-1;
        if (iz==nGrid) iz = nGrid-1;
        force_accumulate(pkd,fft,cid,ix,iy,iz,force, 1.0f);
}
#elif defined(USE_CIC_LIN)
static void cic_addForce(PKD pkd, MDLFFT fft,int cid, int nGrid,
                double x, double y, double z, float* force) {
        double r[] = {x,y,z};
        int    ii[3][2];
        float  H[3][2];
        int i,j,k;

        int           ix, iy, iz, ixp1, iyp1, izp1;
        float         rrx, rry, rrz;
        float         hx, hy, hz;
        float         hx0, hy0, hz0, hxp1, hyp1, hzp1;

        cic_weights(ii,H,r,nGrid);

        /* assign particle according to weights to 8 neighboring nodes */
        for(i=0; i<2; ++i) {
                for(j=0; j<2; ++j) {
                        for(k=0; k<2; ++k) {
                                force_accumulate(pkd,fft,cid,ii[0][i],ii[1][j],ii[2][k],force, H[0][i]*H[1][j]*H[2][k]);
                        }
                }
        }
}

#elif defined(USE_TSC_LIN)
static void tsc_addForce(PKD pkd, MDLFFT fft,int cid, int nGrid,
                double x, double y, double z, float* force) {
        double r[] = {x,y,z};
        int    ii[3][3];
        float  H[3][3];
        int    i,j,k;
        tsc_weights(ii,H,r,nGrid);

        /* assign particle according to weights to 27 neighboring nodes */
        for(i=0; i<3; ++i) {
                for(j=0; j<3; ++j) {
                        for(k=0; k<3; ++k) {
                                force_accumulate(pkd,fft,cid,ii[0][i],ii[1][j],ii[2][k],force,H[0][i]*H[1][j]*H[2][k]);
                        }
                }
        }
}
#else
static void pcs_addForce(PKD pkd, MDLFFT fft, int cid, int nGrid,
                double x, double y, double z, float* force) {
        double r[] = {x,y,z};
        int    ii[3][4];
        float  H[3][4];
        int    i,j,k;

        pcs_weights(ii,H,r,nGrid);

        /* assign particle according to weights to 64 neighboring nodes */
        for(i=0; i<4; ++i) {
                for(j=0; j<4; ++j) {
                        for(k=0; k<4; ++k) {
                                force_accumulate(pkd,fft,cid,ii[0][i],ii[1][j],ii[2][k],force,H[0][i]*H[1][j]*H[2][k]);
                        }
                }
        }
}
#endif

float getLinAcc(PKD pkd, MDLFFT fft,int cid, double r[3]){
    float force = 0;
    int nGrid = fft->rgrid->n1;
    int i;
    /* Recenter, apply periodic boundary and scale to the correct size */
    mdlGridCoord first, last;
    mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&first,&last,1);
    /* Shift particle position by 0.5*fPeriod */
    double r_c[3];
    for(i=0;i<3;i++) {
        r_c[i] = 0.5*pkd->fPeriod[i] + r[i];
        if (r_c[i]>=pkd->fPeriod[i]) r_c[i] -= pkd->fPeriod[i];
        if (r_c[i]< 0 ) r_c[i] += pkd->fPeriod[i];
    }
    assert( r_c[0]>=0.0 && r_c[0]<1.0 && r_c[1]>=0.0 && r_c[1]<1.0 && r_c[2]>=0.0 && r_c[2]<1.0 );
#if defined(USE_NGP_LIN)
    ngp_addForce(pkd, fft, cid, nGrid, r_c[0], r_c[1], r_c[2], &force);
#elif defined(USE_CIC_LIN)
    cic_addForce(pkd, fft, cid, nGrid, r_c[0], r_c[1], r_c[2], &force);
#elif defined(USE_TSC_LIN)
    tsc_addForce(pkd, fft, cid, nGrid, r_c[0], r_c[1], r_c[2], &force);
#else
    pcs_addForce(pkd, fft, cid, nGrid, r_c[0], r_c[1], r_c[2], &force);
#endif
    return force;
}

void pkdLinearKick(PKD pkd, vel_t dtOpen, vel_t dtClose) {
    mdlGridCoord  first, last;
    mdlGridCoordFirstLast(pkd->mdl,pkd->Linfft->rgrid,&first,&last,1);
    FFTW3(real)* forceX = reinterpret_cast<FFTW3(real)*>(mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(real)),pkd->pLite));
    FFTW3(real)* forceY = reinterpret_cast<FFTW3(real)*>(mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(real)),forceX + pkd->Linfft->rgrid->nLocal));
    FFTW3(real)* forceZ = reinterpret_cast<FFTW3(real)*>(mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(real)),forceY + pkd->Linfft->rgrid->nLocal));
    mdlROcache(pkd->mdl,CID_GridLinFx,NULL, forceX, sizeof(FFTW3(real)),last.i );
    mdlROcache(pkd->mdl,CID_GridLinFy,NULL, forceY, sizeof(FFTW3(real)),last.i );
    mdlROcache(pkd->mdl,CID_GridLinFz,NULL, forceZ, sizeof(FFTW3(real)),last.i );
    for ( auto i=0; i<pkd->nLocal; ++i) {
        auto p = pkdParticle(pkd,i);
        auto v = pkdVel(pkd,p);
        double r[3];
        float a[3];
        pkdGetPos1(pkd,p,r);
        a[0] = getLinAcc(pkd,pkd->Linfft,CID_GridLinFx,r);
        a[1] = getLinAcc(pkd,pkd->Linfft,CID_GridLinFy,r);
        a[2] = getLinAcc(pkd,pkd->Linfft,CID_GridLinFz,r);
	v[0] += dtOpen*a[0] + dtClose*a[0];
	v[1] += dtOpen*a[1] + dtClose*a[1];
	v[2] += dtOpen*a[2] + dtClose*a[2];
        }
    mdlFinishCache(pkd->mdl,CID_GridLinFx);
    mdlFinishCache(pkd->mdl,CID_GridLinFy);
    mdlFinishCache(pkd->mdl,CID_GridLinFz);
    }

extern "C"
void pstLinearKick(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inLinearKick *in = reinterpret_cast<struct inLinearKick *>(vin);
    assert( nIn==sizeof(struct inLinearKick) );

    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_LINEARKICK,vin,nIn);
	pstLinearKick(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,rID,vout,pnOut);
	}
    else {
	pkdLinearKick(plcl->pkd,in->dtOpen,in->dtClose);
	}
    if (pnOut) *pnOut = sizeof(struct outMeasureLinPk);
    }

void pkdMeasureLinPk(PKD pkd, int nGrid, double dA, double dBoxSize,
                int nBins,  int iSeed, int bFixed, float fPhase, 
                double *fK, double *fPower, uint64_t *nPower) {
    MDLFFT fft = pkd->Linfft;
    mdlGridCoord first, last, index;
    FFTW3(complex) *fftDataK;
    double ak;
    int i,j,k, idx, ks;
    int iNyquist;
#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("MeasureLinPk");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
    mdlThreadBarrier(pkd->mdl);
    __itt_resume();
#endif

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

    for( i=0; i<nBins; i++ ) {
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
    for( index=first; !mdlGridCoordCompare(&index,&last); mdlGridCoordIncrement(&index) ) {
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
#ifdef USE_ITT
    mdlThreadBarrier(pkd->mdl);
    __itt_pause();
    __itt_task_end(domain);
    mdlThreadBarrier(pkd->mdl);
#endif
    }

extern "C"
void pstMeasureLinPk(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inMeasureLinPk *in = reinterpret_cast<struct inMeasureLinPk *>(vin);
    struct outMeasureLinPk *out = reinterpret_cast<struct outMeasureLinPk *>(vout);
    struct outMeasureLinPk *outUpper;
    int nOut;
    int i;

    assert( nIn==sizeof(struct inMeasureLinPk) );
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl,pst->idUpper,PST_MEASURELINPK,vin,nIn);
	pstMeasureLinPk(pst->pstLower,vin,nIn,vout,pnOut);
	outUpper = reinterpret_cast<struct outMeasureLinPk *>(malloc(sizeof(struct outMeasureLinPk)));
	assert(outUpper != NULL);
	mdlGetReply(pst->mdl,rID,outUpper,&nOut);
	assert(nOut==sizeof(struct outMeasureLinPk));

	for(i=0;i<in->nBins; i++) {
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
    if (pnOut) *pnOut = sizeof(struct outMeasureLinPk);
    }
