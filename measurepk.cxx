#include "pst.h"
#include "master.h"

#include "gridinfo.hpp"
using namespace gridinfo;

static void initPk(void *vpkd, void *g) {
    FFTW3(real) * r = (FFTW3(real) *)g;
    *r = 0.0;
    }
static void combPk(void *vpkd, void *g1, void *g2) {
    FFTW3(real) * r1 = (FFTW3(real) *)g1;
    FFTW3(real) * r2 = (FFTW3(real) *)g2;
    *r1 += *r2;
    }

static void assign_mass(PKD pkd, int iAssignment, double dTotalMass, double fDelta, MDLFFT fft, FFTW3(real) *fftData) {
    int nGrid = fft->rgrid->n1;
    double fftNormalize = 1.0 / (1.0*nGrid*nGrid*nGrid);
    mdlGridCoord first, last;
    int i, j;
    double r[3];

    mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&first,&last,1);
    for( i=first.i; i<last.i; ++i ) fftData[i] = 0.0;
    mdlCOcache(pkd->mdl,CID_PK,NULL,fftData,sizeof(FFTW3(real)),last.i,pkd,initPk,combPk);
    pkdAssignMass(pkd, ROOT, nGrid, fDelta, iAssignment);
    mdlFinishCache(pkd->mdl,CID_PK);

    for( i=first.i; i<last.i; ++i ) {
	assert(fftData[i] >= 0.0);
	}
    double dRhoMean = dTotalMass * fftNormalize;
    double diRhoMean = 1.0 / dRhoMean;

    /*printf( "Calculating density contrast\n" );*/
    for( i=first.i; i<last.i; ++i ) {
	fftData[i] = fftData[i]*diRhoMean - 1.0;
	}

    mdlFFT(pkd->mdl,fft,fftData);
    }

static double deconvolveWindow(int i,int nGrid,int iAssignment) {
    double win = M_PI * i / nGrid;
    if(win>0.1) win = win / sin(win);
    else win=1.0 / (1.0-win*win/6.0*(1.0-win*win/20.0*(1.0-win*win/76.0)));
    double ret = win;
    while(--iAssignment) ret *= win;
    return ret;
    }


extern "C"
void pkdMeasurePk(PKD pkd, double dTotalMass, int iAssignment, int bInterleave,
    int nGrid, int nBins, double *fK, double *fPower, uint64_t *nPower) {
    mdlGridCoord first, last, index;

    auto iNyquist = nGrid / 2;
    auto fft = pkd->fft = mdlFFTInitialize(pkd->mdl,nGrid,nGrid,nGrid,0,0);
    GridInfo G(pkd->mdl,fft);

    mdlGridCoordFirstLast(pkd->mdl,fft->rgrid,&first,&last,1);
    auto fftData1 = reinterpret_cast<FFTW3(real) *>(mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(real)),pkd->pLite));
    auto fftData2 = reinterpret_cast<FFTW3(real) *>(mdlSetArray(pkd->mdl,last.i,sizeof(FFTW3(real)),fftData1 + fft->rgrid->nLocal));

    assign_mass(pkd,iAssignment,dTotalMass,0.0f,fft,fftData1);
    if (bInterleave) {
	//assign_mass(pkd,iAssignment,dTotalMass,0.5f,fft,fftData2);
	}

    complex_array_t K1,K2;
    auto data1 = reinterpret_cast<real_t *>(mdlSetArray(pkd->mdl,0,0,pkd->pLite));
    auto data2 = data1 + fft->rgrid->nLocal;
    G.setupArray(data1,K1);
    G.setupArray(data2,K2);

    for( auto i=0; i<nBins; i++ ) {
	fK[i] = 0.0;
	fPower[i] = 0.0;
	nPower[i] = 0;
	}

    double win_j, win_k;
#ifdef LINEAR_PK
    double scale = nBins * 1.0 / iNyquist;
#else
    double scale = nBins * 1.0 / log(iNyquist+1);
#endif
    int jj, kk, i, j, k;
    i = j = k = -1;
    for( auto index=K1.begin(); index!=K1.end(); ++index ) {
    	auto pos = index.position();
	i = pos[0];
	if ( j != pos[1]) {
	    j = pos[1];
	    jj = j>iNyquist ? nGrid - j : j;
	    win_j = deconvolveWindow(jj,nGrid,iAssignment);
	    }
	if ( k != pos[2] ) {
	    k = pos[2];
	    kk = k>iNyquist ? nGrid - k : k;
            win_k = deconvolveWindow(kk,nGrid,iAssignment);
	    }
	double win = deconvolveWindow(i,nGrid,iAssignment) * win_k * win_j;
	auto ak = sqrt(i*i + jj*jj + kk*kk);
	auto ks = int(ak);
	if ( ks >= 1 && ks <= iNyquist ) {
#ifdef LINEAR_PK
	    ks = floor((ks-1.0) * scale);
#else
	    ks = floor(log(ks) * scale);
#endif
	    assert(ks>=0 && ks <nBins);
	    double delta2 = win*win*std::norm(*index);
	    fK[ks] += ak;
	    fPower[ks] += delta2;
	    nPower[ks] += 1;
	    }
	}
    mdlFFTFinish(pkd->mdl,fft);
    pkd->fft = NULL;
    }
