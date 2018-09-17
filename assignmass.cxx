#include <cstddef>
#include <cstdint>
#include <vector>
#include <numeric>
#include <iostream>
#include "blitz/array.h"
#include "pst.h"
#include "master.h"

typedef blitz::Array<float,3> mass_array_t;
typedef blitz::TinyVector<int,3> shape_t;
typedef blitz::TinyVector<double,3> position_t;
typedef blitz::TinyVector<float,3> float3_t;
typedef blitz::ColumnMajorArray<3> RegularArray;

struct tree_node : public KDN {
    bool is_cell()   { return iLower!=0; }
    bool is_bucket() { return iLower==0; }
    };

template<typename F>
static void ngp_weights(int ii[3],F H[][3],const F r[3]) {
    for( int d=0; d<3; ++d) {
    	ii[d] = r[d];
    	H[0][d] = 1.0;
	}
    }

template<typename F>
static void cic_weights(int ii[3],F H[][3],const F r[3]) {
    for(int d=0; d<3; ++d) {
	F rr = r[d] - 0.5;              /* Center on grid boundary */
	ii[d] = floorf(rr);                 /* index of first grid point [0,size] */
	F h = rr - ii[d];  	    /* distance to first grid point */
	H[0][d] = 1.0 - h;           	    /* calculate CIC weights */
	H[1][d] = h;
	}
    }

template<typename F>
static void tsc_weights(int ii[3],F H[][3],const F r[3]) {
    auto K0 = [](F h) { return 0.75 - h*h; };
    auto K1 = [](F h) { return 0.50 * h*h; };
    for(int d=0; d<3; ++d) {
	F rr = r[d];		    /* */
	int i = rr;      		    /* Middle cell (of three) */
	ii[d] = i - 1;
	F h = rr - i - 0.5;
	H[0][d] = K1(0.5 - h);
	H[1][d] = K0(h);
	H[2][d] = K1(0.5 + h);
	}
    }

template<typename F>
static void pcs_weights(int ii[3],F H[][3],const F r[3]) {
    auto pow3 = [](F x) { return x*x*x; };
    auto K0   = [](F h) { return 1.0/6.0 * ( 4.0 - 6.0*h*h + 3.0*h*h*h); };
    auto K1   = [&pow3](F h) { return 1.0/6.0 * pow3(2.0 - h); };
    for(int d=0; d<3; ++d) {
	F rr = r[d];              	 	    /* coordinates in subcube units [0,size] */
	ii[d] = (floorf(rr-1.5));
	F h = rr - (ii[d]+0.5);
	H[0][d] = K1(h);
	H[1][d] = K0(h-1);
	H[2][d] = K0(2-h);
	H[3][d] = K1(3-h);
	}
    }

template<typename F>
static void assign_mass(mass_array_t &masses, const F r[3], F mass,int iAssignment=3) {
    int    ii[3];
    F  H[4][3];
    switch(iAssignment) { // Get weights/indexes based on choosen assignment scheme
	case 0: ngp_weights(ii,H,r); break;
	case 1: cic_weights(ii,H,r); break;
	case 2: tsc_weights(ii,H,r); break;
	case 3: pcs_weights(ii,H,r); break;
	default: assert(iAssignment>=0 && iAssignment<=3); abort();
	}

    /* assign particle according to weights to neighboring nodes */
    for(int i=0; i<=iAssignment; ++i) {
	for(int j=0; j<=iAssignment; ++j) {
	    for(int k=0; k<=iAssignment; ++k) {
		masses(ii[0]+i,ii[1]+j,ii[2]+k) += H[i][0]*H[j][1]*H[k][2] * mass;
		}
	    }
	}
    }

static void flush_masses(PKD pkd,int nGrid,const mass_array_t &masses, const shape_t &lower) {
    auto wrap = [&nGrid](int i) { if (i>=nGrid) i-=nGrid; else if (i<0) i+=nGrid; return i; };
    for(auto i=masses.begin(); i!=masses.end(); ++i) {
    	if ( *i > 0.0f) {
	    shape_t loc = i.position() + lower;
	    loc[0] = wrap(loc[0]); loc[1] = wrap(loc[1]); loc[2] = wrap(loc[2]);
	    auto id = mdlFFTrId(pkd->mdl,pkd->fft,loc[0],loc[1],loc[2]);
	    auto idx = mdlFFTrIdx(pkd->mdl,pkd->fft,loc[0],loc[1],loc[2]);
	    float *p = reinterpret_cast<float*>(mdlVirtualFetch(pkd->mdl,CID_PK,idx,id));
	    *p += *i;
	    }
	}
    }


extern "C"
void dumpDensity(float *fftData,int nGrid){
    blitz::Array<float,3> raw(fftData,blitz::shape(nGrid+2,nGrid,nGrid),blitz::neverDeleteData,RegularArray());
    blitz::Array<float,3> masses = raw(blitz::Range(0,nGrid-1),blitz::Range(0,nGrid-1),blitz::Range(0,nGrid-1));

    int i, j, k;
    for(k=0; k<nGrid; ++k) {
        for(j=0; j<nGrid; ++j) {
            float m = 0.0;
	    for(i=0; i<nGrid; ++i) {
	    //for(i=10; i<20; ++i) {
	    	if (masses(i,j,k) > m) m = masses(i,j,k);
		}
	    std::clog << k << " " << j << " " << m << std::endl;
	    }
	std::clog << std::endl;
	}
    }

extern "C"
void pkdAssignMass(PKD pkd, uint32_t iLocalRoot, int nGrid, int iAssignment) {
    const std::size_t maxSize = 1000000;
    std::vector<float> data;
    data.reserve(maxSize); // Reserve maximum number
    shape_t index;
    position_t fPeriod(pkd->fPeriod), ifPeriod = 1.0 / fPeriod;

    assert(iAssignment>=0 && iAssignment<=3);

    std::vector<std::uint32_t> stack;
    stack.push_back(iLocalRoot);
    while( !stack.empty()) {
	tree_node *kdn = reinterpret_cast<tree_node *>(pkdTreeNode(pkd,stack.back()));
	stack.pop_back(); // Go to the next node in the tree
	BND bnd = pkdNodeGetBnd(pkd, kdn);
	position_t fCenter(bnd.fCenter), fMax(bnd.fMax);
	const shape_t ilower = shape_t(((fCenter - fMax) * ifPeriod + 0.5) * nGrid) - (iAssignment+1)/2;
	const shape_t iupper = shape_t(((fCenter + fMax) * ifPeriod + 0.5) * nGrid) + (iAssignment+1)/2;
	const shape_t ishape = iupper - ilower + 1;
	const float3_t flower = ilower;
	std::size_t size = blitz::product(ishape);

	if (size > maxSize) { // This cell is too large, so we split it and move on
	    assert(kdn->is_cell()); // At the moment we cannot handle enormous buckets
	    stack.push_back(kdn->iLower+1);
	    stack.push_back(kdn->iLower);
	    }
	else { // Assign the mass for this range of particles
	    data.resize(size); // Hold the right number of masses
	    mass_array_t masses(data.data(),ishape,blitz::neverDeleteData,RegularArray());
	    masses = 0.0f;
	    //masses.dumpStructureInformation(std::cout);
	    for( int i=kdn->pLower; i<=kdn->pUpper; ++i) { // All particles in this tree cell
	    	PARTICLE *p = pkdParticle(pkd,i);
	    	position_t dr; pkdGetPos1(pkd,p,dr.data()); // Centered on 0 with period fPeriod
	    	float3_t r(dr);
		r = (r * ifPeriod + 0.5) * nGrid - flower; // Scale and shift to fit in subcube
		assign_mass(masses, r.data(), pkdMass(pkd,p), iAssignment);
		}
	    flush_masses(pkd,nGrid,masses,ilower);
	    }
	}
    }

extern "C"
void pstAssignMass(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inAssignMass *in = reinterpret_cast<struct inAssignMass *>(vin);
    assert (nIn==sizeof(struct inAssignMass) );
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl, pst->idUpper, PST_ASSIGN_MASS, vin, nIn);
	pstAssignMass(pst->pstLower, vin, nIn, vout, pnOut);
	mdlGetReply(pst->mdl,rID, vout,pnOut);
	}
    else {
    	pkdAssignMass(plcl->pkd,ROOT,in->nGrid,in->iAssignment);
	}
    }

extern "C"
void msrAssignMass(MSR msr,int iAssignment,int nGrid) {
    static const char *schemes[] = {
    	"Nearest Grid Point (NGP)", "Cloud in Cell (CIC)",
        "Triangular Shaped Cloud (TSC)", "Piecewise Cubic Spline (PCS)" };
    struct inAssignMass mass;
    assert(iAssignment>=0 && iAssignment<=3);
    printf("Assigning mass using %s (order %d)\n",schemes[iAssignment],iAssignment);
    mass.nGrid = nGrid;
    mass.iAssignment = iAssignment;
    auto sec = msrTime();
    pstAssignMass(msr->pst, &mass, sizeof(mass), NULL, NULL);
    printf("Mass assignment complete, Wallclock: %f secs\n",msrTime() - sec);
    }
