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
#include <cstddef>
#include <vector>
#include <numeric>
#include <iostream>
#include "blitz/array.h"
#include "pst.h"
#include "master.h"
#include "aweights.hpp"

typedef blitz::Array<float,3> mass_array_t;
typedef blitz::TinyVector<int,3> shape_t;
typedef blitz::TinyVector<double,3> position_t;
typedef blitz::TinyVector<float,3> float3_t;
typedef blitz::ColumnMajorArray<3> RegularArray;

struct tree_node : public KDN {
    bool is_cell()   { return iLower!=0; }
    bool is_bucket() { return iLower==0; }
    };

template<int Order,typename F>
static void assign(mass_array_t &masses, const F r[3], F mass) {
    AssignmentWeights<Order,F> Hx(r[0]),Hy(r[1]),Hz(r[2]);

    for(int i=0; i<Order; ++i) {
	for(int j=0; j<Order; ++j) {
	    for(int k=0; k<Order; ++k) {
		masses(Hx.i+i,Hy.i+j,Hz.i+k) += Hx.H[i]*Hy.H[j]*Hz.H[k] * mass;
		}
	    }
	}
    }

template<typename F>
static void assign_mass(mass_array_t &masses, const F r[3], F mass,int iAssignment=4) {
    switch(iAssignment) {
	case 1: assign<1,F>(masses,r,mass); break;
	case 2: assign<2,F>(masses,r,mass); break;
	case 3: assign<3,F>(masses,r,mass); break;
	case 4: assign<4,F>(masses,r,mass); break;
	default: assert(iAssignment>=1 && iAssignment<=4); abort();
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


#if 0
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
#endif

extern "C"
void pkdAssignMass(PKD pkd, uint32_t iLocalRoot, int nGrid, float fDelta, int iAssignment) {
    const std::size_t maxSize = 100000; // We would like this to remain in L2 cache
    std::vector<float> data;
    data.reserve(maxSize); // Reserve maximum number
    shape_t index;
    position_t fPeriod(pkd->fPeriod), ifPeriod = 1.0 / fPeriod;

    assert(iAssignment>=1 && iAssignment<=4);

    std::vector<std::uint32_t> stack;
    stack.push_back(iLocalRoot);
    while( !stack.empty()) {
	tree_node *kdn = reinterpret_cast<tree_node *>(pkdTreeNode(pkd,stack.back()));
	stack.pop_back(); // Go to the next node in the tree
	BND bnd = pkdNodeGetBnd(pkd, kdn);
	position_t fCenter(bnd.fCenter), fMax(bnd.fMax);
	shape_t ilower = shape_t(floor(((fCenter - fMax) * ifPeriod + 0.5) * nGrid + fDelta)) - iAssignment/2;
	shape_t iupper = shape_t(floor(((fCenter + fMax) * ifPeriod + 0.5) * nGrid + fDelta)) + iAssignment/2;
	shape_t ishape = iupper - ilower + 1;
	float3_t flower = ilower;
	std::size_t size = blitz::product(ishape);

	if (size > maxSize) { // This cell is too large, so we split it and move on
	    if (kdn->is_cell()) { // At the moment we cannot handle enormous buckets
		stack.push_back(kdn->iLower+1);
		stack.push_back(kdn->iLower);
		}
	    // Huge bucket. Do a particle at a time.
	    else for( int i=kdn->pLower; i<=kdn->pUpper; ++i) { // All particles in this tree cell
	        PARTICLE *p = pkdParticle(pkd,i);
	        position_t dr; pkdGetPos1(pkd,p,dr.data()); // Centered on 0 with period fPeriod
	        float3_t r(dr);
	        r = (r * ifPeriod + 0.5) * nGrid + fDelta;
		ilower = shape_t(r) - iAssignment/2;
		iupper = shape_t(r) + iAssignment/2;
		ishape = iupper - ilower + 1;
		flower = ilower;
		size = blitz::product(ishape);
		r = r - flower; // Scale and shift to fit in subcube
		data.resize(size); // Hold the right number of masses
		mass_array_t masses(data.data(),ishape,blitz::neverDeleteData,RegularArray());
		masses = 0.0f;
		assign_mass(masses, r.data(), pkdMass(pkd,p), iAssignment);
		flush_masses(pkd,nGrid,masses,ilower);
		}
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
		r = (r * ifPeriod + 0.5) * nGrid + fDelta - flower; // Scale and shift to fit in subcube
		assign_mass(masses, r.data(), pkdMass(pkd,p), iAssignment);
		}
	    flush_masses(pkd,nGrid,masses,ilower);
	    }
	}
    }

extern "C"
int pstAssignMass(PST pst,void *vin,int nIn,void *vout,int nOut) {
    LCL *plcl = pst->plcl;
    struct inAssignMass *in = reinterpret_cast<struct inAssignMass *>(vin);
    assert (nIn==sizeof(struct inAssignMass) );
    if (pstNotCore(pst)) {
	int rID = mdlReqService(pst->mdl, pst->idUpper, PST_ASSIGN_MASS, vin, nIn);
	pstAssignMass(pst->pstLower, vin, nIn, NULL, 0);
	mdlGetReply(pst->mdl,rID,NULL,NULL);
	}
    else {
    	pkdAssignMass(plcl->pkd,ROOT,in->nGrid,0.0f,in->iAssignment);
	}
    return 0;
    }

void MSR::AssignMass(int iAssignment,int nGrid) {
    static const char *schemes[] = {
    	"Nearest Grid Point (NGP)", "Cloud in Cell (CIC)",
        "Triangular Shaped Cloud (TSC)", "Piecewise Cubic Spline (PCS)" };
    struct inAssignMass mass;
    assert(iAssignment>=1 && iAssignment<=4);
    printf("Assigning mass using %s (order %d)\n",schemes[iAssignment-1],iAssignment);
    mass.nGrid = nGrid;
    mass.iAssignment = iAssignment;
    auto sec = MSR::Time();
    pstAssignMass(pst, &mass, sizeof(mass), NULL, 0);
    printf("Mass assignment complete, Wallclock: %f secs\n",MSR::Time() - sec);
    }
