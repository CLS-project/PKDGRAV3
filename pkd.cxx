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
#include <cinttypes>
#include "io/iomodule.h"
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
    #include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <limits.h>
#ifdef HAVE_UNISTD_H
    #include <unistd.h>
    #ifdef __linux__
        #include <linux/fs.h>
    #else
        #define O_DIRECT 0
    #endif
#endif
#ifdef HAVE_SYS_TIME_H
    #include <sys/time.h>
#endif
#ifdef HAVE_SYS_TYPES_H
    #include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
    #include <sys/stat.h>
#endif
#ifdef __linux__
    #include <sys/resource.h>
#endif
#include <numeric>
#include <algorithm>
#include <boost/range/adaptor/reversed.hpp>
#include <gsl/gsl_spline.h>

#ifdef USE_ITT
    #include "ittnotify.h"
#endif
#include "cuda/cudautil.h"
#include "pkd.h"
#include "gravity/ewald.h"
#include "gravity/walk.h"
#include "gravity/grav.h"
#include "mdl.h"
#include "io/outtype.h"
#include "cosmo.h"
#include "SPH/SPHEOS.h"
#include "SPH/SPHpredict.h"
#include <stack>
extern "C" {
#include "core/healpix.h"
}
#ifdef COOLING
    #include "cooling/cooling.h"
#endif
#if defined(EEOS_JEANS) || defined(EEOS_POLYTROPE)
    #include "eEOS/eEOS.h"
#endif
#ifdef BLACKHOLES
    #include "blackhole/evolve.h"
#endif
#ifdef GRACKLE
    #include "cooling_grackle/cooling_grackle.h"
#endif
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#ifdef FEEDBACK
    #include "starformation/feedback.h"
#endif

#ifdef _MSC_VER
    #define FILE_PROTECTION (_S_IREAD | _S_IWRITE)
    typedef int ssize_t;
    #define open _open
    #define write _write
    #define read _read
    #define close _close
#else
    #define FILE_PROTECTION (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP)
#endif

using blitz::TinyVector;
using blitz::dot;
using blitz::max;
using blitz::min;

static void initLightBallOffsets(PKD pkd,double mrLCP) {
    Bound bnd(TinyVector<double,3>(-0.5),TinyVector<double,3>(0.5));
    double min2;
    int ix,iy,iz,l,nBox,nBoxMax;

    pkd->nLayerMax = ceil(mrLCP);
    pkd->nBoxLC = (int *)malloc(pkd->nLayerMax*sizeof(int));
    assert(pkd->nBoxLC != NULL);

    pkd->lcOffset0 = NULL;
    pkd->lcOffset1 = NULL;
    pkd->lcOffset2 = NULL;
    /*
    ** Set up the light cone offsets such that they proceed from the inner 8
    ** unit boxes outward layer by layer so that we can skip checks of the
    ** outer layers if we want.
    */
    nBox = nBoxMax = 0;
    for (l=0; l<pkd->nLayerMax; ++l) {
        for (ix=-l; ix<=l+1; ++ix) {
            for (iy=-l; iy<=l+1; ++iy) {
                for (iz=-l; iz<=l+1; ++iz) {
                    if ((ix>-l) && (ix<l+1) && (iy>-l) && (iy<l+1) && (iz>-l) && (iz<l+1)) continue;
                    TinyVector<double,3> r(ix - 0.5, iy - 0.5, iz - 0.5);
                    min2 = bnd.mindist(r);
                    if (min2 < mrLCP*mrLCP) {
                        if (pkd->Self()==0) printf("Lightcone replica:%d layer:%d r:<%f %f %f> min2:%f\n",nBox,l,r[0],r[1],r[2],min2);
                        if (nBox == nBoxMax) {
                            nBoxMax += 100;
                            pkd->lcOffset0 = (double *)realloc((void *)pkd->lcOffset0,nBoxMax*sizeof(double));
                            assert(pkd->lcOffset0 != NULL);
                            pkd->lcOffset1 = (double *)realloc((void *)pkd->lcOffset1,nBoxMax*sizeof(double));
                            assert(pkd->lcOffset1 != NULL);
                            pkd->lcOffset2 = (double *)realloc((void *)pkd->lcOffset2,nBoxMax*sizeof(double));
                            assert(pkd->lcOffset2 != NULL);
                        }
                        pkd->lcOffset0[nBox] = r[0];
                        pkd->lcOffset1[nBox] = r[1];
                        pkd->lcOffset2[nBox] = r[2];
                        ++nBox;
                    }
                }
            }
        }
        /*
        ** Can save the nBox here for skipping the deeper parts of the lightcone checks!
        */
        pkd->nBoxLC[l] = nBox;
    }
}

static void initLightConeOffsets(PKD pkd,int bBowtie,blitz::TinyVector<double,3> h,double alpha,double mrLCP) {
    int l,ix,iy,iz,nBox,nBoxMax;
    int *xy;
    int *xz;
    int *yz;
    double hm = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    double min2;

    pkd->nLayerMax = ceil(mrLCP);
    for (int j=0; j<3; ++j) h[j] /= hm;
    if (pkd->Self()==0) printf("mrLCP:%f nLayerMax:%d vector width:%d h:<%f %f %f>\n",mrLCP,pkd->nLayerMax,dvec::width(),h[0],h[1],h[2]);
    /*
    ** Check that all the components of h are positive and that h is not too close to the
    ** walls!
    */
    assert(h[0] > 0);
    assert(h[1] > 0);
    assert(h[2] > 0);

    pkd->nBoxLC = (int *)malloc(pkd->nLayerMax*sizeof(int));
    assert(pkd->nBoxLC != NULL);

    xy = (int *)malloc(pkd->nLayerMax*pkd->nLayerMax*sizeof(int));
    assert(xy != NULL);
    xz = (int *)malloc(pkd->nLayerMax*pkd->nLayerMax*sizeof(int));
    assert(xz != NULL);
    yz = (int *)malloc(pkd->nLayerMax*pkd->nLayerMax*sizeof(int));
    assert(yz != NULL);
    for (int j=0; j<pkd->nLayerMax*pkd->nLayerMax; ++j) xy[j] = xz[j] = yz[j] = 1; // by default all projected cells need to be checked
    double alpha_xy = atan2(tan(alpha),sqrt(h[0]*h[0] + h[1]*h[1])); // angle projected onto the xy plane
    double alpha_xz = atan2(tan(alpha),sqrt(h[0]*h[0] + h[2]*h[2])); // angle projected onto the xz plane
    double alpha_yz = atan2(tan(alpha),sqrt(h[1]*h[1] + h[2]*h[2])); // angle projected onto the xy plane
    double beta_yx = atan2(h[0],h[1]);
    assert(beta_yx > alpha_xy); // make sure it is far enough from the walls
    double beta_zx = atan2(h[0],h[2]);
    assert(beta_zx > alpha_xz); // make sure it is far enough from the walls
    double beta_zy = atan2(h[1],h[2]);
    assert(beta_zy > alpha_yz); // make sure it is far enough from the walls

    double beta_xy = atan2(h[1],h[0]);
    assert(beta_xy > alpha_xy); // make sure it is far enough from the walls
    double beta_xz = atan2(h[2],h[0]);
    assert(beta_xz > alpha_xz); // make sure it is far enough from the walls
    double beta_yz = atan2(h[2],h[1]);
    assert(beta_yz > alpha_yz); // make sure it is far enough from the walls
    if (pkd->Self()==0) printf("alpha:%f beta_xy:%f alpha_xy:%f beta_xz:%f alpha_xz:%f beta_yz:%f alpha_yz:%f\n",alpha,beta_xy,alpha_xy,beta_xz,alpha_xz,beta_yz,alpha_yz);
    for (int j=0; j<pkd->nLayerMax; ++j) {
        for (int i=0; i<pkd->nLayerMax; ++i) {
            double a1 = atan2(j,i+1);
            double a2 = atan2(j+1,i);
            if (beta_xy > a2 + alpha_xy || beta_xy < a1 - alpha_xy) xy[j*pkd->nLayerMax+i] = 0; // doesn't intersect
            if (beta_xz > a2 + alpha_xz || beta_xz < a1 - alpha_xz) xz[j*pkd->nLayerMax+i] = 0; // doesn't intersect
            if (beta_yz > a2 + alpha_yz || beta_yz < a1 - alpha_yz) yz[j*pkd->nLayerMax+i] = 0; // doesn't intersect
        }
    }
    pkd->lcOffset0 = NULL;
    pkd->lcOffset1 = NULL;
    pkd->lcOffset2 = NULL;
    /*
    ** Set up the light cone offsets such that they proceed from the inner 8
    ** unit boxes outward layer by layer so that we can skip checks of the
    ** outer layers if we want.
    */
    nBox = nBoxMax = 0;
    for (l=0; l<pkd->nLayerMax; ++l) {
        for (ix=0; ix<=l; ++ix) {
            for (iy=0; iy<=l; ++iy) {
                for (iz=0; iz<=l; ++iz) {
                    if (ix<l && iy<l && iz<l) continue;
                    if (xy[iy*pkd->nLayerMax+ix] && xz[iz*pkd->nLayerMax+ix] && yz[iz*pkd->nLayerMax+iy]) {
                        double r[3] = {ix + 0.5, iy + 0.5, iz + 0.5};
                        Bound bnd(TinyVector<double,3>(-0.5),TinyVector<double,3>(0.5));
                        min2 = bnd.mindist(r);
                        if (min2 > mrLCP*mrLCP) continue;
                        if (pkd->Self()==0) printf("Lightcone replica:%d layer:%d r:<%f %f %f>\n",nBox,l,r[0],r[1],r[2]);
                        if (nBox == nBoxMax) {
                            nBoxMax += 100;
                            pkd->lcOffset0 = (double *)realloc((void *)pkd->lcOffset0,nBoxMax*sizeof(double));
                            assert(pkd->lcOffset0 != NULL);
                            pkd->lcOffset1 = (double *)realloc((void *)pkd->lcOffset1,nBoxMax*sizeof(double));
                            assert(pkd->lcOffset1 != NULL);
                            pkd->lcOffset2 = (double *)realloc((void *)pkd->lcOffset2,nBoxMax*sizeof(double));
                            assert(pkd->lcOffset2 != NULL);
                        }
                        pkd->lcOffset0[nBox] = r[0];
                        pkd->lcOffset1[nBox] = r[1];
                        pkd->lcOffset2[nBox] = r[2];
                        ++nBox;
                        if (bBowtie) {
                            /*
                            ** Now include also the -,-,- octant for checks, as we will be producing a bowtie.
                            */
                            for (int j=0; j<3; ++j) r[j] = -r[j];
                            if (pkd->Self()==0) printf("Lightcone replica:%d layer:%d r:<%f %f %f>\n",nBox,l,r[0],r[1],r[2]);
                            if (nBox == nBoxMax) {
                                nBoxMax += 100;
                                pkd->lcOffset0 = (double *)realloc((void *)pkd->lcOffset0,nBoxMax*sizeof(double));
                                assert(pkd->lcOffset0 != NULL);
                                pkd->lcOffset1 = (double *)realloc((void *)pkd->lcOffset1,nBoxMax*sizeof(double));
                                assert(pkd->lcOffset1 != NULL);
                                pkd->lcOffset2 = (double *)realloc((void *)pkd->lcOffset2,nBoxMax*sizeof(double));
                                assert(pkd->lcOffset2 != NULL);
                            }
                            pkd->lcOffset0[nBox] = r[0];
                            pkd->lcOffset1[nBox] = r[1];
                            pkd->lcOffset2[nBox] = r[2];
                            ++nBox;
                        }
                    }
                }
            }
        }
        /*
        ** Can save the nBox here for skipping the deeper parts of the lightcone checks!
        */
        pkd->nBoxLC[l] = nBox;
    }
    /*
    ** Add more boxes to get to a multiple of the vector width.
    ** Just make sure the added boxes far enough away to not get
    ** included.
    */
    l = pkd->nLayerMax - 1;
    double r[3] = {2*l + 0.5, 2*l + 0.5, 2*l + 0.5};
    while (nBox%dvec::width()) {
        if (pkd->Self()==0) printf("Lightcone replica:%d layer:%d r:<%f %f %f>\n",nBox,l,r[0],r[1],r[2]);
        if (nBox == nBoxMax) {
            nBoxMax += 100;
            pkd->lcOffset0 = (double *)realloc((void *)pkd->lcOffset0,nBoxMax*sizeof(double));
            assert(pkd->lcOffset0 != NULL);
            pkd->lcOffset1 = (double *)realloc((void *)pkd->lcOffset1,nBoxMax*sizeof(double));
            assert(pkd->lcOffset1 != NULL);
            pkd->lcOffset2 = (double *)realloc((void *)pkd->lcOffset2,nBoxMax*sizeof(double));
            assert(pkd->lcOffset2 != NULL);
        }
        pkd->lcOffset0[nBox] = r[0];
        pkd->lcOffset1[nBox] = r[1];
        pkd->lcOffset2[nBox] = r[2];
        ++nBox;
    }
    pkd->nBoxLC[l] = nBox;
}

pkdContext::pkdContext(mdl::mdlClass *mdl,
                       int nStore,uint64_t nMinTotalStore,uint64_t nMinEphemeral,uint32_t nEphemeralBytes,
                       int nTreeBitsLo, int nTreeBitsHi,
                       int iCacheSize,int iCacheMaxInflight,int iWorkQueueSize,
                       const TinyVector<double,3> &fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,uint64_t nBH,
                       uint64_t mMemoryModel) : mdl(mdl),
    pLightCone(nullptr), pHealpixData(nullptr), csm(nullptr) {
    PARTICLE *p;
    uint32_t pi;
    int j,ism;

    io_init(&afiLightCone,0,0,IO_AIO|IO_LIBAIO);

#define RANDOM_SEED 1
    srand(RANDOM_SEED);

    this->nDark = nDark;
    this->nGas = nGas;
    this->nBH = nBH;
    this->nStar = nStar;
    this->nRejects = 0;
    this->fPeriod = fPeriod;

    this->uMinRungActive  = 0;
    this->uMaxRungActive  = 255;
    for (j=0; j<=IRUNGMAX; ++j) this->nRung[j] = 0;

    this->psGroupTable.nGroups = 0;
    this->psGroupTable.pGroup = NULL;
    this->veryTinyGroupTable = NULL;

#ifdef MDL_FFTW
    this->fft = NULL;
#endif

    /*
    ** Calculate the amount of memory (size) of each particle.  This is the
    ** size of a base particle (PARTICLE), plus any extra fields as defined
    ** by the current memory model.  Fields need to be added in order of
    ** descending size (i.e., doubles & int64 and then float & int32)
    */
    this->bNoParticleOrder = (mMemoryModel&PKD_MODEL_UNORDERED) ? 1 : 0;
    this->bIntegerPosition = (mMemoryModel&PKD_MODEL_INTEGER_POS) ? 1 : 0;

    particles.initialize(this->bIntegerPosition,this->bNoParticleOrder);

    if (!this->bIntegerPosition) particles.add<double[3]>(PKD_FIELD::oPosition,"r");
    if ( mMemoryModel & PKD_MODEL_PARTICLE_ID ) particles.add<int64_t>(PKD_FIELD::oParticleID,"id");
    /*
    ** Add a global group id. This is used when outputing an array of particles with
    ** one group assignment per particle.Usually only used in testing as it adds a
    ** 64 bit integer.
    */
    if ( mMemoryModel & PKD_MODEL_GLOBALGID ) particles.add<int64_t>(PKD_FIELD::oGlobalGid,"gid");
    if ( mMemoryModel & PKD_MODEL_VELOCITY && sizeof(vel_t) == sizeof(double))
        particles.add<double[3]>(PKD_FIELD::oVelocity,"v");
    if (this->bIntegerPosition) particles.add<int32_t[3]>(PKD_FIELD::oPosition,"r");
    if ( mMemoryModel & PKD_MODEL_SPH )
#ifdef OPTIM_UNION_EXTRAFIELDS
        particles.add<meshless::FIELDS,meshless::EXTRAFIELDS>(PKD_FIELD::oSph,"hydro");
#else
        particles.add<meshless::FIELDS>(PKD_FIELD::oSph,"hydro");
#endif

    if ( mMemoryModel & PKD_MODEL_NEW_SPH ) particles.add<sph::FIELDS>(PKD_FIELD::oNewSph,"hydro");
    if ( mMemoryModel & PKD_MODEL_STAR ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
        if (!particles.present(PKD_FIELD::oSph)) particles.add<meshless::FIELDS,meshless::EXTRAFIELDS>(PKD_FIELD::oSph,"hydro");
        // The star field is overlaid on the sph field
        particles.add<meshless::BLACKHOLE>(PKD_FIELD::oStar,"star",particles.offset(PKD_FIELD::oSph));
#else
        particles.add<meshless::STAR>(PKD_FIELD::oStar,"star");
#endif
    }

#ifdef BLACKHOLES
    if ( mMemoryModel & PKD_MODEL_BH ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
        if (!particles.present(PKD_FIELD::oSph)) particles.add<meshless::FIELDS,meshless::EXTRAFIELDS>(PKD_FIELD::oSph,"hydro");
        // The blackhole field is overlaid on the sph field
        particles.add<meshless::BLACKHOLE>(PKD_FIELD::oBH,"bh",particles.offset(PKD_FIELD::oSph));
#else
        particles.add<meshless::BLACKHOLE>(PKD_FIELD::oBH,"bh");
#endif
    }
#endif // BLACKHOLES

    if ( mMemoryModel & PKD_MODEL_VELSMOOTH )
        particles.add<VELSMOOTH>(PKD_FIELD::oVelSmooth,"vsmooth");
    if ( mMemoryModel & PKD_MODEL_VELOCITY ) {
        if (sizeof(vel_t) == sizeof(float)) {
            particles.add<float[3]>(PKD_FIELD::oVelocity,"v");
        }
    }
    if ( mMemoryModel & PKD_MODEL_ACCELERATION )
        particles.add<float[3]>(PKD_FIELD::oAcceleration,"a");

    if ( mMemoryModel & PKD_MODEL_MASS )
        particles.add<float>(PKD_FIELD::oMass,"m");

    if ( mMemoryModel & PKD_MODEL_SOFTENING )
        particles.add<float>(PKD_FIELD::oSoft,"soft");

    if ( mMemoryModel & (PKD_MODEL_SPH|PKD_MODEL_NEW_SPH|PKD_MODEL_BALL) )
        particles.add<float>(PKD_FIELD::oBall,"ball");
    if ( mMemoryModel & (PKD_MODEL_SPH|PKD_MODEL_NEW_SPH|PKD_MODEL_DENSITY) )
        particles.add<float>(PKD_FIELD::oDensity,"rho");

    if ( (mMemoryModel & PKD_MODEL_GROUPS) && !this->bNoParticleOrder) {
        particles.add<int32_t>(PKD_FIELD::oGroup,"group");
    }

    if ( mMemoryModel & PKD_MODEL_POTENTIAL ) {
        particles.add<float>(PKD_FIELD::oPotential,"phi");
    }
    particles.align();

    /*
    ** Tree node memory models
    */
    if (this->bIntegerPosition) tree.add<int32_t[3]>(KDN_FIELD::oNodePosition,"r");
    else tree.add<double[3]>(KDN_FIELD::oNodePosition,"r");
    if ( mMemoryModel & PKD_MODEL_NODE_BND ) {
        if (this->bIntegerPosition) tree.add<IntegerBound>(KDN_FIELD::oNodeBnd,"bnd");
        else tree.add<Bound>(KDN_FIELD::oNodeBnd,"bnd");
    }
    if ( mMemoryModel & PKD_MODEL_NODE_VBND )
        tree.add<Bound>(KDN_FIELD::oNodeVBnd,"vbnd");
    if ( (mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(double))
        tree.add<double[3]>(KDN_FIELD::oNodeVelocity,"v");
    if ( mMemoryModel & (PKD_MODEL_SPH|PKD_MODEL_BH) ) {
#ifdef OPTIM_REORDER_IN_NODES
        tree.add<int32_t>(KDN_FIELD::oNodeNgas,"ngas");
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
        tree.add<int32_t>(KDN_FIELD::oNodeNstar,"nstar");
#endif
        tree.add<int32_t>(KDN_FIELD::oNodeNbh,"nbh");
#endif
    }
    /*
    ** Three extra bounds are required by the fast gas SPH code.
    */
    if ( mMemoryModel & PKD_MODEL_NODE_SPHBNDS )
        tree.add<SPHBNDS>(KDN_FIELD::oNodeSphBounds,"sphbnd");

    if ( mMemoryModel & PKD_MODEL_NODE_BOB )
        tree.add<SPHBOB>(KDN_FIELD::oNodeBOB,"bob");

    if ( mMemoryModel & PKD_MODEL_NODE_MOMENT ) {
        tree.add<FMOMR>(KDN_FIELD::oNodeMom,"mom");
        tree.add<mass_t>(KDN_FIELD::oNodeMass,"mass",tree.offset(KDN_FIELD::oNodeMom) + offsetof(FMOMR,m));
    }
    else tree.add<mass_t>(KDN_FIELD::oNodeMass,"mass");

    /* The acceleration is required for the new time step criteria */
    if ( mMemoryModel & PKD_MODEL_NODE_ACCEL )
        tree.add<float[3]>(KDN_FIELD::oNodeAcceleration,"a");

    if ( (mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(float))
        tree.add<float[3]>(KDN_FIELD::oNodeVelocity,"v");

    assert(tree.ElementSize() > 0);
    if (tree.ElementSize() > pkdContext::MaxNodeSize()) {
        fprintf(stderr, "Node size is too large. Node size=%" PRIu64 ", max node size=%" PRIu64 "\n",
                (uint64_t)tree.ElementSize(), (uint64_t)MaxNodeSize());
    }
    assert(tree.ElementSize()<=MaxNodeSize());
    tree.align();

    /*
    ** We need to allocate one large chunk of memory for:
    **   (1) The particles, and,
    **   (2) The emphemeral storage,
    ** subject to the following constraints:
    **   (a) We must be able to store nStore+1 particles in both stores,
    **   (b) The total per-node ephemeral storage must be at least nMinEphemeral
    **   (c) The total size of the storage block must be at least nMinTotalStore
    **   (d) pStore must be page aligned on each thread
    **
    ** If (b) is not met then we increase the effective size of the ephemeral storage.
    ** If (c) is not met (even after increasing the ephemeral storage), then we increase
    ** the size of the block and use the left over for the first parts of the tree.
    **
    ** STILL TRUE?: We need one EXTRA storage location at the very end to use for
    ** calculating acceleration on arbitrary positions in space, for example
    ** determining the force on the sun. The easiest way to do this is to
    ** allocate one hidden particle, which won't interfere with the rest of
    ** the code (hopefully). this->pStore[this->nStore] is this particle.
    **
    */
    this->nEphemeralBytes = nEphemeralBytes;

    uint64_t nEphemeral = nStore + 1;
    if (nEphemeral * mdl->Cores()*nEphemeralBytes < nMinEphemeral)
        nEphemeral = (nMinEphemeral + mdl->Cores()*nEphemeralBytes - 1) / (mdl->Cores()*nEphemeralBytes);

    uint64_t nElements[3] = {uint64_t(nStore+1),nEphemeral,0};
    uint64_t nBytesPerElement[3] = {particles.ParticleSize(),nEphemeralBytes,sizeof(uint64_t)};
    void *pSegments[3];
    storageSize = mdl->new_shared_array(pSegments,3,nElements,nBytesPerElement,nMinTotalStore);
    if (storageSize == 0) {
        fprintf(stderr, "ERROR: unable to allocate main storage\n");
        abort();
    }

    storageBase = pSegments[0];
    particles.setStore(pSegments[0],nStore);
    this->pLite = pSegments[1];                 // Ephemeral storage
    /*
    ** Now we setup the node storage for the tree.  This storage is no longer
    ** continguous as the MDL now supports non-contiguous arrays.  We allocate
    ** a single "tile" for the tree.  If this is not sufficient, then additional
    ** tiles are allocated dynamically.  The default parameters allow for 2^32
    ** nodes total which is the integer limit anyway. We may use the extra storage
    ** from above if constraint (c) could not otherwise be met.
    */
    tree.initialize(&particles,this->bIntegerPosition,nTreeBitsLo,nTreeBitsHi,nElements[2]*nBytesPerElement[2],pSegments[2]);

    /*
    ** We also allocate a temporary particle used for swapping.  We need to do
    ** this now because the outside world can no longer know the size of a
    ** particle a priori.
    */
    this->pTempPRIVATE = static_cast<PARTICLE *>(malloc(particles.ParticleSize()));
    mdlassert(mdl,this->pTempPRIVATE != NULL);

    if ( iCacheSize > 0 ) mdlSetCacheSize(this->mdl,iCacheSize);
    mdl->SetCacheMaxInflight(iCacheMaxInflight);

    // This is cheeserific - chooses the largest specified
#if defined(USE_CUDA)
    this->cudaClient = new CudaClient(this->mdl->mpi->cuda,this->mdl->gpu);
    mdlSetCudaBufferSize(this->mdl,PP_CUDA_MEMORY_LIMIT,PP_CUDA_MEMORY_LIMIT);
#endif
#if defined(USE_METAL)
    this->metalClient = new MetalClient(*this->mdl);
#endif
    /*
    ** Initialize global group id.
    */
    if (particles.present(PKD_FIELD::oGlobalGid)) {
        for (pi=0; pi<(particles.FreeStore()+1); ++pi) {
            p = Particle(pi);
            particles.global_gid(p) = 0;
        }
    }
    /*
    ** Ewald stuff!
    */
    this->ew.nMaxEwhLoop = 0;
    /*
    ** Allocate Checklist.
    */
    this->cl = new clList(clFreeList);
    this->clNew = new clList(clFreeList);
    /*
    ** Allocate the stack.
    */
    this->nMaxStack = 30;
    this->S = new CSTACK[this->nMaxStack];
    assert(this->S != NULL);
    for (ism=0; ism<this->nMaxStack; ++ism) {
        this->S[ism].cl = new clList(clFreeList);
    }
    this->ga = NULL;

    this->profileBins = NULL;
    this->groupBin = NULL;

    this->tmpHopGroups = NULL;
    this->hopGroups = NULL;
    this->hopRootIndex = NULL;
    this->hopRoots = NULL;

    this->SPHoptions.TuFac = -1.0f;
    assert(NodeSize() > 0);
}

/// @brief Set the number of local particles
/// @param n number of local particles
/// @return number of local particles
int pkdContext::SetLocal(int n) {
    // It is important that if the particle store has changed that we invalidate the tree
    // This is done by setting the root node to all particles
    tree[ROOT]->set_local(0,n - 1);
    return particles.SetLocal(n);
}

/// @brief Increment the number of local particles
/// @param n number of particles to add
/// @return number of local particles
int pkdContext::AddLocal(int n) {
    return SetLocal(Local()+n);
}

pkdContext::~pkdContext() {
    PARTICLE *p;
    char **ppCList;
    uint32_t pi;
    int ism;

    /*
    ** Close caching space and free up nodes.
    */
    if (mdlCacheStatus(mdl,CID_CELL))      mdlFinishCache(mdl,CID_CELL);
    if (mdlCacheStatus(mdl,CID_CELL2))     mdlFinishCache(mdl,CID_CELL2);
    if (mdlCacheStatus(mdl,CID_PARTICLE2)) mdlFinishCache(mdl,CID_PARTICLE2);
    /*
    ** Free checklist.
    */
    delete cl;
    delete clNew;
    /*
    ** Free Stack.
    */
    for (ism=0; ism<nMaxStack; ++ism) {
        delete S[ism].cl;
    }
    delete [] S;
    if (ew.nMaxEwhLoop) {
        delete [] ewt.hx.f;
        delete [] ewt.hy.f;
        delete [] ewt.hz.f;
        delete [] ewt.hCfac.f;
        delete [] ewt.hSfac.f;
    }

    /* Only thread zero allocated this memory block  */
    mdlThreadBarrier(mdl);
    if (mdlCore(mdl) == 0) {
        mdl->delete_shared_array(storageBase,storageSize);
    }
    free(pTempPRIVATE);
    if (pLightCone) {
#ifdef _MSC_VER
        _aligned_free(pLightCone);
#else
        free(pLightCone);
#endif
    }
    if (pHealpixData) free(pHealpixData);
    io_free(&afiLightCone);
    if (csm) {
        csmFinish(csm);
        csm = NULL;
    }
#ifdef COOLING
    cooling_clean(cooling);
#endif
#ifdef STELLAR_EVOLUTION
    free(StelEvolData);
#endif
}

size_t pkdClCount(PKD pkd) {
    size_t nCount = pkd->cl->count();
    int i;
    for (i=0; i<pkd->nMaxStack; ++i)
        nCount += pkd->S[i].cl->count();
    return nCount;
}

size_t pkdClMemory(PKD pkd) {
    return pkd->cl->memory();
}

size_t pkdIlpMemory(PKD pkd) {
    return pkd->ilp.memory();
}

size_t pkdIlcMemory(PKD pkd) {
    return pkd->ilc.memory();
}

size_t pkdIllMemory(PKD pkd) {
    return pkd->ill.memory();
}

void pkdReadFIO(PKD pkd,FIO fio,uint64_t iFirst,int nLocal,double dvFac, double dTuFac) {
    float dummypot;
    TinyVector<double,3> r, vel;
    TinyVector<float,ELEMENT_COUNT> metals;
    float fMass,fSoft,fDensity,u,fTimer;
    FIO_SPECIES eSpecies;
    uint64_t iParticleID;

    mdlassert(pkd->mdl,fio != NULL);
#ifdef USE_ITT
    __itt_domain *domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle *shMyTask = __itt_string_handle_create("Read");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
    pkd->particles.clearClasses();
    if (pkd->particles.present(PKD_FIELD::oStar)) {
        /* Make sure star class established -- how do all procs know of these classes? How do we ensure they agree on the class identifiers? */
        auto p = pkd->particles[pkd->Local()];
        pkd->particles.setClass(0,0,0,FIO_SPECIES_STAR,&p);
    }
#ifdef BLACKHOLES
    assert(pkd->particles.present(PKD_FIELD::oMass));
    auto p = pkd->particles[pkd->Local()];
    pkd->particles.setClass(0,0,0,FIO_SPECIES_BH,&p);
#endif

    // Protect against uninitialized values
    fMass = 0.0f;
    fSoft = 0.0f;

    fioSeek(fio,iFirst,FIO_SPECIES_ALL);
    for (auto i = 0; i < nLocal; ++i) {
        auto p = pkd->particles[pkd->Local()+i];
        /*
        ** General initialization.
        */
        p.set_rung(0);
        if (!pkd->bNoParticleOrder) p.set_new_rung(0);
        p.set_marked(true);
        p.set_density(0.0);
        if (p.have_ball()) p.set_ball(0.0);
        /*
        ** Clear the accelerations so that the timestepping calculations do not
        ** get funny uninitialized values!
        */
        if ( p.have_acceleration() ) p.acceleration() = 0;
        float *pPot = p.have_potential() ? &p.potential() : &dummypot;
        p.set_group(0);

        /* Initialize New SPH fields if present */
        if (p.have_newsph()) {
            auto &NewSph = p.newsph();
            NewSph.u = NewSph.uDot = NewSph.divv = NewSph.Omega = 0.0;
        }

        /* Initialize Star fields if present */
        if (p.have_star()) {
            auto &Star = p.star();
            Star.fTimer = 0;
            /*      Star.iGasOrder = IORDERMAX;*/
        }

        eSpecies = fioSpecies(fio);
        switch (eSpecies) {
        case FIO_SPECIES_SPH:
            ;
            float afSphOtherData[2];
            fioReadSph(fio,&iParticleID,r.data(),vel.data(),&fMass,&fSoft,pPot,
                       &fDensity,&u,metals.data(),afSphOtherData);
            if (p.have_newsph()) {
                fSoft = 1.0f; // Dummy value, because the field in the file is used as hSmooth
                pkd->particles.setClass(fMass,fSoft,metals[0],eSpecies,&p);
            }
            else {
                pkd->particles.setClass(fMass,fSoft,0,eSpecies,&p);
            }
            p.set_density(fDensity);
            if (p.have_newsph()) {
                auto &NewSph = p.newsph();
                NewSph.u = -u; /* Can't do conversion until density known */
            }
            else {
                assert(dTuFac>0.0);
                p.set_ball(2.*afSphOtherData[0]);
                if (p.have_sph()) {
                    auto &Sph = p.sph();
                    Sph.ElemMass = metals * fMass;
#ifdef HAVE_METALLICITY
                    Sph.fMetalMass = afSphOtherData[1] * fMass;
#endif
                    // If the value is negative, means that it is a temperature
                    u = (u<0.0) ? -u*dTuFac : u;
                    Sph.Frho = 0.0;
                    Sph.Fmom = 0.0;
                    Sph.Fene = 0.0;
                    Sph.E = (u + 0.5*dvFac*dot(vel,vel)) * fMass;
                    Sph.Uint = u * fMass;
                    assert(Sph.E > 0.);
                    Sph.mom = fMass * vel * sqrt(dvFac);
                    Sph.lastMom = 0.;
                    Sph.lastE = Sph.E;
#ifdef ENTROPY_SWITCH
                    Sph.S = 0.0;
                    Sph.lastS = 0.0;
                    Sph.maxEkin = 0.0;
#endif
                    Sph.lastUint = Sph.Uint;
                    Sph.lastHubble = 0.0;
                    Sph.lastMass = fMass;
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
#ifdef BLACKHOLES
                    Sph.BHAccretor.iIndex = NOT_ACCRETED;
                    Sph.BHAccretor.iPid   = NOT_ACCRETED;
#endif
                    Sph.uWake = 0;
                }
            }
            break;
        case FIO_SPECIES_DARK:
            fioReadDark(fio,&iParticleID,r.data(),vel.data(),&fMass,&fSoft,pPot,&fDensity);
            pkd->particles.setClass(fMass,fSoft,0,eSpecies,&p);
            p.set_density(fDensity);
            break;
        case FIO_SPECIES_STAR:
            ;
            float afStarOtherData[4];
            fioReadStar(fio,&iParticleID,r.data(),vel.data(),&fMass,&fSoft,pPot,&fDensity,
                        metals.data(),&fTimer,afStarOtherData);
            pkd->particles.setClass(fMass,fSoft,0,eSpecies,&p);
            p.set_density(fDensity);
            if (p.have_star()) {
                auto &Star = p.star();
                Star.fTimer = fTimer;
                Star.omega  = 0.;
#ifdef FEEDBACK
                // We avoid that star in the IC could explode
                Star.bCCSNFBDone = 1;
                Star.bSNIaFBDone = 1;
                Star.fSNEfficiency = afStarOtherData[3];
#endif
#ifdef STELLAR_EVOLUTION
                Star.ElemAbun = metals;
                Star.fMetalAbun = afStarOtherData[0];
                Star.fInitialMass = afStarOtherData[1];
                Star.fLastEnrichTime = afStarOtherData[2];
#endif
            }
            break;
        case FIO_SPECIES_BH:
            p.set_ball(p.soft());
            float otherData[3];
            fioReadBH(fio,&iParticleID,r.data(),vel.data(),&fMass,&fSoft,pPot,
                      &fDensity,otherData,&fTimer);
            pkd->particles.setClass(fMass,fSoft,0,eSpecies,&p);
            if (p.have_bh()) {
                auto &BH = p.BH();
                BH.omega  = 0.;
                BH.dInternalMass = otherData[0];
                BH.lastUpdateTime = -1.;
                BH.dAccretionRate = otherData[1];
                BH.dFeedbackRate = 0.0;
                BH.dAccEnergy = otherData[2];
                BH.fTimer = fTimer;
                BH.bForceReposition = false;
            }
            break;
        default:
            fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
            assert(0);
        }
        p.set_position(r);
        if (!pkd->bNoParticleOrder) p.set_order(iParticleID);
        if (p.have_particle_id()) p.ParticleID() = iParticleID;

        if (p.have_velocity()) {
            auto &v = p.velocity();
            if (!p.is_gas()) {
                // IA: dvFac = a*a, and for the gas we already provide
                // the peculiar velocity in the IC
                v = vel * dvFac;
            }
            else {
                v = vel * sqrt(dvFac);
            }
        }
    }

    pkd->AddLocal(nLocal);
    pkd->nActive += nLocal;

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
}

void pkdEnforcePeriodic(PKD pkd,Bound bnd) {
    int i;
#if defined(USE_SIMD) && defined(__SSE2__)
    if (pkd->bIntegerPosition) {
        __m128i period = _mm_set1_epi32 (INTEGER_FACTOR);
        __m128i top = _mm_setr_epi32 ( (INTEGER_FACTOR/2)-1,(INTEGER_FACTOR/2)-1,(INTEGER_FACTOR/2)-1,0x7fffffff );
        __m128i bot = _mm_setr_epi32 (-(INTEGER_FACTOR/2),-(INTEGER_FACTOR/2),-(INTEGER_FACTOR/2),-0x80000000 );

        auto pPos = &pkd->particles.get<char>(pkd->Particle(0),PKD_FIELD::oPosition);
        const int iSize = pkd->particles.ParticleSize();
        for (i=0; i<pkd->Local(); ++i) {
            __m128i v = _mm_loadu_si128((__m128i *)pPos);
            __m128i *r = (__m128i *)pPos;
            pPos += iSize;
            _mm_prefetch(pPos,_MM_HINT_T0);
            v = _mm_sub_epi32(_mm_add_epi32(v,_mm_and_si128(_mm_cmplt_epi32(v,bot),period)),
                              _mm_and_si128(_mm_cmpgt_epi32(v,top),period));
            /* The fourth field (not part of position) is not modified because of bot/top */
            _mm_storeu_si128(r,v);
            //  r[0] = _mm_extract_epi32 (v,0);
            //  r[1] = _mm_extract_epi32 (v,1);
            //  r[2] = _mm_extract_epi32 (v,2);
        }
    }
    else
#endif
    {
        for (auto &p : pkd->particles) {
            p.set_position(bnd.wrap(p.position()));
        }
    }
}

/*
** x and y must have range [1,2) !
** returns key in the range [0,2^63-1]
*/
uint64_t hilbert2d(float x,float y) {
    uint64_t s = 0;
    uint32_t m,ux,uy,ut;

    union {
        float    f;
        uint32_t u;
    } punner;

    punner.f = x;
    ux = punner.u >> 2;
    punner.f = y;
    uy = punner.u >> 2;

    m = 0x00100000;

    while (m) {
        s = s << 2;
        if (ux&m) {
            if (uy&m) {
                s |= 2;
            }
            else {
                ut = ux;
                ux = ~uy;
                uy = ~ut;
                s |= 3;
            }
        }
        else {
            if (uy&m) {
                s |= 1;
            }
            else {
                ut = ux;
                ux = uy;
                uy = ut;
            }
        }
        m = m >> 1;
    }
    return s;
}

/*
** x, y and z must have range [1,2) !
*/
uint64_t hilbert3d(float x,float y,float z) {
    uint64_t s = 0;
    uint32_t m,ux,uy,uz,ut;

    union {
        float    f;
        uint32_t u;
    } punner;

    punner.f = x;
    ux = punner.u >> 2;
    punner.f = y;
    uy = punner.u >> 2;
    punner.f = z;
    uz = punner.u >> 2;
    /* Or: ux = (uint32_t)((x-1.0f) * 0x00200000)*/

    m = 0x00100000;
    while (m) {
        s = s << 3;

        if (ux&m) {
            if (uy&m) {
                if (uz&m) {
                    ut = ux;
                    ux = uy;
                    uy = ~uz;
                    uz = ~ut;
                    s |= 5;
                }
                else {
                    ut = uz;
                    uz = ux;
                    ux = uy;
                    uy = ut;
                    s |= 2;
                }
            }
            else {
                ux = ~ux;
                uy = ~uy;
                if (uz&m) {
                    s |= 4;
                }
                else {
                    s |= 3;
                }
            }
        }
        else {
            if (uy&m) {
                if (uz&m) {
                    ut = ux;
                    ux = uy;
                    uy = ~uz;
                    uz = ~ut;
                    s |= 6;
                }
                else {
                    ut = uz;
                    uz = ux;
                    ux = uy;
                    uy = ut;
                    s |= 1;
                }
            }
            else {
                if (uz&m) {
                    ut = uy;
                    uy = ux;
                    ux = ~uz;
                    uz = ~ut;
                    s |= 7;
                }
                else {
                    ut = uy;
                    uy = ux;
                    ux = uz;
                    uz = ut;
                    s |= 0;
                }
            }
        }
        m = m >> 1;
    }
    return s;
}

/*
** Partition particles between iFrom and iTo into those < fSplit and
** those >= to fSplit.  Find number and weight in each partition.
*/
int pkdWeight(PKD pkd,int d,double fSplit,int iSplitSide,int iFrom,int iTo,
              int *pnLow,int *pnHigh,double *pfLow,double *pfHigh) {
    int iPart;
    double fLower,fUpper;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
        iPart = pkdLowerPart(pkd,d,fSplit,iFrom,iTo);
        *pnLow = pkd->Local()-iPart;
        *pnHigh = iPart;
    }
    else {
        iPart = pkdUpperPart(pkd,d,fSplit,iFrom,iTo);
        *pnLow = iPart;
        *pnHigh = pkd->Local()-iPart;
    }
    /*
    ** Calculate the lower weight and upper weight BETWEEN the particles
    ** iFrom to iTo!
    */
    fLower = iPart - iFrom;
    fUpper = iTo - iPart + 1;
    if (iSplitSide) {
        *pfLow = fUpper;
        *pfHigh = fLower;
    }
    else {
        *pfLow = fLower;
        *pfHigh = fUpper;
    }
    return (iPart);
}

/*
** Partition particles between iFrom and iTo into those < fSplit and
** those >= to fSplit.  Find number and weight in each partition.
*/
int pkdWeightWrap(PKD pkd,int d,double fSplit,double fSplit2,int iSplitSide,
                  int iFrom,int iTo,int *pnLow,int *pnHigh) {
    int iPart;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (!iSplitSide) {
        iPart = pkdLowerPartWrap(pkd,d,fSplit,fSplit2,iFrom,iTo);
        *pnLow = iPart;
        *pnHigh = pkd->Local()-iPart;
    }
    else {
        iPart = pkdUpperPartWrap(pkd,d,fSplit,fSplit2,iFrom,iTo);
        *pnHigh = iPart;
        *pnLow = pkd->Local()-iPart;
    }
    return (iPart);
}

int pkdOrdWeight(PKD pkd,uint64_t iOrdSplit,int iSplitSide,int iFrom,int iTo,
                 int *pnLow,int *pnHigh) {
    int iPart;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
        iPart = pkdLowerOrdPart(pkd,iOrdSplit,iFrom,iTo);
        *pnLow = pkd->Local()-iPart;
        *pnHigh = iPart;
    }
    else {
        iPart = pkdUpperOrdPart(pkd,iOrdSplit,iFrom,iTo);
        *pnLow = iPart;
        *pnHigh = pkd->Local()-iPart;
    }
    return (iPart);
}

int pkdLowerPart(PKD pkd,int d,double fSplit,int i,int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto ii = std::partition(pi,pj+1,[d,fSplit](auto &p) {return p.position(d)>=fSplit;});
    return ii - pkd->particles.begin();
}

int pkdUpperPart(PKD pkd,int d,double fSplit,int i,int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto ii = std::partition(pi,pj+1,[d,fSplit](auto &p) {return p.position(d)<fSplit;});
    return ii - pkd->particles.begin();
}

int pkdLowerPartWrap(PKD pkd,int d,double fSplit1,double fSplit2,int i,int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    if (fSplit1 > fSplit2) {
        auto ii = std::partition(pi,pj+1,
        [d,fSplit1,fSplit2](auto &p) {
            auto r = p.position(d);
            return r < fSplit2 || r >= fSplit1;
        });
        return ii - pkd->particles.begin();
    }
    else {
        auto ii = std::partition(pi,pj+1,
        [d,fSplit1,fSplit2](auto &p) {
            auto r = p.position(d);
            return r < fSplit2 && r >= fSplit1;
        });
        return ii - pkd->particles.begin();
    }
}

int pkdUpperPartWrap(PKD pkd,int d,double fSplit1,double fSplit2,int i,int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;

    if (fSplit1 > fSplit2) {
        auto ii = std::partition(pi,pj+1,
        [d,fSplit1,fSplit2](auto &p) {
            auto r = p.position(d);
            return r >= fSplit2 && r < fSplit1;
        });
        return ii - pkd->particles.begin();
    }
    else {
        auto ii = std::partition(pi,pj+1,
        [d,fSplit1,fSplit2](auto &p) {
            auto r = p.position(d);
            return r >= fSplit2 || r < fSplit1;
        });
        return ii - pkd->particles.begin();
    }
}

int pkdLowerOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto split = [nOrdSplit](auto &p) {return p.order()>=nOrdSplit;};
    auto ii = std::partition(pi,pj+1,split);
    //assert (std::all_of(pi,ii,split) && std::none_of(ii,pj+1,split));
    return ii - pkd->particles.begin();
}

int pkdUpperOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto split = [nOrdSplit](auto &p) {return p.order()<nOrdSplit;};
    auto ii = std::partition(pi,pj+1,split);
    //assert (std::all_of(pi,ii,split) && std::none_of(ii,pj+1,split));
    return ii - pkd->particles.begin();
}

int pkdActiveOrder(PKD pkd) {
    auto i = std::partition(pkd->particles.begin(),pkd->particles.end(),[](auto &p) {return p.is_active();});
    return (pkd->nActive = i - pkd->particles.begin());
}

int pkdColRejects(PKD pkd,int nSplit) {
    int iRejects,i;

    mdlassert(pkd->mdl,pkd->nRejects == 0);

    pkd->nRejects = pkd->Local() - nSplit;
    iRejects = pkd->FreeStore() - pkd->nRejects;
    /*
    ** Move rejects to High memory.
    */
    if (pkd->Local() != pkd->FreeStore()) {
        for (i=pkd->nRejects-1; i>=0; --i)
            pkd->particles[iRejects+i] = pkd->particles[nSplit+i];
    }
    pkd->SetLocal(nSplit);
    return (pkd->nRejects);
}

int pkdSwapRejects(PKD pkd,int idSwap) {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;

    if (idSwap != -1) {
        nBuf = (pkdSwapSpace(pkd))*pkd->particles.ParticleSize();
        nOutBytes = pkd->nRejects*pkd->particles.ParticleSize();
        mdlassert(pkd->mdl,pkd->Local() + pkd->nRejects <= pkd->FreeStore());
        mdlSwap(pkd->mdl,idSwap,nBuf,pkd->Particle(pkd->Local()),
                nOutBytes,&nSndBytes,&nRcvBytes);
        pkd->AddLocal(nRcvBytes/pkd->particles.ParticleSize());
        pkd->nRejects -= nSndBytes/pkd->particles.ParticleSize();
    }
    return (pkd->nRejects);
}

void pkdSwapAll(PKD pkd, int idSwap) {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;
    int i;
    int iBuf;

    /*
    ** Move particles to High memory.
    */
    iBuf = pkdSwapSpace(pkd);
    for (i=pkd->Local()-1; i>=0; --i)
        pkd->CopyParticle(pkd->Particle(iBuf+i),pkd->Particle(i));
    nBuf = pkd->FreeStore()*pkd->particles.ParticleSize();
    nOutBytes = pkd->Local()*pkd->particles.ParticleSize();
    mdlSwap(pkd->mdl,idSwap,nBuf,pkd->particles, nOutBytes,
            &nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl,nSndBytes/pkd->particles.ParticleSize() == pkd->Local());
    pkd->SetLocal(nRcvBytes/pkd->particles.ParticleSize());
}

int pkdSwapSpace(PKD pkd) {
    return (pkd->FreeStore() - pkd->Local());
}

int pkdActive(PKD pkd) {
    return (pkd->nActive);
}

int pkdInactive(PKD pkd) {
    return (pkd->Local() - pkd->nActive);
}

/*
** Returns a pointer to the i'th KDN in the tree.  Used for fetching
** cache element.  Normal code should use pkd->tree[].
*/
void *pkdTreeNodeGetElement(void *vData,int i,int iDataSize) {
    auto pkd = static_cast<PKD>(vData);
    return pkd->tree[i];
}

int pkdColOrdRejects(PKD pkd,uint64_t nOrdSplit,int iSplitSide) {
    int nSplit;
    if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkd->Local()-1);
    else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkd->Local()-1);
    return pkdColRejects(pkd,nSplit);
}

#ifndef NEW_REORDER
void pkdLocalOrder(PKD pkd,uint64_t iMinOrder, uint64_t iMaxOrder) {
    int i;
    assert(pkd->Local() == iMaxOrder - iMinOrder + 1);
    for (i=0; i<pkd->Local(); ++i) {
        auto p1 = pkd->particles[i];

        assert(p1.order() >= iMinOrder && p1.order() <= iMaxOrder);
        while (p1.order() - iMinOrder !=  i) {
            auto p2 = pkd->particles[p1.order()-iMinOrder];
            swap(p1,p2);
        }
    }
    /* Above replaces: qsort(pkd->ParticleBase(),pkd->Local(),pkd->particles.ParticleSize(),cmpParticles); */
}
#endif

#define MAX_IO_BUFFER_SIZE (8*1024*1024)

void pkdCheckpoint(PKD pkd,const char *fname) {
    asyncFileInfo info;
    size_t nFileSize;
    int fd;
    io_init(&info, IO_MAX_ASYNC_COUNT, 0, IO_AIO|IO_LIBAIO);
    fd = io_create(&info, fname);
    if (fd<0) {
        perror(fname);
        abort();
    }
    nFileSize = pkd->particles.ParticleSize() * pkd->Local();
    char *pBuffer = (char *)pkd->particles.Element(0);
    while (nFileSize) {
        size_t count = nFileSize > MAX_IO_BUFFER_SIZE ? MAX_IO_BUFFER_SIZE : nFileSize;
        io_write(&info, pBuffer, count);
        pBuffer += count;
        nFileSize -= count;
    }
    io_close(&info);
}

/*****************************************************************************\
* Write particles received from another node
\*****************************************************************************/

static void writeParticle(PKD pkd,FIO fio,double dvFac,double dvFacGas,Bound bnd,particleStore::Particle &p) {
    TinyVector<double,3> v,r;
    TinyVector<float,ELEMENT_COUNT> metals;
    float fTimer;
    uint64_t iParticleID;

    float fPot = p.have_potential() ? p.potential() :0;
    if (p.have_velocity()) {
        /* IA: the gas velocity in the code is v = a \dot x
         *  and the dm/star velocity v = a^2 \dot x
         *
         *  However, we save both as \dot x such that
         *  they can be be directly compared at the output
         */
        v = p.velocity() * (p.is_gas() ? dvFacGas : dvFac);
    }
    else v = 0.0;

    auto fMass = p.mass();
    auto fSoft = p.soft();
    if (pkd->particles.fixedsoft() >= 0.0) fSoft = 0.0;
    if (p.have_particle_id()) iParticleID = p.ParticleID();
    else if (!pkd->bNoParticleOrder) iParticleID = p.order();
    else iParticleID = 0;
    float fDensity = p.have_density() ? p.density() : 0;

    r = bnd.wrap(p.position()); // Enforce periodic boundaries */
    // If it still doesn't lie in the "unit" cell then something has gone quite wrong with the
    // simulation. Either we have a super fast particle or the initial condition is somehow not conforming
    // to the specified periodic box in a gross way.
    assert(all(r>=bnd.lower() && r<bnd.upper()));

    /* IA: In the case of cosmological boxes, it is typical to have a box not centered on the origin.
     *  Such boxes are defined in the (+,+,+) octant. To convert to that system of reference,
     *  we just add half the box size to all positions.
     *
     *  Another option would be to change fCenter, but I do not if that could break other parts of the code...
     *  15/06/19: I have stoped using this as this can lead problems with mpi because we do not
     *  have available pkd->csm. This should be easy to fix for all halo-finder and post-process routines.
     */
    //if (bComove){
    //   for (j=0;j<3;++j){
    //     r[j] += bnd->fMax[j];
    //   }
    //}
    switch (p.species()) {
    case FIO_SPECIES_SPH:
        if (p.have_newsph()) {
            const auto &NewSph = p.newsph();
            assert(pkd->SPHoptions.TuFac > 0.0f);
            double T;
            float otherData[3];
            otherData[0] = otherData[1] = otherData[2] = 0.0f;
            T = SPHEOSTofRhoU(pkd,fDensity, NewSph.u, p.imaterial(), &pkd->SPHoptions);
            metals = 0.0f;
            metals[0] = p.imaterial();
            fSoft = p.ball() / 2.0f;
            fioWriteSph(fio,iParticleID,r.data(),v.data(),fMass,fSoft,fPot,
                        fDensity,T,metals.data(),0.0f,T,otherData);
        }
        else {
            assert(p.have_sph());
            auto &Sph = p.sph();
            {
#if defined(COOLING)
                const double dRedshift = dvFacGas - 1.;
                float temperature = cooling_get_temperature(pkd, dRedshift, pkd->cooling, p, &Sph);
#elif defined(GRACKLE)
                gr_float fMetalDensity = Sph.fMetalMass * Sph.omega;
                gr_float fSpecificUint = Sph.Uint / fMass;

                // Set field arrays.
                pkd->grackle_field->density[0]         = fDensity;
                pkd->grackle_field->internal_energy[0] = fSpecificUint;
                pkd->grackle_field->x_velocity[0]      = 1.; // Velocity input is not used
                pkd->grackle_field->y_velocity[0]      = 1.;
                pkd->grackle_field->z_velocity[0]      = 1.;
                // for metal_cooling = 1
                pkd->grackle_field->metal_density[0]   = fMetalDensity;

                int err;
                gr_float temperature;
                err = local_calculate_temperature(pkd->grackle_data, pkd->grackle_rates, pkd->grackle_units, pkd->grackle_field,
                                                  &temperature);
                if (err == 0) fprintf(stderr, "Error in calculate_temperature.\n");
#else
                float temperature = 0;
#endif

                metals = Sph.ElemMass / fMass;

#ifdef STAR_FORMATION
                float SFR = Sph.SFR;
#else
                float SFR=0.;
#endif

                float ph = 0.5 * p.ball();
                float otherData[3];
                otherData[0] = SFR;
                // Casting integers to floats will become a problem if the number of groups
                // reaches 2^24, but for the moment we have to live with it
                otherData[1] = p.have_global_gid() ? p.global_gid() : -1.;
#ifdef HAVE_METALLICITY
                otherData[2] = Sph.fMetalMass / fMass;
#endif

                fioWriteSph(fio,iParticleID,r.data(),v.data(),fMass,fSoft,fPot,fDensity,
                            Sph.Uint/fMass,metals.data(),ph,temperature,&otherData[0]);
            }
        }
        break;
    case FIO_SPECIES_DARK: {
        float otherData[2];
        otherData[0] = p.have_global_gid() ? p.global_gid() : -1.;
        fioWriteDark(fio,iParticleID,r.data(),v.data(),fMass,fSoft,fPot,fDensity, &otherData[0]);
    }
    break;
    case FIO_SPECIES_STAR: {
        auto &Star = p.star();
#ifdef STELLAR_EVOLUTION
        metals = Star.ElemAbun;
#endif
        float otherData[6];
        otherData[0] = Star.fTimer;
        otherData[1] = p.have_global_gid() ? p.global_gid() : -1.;
#ifdef STELLAR_EVOLUTION
        otherData[2] = Star.fMetalAbun;
        otherData[3] = Star.fInitialMass;
        otherData[4] = Star.fLastEnrichTime;
#endif
#ifdef FEEDBACK
        otherData[5] = Star.fSNEfficiency;
#endif
        fioWriteStar(fio,iParticleID,r.data(),v.data(),fMass,fSoft,fPot,fDensity,
                     metals.data(),&otherData[0]);
    }
    break;
    case FIO_SPECIES_BH: {
        const auto &BH = p.BH();
        fTimer = BH.fTimer;
        float otherData[6];
        otherData[0] = BH.dInternalMass;
        otherData[1] = BH.dAccretionRate;
        otherData[2] = BH.dEddingtonRatio;
        otherData[3] = BH.dFeedbackRate;
        otherData[4] = BH.dAccEnergy;
        otherData[5] = p.have_global_gid() ? p.global_gid() : -1.;
        fioWriteBH(fio,iParticleID,r.data(),v.data(),fMass,fSoft,fPot,fDensity,
                   otherData,fTimer);
    }
    break;
    case FIO_SPECIES_UNKNOWN:
        break;
    default:
        fprintf(stderr,"Unsupported particle type: %d\n",p.species());
        assert(0);
    }
}

struct packWriteCtx {
    PKD pkd;
    FIO fio;
    Bound bnd;
    double dvFac;
    double dTuFac;
    int iIndex;
};

static int unpackWrite(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packWriteCtx *ctx = (struct packWriteCtx *)vctx;
    PKD pkd = ctx->pkd;
    auto p = &pkd->particles[static_cast<PARTICLE *>(vBuff)];
    int n = nSize / pkd->particles.ParticleSize();
    int i;
    double dvFacGas = sqrt(ctx->dvFac);
    assert( n*pkd->particles.ParticleSize() == nSize);
    for (i=0; i<n; ++i,++p) {
        writeParticle(pkd,ctx->fio,ctx->dvFac,dvFacGas,ctx->bnd,*p);
    }
    return 1;
}

void pkdWriteFromNode(PKD pkd,int iNode, FIO fio,double dvFac,double dTuFac,Bound bnd) {
    struct packWriteCtx ctx;
    ctx.pkd = pkd;
    ctx.fio = fio;
    ctx.bnd = bnd;
    ctx.dvFac = dvFac;
    ctx.dTuFac = dTuFac;
    ctx.iIndex = 0;
#ifdef MPI_VERSION
    mdlRecv(pkd->mdl,iNode,unpackWrite,&ctx);
#endif
}

/*****************************************************************************\
* Send particles to be written
\*****************************************************************************/

static int packWrite(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packWriteCtx *ctx = (struct packWriteCtx *)vctx;
    PKD pkd = ctx->pkd;
    int nLeft = pkd->Local() - ctx->iIndex;
    int n = nSize / pkd->particles.ParticleSize();
    if ( n > nLeft ) n = nLeft;
    nSize = n*pkd->particles.ParticleSize();
    memcpy(vBuff,pkd->Particle(ctx->iIndex), nSize );
    ctx->iIndex += n;
    return nSize;
}

/* Send all particled data to the specified node for writing */
void pkdWriteViaNode(PKD pkd, int iNode) {
    struct packWriteCtx ctx;
    ctx.pkd = pkd;
    ctx.fio = NULL;
    ctx.dvFac = 1.0;
    ctx.iIndex = 0;
#ifdef MPI_VERSION
    mdlSend(pkd->mdl,iNode,packWrite, &ctx);
#endif
}

void pkdWriteHeaderFIO(PKD pkd, FIO fio, double dScaleFactor, double dTime,
                       uint64_t nDark, uint64_t nGas, uint64_t nStar, uint64_t nBH,
                       double dBoxSize, double h, int nProcessors, UNITS units) {
    char version[] = PACKAGE_VERSION;
    fioSetAttr(fio, HDF5_HEADER_G, "PKDGRAV version", FIO_TYPE_STRING, 1, version);
    fioSetAttr(fio, HDF5_HEADER_G, "Time", FIO_TYPE_DOUBLE, 1, &dTime);
    if (pkd->csm->val.bComove) {
        double z = 1./dScaleFactor - 1.;
        fioSetAttr(fio, 0, "Redshift", FIO_TYPE_DOUBLE, 1, &z);
    }
    int flag;
#ifdef STAR_FORMATION
    flag=1;
#else
    flag=0;
#endif
    fioSetAttr(fio, HDF5_HEADER_G, "Flag_Sfr", FIO_TYPE_INT, 1, &flag);

#ifdef FEEDBACK
    flag=1;
#else
    flag=0;
#endif
    fioSetAttr(fio, HDF5_HEADER_G, "Flag_Feedback", FIO_TYPE_INT, 1, &flag);

#ifdef COOLING
    flag=1;
#else
    flag=0;
#endif
    fioSetAttr(fio, HDF5_HEADER_G, "Flag_Cooling", FIO_TYPE_INT, 1, &flag);

#ifdef STAR_FORMATION
    flag=1;
#else
    flag=0;
#endif
    fioSetAttr(fio, HDF5_HEADER_G, "Flag_StellarAge", FIO_TYPE_INT, 1, &flag);

#ifdef HAVE_METALLICITY
    flag=1;
#else
    flag=0;
#endif
    fioSetAttr(fio, HDF5_HEADER_G, "Flag_Metals", FIO_TYPE_INT, 1, &flag);

    // In HDF5 format the position and velocities are always stored as doubles
    flag=1;
    fioSetAttr(fio, HDF5_HEADER_G, "Flag_DoublePrecision", FIO_TYPE_INT, 1, &flag);

    fioSetAttr(fio, HDF5_HEADER_G, "BoxSize", FIO_TYPE_DOUBLE, 1, &dBoxSize);
    fioSetAttr(fio, HDF5_HEADER_G, "NumFilesPerSnapshot", FIO_TYPE_INT, 1, &nProcessors);

    /* Prepare the particle information in tables
     * For now, we only support one file per snapshot,
     * this bParaWrite=0
     *
     * TODO: Check and debug parallel HDF5
     */
    unsigned int numPart_file[6] = {0,0,0,0,0,0};
    fioSetAttr(fio, HDF5_HEADER_G, "NumPart_Total_HighWord", FIO_TYPE_UINT32, 6, &numPart_file[0]);
    //int numPart_all[6];

    numPart_file[0] = nGas;
    numPart_file[1] = nDark;
    numPart_file[2] = 0;
    numPart_file[3] = 0;
    numPart_file[4] = nStar;
    numPart_file[5] = nBH;

    fioSetAttr(fio, HDF5_HEADER_G, "NumPart_ThisFile", FIO_TYPE_UINT32, 6, &numPart_file[0]);
    fioSetAttr(fio, HDF5_HEADER_G, "NumPart_Total", FIO_TYPE_UINT32, 6, &numPart_file[0]);

    double massTable[6] = {0,0,0,0,0,0};
    // This is not yet fully supported, as the classes do not have to match the
    //  six available particle types.
    // However, we add this in the header so it can be parsed by other tools
    fioSetAttr(fio, HDF5_HEADER_G, "MassTable", FIO_TYPE_DOUBLE, 6, &massTable[0]);

    float fSoft = pkdSoft(pkd,pkd->Particle(0)); // we take any particle
    fioSetAttr(fio, HDF5_HEADER_G, "Softening", FIO_TYPE_FLOAT, 1, &fSoft);

    /*
     * Cosmology header
     */
    if (pkd->csm->val.bComove) {
        flag = 1;
        fioSetAttr(fio, HDF5_COSMO_G, "Omega_m", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dOmega0);
        fioSetAttr(fio, HDF5_COSMO_G, "Omega_lambda", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dLambda);
        fioSetAttr(fio, HDF5_COSMO_G, "Omega_b", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dOmegab);
        fioSetAttr(fio, HDF5_COSMO_G, "Hubble0", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dHubble0);
        fioSetAttr(fio, HDF5_COSMO_G, "Cosmological run", FIO_TYPE_INT, 1, &flag);
        fioSetAttr(fio, HDF5_COSMO_G, "HubbleParam", FIO_TYPE_DOUBLE, 1, &h);

        // Keep a copy also in the Header for increased compatibility
        fioSetAttr(fio, HDF5_HEADER_G, "Omega0", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dOmega0);
        fioSetAttr(fio, HDF5_HEADER_G, "OmegaLambda", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dLambda);
        fioSetAttr(fio, HDF5_HEADER_G, "OmegaB", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dOmegab);
        fioSetAttr(fio, HDF5_HEADER_G, "Hubble0", FIO_TYPE_DOUBLE, 1, &pkd->csm->val.dHubble0);
        fioSetAttr(fio, HDF5_HEADER_G, "Cosmological run", FIO_TYPE_INT, 1, &flag);
        fioSetAttr(fio, HDF5_HEADER_G, "HubbleParam", FIO_TYPE_DOUBLE, 1, &h);
    }
    else {
        flag = 0;
        fioSetAttr(fio, HDF5_COSMO_G, "Cosmological run", FIO_TYPE_INT, 1, &flag);
    }

    /*
     * Units header
     */
    fioSetAttr(fio, HDF5_UNITS_G, "MsolUnit", FIO_TYPE_DOUBLE, 1, &units.dMsolUnit);
    fioSetAttr(fio, HDF5_UNITS_G, "KpcUnit", FIO_TYPE_DOUBLE, 1, &units.dKpcUnit);
    fioSetAttr(fio, HDF5_UNITS_G, "SecUnit", FIO_TYPE_DOUBLE, 1, &units.dSecUnit);
    fioSetAttr(fio, HDF5_UNITS_G, "KmPerSecUnit", FIO_TYPE_DOUBLE, 1, &units.dKmPerSecUnit);
    fioSetAttr(fio, HDF5_UNITS_G, "GmPerCcUnit", FIO_TYPE_DOUBLE, 1, &units.dGmPerCcUnit);
    fioSetAttr(fio, HDF5_UNITS_G, "ErgPerGmUnit", FIO_TYPE_DOUBLE, 1, &units.dErgPerGmUnit);
    fioSetAttr(fio, HDF5_UNITS_G, "ErgUnit", FIO_TYPE_DOUBLE, 1, &units.dErgUnit);
    fioSetAttr(fio, HDF5_UNITS_G, "GasConst", FIO_TYPE_DOUBLE, 1, &units.dGasConst);

}

/*****************************************************************************\
* Send an array/vector to the specified node
\*****************************************************************************/
struct packArrayCtx {
    PKD pkd;
    double dvFac;
    PKD_FIELD field;
    int iIndex;
    int iUnitSize;
    int bMarked;
};

char *pkdPackArray(PKD pkd,int iSize,void *vBuff,int *piIndex,int n,PKD_FIELD field,int iUnitSize,double dvFac,int bMarked) {
    auto pBuff = static_cast<char *>(vBuff);
    int iIndex = *piIndex;

    while (iIndex<n && iSize>=iUnitSize) {
        PARTICLE *p = pkd->Particle(iIndex++);
        if (bMarked && !p->bMarked) continue;
        if (field==PKD_FIELD::oPosition) {
            auto &d = * reinterpret_cast<TinyVector<double,3>*>(pBuff);
            d = pkd->particles.position(p);
        }
        else if (field==PKD_FIELD::oVelocity) {
            auto &V = * reinterpret_cast<TinyVector<float,3>*>(pBuff);
            V = pkd->particles.velocity(p) * dvFac;
        }
        else {
            const char &src = pkd->particles.get<char>(p,field);
            memcpy(pBuff,&src,iUnitSize);
        }
        pBuff += iUnitSize;
        iSize -= iUnitSize;
    }
    *piIndex = iIndex;
    return pBuff;
}

static int packArray(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packArrayCtx *ctx = (struct packArrayCtx *)vctx;
    PKD pkd = ctx->pkd;
    char *pBuff = (char *)vBuff;
    char *pEnd = pkdPackArray(pkd,nSize,pBuff,&ctx->iIndex,pkd->Local(),ctx->field,ctx->iUnitSize,ctx->dvFac,ctx->bMarked);
    return pEnd-pBuff;
}

/* Send all particled data to the specified node for writing */
void pkdSendArray(PKD pkd, int iNode, PKD_FIELD field, int iUnitSize,double dvFac,int bMarked) {
    struct packArrayCtx ctx;
    ctx.pkd = pkd;
    ctx.dvFac = dvFac;
    ctx.field = field;
    ctx.iUnitSize = iUnitSize;
    ctx.iIndex = 0;
    ctx.bMarked = bMarked;
#ifdef MPI_VERSION
    mdlSend(pkd->mdl,iNode,packArray, &ctx);
#endif
}

/*****************************************************************************\
* Receive an array/vector from a specified node
\*****************************************************************************/

struct unpackArrayCtx {
    char *pDest;
};

static int unpackArray(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct unpackArrayCtx *ctx = (struct unpackArrayCtx *)vctx;
    memcpy(ctx->pDest,vBuff,nSize);
    ctx->pDest += nSize;
    return 1;
}

void *pkdRecvArray(PKD pkd,int iNode, void *pDest, int iUnitSize) {
    struct unpackArrayCtx ctx;
    ctx.pDest = static_cast<char *>(pDest);
#ifdef MPI_VERSION
    mdlRecv(pkd->mdl,iNode,unpackArray,&ctx);
#endif
    return ctx.pDest;
}

/*****************************************************************************\
*
\*****************************************************************************/

uint32_t pkdWriteFIO(PKD pkd,FIO fio,double dvFac,double dTuFac,Bound bnd) {
    double dvFacGas = sqrt(dvFac);
    for (auto &p : pkd->particles) {
        writeParticle(pkd,fio,dvFac,dvFacGas,bnd,p);
    }
    return pkd->particles.Local();
}

void pkdSetSoft(PKD pkd,double dSoft) {
    pkd->particles.SetSoft(dSoft);
}

void pkdSetSmooth(PKD pkd,double dSmooth) {
    for (auto &p : pkd->particles) {
        if ((p.is_gas()||p.is_star()||p.is_bh()) && p.ball()==0.0)
            p.set_ball(dSmooth);
    }
}

void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul) {
    pkd->particles.PhysicalSoft(dSoftMax,dFac,bSoftMaxMul);
}

static void initSetMarked(void *vpkd, void *v) {}
static void combSetMarked(void *vpkd, void *v1, const void *v2) {
    PARTICLE *p1 = (PARTICLE *)v1;
    const PARTICLE *p2 = (const PARTICLE *)v2;
    if (p2->bMarked) p1->bMarked = 1;
#ifdef NN_FLAG_IN_PARTICLE
    if (p2->bNNflag) p1->bNNflag = 1;
#endif
}

void extensiveMarkerTest(PKD pkd, struct pkdTimestepParameters *ts, SPHOptions *SPHoptions) {
    std::stack<std::pair<int,int>> cellStack;

    // Add the toptree cell corresponding to the ROOT cell to the stack
    cellStack.push(std::make_pair(pkd->iTopTree[ROOT],pkd->Self()));
    int nParticles = 0;

    while (!cellStack.empty()) {
        auto [iCell, id] = cellStack.top();
        cellStack.pop();
        auto c = (id == pkd->Self()) ? pkd->tree[iCell] : pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,CID_CELL,iCell,id))];

        // Walking the top tree
        if (c->is_top_tree()) {
            auto [iLower, idLower, iUpper, idUpper] = c->get_child_cells(id);
            cellStack.push(std::make_pair(iUpper, idUpper));
            cellStack.push(std::make_pair(iLower, idLower));
            continue;
        }

        /* We are either on
        ** the toptree cell that corresponds to my local cell
        ** or on one of the 2 children for cell corresponding to remote threads id != pkd->Self()
        ** so we can loop through the particles between c->lower() and c->upper()
        ** and mdlFetch them from that thread
        */
        for (auto pj=c->lower(); pj<=c->upper(); ++pj) {
            auto p = (id == pkd->Self()) ? pkd->particles[pj] : pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,CID_PARTICLE,pj,id))];
            ++nParticles;

            // Shortcut if necessary flag already set, as the test will pass in the end
            if ((SPHoptions->doSetDensityFlags && p.marked()) || (SPHoptions->doSetNNflags && p.NN_flag())) {
                if (id != pkd->Self()) mdlRelease(pkd->mdl,CID_PARTICLE,&p);
                continue;
            }

            // Get position and ball size of particle p
            auto pr = p.position();
            float fBallFactor = (SPHoptions->dofBallFactor) ? SPHoptions->fBallFactor : 1.0f;
            float pBall2 = std::min(SPHoptions->ballSizeLimit, p.ball() * fBallFactor);
            pBall2 *= pBall2;

            int needsCheck = 0;
            // Loop over all particles in my root
            auto rootc = pkd->tree[ROOT];
            for (auto qj=rootc->lower(); qj<=rootc->upper(); ++qj) {
                auto q = pkd->particles[qj];

                // Skip those that are not active
                if (!q.is_rung_range(ts->uRungLo,ts->uRungHi) && !(SPHoptions->useDensityFlags && q.marked()) && !(SPHoptions->useNNflags && q.NN_flag())) continue;

                // Get position and ball size of particle q
                auto qr = q.position();
                float qBall2 = std::min(SPHoptions->ballSizeLimit, q.ball() * fBallFactor);
                qBall2 *= qBall2;

                // Calculate distance squared between particle p and q
                auto dist = pr - qr;
                float dist2 = dot(dist,dist);

                // Do check for gather
                if (dist2 < qBall2) {
                    needsCheck = 1;
                    break;
                }
                // Do check for scatter
                if (dist2 < pBall2) {
                    needsCheck = 1;
                    break;
                }
            }

            // Now we need to check
            if (needsCheck) {
                // Check that NN flag is set if applicable
                if (SPHoptions->doSetNNflags) assert(p.NN_flag());
                // Check that density flag is set if applicable
                if (SPHoptions->doSetDensityFlags) assert(p.marked());
            }
            if (id != pkd->Self()) mdlRelease(pkd->mdl,CID_PARTICLE,&p);
        }
    }
    assert(nParticles == pkd->nGas);
}

void pkdGravAll(PKD pkd,
                struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
                double dTime,int nReps,int bPeriodic,int bGPU,
                int bEwald,int iRoot1, int iRoot2,
                double fEwCut,double fEwhCut,double dThetaMin,SPHOptions *SPHoptions,
                uint64_t *pnActive,
                double *pdPart,double *pdPartNumAccess,double *pdPartMissRatio,
                double *pdCell,double *pdCellNumAccess,double *pdCellMissRatio,
                double *pdFlop,uint64_t *pnRung) {

    double dActive;
    double dPartSum;
    double dCellSum;
    int i;

#ifdef USE_ITT
    __itt_domain *domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle *shMyTask = __itt_string_handle_create("Gravity");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif

    /*
    ** Clear all the rung counters to be safe.
    */
    for (i=0; i<=IRUNGMAX; ++i) pkd->nRung[i] = 0;

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    mdlTimeReset(pkd->mdl);
#endif

    /*
    ** Set up Ewald tables and stuff.
    */
    if (bPeriodic && bEwald && SPHoptions->doGravity) {
        pkdEwaldInit(pkd,nReps,fEwCut,fEwhCut,bGPU); /* ignored in Flop count! */
    }
    pkdCopySPHOptionsToDevice(pkd, SPHoptions, bGPU);
    /*
    ** Start particle caching space (cell cache already active).
    */
    if (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) {
        mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,pkd->particles,pkd->particles.ParticleSize(),
                   pkd->Local(),NULL,initSetMarked,combSetMarked);
    }
    else {
        mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkd->particles,pkd->particles.ParticleSize(),
                   pkd->Local());
    }

    /*
    ** Calculate newtonian gravity, including replicas if any.
    */
    *pdFlop = 0.0;
    dPartSum = 0.0;
    dCellSum = 0.0;
    pkd->dFlopSingleCPU = pkd->dFlopDoubleCPU = 0.0;
    pkd->dFlopSingleGPU = pkd->dFlopDoubleGPU = 0.0;
    pkd->nWpPending = 0;
    pkd->nTilesTotal = 0;
    pkd->nTilesCPU = 0;

    *pnActive = pkdGravWalk(pkd,kick,lc,ts,
                            dTime,nReps,bPeriodic && bEwald,bGPU,
                            iRoot1,iRoot2,0,dThetaMin,pdFlop,&dPartSum,&dCellSum,SPHoptions);

    assert(pkd->nWpPending == 0);

    if (SPHoptions->doExtensiveILPTest && (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags)) {
        mdlFlushCache(pkd->mdl,CID_PARTICLE);
        mdlCacheBarrier(pkd->mdl,CID_PARTICLE);
        extensiveMarkerTest(pkd, ts, SPHoptions);
    }

    dActive = (double)(*pnActive);
    if (*pnActive) {
        *pdPart = dPartSum/dActive;
        *pdCell = dCellSum/dActive;
    }
    else {
        assert(dPartSum == 0 && dCellSum == 0);
        *pdPart = 0;  /* for the statistics we don't count this processor, see pstGravity(). */
        *pdCell = 0;
    }
    /*
    ** Get caching statistics.
    */
    if (*pnActive) {
        *pdCellNumAccess = mdlNumAccess(pkd->mdl,CID_CELL)/dActive;
        *pdPartNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE)/dActive;
    }
    else {
        *pdCellNumAccess = 0;
        *pdPartNumAccess = 0;
    }
    *pdCellMissRatio = 100.0*mdlMissRatio(pkd->mdl,CID_CELL);      /* as a percentage */
    *pdPartMissRatio = 100.0*mdlMissRatio(pkd->mdl,CID_PARTICLE);  /* as a percentage */
    /*
    ** Output flops count in GFlops!
    */
    *pdFlop *= 1e-9;
    /*
    ** Stop particle caching space.
    */
    mdlFinishCache(pkd->mdl,CID_PARTICLE);

    for (i=0; i<=IRUNGMAX; ++i) pnRung[i] = pkd->nRung[i];

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
}

/* This became easier; we already calculated these values on the fly */
void pkdCalcEandL(PKD pkd,double &T,double &U,double &Eth,TinyVector<double,3> &L,TinyVector<double,3> &F,double &W) {
    T = pkd->dEnergyT;
    U = pkd->dEnergyU;
    W = pkd->dEnergyW;
    L = pkd->dEnergyL;
    F = pkd->dEnergyF;
    Eth = 0.0;
    if (pkd->particles.present(PKD_FIELD::oSph)) {
        for (auto &p : pkd->particles) {
            if (p.is_gas()) Eth += p.sph().E;
        }
    }
}

static void flushLightCone(PKD pkd) {
    size_t count = pkd->nLightCone * sizeof(LIGHTCONEP);
    io_write(&pkd->afiLightCone,pkd->pLightCone,count);
    pkd->nLightCone = 0;
}

static void initHealpix(void *vpkd, void *v) {
    healpixData *m = (healpixData *)v;
    m->nGrouped = 0;
    m->nUngrouped = 0;
    m->fPotential = 0;
}

static uint32_t SumWithSaturate(uint32_t a,uint32_t b) {
    uint64_t sum = (uint64_t)a + b;
    if (sum > 0xffffffffu) return 0xffffffffu;
    else return sum;
}

static void combHealpix(void *vctx, void *v1, const void *v2) {
    healpixData *m1 = (healpixData *)v1;
    const healpixData *m2 = (const healpixData *)v2;
    m1->nGrouped = SumWithSaturate(m1->nGrouped,m2->nGrouped);
    m1->nUngrouped = SumWithSaturate(m1->nUngrouped,m2->nUngrouped);
    m1->fPotential += m2->fPotential;
}

void pkdLightConeClose(PKD pkd,const char *healpixname) {
    int i;
    size_t nWrite;
    if (pkd->afiLightCone.fd > 0) {
        flushLightCone(pkd);
        io_close(&pkd->afiLightCone);
    }
    if (pkd->nSideHealpix) {
        assert(healpixname && healpixname[0]);
        mdlFinishCache(pkd->mdl,CID_HEALPIX);
        int fd = open(healpixname,O_CREAT|O_WRONLY|O_TRUNC,FILE_PROTECTION);
        if (fd<0) {
            perror(healpixname);
            abort();
        }
        for (i=0; i<pkd->nHealpixPerDomain; ++i) {
            uint64_t sum = pkd->pHealpixData[i].nGrouped;
            sum += pkd->pHealpixData[i].nUngrouped;
            if (sum) pkd->pHealpixData[i].fPotential /= sum;
        }
        nWrite = pkd->nHealpixPerDomain * sizeof(*pkd->pHealpixData);
        if (write(fd,pkd->pHealpixData,nWrite) != nWrite) {
            perror("Wrong size writing Healpix");
        }
        close(fd);
    }
}

void pkdLightConeOpen(PKD pkd,const char *fname,int nSideHealpix) {
    int i;
    if (fname[0]) {
        if (io_create(&pkd->afiLightCone,fname) < 0) {
            perror(fname);
            abort();
        }
    }
    else pkd->afiLightCone.fd = -1;

    pkd->nSideHealpix = nSideHealpix;
    if (pkd->nSideHealpix) {
        pkd->nHealpixPerDomain = ( nside2npix64(pkd->nSideHealpix) + mdlThreads(pkd->mdl) - 1) / mdlThreads(pkd->mdl);
        pkd->nHealpixPerDomain = (pkd->nHealpixPerDomain+15) & ~15;
        if (pkd->pHealpixData==NULL) {
            pkd->pHealpixData = static_cast<healpixData *>(malloc(pkd->nHealpixPerDomain * sizeof(*pkd->pHealpixData)));
            assert(pkd->pHealpixData!=NULL);
        }
        for (i=0; i<pkd->nHealpixPerDomain; ++i) {
            pkd->pHealpixData[i].nGrouped = 0;
            pkd->pHealpixData[i].nUngrouped = 0;
            pkd->pHealpixData[i].fPotential = 0;
        }
        mdlCOcache(pkd->mdl,CID_HEALPIX,NULL,
                   pkd->pHealpixData,sizeof(*pkd->pHealpixData),
                   pkd->nHealpixPerDomain,pkd,initHealpix,combHealpix);
    }
}

void addToLightCone(PKD pkd,double dvFac,double *r,float fPot,PARTICLE *p,int bParticleOutput) {
    const auto &v = pkd->particles.velocity(p);
    if (pkd->afiLightCone.fd>0 && bParticleOutput) {
        LIGHTCONEP *pLC = pkd->pLightCone;
        pLC[pkd->nLightCone].id = p->iOrder;
        pLC[pkd->nLightCone].pos[0] = r[0];
        pLC[pkd->nLightCone].pos[1] = r[1];
        pLC[pkd->nLightCone].pos[2] = r[2];
        /*
        ** Convert momentum (a^2 * x_dot to physical peculiar velocity (a * x_dot in sim units)
        ** This is easy to do now, and somewhat more involved to do later.
        */
        pLC[pkd->nLightCone].vel[0] = v[0]*dvFac;
        pLC[pkd->nLightCone].vel[1] = v[1]*dvFac;
        pLC[pkd->nLightCone].vel[2] = v[2]*dvFac;
        pLC[pkd->nLightCone].pot    = fPot;
        if (++pkd->nLightCone == pkd->nLightConeMax) flushLightCone(pkd);
    }
    if (pkd->nSideHealpix) {
        int64_t iPixel = vec2pix_ring64(pkd->nSideHealpix, r);
        assert(iPixel >= 0);
        int id  = iPixel / pkd->nHealpixPerDomain;
        int idx = iPixel - id*pkd->nHealpixPerDomain;
        assert(id<mdlThreads(pkd->mdl));
        assert(idx < pkd->nHealpixPerDomain);
        auto m = static_cast<healpixData *>(mdlVirtualFetch(pkd->mdl,CID_HEALPIX,idx,id));
        if (pkdGetGroup(pkd,p)) {
            if (m->nGrouped < 0xffffffffu) ++m->nGrouped; /* Increment with saturate */
        }
        else {
            if (m->nUngrouped < 0xffffffffu) ++m->nUngrouped; /* Increment with saturate */
        }
        m->fPotential += fPot;
    }
}

/*
** Drift particles whose Rung falls between uRungLo (large step) and uRungHi (small step) inclusive,
** and those whose destination activity flag is set.
**
** Note that the drift funtion no longer wraps the particles around the periodic "unit" cell. This is
** now done by Domain Decomposition only.
*/
void pkdDrift(PKD pkd,int iRoot,double dTime,double dDelta,double dDeltaVPred,double dDeltaTime,int bDoGas) {
    int i;
    TinyVector<double,3> rfinal,r0,dr;
    int pLower, pUpper;

    if (iRoot>=0) {
        auto pRoot = pkd->tree[iRoot];
        pLower = pRoot->lower();
        pUpper = pRoot->upper();
    }
    else {
        pLower = 0;
        pUpper = pkd->Local();
    }

    mdlDiag(pkd->mdl, "Into pkdDrift\n");
    assert(pkd->particles.present(PKD_FIELD::oVelocity));

    auto dMin = pkd->bnd.lower();
    auto dMax = pkd->bnd.upper();
    /*
    ** Update particle positions
    */
    if (bDoGas) {
        if (pkd->particles.present(PKD_FIELD::oSph) /*pkd->param.bMeshlessHydro*/) {
            assert(pkd->particles.present(PKD_FIELD::oSph));
            assert(pkd->particles.present(PKD_FIELD::oAcceleration));
            for (i=pLower; i<=pUpper; ++i) {
                auto p = pkd->particles[i];
                const auto &v = p.velocity();

                if (p.is_gas()) {
                    // As for gas particles we use dx/dt = v/a, we must use the
                    // Kick factor provided by pkdgrav.
                    // See Stadel 2001, Appendix 3.8
                    dr = v * dDeltaVPred;
                }
                else {
                    dr = v * dDelta;
                }

#ifdef FORCE_1D
                dr[2] = 0.0;
                dr[1] = 0.0;
#endif
#ifdef FORCE_2D
                dr[2] = 0.0;
#endif
                r0 = p.position();
                p.set_position(rfinal = r0 + dr);
                assert(isfinite(rfinal[0]));
                assert(isfinite(rfinal[1]));
                assert(isfinite(rfinal[2]));
                dMin = min(dMin,rfinal);
                dMax = max(dMax,rfinal);
            }
        }
        else {
            assert(pkd->particles.present(PKD_FIELD::oNewSph));
            for (i=pLower; i<=pUpper; ++i) {
                auto p = pkd->particles[i];
                const auto &v = p.velocity();
                // if (p.is_gas()) {
                // auto &NewSph = p.newsph();
                // float dfBalldt = 1.0f / 3.0f * p.ball() * p.density() * NewSph.divv;
                // p.set_ball(p.ball() + dDelta * dfBalldt);
                // }
                r0 = p.position();
                p.set_position(rfinal = r0 + dDelta*v);
                dMin = min(dMin,rfinal);
                dMax = max(dMax,rfinal);
            }
        }
    }
    else {
        for (i=pLower; i<=pUpper; ++i) {
            auto p = pkd->particles[i];
            const auto &v = p.velocity();
            r0 = p.position();
            p.set_position(rfinal = r0 + dDelta*v);
            dMin = min(dMin,rfinal);
            dMax = max(dMax,rfinal);
        }
    }
    pkd->bnd = Bound(dMin,dMax);
    mdlDiag(pkd->mdl, "Out of pkdDrift\n");
}

#ifdef OPTIM_REORDER_IN_NODES
/// @brief Order particles in buckets by type
/// @param pkd
/// After this function is complete, particles in buckets will be in the following order:
///     Gas, Star, Blackhole, Dark
void pkdReorderWithinNodes(PKD pkd) {
    for (int i=NRESERVED_NODES; i<pkd->Nodes(); i++) {
        auto node = pkd->tree[i];
        if (node->is_bucket()) { // We are in a bucket
            auto start = node->begin();
            auto i = std::partition(start,node->end(),[](auto &p) {return p.is_gas();});
            node->Ngas() = i - start;
            start = i;
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
            // We perform another swap, just to have nGas->nStar->DM
            i = std::partition(start,node->end(),[](auto &p) {return p.is_star();});
            node->Nstar() = i - start;
            start = i;
#endif
            i = std::partition(start,node->end(),[](auto &p) {return p.is_bh();});
            node->Nbh() = i - start;
        }
    }
    if (mdlCacheStatus(pkd->mdl,CID_CELL)) mdlFinishCache(pkd->mdl,CID_CELL);
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,pkd->NodeSize(),pkd->Nodes());
}
#endif

void pkdEndTimestepIntegration(PKD pkd, struct inEndTimestep in) {
    double pDelta, dScaleFactor, dHubble;
    TinyVector<double,3> pa;

    int bComove = pkd->csm->val.bComove;

    mdlDiag(pkd->mdl, "Into pkdComputePrimiteVars\n");
    assert(pkd->particles.present(PKD_FIELD::oVelocity));
#ifndef USE_MFM
    assert(pkd->particles.present(PKD_FIELD::oMass));
#endif

    if (bComove) {
        dScaleFactor = csmTime2Exp(pkd->csm,in.dTime);
        dHubble = csmTime2Hub(pkd->csm,in.dTime);
    }
    else {
        dScaleFactor = 1.0;
        dHubble = 0.0;
    }
#ifdef GRACKLE
    pkdGrackleUpdate(pkd, dScaleFactor, in.achCoolingTable, in.units);
#endif
    for (auto &p : pkd->particles) {
        if (p.is_gas() && p.is_active()) {
            auto &sph = p.sph();

            // ##### Add ejecta from stellar evolution
#ifdef STELLAR_EVOLUTION
            if (in.bChemEnrich && sph.fReceivedMass > 0.0f) {
                pkdAddStellarEjecta(pkd, p, sph, in.dConstGamma);
            }
#endif

            if (in.dDelta > 0) {
                pDelta = in.dTime - sph.lastUpdateTime;
            }
            else {
                pDelta = 0.0;
            }

            // ##### Gravity
            hydroSourceGravity(pkd, p, &sph, pDelta, pa, dScaleFactor, bComove);

            // ##### Expansion effects
            hydroSourceExpansion(pkd, p, &sph,
                                 pDelta, dScaleFactor, dHubble, bComove, in.dConstGamma);

            // ##### Synchronize Uint, Etot (and possibly S)
            hydroSyncEnergies(pkd, p, &sph, pa, in.dConstGamma);

            // ##### Cooling
#ifdef COOLING
            double dRedshift = 1./dScaleFactor - 1.;
            const float delta_redshift = -pDelta * dHubble * (dRedshift + 1.);
            cooling_cool_part(pkd, pkd->cooling, p, &sph, pDelta, in.dTime, delta_redshift, dRedshift);
#endif
#ifdef GRACKLE
            pkdGrackleCooling(pkd, p, pDelta, in.dTuFac);
#endif

#if defined(EEOS_JEANS) || defined(EEOS_POLYTROPE)
            // ##### Effective Equation Of State
            const double a_inv3 = 1. / (dScaleFactor * dScaleFactor * dScaleFactor);
            const double dFlooru = eEOSEnergyFloor<vec<double,double>,mmask<bool>>(a_inv3, p.density(), p.ball(),
                                   in.dConstGamma, in.eEOS);
            if (dFlooru != NOT_IN_EEOS) {
                const double dEOSUint = p.mass() * dFlooru;
                if (sph.Uint < dEOSUint) {
                    sph.E = sph.E - sph.Uint;
                    sph.Uint = dEOSUint;
                    sph.E = sph.E + sph.Uint;
                }
            }
#endif

#if defined(FEEDBACK) || defined(BLACKHOLES)
            // ##### Apply feedback
            pkdAddFBEnergy(pkd, p, &sph, in.dConstGamma);
#endif

            // Actually set the primitive variables
            hydroSetPrimitives(pkd, p, &sph, in.dTuFac, in.dConstGamma);

            // Set 'last*' variables for next timestep
            hydroSetLastVars(pkd, p, &sph, pa, dScaleFactor, in.dTime, in.dDelta, in.dConstGamma);

            hydroResetFluxes(&sph);
        }
        else if (p.is_bh() && p.is_active()) {
#ifdef BLACKHOLES
            pkdBHIntegrate(pkd, p, in.dTime, in.dDelta, in.dBHRadiativeEff);
#endif
        }
    }
}

void pkdSetupInterpScale(PKD pkd,double dBoxSize,double mrMax) {
    const int nTable=1000;
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    double dr,rt[nTable],at_inv[nTable];
    /*
    ** Setup lookup table.
    */
    dr = mrMax/(nTable-1);
    for (int i=0; i<nTable; ++i) {
        rt[i] = i*dr;
        at_inv[i] = 1.0/csmComoveLookbackTime2Exp(pkd->csm,rt[i]/dLightSpeed);
    }
    pkd->interp_scale = gsl_spline_alloc(gsl_interp_cspline,nTable);
    gsl_spline_init(pkd->interp_scale,rt,at_inv,nTable);
}

void pkdLightConeVel(PKD pkd,double dBoxSize) {
    const int nTable=1000;
    const double rMax=3.0;
    gsl_spline *scale;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    double dr,rt[nTable],at_inv[nTable];
    double dvFac;
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    int i;

    assert(pkd->particles.present(PKD_FIELD::oVelocity));
    /*
    ** Setup lookup table.
    */
    dr = rMax/(nTable-1);
    for (i=0; i<nTable; ++i) {
        rt[i] = i*dr;
        at_inv[i] = 1.0/csmComoveLookbackTime2Exp(pkd->csm,rt[i]/dLightSpeed);
    }
    scale = gsl_spline_alloc(gsl_interp_cspline,nTable);
    gsl_spline_init(scale,rt,at_inv,nTable);
    /*
    ** Now loop over all particles using this table.
    */
    for (auto &p : pkd->particles) {
        auto &v = p.velocity();
        auto r = p.position();
        auto r2 = dot(r,r);
        /*
        ** Use r -> 1/a spline table.
        */
        dvFac = gsl_spline_eval(scale,sqrt(r2),acc);
        /* input velocities are momenta p = a^2*x_dot and we want v_pec = a*x_dot */
        v *= dvFac;
    }
    gsl_spline_free(scale);
    gsl_interp_accel_free(acc);
}

/*
 * Stripped down versions of routines from master.c
 */
void pkdKickKDKOpen(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    if (pkd->csm->val.bComove) {
        dDelta = csmComoveKickFac(pkd->csm,dTime,dDelta);
    }
    pkdKick(pkd,dTime,dDelta,0,0,0,0,uRungLo,uRungHi);
}

void pkdKickKDKClose(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    if (pkd->csm->val.bComove) {
        dDelta = csmComoveKickFac(pkd->csm,dTime,dDelta);
    }
    pkdKick(pkd,dTime,dDelta,0,0,0,0,uRungLo,uRungHi);
}

void pkdKick(PKD pkd,double dTime,double dDelta,int bDoGas,double dDeltaVPred,double dDeltaU,double dDeltaUPred,uint8_t uRungLo,uint8_t uRungHi) {
    assert(pkd->particles.present(PKD_FIELD::oVelocity));
    assert(pkd->particles.present(PKD_FIELD::oAcceleration));

    if (bDoGas) { // pkd->param.bMeshlessHydro
        assert(pkd->particles.present(PKD_FIELD::oSph));
        for (auto &p : pkd->particles) {
            if (!p.is_gas()) {
                if (p.is_rung_range(uRungLo,uRungHi)) {
                    auto &a = p.acceleration();
                    auto &v = p.velocity();
                    v += a*dDelta;
                }
            }
        }
    }
    else {
        for (auto &p : pkd->particles) {
            if (p.is_rung_range(uRungLo,uRungHi)) {
                auto &a = p.acceleration();
                auto &v = p.velocity();
                v += a*dDelta;
            }
        }
    }

    mdlDiag(pkd->mdl, "Done pkdkick\n");
}

/* Kick the tree at iRoot. */
void pkdKickTree(PKD pkd,double dTime,double dDelta,double dDeltaVPred,double dDeltaU,double dDeltaUPred,int iRoot) {
    /* Skip to local tree */
    auto c = pkd->tree[iRoot];
    while (c->is_remote()) c = pkd->tree[iRoot = c->lchild()];

    /* Now just kick all of the particles in the tree */
    for (auto &p : *c) {
        auto &v = p.velocity();
        auto &a = p.acceleration();
        v += a*dDelta;
        a = 0.0;
    }
}

void pkdInitCosmology(PKD pkd, struct csmVariables *cosmo) {
    /*
    ** Need to be careful to correctly copy the cosmo
    ** parameters. This is a bit ugly.
    */
    if (pkd->csm) csmFinish(pkd->csm);
    csmInitialize(&pkd->csm);
    pkd->csm->val = *cosmo;
    if (pkd->csm->val.classData.bClass) {
        csmClassGslInitialize(pkd->csm);
    }
}

/*
** Initialize Lightcone stuff.
*/
void pkdInitLightcone(PKD pkd,int bBowtie,int bLightConeParticles,double dBoxSize,double dRedshiftLCP,double alphaLCP,blitz::TinyVector<double,3> hLCP) {
#ifdef __linux__
    uint64_t nPageSize = sysconf(_SC_PAGESIZE);
#else
    uint64_t nPageSize = 512;
#endif

    /*
    ** Initialize Lookup table for converting lightcone velocities to physical (sim units)
    */
    double dTimeLCP = csmExp2Time(pkd->csm,1.0/(1.0+dRedshiftLCP));
    double mrLCP = dLightSpeedSim(dBoxSize)*csmComoveKickFac(pkd->csm,dTimeLCP,(csmExp2Time(pkd->csm,1.0) - dTimeLCP));
    pkdSetupInterpScale(pkd,dBoxSize,mrLCP);
    /*
    ** Initialize light cone offsets.
    */
    if (pkd->Self() == 0) printf("alphaLCP = %f\n",alphaLCP);
    if (alphaLCP < 0) {
        initLightBallOffsets(pkd,mrLCP);
    }
    else {
        initLightConeOffsets(pkd,bBowtie,hLCP,alphaLCP,mrLCP);
    }
    /*
    ** allocate enough space for light cone particle output
    */
    uint64_t nLightConeBytes = (1024*1024*16);
    pkd->nLightConeMax = nLightConeBytes / sizeof(LIGHTCONEP);
    pkd->nLightCone = 0;
    if (bLightConeParticles) {
        void *v;
#ifdef _MSC_VER
        pkd->pLightCone = _aligned_malloc(nLightConeBytes, nPageSize);
#else
        if (posix_memalign(&v, nPageSize, nLightConeBytes)) pkd->pLightCone = NULL;
        else pkd->pLightCone = static_cast<LIGHTCONEP *>(v);
#endif
        mdlassert(pkd->mdl,pkd->pLightCone != NULL);
        io_init(&pkd->afiLightCone,8,2*1024*1024,IO_AIO|IO_LIBAIO);
    }
    else {
        pkd->afiLightCone.nBuffers = 0;
        pkd->pLightCone = NULL;
    }
    pkd->afiLightCone.fd = -1;
    pkd->pHealpixData = NULL;
}

void pkdZeroNewRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung) {  /* JW: Ugly -- need to clean up */
    if (!pkd->bNoParticleOrder) {
        for (auto &p : pkd->particles) {
            if (p.is_active()) p.set_new_rung(0);
        }
    }
}

void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
                  double dDelta, int iMaxRung,
                  double dEta,double dVelFac,double dAccFac,
                  int bDoGravity,int bEpsAcc,double dhMinOverSoft) {
    assert(pkd->particles.present(PKD_FIELD::oVelocity));
    assert(pkd->particles.present(PKD_FIELD::oAcceleration));
    assert(!pkd->bNoParticleOrder);

    for (auto &p : pkd->particles) {
        if (p.is_active()) {
            const auto &v = p.velocity();
            const auto &a = p.acceleration();
            double fSoft = p.soft();
            double vel = dot(v,v);
            double acc = dot(a,a);
            mdlassert(pkd->mdl,vel >= 0);
            vel = sqrt(vel)*dVelFac;
            mdlassert(pkd->mdl,acc >= 0);
            acc = sqrt(acc)*dAccFac;
            double dT = (acc>0 && bEpsAcc) ? dEta*sqrt(fSoft/acc) : FLOAT_MAXVAL;
            int uNewRung = pkdDtToRung(dT,dDelta,iMaxRung);
            if (uNewRung > p.new_rung()) p.set_new_rung(uNewRung);
        }
    }
}

void pkdChemCompInit(PKD pkd, struct inChemCompInit in) {

    for (auto &p : pkd->particles) {
        if (p.is_gas()) {
            auto &Sph = p.sph();
            float fMass = p.mass();

            if (Sph.ElemMass[ELEMENT_H] < 0.0f) {
                Sph.ElemMass[ELEMENT_H]  = in.dInitialH  * fMass;
#ifdef HAVE_HELIUM
                Sph.ElemMass[ELEMENT_He] = in.dInitialHe * fMass;
#endif
#ifdef HAVE_CARBON
                Sph.ElemMass[ELEMENT_C]  = in.dInitialC  * fMass;
#endif
#ifdef HAVE_NITROGEN
                Sph.ElemMass[ELEMENT_N]  = in.dInitialN  * fMass;
#endif
#ifdef HAVE_OXYGEN
                Sph.ElemMass[ELEMENT_O]  = in.dInitialO  * fMass;
#endif
#ifdef HAVE_NEON
                Sph.ElemMass[ELEMENT_Ne] = in.dInitialNe * fMass;
#endif
#ifdef HAVE_MAGNESIUM
                Sph.ElemMass[ELEMENT_Mg] = in.dInitialMg * fMass;
#endif
#ifdef HAVE_SILICON
                Sph.ElemMass[ELEMENT_Si] = in.dInitialSi * fMass;
#endif
#ifdef HAVE_IRON
                Sph.ElemMass[ELEMENT_Fe] = in.dInitialFe * fMass;
#endif
            }
#ifdef HAVE_METALLICITY
            if (Sph.fMetalMass < 0.0f)
                Sph.fMetalMass = in.dInitialMetallicity * fMass;
#endif
        }

#ifdef STELLAR_EVOLUTION
        else if (p.is_star()) {
            auto &Star = p.star();

            if (Star.ElemAbun[ELEMENT_H] < 0.0f) {
                Star.ElemAbun[ELEMENT_H]  = in.dInitialH;
#ifdef HAVE_HELIUM
                Star.ElemAbun[ELEMENT_He] = in.dInitialHe;
#endif
#ifdef HAVE_CARBON
                Star.ElemAbun[ELEMENT_C]  = in.dInitialC;
#endif
#ifdef HAVE_NITROGEN
                Star.ElemAbun[ELEMENT_N]  = in.dInitialN;
#endif
#ifdef HAVE_OXYGEN
                Star.ElemAbun[ELEMENT_O]  = in.dInitialO;
#endif
#ifdef HAVE_NEON
                Star.ElemAbun[ELEMENT_Ne] = in.dInitialNe;
#endif
#ifdef HAVE_MAGNESIUM
                Star.ElemAbun[ELEMENT_Mg] = in.dInitialMg;
#endif
#ifdef HAVE_SILICON
                Star.ElemAbun[ELEMENT_Si] = in.dInitialSi;
#endif
#ifdef HAVE_IRON
                Star.ElemAbun[ELEMENT_Fe] = in.dInitialFe;
#endif
            }
            if (Star.fMetalAbun < 0.0f)
                Star.fMetalAbun = in.dInitialMetallicity;
        }
#endif //STELLAR_EVOLUTION
    }
}

void pkdCorrectEnergy(PKD pkd, double dTuFac, double z, double dTime, int iDirection ) {
    /*PARTICLE *p;
    meshless::FIELDS *sph;
    int i;
    double T,E;*/
    switch (iDirection)  {
    case CORRECTENERGY_IN:
        break;
    /* Careful using this -- it permanenty converts the thermal energy */
    case CORRECTENERGY_OUT:
        break;
    case CORRECTENERGY_SPECIAL:
        break;
    default:
        assert(0);
        break;
    }
}

void pkdDensityStep(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, int iMaxRung, double dDelta, double dEta, double dRhoFac) {
    assert(!pkd->bNoParticleOrder);
    for (auto &p : pkd->particles) {
        if (p.is_active()) {
            double dT = dEta/sqrt(p.density()*dRhoFac);
            p.set_new_rung(pkdDtToRung(dT,dDelta,iMaxRung));
        }
    }
}

uint8_t pkdDtToRung(double dT, double dDelta, uint8_t uMaxRung) {
    union {
        double d;
        struct {
            uint64_t mantisa : 52;
            uint64_t exponent : 11;
            uint64_t sign : 1;
        } ieee;
    } T;
    int iRung;
    T.d = fabs(dDelta)/dT;
    if (T.d<=1.0) return 0;
    iRung = T.ieee.exponent - 1023; /* log2(d) */
    if (iRung > uMaxRung) return uMaxRung;
    else return iRung;
}

int pkdUpdateRung(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
                  uint8_t uRung,int iMaxRung,uint64_t *nRungCount) {
    int iTempRung;
    assert(!pkd->bNoParticleOrder);
    for (auto i=0; i<iMaxRung; ++i) nRungCount[i] = 0;
    for (auto &p : pkd->particles) {
        if ( p.is_active() ) {
            if ( p.new_rung() > iMaxRung ) p.set_new_rung(iMaxRung);
            if ( p.new_rung() >= uRung ) p.set_rung(p.new_rung());
            else if ( p.rung() > uRung) p.set_rung(uRung);
        }
        /*
        ** Now produce a count of particles in rungs.
        */
        nRungCount[p.rung()] += 1;
    }
    iTempRung = iMaxRung;
    while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    return iTempRung;
}

void pkdDeleteParticle(PKD pkd, particleStore::ParticleReference &p) {
    /* p.iOrder = -2 - p.iOrder; JW: Not needed -- just preserve iOrder */
    int pSpecies = p.species();
    p.set_class(0.0,0.0,0,FIO_SPECIES_UNKNOWN); /* Special "DELETED" class == FIO_SPECIES_UNKNOWN */

    // IA: We copy the last particle into this position, the tree will no longer be valid!!!
    //
    // If we do something else with the tree before a reconstruction (such as finishing the particle loop), we need to be extra careful and:
    //   * check that the particle is marked
    //   * check that the particle type correspond to its placing inside the node
    //
    // Even with this, some weird bugs may appear, so seriously, be EXTRA careful!!
    //PARTICLE* lastp = pkd->Particle( pkd->Local() - 1);

    //
    //assert(!pkdIsDeleted(pkd,lastp));
    // IA: We can encounter a case where the last particle is the one being deleted; workaround:
    //while (pkdIsDeleted(pkd,lastp) || lastp==p){
    //   lastp--;
    //}

    //pkd->CopyParticle(p, lastp);
    //pkd->Local() -= 1;
    //pkd->TreeNode(ROOT)->pUpper -= 1;
    switch (pSpecies) {
    case FIO_SPECIES_DARK:
        pkd->nDark -= 1;
        break;
    case FIO_SPECIES_SPH:
        pkd->nGas -= 1;
        break;
    case FIO_SPECIES_STAR:
        pkd->nStar -= 1;
        break;
    case FIO_SPECIES_BH:
        pkd->nBH -= 1;
        break;
    default:
        printf("Deleting particle with unknown type");
        abort();
    }
    p.set_marked(false);
    p.set_rung(0);
    p.set_new_rung(0);
}

/* IA: We replace the deleted particles with those at the end of the particle array (which are still valid), making the tree
 *  no longer usable, unless *extreme* care is taken.
 *
 *  We also update the number of particles in the master level, assuming that the values in PKD are correct (i.e., pkdDeleteParticle was called correctly)
 */
void pkdMoveDeletedParticles(PKD pkd, total_t *n, total_t *nGas, total_t *nDark, total_t *nStar, total_t *nBH) {
    auto i = std::partition(pkd->particles.begin(),pkd->particles.end(),[](auto &p) {return !p.is_deleted();});
    pkd->SetLocal(i - pkd->particles.begin());
    *n  = pkd->Local();
    *nGas = pkd->nGas;
    *nDark = pkd->nDark;
    *nStar = pkd->nStar;
    *nBH = pkd->nBH;
}

void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
                  int *nDeltaStar) {
    int pi, pj;
    int nNew;
    int ndGas;
    int ndDark;
    int ndStar;
    int newnLocal;
    PARTICLE *p;

    nNew = 0;
    ndGas = 0;
    ndDark = 0;
    ndStar = 0;
    newnLocal = pkd->Local();
    for (pi = 0, pj = 0; pi < pkd->Local(); pi++) {
        p = pkd->Particle(pi);
        if (pj < pi)
            pkd->CopyParticle(pkd->Particle(pj),p);
        if (pkdIsNew(pkd,p)) {
            ++pj;
            ++nNew;
            if (pkdIsGas(pkd, p))
                ++ndGas;
            else if (pkdIsDark(pkd, p))
                ++ndDark;
            else if (pkdIsStar(pkd, p))
                ++ndStar;
            else
                mdlassert(pkd->mdl,0);
            if (pkdIsActive(pkd,p))
                ++pkd->nActive;
            continue;
        }
        else if (pkdIsDeleted(pkd,p)) {
            --newnLocal; /* no idea about type now -- type info lost */
            --ndGas; /* JW: Hack: assume only gas deleted fix this! */
            /*      if (pkdIsGas(pkd, p))
            --ndGas;
                else if (pkdIsDark(pkd, p))
              --ndDark;
                else if (pkdIsStar(pkd, p))
              --ndStar;
                else
                mdlassert(pkd->mdl,0);*/
            if (pkdIsActive(pkd,p))
                --pkd->nActive;
        }
        else {
            ++pj;
        }
    }

    *pnNew = nNew;
    *nDeltaGas = ndGas;
    *nDeltaDark = ndDark;
    *nDeltaStar = ndStar;
    pkd->SetLocal(newnLocal);
}

void pkdNewOrder(PKD pkd,int nStart) {
    PARTICLE *p;
    int pi;

    for (pi=0; pi<pkd->Local(); pi++) {
        p = pkd->Particle(pi);
        if (p->iOrder == IORDERMAX) {
            p->iOrder = nStart++;
        }
    }
}

/// @brief Count the number of each species of particle
void pkdGetNParts(PKD pkd, struct outGetNParts *out ) {
    TinyVector<int,FIO_SPECIES_LAST> counts = 0;
    total_t iMaxOrder = 0;
    // Loop over all particles and accumulate counts for each species.
    // Also keep track of the maximum iOrder found.
    std::for_each(pkd->particles.begin(),pkd->particles.end(),
    [&counts,&iMaxOrder](auto &p) {
        iMaxOrder = std::max(iMaxOrder,p.order());
        ++counts[p.species()];
    });

    out->n     = blitz::sum(counts);
    out->nGas  = counts[FIO_SPECIES_SPH];
    out->nDark = counts[FIO_SPECIES_DARK];
    out->nStar = counts[FIO_SPECIES_STAR];
    out->nBH   = counts[FIO_SPECIES_BH];
    out->nMaxOrder = iMaxOrder;

    pkdSetNParts(pkd, out->nGas, out->nDark, out->nStar, out->nBH);
}

void pkdSetNParts(PKD pkd,int nGas,int nDark,int nStar, int nBH) {
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
    pkd->nBH = nBH;
}

int pkdIsGas(PKD pkd,PARTICLE *p) {
    return pkdSpecies(pkd,p) == FIO_SPECIES_SPH;
}

int pkdIsDark(PKD pkd,PARTICLE *p) {
    return pkdSpecies(pkd,p) == FIO_SPECIES_DARK;
}

int pkdIsStar(PKD pkd,PARTICLE *p) {
    return pkdSpecies(pkd,p) == FIO_SPECIES_STAR;
}

int pkdIsBH(PKD pkd,PARTICLE *p) {
    return pkdSpecies(pkd,p) == FIO_SPECIES_BH;
}

double pkdTotalMass(PKD pkd) {
    return std::accumulate(pkd->particles.begin(),pkd->particles.end(),0.0,
    [](double a,auto &p) {return a + p.mass(); });
}

uint8_t pkdGetMinDt(PKD pkd) {
    return std::accumulate(pkd->particles.begin(),pkd->particles.end(),0,
    [](uint8_t a,auto &p) {return std::max(a,p.new_rung()); });
}

void pkdSetGlobalDt(PKD pkd, uint8_t minDt) {
    for (auto &p : pkd->particles) p.set_new_rung(minDt);
}

void pkdOutPsGroup(PKD pkd,char *pszFileName,int iType) {
    FILE *fp;
    int i;

    if (iType == OUT_PSGROUP_STATS) {
        fp = fopen(pszFileName,"a+");
        assert(fp != NULL);
        struct psGroup *gd = pkd->psGroupTable.pGroup;

        for (i=1; i<pkd->psGroupTable.nGroups; ++i) {
            if (gd[i].iPid != pkd->Self()) continue;
            fprintf(fp,"%d",gd[i].iGlobalId);
            fprintf(fp," %10" PRIu64 "",gd[i].nTotal);
            fprintf(fp," %12.8e",gd[i].fMass);
            fprintf(fp," %12.8e",gd[i].fRMSRadius);
            fprintf(fp," %12.8e",gd[i].r[0]);
            fprintf(fp," %12.8e",gd[i].r[1]);
            fprintf(fp," %12.8e",gd[i].r[2]);
            fprintf(fp," %12.8e",gd[i].v[0]);
            fprintf(fp," %12.8e",gd[i].v[1]);
            fprintf(fp," %12.8e",gd[i].v[2]);
#if 0
            fprintf(fp,"%.11g ",pkd->groupData[i].rcom[0]);
            fprintf(fp,"%.11g ",pkd->groupData[i].rcom[1]);
            fprintf(fp,"%.11g ",pkd->groupData[i].rcom[2]);
            fprintf(fp,"%.11g ",pkd->groupData[i].r[0]);
            fprintf(fp,"%.11g ",pkd->groupData[i].r[1]);
            fprintf(fp,"%.11g ",pkd->groupData[i].r[2]);
            fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[0]);
            fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[1]);
            fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[2]);
#endif
            fprintf(fp,"\n");
        }
        if (fclose(fp) == EOF) {
            perror("pkdOutPsGroup: could not close file");
            exit(1);
        }
    }
    else
        assert(0);
}

int pkdGetParticles(PKD pkd, int nIn, uint64_t *ID, struct outGetParticles *out) {
    assert(!pkd->bNoParticleOrder); /* We need particle IDs */
    int nOut = 0;
    for (auto &p : pkd->particles) {
        for (auto j=0; j<nIn; ++j) {
            if (ID[j] == p.order()) {
                out[nOut].id = ID[j];
                out[nOut].mass = p.mass();
                out[nOut].phi = p.have_potential() ? p.potential() : 0.0;
                out[nOut].r = p.position();
                out[nOut].v = p.velocity();
                ++nOut;
                break;
            }
        }
    }
    return nOut;
}

void pkdUpdateGasValues(PKD pkd, struct pkdKickParameters *kick, SPHOptions *SPHoptions) {
    int doUConversion = SPHoptions->doUConversion;
    if (SPHoptions->doUConversion) {
        for (auto &p : pkd->particles) {
            auto &NewSph = p.newsph();
            NewSph.u = SPHEOSUofRhoT(pkd,p.density(),NewSph.u,p.imaterial(),SPHoptions);
            NewSph.oldRho = p.density();
        }
    }
    if (doUConversion) SPHoptions->doUConversion = 0;
    for (auto &p : pkd->particles) {
        if (SPHoptions->useDensityFlags && p.rung() < SPHoptions->nPredictRung && !p.marked()) continue;
        auto &NewSph = p.newsph();
        SPHpredictInDensity(pkd, p, kick, SPHoptions->nPredictRung, &NewSph.P, &NewSph.cs, &NewSph.T, SPHoptions);
    }
    if (doUConversion) SPHoptions->doUConversion = 1;
}

/*
** Initialize the EOS tables
*/
void pkdInitializeEOS(PKD pkd) {
    auto materials = pkd->particles.getMaterials();
    for (auto iMat : materials) {
        if (iMat == 0 && pkd->SPHoptions.useBuiltinIdeal) {
            // Nothing to do
        }
        else {
#ifdef HAVE_EOSLIB_H
            if (pkd->materials[iMat] == NULL) {
                if (iMat == MAT_IDEALGAS) {
                    struct igeosParam param;
                    param.dConstGamma = pkd->SPHoptions.gamma;
                    param.dMeanMolMass = pkd->SPHoptions.dMeanMolWeight;
                    pkd->materials[iMat] = EOSinitMaterial(iMat, pkd->SPHoptions.dKpcUnit, pkd->SPHoptions.dMsolUnit, &param);
                    if (pkd->SPHoptions.useIsentropic) {
                        EOSinitIsentropicLookup(pkd->materials[iMat],NULL);
                    }
                }
                else {
                    pkd->materials[iMat] = EOSinitMaterial(iMat, pkd->SPHoptions.dKpcUnit, pkd->SPHoptions.dMsolUnit, NULL);
                    if (pkd->SPHoptions.useIsentropic) {
                        EOSinitIsentropicLookup(pkd->materials[iMat],NULL);
                    }
                }
            }
#else
            printf("Trying to initialize an EOSlib material, but EOSlib was not compiled in!\n");
            assert(0);
#endif
        }
    }
}

void pkdCopySPHOptionsToDevice(PKD pkd, SPHOptions *SPHoptions, int bGPU) {
    if (bGPU) {
#ifdef USE_CUDA
        auto cuda = reinterpret_cast<CudaClient *>(pkd->cudaClient);
        // Only one thread needs to transfer the SPHoptions to the GPU
        if (pkd->mdl->Core()==0) {
            SPHOptionsGPU SPHoptionsGPU;
            copySPHOptionsGPU(SPHoptions, &SPHoptionsGPU);
            cuda->setupSPHOptions(&SPHoptionsGPU);
        }
        pkd->mdl->ThreadBarrier();
#endif
    }
}
