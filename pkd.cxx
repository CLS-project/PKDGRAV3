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

#include "core/decorator.h"
#include "core/decorator_fields.h"
#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <cinttypes>
#include "io/iomodule.h"
#include "io/iostaggered.h"
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
#include "fmt/format.h"        // This will be part of c++20
using namespace fmt::literals; // Gives us ""_a and ""_format literals

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
#include "cooling/cooling.h"
#include "blackhole/evolve.h"
#ifdef GRACKLE
    #include "cooling_grackle/cooling_grackle.h"
#endif
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
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

using blitz::dot;
using blitz::max;
using blitz::min;
using blitz::TinyVector;

static void initLightBallOffsets(PKD pkd, double mrLCP) {
    Bound bnd(TinyVector<double, 3>(-0.5), TinyVector<double, 3>(0.5));
    double min2;
    int ix, iy, iz, l;

    pkd->nLayerMax = ceil(mrLCP);
    pkd->nBoxLC.clear();

    pkd->lcOffset0.clear();
    pkd->lcOffset1.clear();
    pkd->lcOffset2.clear();
    /*
    ** Set up the light cone offsets such that they proceed from the inner 8
    ** unit boxes outward layer by layer so that we can skip checks of the
    ** outer layers if we want.
    */
    for (l = 0; l < pkd->nLayerMax; ++l) {
        for (ix = -l; ix <= l + 1; ++ix) {
            for (iy = -l; iy <= l + 1; ++iy) {
                for (iz = -l; iz <= l + 1; ++iz) {
                    if ((ix > -l) && (ix < l + 1) && (iy > -l) && (iy < l + 1) && (iz > -l) && (iz < l + 1)) continue;
                    TinyVector<double, 3> r(ix - 0.5, iy - 0.5, iz - 0.5);
                    min2 = bnd.mindist(r);
                    if (min2 < mrLCP * mrLCP) {
                        if (pkd->Self() == 0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f> min2:%f\n", pkd->lcOffset0.size(), l, r[0], r[1], r[2], min2);
                        pkd->lcOffset0.push_back(r[0]);
                        pkd->lcOffset1.push_back(r[1]);
                        pkd->lcOffset2.push_back(r[2]);
                    }
                }
            }
        }
        /*
        ** Can save the nBox here for skipping the deeper parts of the lightcone checks!
        */
        pkd->nBoxLC.push_back(pkd->lcOffset0.size());
    }
}

static void initLightConeOffsets(PKD pkd, int bBowtie, blitz::TinyVector<double, 3> h, double alpha, double mrLCP) {
    int l, ix, iy, iz;
    double hm = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
    double min2;

    pkd->nLayerMax = ceil(mrLCP);
    for (int j = 0; j < 3; ++j) h[j] /= hm;
    if (pkd->Self() == 0) printf("mrLCP:%f nLayerMax:%d vector width:%d h:<%f %f %f>\n", mrLCP, pkd->nLayerMax, dvec::width(), h[0], h[1], h[2]);
    /*
    ** Check that all the components of h are positive and that h is not too close to the
    ** walls!
    */
    assert(h[0] > 0);
    assert(h[1] > 0);
    assert(h[2] > 0);

    pkd->nBoxLC.clear();

    // by default all projected cells need to be checked
    std::vector<char> xy(pkd->nLayerMax * pkd->nLayerMax, true);
    std::vector<char> xz(pkd->nLayerMax * pkd->nLayerMax, true);
    std::vector<char> yz(pkd->nLayerMax * pkd->nLayerMax, true);

    double alpha_xy = atan2(tan(alpha), sqrt(h[0] * h[0] + h[1] * h[1])); // angle projected onto the xy plane
    double alpha_xz = atan2(tan(alpha), sqrt(h[0] * h[0] + h[2] * h[2])); // angle projected onto the xz plane
    double alpha_yz = atan2(tan(alpha), sqrt(h[1] * h[1] + h[2] * h[2])); // angle projected onto the xy plane
    double beta_yx = atan2(h[0], h[1]);
    assert(beta_yx > alpha_xy); // make sure it is far enough from the walls
    double beta_zx = atan2(h[0], h[2]);
    assert(beta_zx > alpha_xz); // make sure it is far enough from the walls
    double beta_zy = atan2(h[1], h[2]);
    assert(beta_zy > alpha_yz); // make sure it is far enough from the walls

    double beta_xy = atan2(h[1], h[0]);
    assert(beta_xy > alpha_xy); // make sure it is far enough from the walls
    double beta_xz = atan2(h[2], h[0]);
    assert(beta_xz > alpha_xz); // make sure it is far enough from the walls
    double beta_yz = atan2(h[2], h[1]);
    assert(beta_yz > alpha_yz); // make sure it is far enough from the walls
    if (pkd->Self() == 0)
        printf("alpha:%f beta_xy:%f alpha_xy:%f beta_xz:%f alpha_xz:%f beta_yz:%f alpha_yz:%f\n", alpha, beta_xy, alpha_xy, beta_xz, alpha_xz, beta_yz,
               alpha_yz);
    for (int j = 0; j < pkd->nLayerMax; ++j) {
        for (int i = 0; i < pkd->nLayerMax; ++i) {
            double a1 = atan2(j, i + 1);
            double a2 = atan2(j + 1, i);
            if (beta_xy > a2 + alpha_xy || beta_xy < a1 - alpha_xy) xy[j * pkd->nLayerMax + i] = false; // doesn't intersect
            if (beta_xz > a2 + alpha_xz || beta_xz < a1 - alpha_xz) xz[j * pkd->nLayerMax + i] = false; // doesn't intersect
            if (beta_yz > a2 + alpha_yz || beta_yz < a1 - alpha_yz) yz[j * pkd->nLayerMax + i] = false; // doesn't intersect
        }
    }
    pkd->lcOffset0.clear();
    pkd->lcOffset1.clear();
    pkd->lcOffset2.clear();
    /*
    ** Set up the light cone offsets such that they proceed from the inner 8
    ** unit boxes outward layer by layer so that we can skip checks of the
    ** outer layers if we want.
    */
    for (l = 0; l < pkd->nLayerMax; ++l) {
        for (ix = 0; ix <= l; ++ix) {
            for (iy = 0; iy <= l; ++iy) {
                for (iz = 0; iz <= l; ++iz) {
                    if (ix < l && iy < l && iz < l) continue;
                    if (xy[iy * pkd->nLayerMax + ix] && xz[iz * pkd->nLayerMax + ix] && yz[iz * pkd->nLayerMax + iy]) {
                        double r[3] = {ix + 0.5, iy + 0.5, iz + 0.5};
                        Bound bnd(TinyVector<double, 3>(-0.5), TinyVector<double, 3>(0.5));
                        min2 = bnd.mindist(r);
                        if (min2 > mrLCP * mrLCP) continue;
                        if (pkd->Self() == 0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f>\n", pkd->lcOffset0.size(), l, r[0], r[1], r[2]);
                        pkd->lcOffset0.push_back(r[0]);
                        pkd->lcOffset1.push_back(r[1]);
                        pkd->lcOffset2.push_back(r[2]);
                        if (bBowtie) {
                            /*
                            ** Now include also the -,-,- octant for checks, as we will be producing a bowtie.
                            */
                            for (int j = 0; j < 3; ++j) r[j] = -r[j];
                            if (pkd->Self() == 0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f>\n", pkd->lcOffset0.size(), l, r[0], r[1], r[2]);
                            pkd->lcOffset0.push_back(r[0]);
                            pkd->lcOffset1.push_back(r[1]);
                            pkd->lcOffset2.push_back(r[2]);
                        }
                    }
                }
            }
        }
        /*
        ** Can save the nBox here for skipping the deeper parts of the lightcone checks!
        */
        pkd->nBoxLC.push_back(pkd->lcOffset0.size());
    }
    /*
    ** Add more boxes to get to a multiple of the vector width.
    ** Just make sure the added boxes far enough away to not get
    ** included.
    */
    l = pkd->nLayerMax - 1;
    double r[3] = {2 * l + 0.5, 2 * l + 0.5, 2 * l + 0.5};
    while (pkd->lcOffset0.size() % dvec::width()) {
        if (pkd->Self() == 0) printf("Lightcone replica:%zu layer:%d r:<%f %f %f>\n", pkd->lcOffset0.size(), l, r[0], r[1], r[2]);
        pkd->lcOffset0.push_back(r[0]);
        pkd->lcOffset1.push_back(r[1]);
        pkd->lcOffset2.push_back(r[2]);
    }
    pkd->nBoxLC.push_back(pkd->lcOffset0.size());
}

void pkdContext::set_factor(std::uint32_t factor) {
    assert(particles.size() == 0); // Must be called before any particles are allocated
    assert(tree.size() == 0);      // Must be called before any nodes are allocated
    Integerize::set_factor(factor);
    particles.set_factor(factor);
    tree.set_factor(factor);
}

pkdContext::pkdContext(mdl::mdlClass* mdl, int nStore, int nDecoratorStore, int nBucket, uint64_t nMinTotalStore, uint64_t nMinEphemeral,
                       uint32_t nEphemeralBytes, int nTreeBitsLo, int nTreeBitsHi, int iCacheSize, int iCacheMaxInflight, int iWorkQueueSize,
                       mdl::MPI_RECV_MODE mpi_recv_mode, int mpi_recv_count, bool mpi_RDMA, const TinyVector<double, 3>& fPeriod, uint64_t nDark, uint64_t nGas,
                       uint64_t nStar, uint64_t nBH, uint64_t mMemoryModel, uint32_t nIntegerFactor)
    : mdl(mdl), pLightCone((LIGHTCONEP*)nullptr), pHealpixData(nullptr), csm(nullptr) {
    PARTICLE* p;
    uint32_t pi;
    int j, ism;

    set_factor(nIntegerFactor);
    for (auto o = 0; o < 8; ++o) io_init(&afiLightCone[o], 0, 0, IO_AIO | IO_LIBAIO);

#define RANDOM_SEED 1
    srand(RANDOM_SEED);

    this->nDark = nDark;
    this->nGas = nGas;
    this->nBH = nBH;
    this->nStar = nStar;
    this->nRejects = 0;
    this->nDecoratorRejects = 0;
    this->fPeriod = fPeriod;

    this->uMinRungActive = 0;
    this->uMaxRungActive = 255;
    for (j = 0; j <= IRUNGMAX; ++j) this->nRung[j] = 0;

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
    this->bNoParticleOrder = (mMemoryModel & PKD_MODEL_UNORDERED) ? 1 : 0;
    this->bIntegerPosition = (mMemoryModel & PKD_MODEL_INTEGER_POS) ? 1 : 0;

    if (mMemoryModel & (PKD_MODEL_MFM | PKD_MODEL_MFV | PKD_MODEL_STAR | PKD_MODEL_BH)) {
        particles.add<decorators::FIELDS*>(PKD_FIELD::oDecorator, "decorator");
    }

    if (mMemoryModel & (PKD_MODEL_MFM | PKD_MODEL_MFV)) {
        auto gas = decorators.overlay(DECORATOR_FIELD::oMFM);
        gas.add<decorators::MFM>(DECORATOR_FIELD::oMFM, "hydro");
        if (mMemoryModel & PKD_MODEL_MFV) gas.add<decorators::MFV>(DECORATOR_FIELD::oMFV, "mfv");
    }

    this->bStars = (mMemoryModel & PKD_MODEL_STAR) ? 1 : 0;
    if (mMemoryModel & PKD_MODEL_STAR) {
        auto star = decorators.overlay(DECORATOR_FIELD::oMFM); // Overlay the star field on the MFM field
        star.add<decorators::STAR>(DECORATOR_FIELD::oStar, "star");
    }

    this->bBlackholes = (mMemoryModel & PKD_MODEL_BH) ? 1 : 0;
    if (mMemoryModel & PKD_MODEL_BH) {
        auto bh = decorators.overlay(DECORATOR_FIELD::oMFM); // Overlay the blackhole field on the MFM field
        bh.add<decorators::BLACKHOLE>(DECORATOR_FIELD::oBH, "bh");
    }

    decorators.commit();

    particles.initialize(this->bIntegerPosition, this->bNoParticleOrder, &decorators);

    if (this->bIntegerPosition) particles.add<int32_t[3]>(PKD_FIELD::oPosition, "r");
    else particles.add<double[3]>(PKD_FIELD::oPosition, "r");
    if (mMemoryModel & PKD_MODEL_VELOCITY) {
        if (sizeof(vel_t) == sizeof(double)) particles.add<double[3]>(PKD_FIELD::oVelocity, "v");
        else particles.add<float[3]>(PKD_FIELD::oVelocity, "v");
    }

    if (mMemoryModel & PKD_MODEL_PARTICLE_ID) particles.add<int64_t>(PKD_FIELD::oParticleID, "id");
    /*
    ** Add a global group id. This is used when outputing an array of particles with
    ** one group assignment per particle.Usually only used in testing as it adds a
    ** 64 bit integer.
    */
    if (mMemoryModel & PKD_MODEL_GLOBALGID) particles.add<int64_t>(PKD_FIELD::oGlobalGid, "gid");

    if (mMemoryModel & PKD_MODEL_NEW_SPH) particles.add<sph::FIELDS>(PKD_FIELD::oNewSph, "hydro");
    if (mMemoryModel & PKD_MODEL_STRENGTH) particles.add<sph::STRENGTHFIELDS>(PKD_FIELD::oNewSphStr, "strength");
    if (mMemoryModel & PKD_MODEL_PEAKPRESSURE) particles.add<sph::PEAKPRESSUREFIELDS>(PKD_FIELD::oNewSphPeakPressure, "peakpressure");

    if (mMemoryModel & PKD_MODEL_VELSMOOTH) particles.add<VELSMOOTH>(PKD_FIELD::oVelSmooth, "vsmooth");
    if (mMemoryModel & PKD_MODEL_ACCELERATION) particles.add<float[3]>(PKD_FIELD::oAcceleration, "a");

    if (mMemoryModel & PKD_MODEL_MASS) particles.add<float>(PKD_FIELD::oMass, "m");

    if (mMemoryModel & PKD_MODEL_SOFTENING) particles.add<float>(PKD_FIELD::oSoft, "soft");

    if (mMemoryModel & (PKD_MODEL_MFM | PKD_MODEL_NEW_SPH | PKD_MODEL_BALL)) particles.add<float>(PKD_FIELD::oBall, "ball");
    if (mMemoryModel & (PKD_MODEL_MFM | PKD_MODEL_NEW_SPH | PKD_MODEL_DENSITY)) particles.add<float>(PKD_FIELD::oDensity, "rho");

    if ((mMemoryModel & PKD_MODEL_GROUPS) && !this->bNoParticleOrder) {
        particles.add<int32_t>(PKD_FIELD::oGroup, "group");
    }
    if ((mMemoryModel & PKD_MODEL_NEWRUNG) && this->bNoParticleOrder) {
        // If not bNoParticleOrder (i.e., ordered) the field
        // is already present
        particles.add<uint8_t>(PKD_FIELD::oNewRung, "newrung");
    }

    if (mMemoryModel & PKD_MODEL_POTENTIAL) {
        particles.add<float>(PKD_FIELD::oPotential, "phi");
    }
    particles.commit();

    /*
    ** Tree node memory models
    */
    if (this->bIntegerPosition) tree.add<int32_t[3]>(KDN_FIELD::oNodePosition, "r");
    else tree.add<double[3]>(KDN_FIELD::oNodePosition, "r");
    if (mMemoryModel & PKD_MODEL_NODE_BND) {
        if (this->bIntegerPosition) tree.add<IntegerBound>(KDN_FIELD::oNodeBnd, "bnd");
        else tree.add<Bound>(KDN_FIELD::oNodeBnd, "bnd");
    }
    if (mMemoryModel & PKD_MODEL_NODE_VBND) tree.add<Bound>(KDN_FIELD::oNodeVBnd, "vbnd");
    if ((mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(double)) tree.add<double[3]>(KDN_FIELD::oNodeVelocity, "v");
#ifdef OPTIM_REORDER_IN_NODES
    if (mMemoryModel & (PKD_MODEL_MFM | PKD_MODEL_BH)) {
        tree.add<int32_t>(KDN_FIELD::oNodeNgas, "ngas");
        tree.add<int32_t>(KDN_FIELD::oNodeNstar, "nstar");
        tree.add<int32_t>(KDN_FIELD::oNodeNbh, "nbh");
#endif
    }
    /*
    ** Three extra bounds are required by the fast gas SPH code.
    */
    if (mMemoryModel & PKD_MODEL_NODE_SPHBNDS) tree.add<SPHBNDS>(KDN_FIELD::oNodeSphBounds, "sphbnd");

    if (mMemoryModel & PKD_MODEL_NODE_BOB) tree.add<SPHBOB>(KDN_FIELD::oNodeBOB, "bob");

    if (mMemoryModel & PKD_MODEL_NODE_MOMENT) {
        tree.add<FMOMR>(KDN_FIELD::oNodeMom, "mom");
        auto moment = tree.overlay(KDN_FIELD::oNodeMom);
        moment.add<mass_t>(KDN_FIELD::oNodeMass, "mass", offsetof(FMOMR, m));
        // tree.add<mass_t>(KDN_FIELD::oNodeMass,"mass",tree.offset(KDN_FIELD::oNodeMom) + offsetof(FMOMR,m));
    }
    else tree.add<mass_t>(KDN_FIELD::oNodeMass, "mass");

    /* The acceleration is required for the new time step criteria */
    if (mMemoryModel & PKD_MODEL_NODE_ACCEL) tree.add<float[3]>(KDN_FIELD::oNodeAcceleration, "a");

    if ((mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(float)) tree.add<float[3]>(KDN_FIELD::oNodeVelocity, "v");

    assert(tree.ElementSize() > 0);
    if (tree.ElementSize() > pkdContext::MaxNodeSize()) {
        fprintf(stderr, "Node size is too large. Node size=%" PRIu64 ", max node size=%" PRIu64 "\n", (uint64_t)tree.ElementSize(), (uint64_t)MaxNodeSize());
    }
    assert(tree.ElementSize() <= MaxNodeSize());
    tree.commit();

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
    if (nEphemeral * mdl->Cores() * nEphemeralBytes < nMinEphemeral)
        nEphemeral = (nMinEphemeral + mdl->Cores() * nEphemeralBytes - 1) / (mdl->Cores() * nEphemeralBytes);

    // For the fourth storage we use a rather conservative minimum number of tree nodes.
    // Ideally we would have nStore/nBucket leaf nodes, and thus twice that total tree nodes.
    // The buckets are not always full, and the tree is not completely balanced, so it will be more in practice.
    ephemeral.reset(); // The monotonic buffer resource will be invalidated by the new allocation
#ifdef FIXED_SIZE_TREES
    uint64_t nTreeNodes = (nStore + 1) / uint64_t(std::max(1.0, std::ceil(nBucket - std::sqrt(nBucket)))) * 3;
#else
uint64_t nTreeNodes = (nStore + 1) / nBucket * 2;
#endif
    static_assert(nSegments == 4);
    TinyVector<size_t, nSegments> nElements{uint64_t(nStore + 1), nEphemeral, uint64_t(nDecoratorStore + 1), nTreeNodes};
    TinyVector<size_t, nSegments> nBytesPerElement{particles.ParticleSize(), nEphemeralBytes, decorators.ElementSize(), tree.ElementSize()};
    storageSize = mdl->new_shared_array(pSegments.data(), 4, nElements.data(), nBytesPerElement.data(), nMinTotalStore);
    if (storageSize == 0) {
        fprintf(stderr, "ERROR: unable to allocate main storage\n");
        abort();
    }
    nSegmentSizeBytes = nElements * nBytesPerElement;

    assert(nElements[0] >= nStore);
    particles.setStore(pSegments[0], nStore);
    decorators.setStore(pSegments[2], nDecoratorStore);

    // Setup the ephemeral storage. If we run out of space, it will throw an exception.
    ephemeral.emplace(pSegments[1], nElements[1] * nBytesPerElement[1], std::pmr::null_memory_resource());

    /*
    ** Now we setup the node storage for the tree.  This storage is no longer
    ** continguous as the MDL now supports non-contiguous arrays.  We allocate
    ** a single "tile" for the tree.  If this is not sufficient, then additional
    ** tiles are allocated dynamically.  The default parameters allow for 2^32
    ** nodes total which is the integer limit anyway. We may use the extra storage
    ** from above if constraint (c) could not otherwise be met.
    */
    tree.initialize(&particles, this->bIntegerPosition, nTreeBitsLo, nTreeBitsHi, nElements[3] * nBytesPerElement[3], pSegments[3]);

    /*
    ** We also allocate a temporary particle used for swapping.  We need to do
    ** this now because the outside world can no longer know the size of a
    ** particle a priori.
    */
    this->pTempPRIVATE = static_cast<PARTICLE*>(malloc(particles.ParticleSize()));
    mdlassert(mdl, this->pTempPRIVATE != NULL);

    if (iCacheSize > 0) mdl->SetCacheSize(iCacheSize);
    mdl->SetCacheMaxInflight(iCacheMaxInflight);
    mdl->SetMpiRecvCount(mpi_recv_count, mpi_recv_mode);
    mdl->SetMpiRDMA(mpi_RDMA);

    // This is cheeserific - chooses the largest specified
#if defined(USE_CUDA)
    this->cudaClient = std::make_unique<CudaClient>(this->mdl->mpi->cuda, this->mdl->gpu);
#endif
#if defined(USE_METAL)
    this->metalClient = std::make_unique<MetalClient>(*this->mdl);
#endif
    /*
    ** Initialize global group id.
    */
    if (particles.present(PKD_FIELD::oGlobalGid)) {
        for (pi = 0; pi < (particles.capacity() + 1); ++pi) {
            p = Particle(pi);
            particles.global_gid(p) = 0;
        }
    }
    /*
    ** Ewald stuff!
    */
    this->ew.nMaxEwhLoop = 0;
    /*
    ** Allocate the stack.
    */
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
    tree[ROOT]->set_local(0, n - 1);
    particles.resize(n);
    return n;
}

/// @brief Increment the number of local particles
/// @param n number of particles to add
/// @return number of local particles
int pkdContext::AddLocal(int n) {
    return SetLocal(Local() + n);
}

pkdContext::~pkdContext() {
    int ism;

    /*
    ** Close caching space and free up nodes.
    */
    if (mdl->CacheStatus(CID_CELL)) mdl->FinishCache(CID_CELL);
    if (mdl->CacheStatus(CID_CELL2)) mdl->FinishCache(CID_CELL2);
    if (mdl->CacheStatus(CID_PARTICLE2)) mdl->FinishCache(CID_PARTICLE2);
    /*
    ** Free checklist.
    */
    if (ew.nMaxEwhLoop) {
        delete[] ewt.hx.f;
        delete[] ewt.hy.f;
        delete[] ewt.hz.f;
        delete[] ewt.hCfac.f;
        delete[] ewt.hSfac.f;
    }

    /* Only thread zero allocated this memory block  */
    mdl->ThreadBarrier();
    mdl->delete_shared_array(pSegments[0], storageSize);
    free(pTempPRIVATE);
    for (auto octant = 0; octant < 8; ++octant)
        if (pLightCone[octant]) free(pLightCone[octant]);
    if (pHealpixData) free(pHealpixData);
    for (auto o = 0; o < 8; ++o) io_free(&afiLightCone[o]);
    if (csm) {
        csmFinish(csm);
        csm = NULL;
    }
    if (interp_scale) gsl_spline_free(interp_scale);
    if (interp_accel) gsl_interp_accel_free(interp_accel);
    cooling_clean(cooling);
#ifdef STELLAR_EVOLUTION
    free(StelEvolData);
#endif
}

size_t pkdClCount(PKD pkd) {
    return pkd->cl.count() + pkd->clFreeList.count();
}

size_t pkdClMemory(PKD pkd) {
    return pkd->cl.memory() + pkd->clFreeList.memory();
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

void pkdReadFIO(PKD pkd, FIO fio, uint64_t iFirst, int nLocal, double dvFac, double dTuFac, enum FIO_SPECIES species) {
    float dummypot;
    TinyVector<double, 3> r, vel;
    TinyVector<float, ELEMENT_COUNT> metals;
    float fMass, fSoft, fDensity, u, fTimer;
    FIO_SPECIES eSpecies;
    uint64_t iParticleID;
    int i;

    mdlassert(pkd->mdl, fio != NULL);
#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Read");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif

    if (species == FIO_SPECIES_ALL) pkd->particles.clearClasses();

    if (pkd->decorators.present(DECORATOR_FIELD::oStar)) {
        /* Make sure star class established -- how do all procs know of these classes? How do we ensure they agree on the class identifiers? */
        auto p = pkd->particles[pkd->Local()];
        pkd->particles.setClass(0, 0, 0, 1.0, FIO_SPECIES_STAR, &p);
    }
    if (pkd->decorators.present(DECORATOR_FIELD::oBH)) {
        auto p = pkd->particles[pkd->Local()];
        pkd->particles.setClass(0, 0, 0, 1.0, FIO_SPECIES_BH, &p);
    }

    // Protect against uninitialized values
    fMass = 0.0f;
    fSoft = 0.0f;

    fioSeek(fio, iFirst, species);
    for (i = 0; i < nLocal; ++i) {
        auto p = pkd->particles[pkd->Local() + i];
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
        if (p.have_acceleration()) p.acceleration() = 0;
        float* pPot = p.have_potential() ? &p.potential() : &dummypot;
        p.set_group(0);

        /* Initialize New SPH fields if present */
        if (p.have_newsph()) {
            auto& NewSph = p.newsph();
            NewSph.u = NewSph.uDot = NewSph.divv = NewSph.Omega = 0.0;
            NewSph.cs = NewSph.P = NewSph.oldRho = NewSph.expImb2 = NewSph.T = 0.0;
            NewSph.vpred[0] = NewSph.vpred[1] = NewSph.vpred[2] = 0.0;
            NewSph.maxGroupSize = 32768;
        }

        if (p.have_newsphstr()) {
            auto& NewSphStr = p.newsphstr();
            NewSphStr.Sxx = NewSphStr.Syy = NewSphStr.Sxy = NewSphStr.Sxz = NewSphStr.Syz = 0.0f;
            NewSphStr.SDotxx = NewSphStr.SDotyy = NewSphStr.SDotxy = NewSphStr.SDotxz = NewSphStr.SDotyz = 0.0f;
            NewSphStr.Spredxx = NewSphStr.Spredyy = NewSphStr.Spredxy = NewSphStr.Spredxz = NewSphStr.Spredyz = 0.0f;
            NewSphStr.strainJ2 = 0.0f;
        }

        eSpecies = fioSpecies(fio);
        if (species != FIO_SPECIES_ALL) {
            assert(eSpecies == species);
        }
        switch (eSpecies) {
        case FIO_SPECIES_SPH:;
            float afSphOtherData[2];
            fioReadSph(fio, &iParticleID, r.data(), vel.data(), &fMass, &fSoft, pPot, &fDensity, &u, metals.data(), afSphOtherData);
            if (p.have_newsph()) {
                fSoft = 1.0f; // Dummy value, because the field in the file is used as hSmooth
                pkd->particles.setClass(fMass, fSoft, metals[0] >= 0.0 ? metals[0] : -metals[0], metals[0] >= 0.0 ? 1.0 : 0.0, eSpecies, &p);
            }
            else {
                pkd->particles.setClass(fMass, fSoft, 0, 1.0, eSpecies, &p);
            }
            p.set_density(fDensity);
            if (p.have_newsph()) {
                auto& NewSph = p.newsph();
                NewSph.u = -u; /* Can't do conversion until density known */
            }
            else {
                assert(dTuFac > 0.0);
                p.set_ball(2. * afSphOtherData[0]);
                if (p.have_mfm()) {
                    uint64_t decorator_index = pkd->decorators.append();
                    p.set_decorator_index(decorator_index);
                    auto& Sph = p.sph();
                    Sph.ElemMass = metals * fMass;
#ifdef HAVE_METALLICITY
                    Sph.fMetalMass = afSphOtherData[1] * fMass;
#endif
                    // If the value is negative, means that it is a temperature
                    u = (u < 0.0) ? -u * dTuFac : u;
                    Sph.Frho = 0.0;
                    Sph.Fmom = 0.0;
                    Sph.Fene = 0.0;
                    Sph.E = (u + 0.5 * dvFac * dot(vel, vel)) * fMass;
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
                    Sph.omega = 0.;
                    Sph.spheom = 0.;
                    if (p.have_mfv()) {
                        auto& mfv = p.mfv();
                        mfv.lastDrDotFrho = 0.;
                        mfv.drDotFrho = 0.;
                    }
                    // Sph.fLastBall = 0.0;
                    Sph.lastUpdateTime = -1.;
                    // Sph.nLastNeighs = 100;
                    Sph.SFR = -1.;
                    Sph.fAccFBEnergy = 0.;
                    Sph.BHAccretor.iIndex = NOT_ACCRETED;
                    Sph.BHAccretor.iPid = NOT_ACCRETED;
                    Sph.uWake = 0;
                }
            }
            break;
        case FIO_SPECIES_DARK:
            fioReadDark(fio, &iParticleID, r.data(), vel.data(), &fMass, &fSoft, pPot, &fDensity);
            pkd->particles.setClass(fMass, fSoft, 0, 1.0, eSpecies, &p);
            p.set_density(fDensity);
            break;
        case FIO_SPECIES_STAR:;
            float afStarOtherData[4];
            fioReadStar(fio, &iParticleID, r.data(), vel.data(), &fMass, &fSoft, pPot, &fDensity, metals.data(), &fTimer, afStarOtherData);
            pkd->particles.setClass(fMass, fSoft, 0, 1.0, eSpecies, &p);
            p.set_density(fDensity);
            if (p.have_star()) {
                uint64_t decorator_index = pkd->decorators.append();
                p.set_decorator_index(decorator_index);
                auto& Star = p.star();
                Star.fTimer = fTimer;
                Star.omega = 0.;
                // We avoid that star in the IC could explode
                Star.bCCSNFBDone = 1;
                Star.bSNIaFBDone = 1;
                Star.fSNEfficiency = afStarOtherData[3];
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
            fioReadBH(fio, &iParticleID, r.data(), vel.data(), &fMass, &fSoft, pPot, &fDensity, otherData, &fTimer);
            pkd->particles.setClass(fMass, fSoft, 0, 1.0, eSpecies, &p);
            if (p.have_bh()) {
                uint64_t decorator_index = pkd->decorators.append();
                p.set_decorator_index(decorator_index);
                auto& BH = p.BH();
                BH.omega = 0.;
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
            fprintf(stderr, "Unsupported particle type: %d\n", eSpecies);
            assert(0);
        }
        p.set_position(r);
        if (!pkd->bNoParticleOrder) p.set_order(iParticleID);
        if (p.have_particle_id()) p.ParticleID() = iParticleID;

        if (p.have_velocity()) {
            auto& v = p.velocity();
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

void pkdEnforcePeriodic(PKD pkd, Bound bnd) {
    for (auto& p : pkd->particles) {
        p.set_position(bnd.wrap(p.position()));
    }
}

/*
** x and y must have range [1,2) !
** returns key in the range [0,2^63-1]
*/
uint64_t hilbert2d(float x, float y) {
    uint64_t s = 0;
    uint32_t m, ux, uy, ut;

    union {
        float f;
        uint32_t u;
    } punner;

    punner.f = x;
    ux = punner.u >> 2;
    punner.f = y;
    uy = punner.u >> 2;

    m = 0x00100000;

    while (m) {
        s = s << 2;
        if (ux & m) {
            if (uy & m) {
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
            if (uy & m) {
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
uint64_t hilbert3d(float x, float y, float z) {
    uint64_t s = 0;
    uint32_t m, ux, uy, uz, ut;

    union {
        float f;
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

        if (ux & m) {
            if (uy & m) {
                if (uz & m) {
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
                if (uz & m) {
                    s |= 4;
                }
                else {
                    s |= 3;
                }
            }
        }
        else {
            if (uy & m) {
                if (uz & m) {
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
                if (uz & m) {
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
int pkdWeight(PKD pkd, int d, double fSplit, int iSplitSide, int iFrom, int iTo, int* pnLow, int* pnHigh, double* pfLow, double* pfHigh) {
    int iPart;
    double fLower, fUpper;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
        iPart = pkdLowerPart(pkd, d, fSplit, iFrom, iTo);
        *pnLow = pkd->Local() - iPart;
        *pnHigh = iPart;
    }
    else {
        iPart = pkdUpperPart(pkd, d, fSplit, iFrom, iTo);
        *pnLow = iPart;
        *pnHigh = pkd->Local() - iPart;
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
int pkdWeightWrap(PKD pkd, int d, double fSplit, double fSplit2, int iSplitSide, int iFrom, int iTo, int* pnLow, int* pnHigh) {
    int iPart;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (!iSplitSide) {
        iPart = pkdLowerPartWrap(pkd, d, fSplit, fSplit2, iFrom, iTo);
        *pnLow = iPart;
        *pnHigh = pkd->Local() - iPart;
    }
    else {
        iPart = pkdUpperPartWrap(pkd, d, fSplit, fSplit2, iFrom, iTo);
        *pnHigh = iPart;
        *pnLow = pkd->Local() - iPart;
    }
    return (iPart);
}

int pkdOrdWeight(PKD pkd, uint64_t iOrdSplit, int iSplitSide, int iFrom, int iTo, int* pnLow, int* pnHigh) {
    int iPart;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
        iPart = pkdLowerOrdPart(pkd, iOrdSplit, iFrom, iTo);
        *pnLow = pkd->Local() - iPart;
        *pnHigh = iPart;
    }
    else {
        iPart = pkdUpperOrdPart(pkd, iOrdSplit, iFrom, iTo);
        *pnLow = iPart;
        *pnHigh = pkd->Local() - iPart;
    }
    return (iPart);
}

int pkdLowerPart(PKD pkd, int d, double fSplit, int i, int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto ii = std::partition(pi, pj + 1, [d, fSplit](auto& p) { return p.position(d) >= fSplit; });
    return ii - pkd->particles.begin();
}

int pkdUpperPart(PKD pkd, int d, double fSplit, int i, int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto ii = std::partition(pi, pj + 1, [d, fSplit](auto& p) { return p.position(d) < fSplit; });
    return ii - pkd->particles.begin();
}

int pkdLowerPartWrap(PKD pkd, int d, double fSplit1, double fSplit2, int i, int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    if (fSplit1 > fSplit2) {
        auto ii = std::partition(pi, pj + 1, [d, fSplit1, fSplit2](auto& p) {
            auto r = p.position(d);
            return r < fSplit2 || r >= fSplit1;
        });
        return ii - pkd->particles.begin();
    }
    else {
        auto ii = std::partition(pi, pj + 1, [d, fSplit1, fSplit2](auto& p) {
            auto r = p.position(d);
            return r < fSplit2 && r >= fSplit1;
        });
        return ii - pkd->particles.begin();
    }
}

int pkdUpperPartWrap(PKD pkd, int d, double fSplit1, double fSplit2, int i, int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;

    if (fSplit1 > fSplit2) {
        auto ii = std::partition(pi, pj + 1, [d, fSplit1, fSplit2](auto& p) {
            auto r = p.position(d);
            return r >= fSplit2 && r < fSplit1;
        });
        return ii - pkd->particles.begin();
    }
    else {
        auto ii = std::partition(pi, pj + 1, [d, fSplit1, fSplit2](auto& p) {
            auto r = p.position(d);
            return r >= fSplit2 || r < fSplit1;
        });
        return ii - pkd->particles.begin();
    }
}

int pkdLowerOrdPart(PKD pkd, uint64_t nOrdSplit, int i, int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto split = [nOrdSplit](auto& p) { return p.order() >= nOrdSplit; };
    auto ii = std::partition(pi, pj + 1, split);
    // assert (std::all_of(pi,ii,split) && std::none_of(ii,pj+1,split));
    return ii - pkd->particles.begin();
}

int pkdUpperOrdPart(PKD pkd, uint64_t nOrdSplit, int i, int j) {
    auto pi = pkd->particles.begin() + i;
    auto pj = pkd->particles.begin() + j;
    auto split = [nOrdSplit](auto& p) { return p.order() < nOrdSplit; };
    auto ii = std::partition(pi, pj + 1, split);
    // assert (std::all_of(pi,ii,split) && std::none_of(ii,pj+1,split));
    return ii - pkd->particles.begin();
}

std::pair<int, int> pkdColRejects(PKD pkd, int nSplit) {
    int iRejects, iDecoratorRejects, i;

    mdlassert(pkd->mdl, pkd->nRejects == 0);
    mdlassert(pkd->mdl, pkd->nDecoratorRejects == 0);

    pkd->nRejects = pkd->Local() - nSplit;
    iRejects = pkd->FreeStore() - pkd->nRejects;

    pkdReorderDecorator(pkd);

    // Only for gas only sims
    // assert(nSplit==nDecoratorSplit);

    // We need to keep track of the decorator indeces of the particles to be swapped,
    // for this we use the ephemeral storage
    auto ephem = pkd->Ephemeral<bool>();
    int nDecoratorRejects = 0;

    /*
    ** Move rejects to High memory.
    */
    if (pkd->Local() != pkd->FreeStore()) {
        for (i = pkd->nRejects - 1; i >= 0; --i) {
            pkd->particles[iRejects + i] = pkd->particles[nSplit + i];
            ephem[pkd->nRejects - i - 1] = pkd->particles[iRejects + i].have_decorator();
            if (ephem[pkd->nRejects - i - 1]) nDecoratorRejects++;
            assert(pkd->particles[iRejects + i].mass() == pkd->particles[nSplit + i].mass());
            assert(pkd->particles[iRejects + i].position(0) == pkd->particles[nSplit + i].position(0));
        }
    }
    else {
        // Even if there is no particle movements, we should still fill the
        // decorator flag in the ephemeral storage
        for (i = pkd->nRejects - 1; i >= 0; --i) {
            ephem[pkd->nRejects - i - 1] = pkd->particles[nSplit + i].have_decorator();
            if (ephem[pkd->nRejects - i - 1]) nDecoratorRejects++;
        }
    }

    uint64_t nDecoratorSplit = nSplit == pkd->Local() ? pkd->DecoratorLocal() : pkd->DecoratorLocal() - nDecoratorRejects;
    pkd->nDecoratorRejects = pkd->DecoratorLocal() - nDecoratorSplit;
    iDecoratorRejects = pkd->FreeDecoratorStore() - pkd->nDecoratorRejects;

    // Same for decorators
    if (pkd->DecoratorLocal() != pkd->FreeDecoratorStore()) {
        for (i = pkd->nDecoratorRejects - 1; i >= 0; --i) {
            pkd->decorators[iDecoratorRejects + i] = pkd->decorators[nDecoratorSplit + i];
            assert(pkd->decorators[nDecoratorSplit + i].mfm().Uint == pkd->decorators[iDecoratorRejects + i].mfm().Uint);
        }
    }

    pkd->SetLocal(nSplit);
    pkd->decorators.SetLocal(nDecoratorSplit);

    return std::pair<int, int>(pkd->nRejects, pkd->nDecoratorRejects);
}

std::pair<int, int> pkdSwapRejects(PKD pkd, int idSwap) {
    size_t nBuf;
    size_t nOutBytes, nSndBytes, nRcvBytes;
    size_t nOutDecoratorBytes, nSndDecoratorBytes, nRcvDecoratorBytes, nSkipBytes;

    auto oldLocal = pkd->Local();
    auto oldDecoratorLocal = pkd->DecoratorLocal();

    if (idSwap != -1) {
        // printf("%d have %d(%d) rejects\n", pkd->Self(), pkd->nRejects, pkd->nDecoratorRejects);
        nBuf = (pkdSwapSpace(pkd)) * pkd->particles.ParticleSize();

        nOutBytes = pkd->nRejects * pkd->particles.ParticleSize();

        mdlassert(pkd->mdl, pkd->Local() + pkd->nRejects <= pkd->FreeStore());
        // Swap the particles
        auto outPart = pkd->mdl->Swap(idSwap, nBuf, pkd->Particle(pkd->Local()), nOutBytes, &nSndBytes, &nRcvBytes);

        auto nNewPart = nRcvBytes / pkd->particles.ParticleSize();
        auto nSentPart = nSndBytes / pkd->particles.ParticleSize();
        // if (!outPart)
        //     printf("%d tried to swap %d with %d, but only sent %d and received %d\n",
        //            pkd->Self(), pkd->nRejects, idSwap, nSentPart, nNewPart);

        pkd->AddLocal(nNewPart);

        if (pkd->decorators.DecoratorSize() > 0) {

            // Make sure that we have enough space for receiving the associated
            // decorators to the newly arrived particles
            int nNewDecorator = 0;
            for (auto i = oldLocal; i < pkd->Local(); i++) {
                if (pkd->particles[i].have_decorator()) {
                    nNewDecorator++;
                }
            }
            assert(nNewDecorator <= pkdSwapDecoratorSpace(pkd));
            // Also only share now the associated decorators to the particles that
            // we have just sent

            int nToBeSendDecorator = 0;
            auto ephem = pkd->Ephemeral<bool>();
            for (auto i = 0; i < nSentPart; i++) {
                // NB: We store this in reverse order so that if we need multiple
                // swaps we no not need to know how many elements were copied before.
                if (ephem[pkd->nRejects - i - 1]) {
                    nToBeSendDecorator++;
                }
            }

            // i.e., sanity check, we are not miscounting the decorators
            //        printf("%d %d \t %d %d\n", pkd->Self(), nSentPart,nToBeSendDecorator, pkd->nDecoratorRejects);
            assert(nToBeSendDecorator <= pkd->nDecoratorRejects);

            // i.e., otherwise we can not be sure that we can get all the decorators
            // assert(pkdSwapSpace(pkd) <= pkdSwapDecoratorSpace(pkd));

            nOutDecoratorBytes = pkd->nDecoratorRejects * pkd->decorators.DecoratorSize();
            nSkipBytes = (pkd->nDecoratorRejects - nToBeSendDecorator) * pkd->decorators.DecoratorSize();
            nBuf = pkdSwapDecoratorSpace(pkd) * pkd->decorators.DecoratorSize();

            auto outDecorator =
                pkd->mdl->Swap(idSwap, nBuf, pkd->Decorator(pkd->DecoratorLocal()), nOutDecoratorBytes, &nSndDecoratorBytes, &nRcvDecoratorBytes, nSkipBytes);

            auto nSentDecorator = nSndDecoratorBytes / pkd->decorators.DecoratorSize();

            // i.e., the swap was completed
            assert(nSndDecoratorBytes == nOutDecoratorBytes - nSkipBytes);
            // which should be the same as
            assert(outDecorator == 1);
            // i.e., we received the correct amount of decorators
            assert(nRcvDecoratorBytes == nNewDecorator * pkd->decorators.DecoratorSize());
            // i.e., same amount of particles and decorators was exchanged (only for gas-only sims)
            // assert(nRcvDecoratorBytes/pkd->decorators.DecoratorSize() == nRcvBytes/pkd->particles.ParticleSize());

            // printf("%d %d \t %d \t %" PRIu64 "\t %" PRIu64 "\t %" PRIu64 "\n",
            //         pkd->Self(), outPart, pkd->nRejects, nOutBytes, nSndBytes, nRcvBytes);

            pkd->decorators.Add(nRcvDecoratorBytes / pkd->decorators.DecoratorSize());

            // printf("%d (%d/%d) swapped with %d a total of %d(%d) elements, changing the count from %d(%d) to %d(%d)\n",
            //        pkd->Self(), outPart, outDecorator, idSwap, nSentPart, nSentDecorator,
            //        oldLocal, oldDecoratorLocal, pkd->Local(), pkd->DecoratorLocal());
            // assert (nSentPart == nSentPart); // only for gas-only simulations

            pkd->nDecoratorRejects -= nSentDecorator;
            // printf("%d have %d(%d) remaining rejects\n", pkd->Self(), pkd->nRejects, pkd->nDecoratorRejects);

            // Bookkeep so that newly arrived particles have the correct decorator index
            uint64_t d = oldDecoratorLocal;
            for (auto i = oldLocal; i < pkd->Local(); i++) {
                auto P = pkd->particles[i];
                if (P.have_decorator()) {
                    P.set_decorator_index(d);
                    d++;
                }
            }
            assert(d == pkd->DecoratorLocal());
        }
        pkd->nRejects -= nSentPart;
    }
    return std::pair<int, int>(pkd->nRejects, pkd->nDecoratorRejects);
}

void pkdSwapAll(PKD pkd, int idSwap) {
    size_t nBuf;
    size_t nOutBytes, nSndBytes, nRcvBytes;
    int i;
    int iBuf;

    /*
    ** Move particles to High memory.
    */
    iBuf = pkdSwapSpace(pkd);
    for (i = pkd->Local() - 1; i >= 0; --i) pkd->CopyParticle(pkd->Particle(iBuf + i), pkd->Particle(i));

    nBuf = pkd->FreeStore() * pkd->particles.ParticleSize();
    nOutBytes = pkd->Local() * pkd->particles.ParticleSize();
    pkd->mdl->Swap(idSwap, nBuf, pkd->particles, nOutBytes, &nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl, nSndBytes / pkd->particles.ParticleSize() == pkd->Local());
    pkd->SetLocal(nRcvBytes / pkd->particles.ParticleSize());
    /*
    ** Repeat for decorator.
    */
    if (pkd->decorators.DecoratorSize() == 0) return;

    iBuf = pkdSwapDecoratorSpace(pkd);
    for (i = pkd->DecoratorLocal() - 1; i >= 0; --i) pkd->CopyDecorator(pkd->Decorator(iBuf + i), pkd->Decorator(i));
    nBuf = pkd->FreeDecoratorStore() * pkd->decorators.DecoratorSize();
    nOutBytes = pkd->DecoratorLocal() * pkd->decorators.DecoratorSize();
    pkd->mdl->Swap(idSwap, nBuf, pkd->decorators, nOutBytes, &nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl, nSndBytes / pkd->decorators.DecoratorSize() == pkd->DecoratorLocal());
    pkd->decorators.SetLocal(nRcvBytes / pkd->decorators.DecoratorSize());

    // Handle the newly arrived pointers to the decorators
    // Assume that they were ordered
    uint64_t nDecorator = 0;
    for (auto i = 0; i < pkd->Local(); ++i) {
        auto p = pkd->particles[i];
        if (p.have_decorator()) {
            p.set_decorator_index(nDecorator);
            nDecorator++;
        }
    }
    assert(nDecorator == pkd->DecoratorLocal());
}

int pkdSwapSpace(PKD pkd) {
    return (pkd->FreeStore() - pkd->Local());
}
int pkdSwapDecoratorSpace(PKD pkd) {
    return (pkd->FreeDecoratorStore() - pkd->DecoratorLocal());
}

int pkdActive(PKD pkd) {
    return (pkd->nActive);
}

int pkdInactive(PKD pkd) {
    return (pkd->Local() - pkd->nActive);
}

std::pair<int, int> pkdColOrdRejects(PKD pkd, uint64_t nOrdSplit, int iSplitSide) {
    int nSplit;
    if (iSplitSide) nSplit = pkdLowerOrdPart(pkd, nOrdSplit, 0, pkd->Local() - 1);
    else nSplit = pkdUpperOrdPart(pkd, nOrdSplit, 0, pkd->Local() - 1);
    return pkdColRejects(pkd, nSplit);
}

#ifndef NEW_REORDER
void pkdLocalOrder(PKD pkd, uint64_t iMinOrder, uint64_t iMaxOrder) {
    int i;
    assert(pkd->Local() == iMaxOrder - iMinOrder + 1);
    for (i = 0; i < pkd->Local(); ++i) {
        auto p1 = pkd->particles[i];

        assert(p1.order() >= iMinOrder && p1.order() <= iMaxOrder);
        while (p1.order() - iMinOrder != i) {
            auto p2 = pkd->particles[p1.order() - iMinOrder];
            swap(p1, p2);
        }
    }
    /* Above replaces: qsort(pkd->ParticleBase(),pkd->Local(),pkd->particles.ParticleSize(),cmpParticles); */
}
#endif

#define MAX_IO_BUFFER_SIZE (8 * 1024 * 1024)

void pkdCheckpoint(PKD pkd, const char* fname, const char* fname_decorator) {
    // Order particles by type to prepare for checkpoint
    auto ii = std::partition(pkd->particles.begin(), pkd->particles.end(), [](auto& p) { return !p.have_decorator(); });
    auto split = ii - pkd->particles.begin();
    pkdReorderDecorator(pkd);

    asyncFileInfo info;
    size_t nFileSize;
    int fd;
    io_init(&info, IO_MAX_ASYNC_COUNT, 0, IO_AIO | IO_LIBAIO);
    fd = io_create(&info, fname);
    if (fd < 0) {
        perror(fname);
        abort();
    }
    nFileSize = pkd->particles.ParticleSize() * split;
    char* pBuffer = (char*)pkd->particles.Element(0);
    while (nFileSize) {
        size_t count = nFileSize > MAX_IO_BUFFER_SIZE ? MAX_IO_BUFFER_SIZE : nFileSize;
        io_write(&info, pBuffer, count);
        pBuffer += count;
        nFileSize -= count;
    }
    io_close(&info);

    // Repeat for particles with decorators
    auto nElements = pkd->Local() - split;
    for (auto i = split; i < pkd->Local(); i++) assert(pkd->particles[i].have_decorator());
    auto pDecorators = static_cast<void*>(pkd->decorators.Element(0));
    auto pParticles = static_cast<void*>(pkd->particles.Element(split));
    io_staggered_write(fname_decorator, pParticles, pDecorators, nElements, pkd->particles.ParticleSize(), pkd->decorators.DecoratorSize());
}

/*****************************************************************************\
* Write particles received from another node
\*****************************************************************************/

static void writeParticle(PKD pkd, FIO fio, double dvFac, double dvFacGas, Bound bnd, particleStore::Particle& p) {
    TinyVector<double, 3> v, r;
    TinyVector<float, ELEMENT_COUNT> metals;
    float fTimer;
    uint64_t iParticleID;

    float fPot = p.have_potential() ? p.potential() : 0;
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
    assert(all(r >= bnd.lower() && r < bnd.upper()));

    /* IA: In the case of cosmological boxes, it is typical to have a box not centered on the origin.
     *  Such boxes are defined in the (+,+,+) octant. To convert to that system of reference,
     *  we just add half the box size to all positions.
     *
     *  Another option would be to change fCenter, but I do not if that could break other parts of the code...
     *  15/06/19: I have stoped using this as this can lead problems with mpi because we do not
     *  have available pkd->csm. This should be easy to fix for all halo-finder and post-process routines.
     */
    // if (bComove){
    //    for (j=0;j<3;++j){
    //      r[j] += bnd->fMax[j];
    //    }
    // }
    switch (p.species()) {
    case FIO_SPECIES_SPH:
        if (p.have_newsph()) {
            const auto& NewSph = p.newsph();
            assert(pkd->SPHoptions.TuFac > 0.0f);
            double T;
            float otherData[3];
            otherData[0] = otherData[1] = otherData[2] = 0.0f;
            T = SPHEOSTofRhoU(pkd, fDensity, NewSph.u, p.imaterial(), &pkd->SPHoptions);
            metals = 0.0f;
            metals[0] = p.imaterial();
            fSoft = p.ball() / 2.0f;
            fioWriteSph(fio, iParticleID, r.data(), v.data(), fMass, fSoft, fPot, fDensity, T, metals.data(), 0.0f, T, otherData);
        }
        else {
            assert(p.have_mfm());
            auto& Sph = p.sph();
            {
                float temperature = 0;
                if (pkd->cooling) {
                    const double dRedshift = dvFacGas - 1.;
                    temperature = cooling_get_temperature(pkd, dRedshift, pkd->cooling, p, &Sph);
                }
#ifdef GRACKLE
                gr_float fMetalDensity = Sph.fMetalMass * Sph.omega;
                gr_float fSpecificUint = Sph.Uint / fMass;

                // Set field arrays.
                pkd->grackle_field->density[0] = fDensity;
                pkd->grackle_field->internal_energy[0] = fSpecificUint;
                pkd->grackle_field->x_velocity[0] = 1.; // Velocity input is not used
                pkd->grackle_field->y_velocity[0] = 1.;
                pkd->grackle_field->z_velocity[0] = 1.;
                // for metal_cooling = 1
                pkd->grackle_field->metal_density[0] = fMetalDensity;

                int err;
                gr_float temperature;
                err = local_calculate_temperature(pkd->grackle_data, pkd->grackle_rates, pkd->grackle_units, pkd->grackle_field, &temperature);
                if (err == 0) fprintf(stderr, "Error in calculate_temperature.\n");
#endif

                metals = Sph.ElemMass / fMass;

                // TODO: if is present, then use it, otherwise -1 to signal
                // the IO to not save it
                float SFR = Sph.SFR;

                float ph = 0.5 * p.ball();
                // TODO: do something prettier
                float otherData[3];
                otherData[0] = SFR;
                // Casting integers to floats will become a problem if the number of groups
                // reaches 2^24, but for the moment we have to live with it
                otherData[1] = p.have_global_gid() ? p.global_gid() : -1.;
#ifdef HAVE_METALLICITY
                otherData[2] = Sph.fMetalMass / fMass;
#endif
                // printf("writeParticle %e\n", Sph.Uint);

                fioWriteSph(fio, iParticleID, r.data(), v.data(), fMass, fSoft, fPot, fDensity, Sph.Uint / fMass, metals.data(), ph, temperature,
                            &otherData[0]);
            }
        }
        break;
    case FIO_SPECIES_DARK: {
        float otherData[2];
        otherData[0] = p.have_global_gid() ? p.global_gid() : -1.;
        fioWriteDark(fio, iParticleID, r.data(), v.data(), fMass, fSoft, fPot, fDensity, &otherData[0]);
    } break;
    case FIO_SPECIES_STAR: {
        auto& Star = p.star();
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
        otherData[5] = Star.fSNEfficiency;
        fioWriteStar(fio, iParticleID, r.data(), v.data(), fMass, fSoft, fPot, fDensity, metals.data(), &otherData[0]);
    } break;
    case FIO_SPECIES_BH: {
        const auto& BH = p.BH();
        fTimer = BH.fTimer;
        float otherData[6];
        otherData[0] = BH.dInternalMass;
        otherData[1] = BH.dAccretionRate;
        otherData[2] = BH.dEddingtonRatio;
        otherData[3] = BH.dFeedbackRate;
        otherData[4] = BH.dAccEnergy;
        otherData[5] = p.have_global_gid() ? p.global_gid() : -1.;
        fioWriteBH(fio, iParticleID, r.data(), v.data(), fMass, fSoft, fPot, fDensity, otherData, fTimer);
    } break;
    case FIO_SPECIES_UNKNOWN:
        break;
    default:
        fprintf(stderr, "Unsupported particle type: %d\n", p.species());
        assert(0);
    }
}

void pkdWriteFromNode(PKD pkd, int iNode, FIO fio, double dvFac, double dTuFac, Bound bnd, fioSpeciesList nWritten) {
#ifdef MPI_VERSION
    double dvFacGas = std::sqrt(dvFac);

    pkd->mdl->Recv(iNode, [&](std::span<std::byte> buffer) -> size_t {
        if (buffer.empty()) return 0;

        auto* header = reinterpret_cast<uint64_t*>(buffer.data());

        uint64_t nPart = header[0];
        uint64_t nDecorator = header[1];

        std::byte* pStart = reinterpret_cast<std::byte*>(header + 2);
        std::byte* dStart = pStart + nPart * pkd->particles.ParticleSize();
        assert(2 * sizeof(uint64_t) + nPart * pkd->particles.ParticleSize() + nDecorator * pkd->decorators.DecoratorSize() == buffer.size());

        auto p = &pkd->particles[reinterpret_cast<PARTICLE*>(pStart)];
        auto d = &pkd->decorators[reinterpret_cast<decorators::FIELDS*>(dStart)];

        for (uint64_t i = 0; i < nPart; ++i, ++p) {
            particleStore::Particle& P = *p;

            if (P.have_decorator()) {
                P.set_decorator(d);
                ++d;
            }

            writeParticle(pkd, fio, dvFac, dvFacGas, bnd, P);

            ++nWritten[p->species()];
        }

        return 1;
    });
#endif
}

void pkdWriteHeaderFIO(PKD pkd, FIO fio, double dScaleFactor, double dTime, fioSpeciesList nTotal, fioSpeciesList nThisFile, double dBoxSize, double h,
                       int nProcessors, UNITS units) {
    char version[] = PACKAGE_VERSION;
    fioSetAttr(fio, HDF5_HEADER_G, "PKDGRAV version", FIO_TYPE_STRING, 1, version);
    fioSetAttr(fio, HDF5_HEADER_G, "Time", FIO_TYPE_DOUBLE, 1, &dTime);
    if (pkd->csm->val.bComove) {
        double z = 1. / dScaleFactor - 1.;
        fioSetAttr(fio, 0, "Redshift", FIO_TYPE_DOUBLE, 1, &z);
    }
    int flag;

    flag = (fio->mFlags & FIO_FLAG_DOUBLE_POS) && (fio->mFlags & FIO_FLAG_DOUBLE_VEL);
    fioSetAttr(fio, HDF5_HEADER_G, "Flag_DoublePrecision", FIO_TYPE_INT, 1, &flag);

    fioSetAttr(fio, HDF5_HEADER_G, "BoxSize", FIO_TYPE_DOUBLE, 1, &dBoxSize);
    fioSetAttr(fio, HDF5_HEADER_G, "NumFilesPerSnapshot", FIO_TYPE_INT, 1, &nProcessors);

    // Prepare the particle information in tables
    // Header compatible with Gadget-4, with NumPart_ThisFile and NumPart_Total as uint64,
    // no NumPart_Total_HighWord. This limits the number of particles at 2**64-1
    fioSpeciesList nPart;

    nPart[0] = nThisFile[FIO_SPECIES_SPH];
    nPart[1] = nThisFile[FIO_SPECIES_DARK];
    nPart[2] = 0;
    nPart[3] = 0;
    nPart[4] = nThisFile[FIO_SPECIES_STAR];
    nPart[5] = nThisFile[FIO_SPECIES_BH];
    fioSetAttr(fio, HDF5_HEADER_G, "NumPart_ThisFile", FIO_TYPE_UINT64, 6, &nPart[0]);

    nPart[0] = nTotal[FIO_SPECIES_SPH];
    nPart[1] = nTotal[FIO_SPECIES_DARK];
    nPart[2] = 0;
    nPart[3] = 0;
    nPart[4] = nTotal[FIO_SPECIES_STAR];
    nPart[5] = nTotal[FIO_SPECIES_BH];
    fioSetAttr(fio, HDF5_HEADER_G, "NumPart_Total", FIO_TYPE_UINT64, 6, &nPart[0]);

    double massTable[6] = {0, 0, 0, 0, 0, 0};
    // This is not yet fully supported, as the classes do not have to match the
    //  six available particle types.
    // However, we add this in the header so it can be parsed by other tools
    fioSetAttr(fio, HDF5_HEADER_G, "MassTable", FIO_TYPE_DOUBLE, 6, &massTable[0]);

    float fSoft = pkd->particles[0].soft(); // we take any particle
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
auto pkdPackArray(PKD pkd, int iSize, void* vBuff, int* piIndex, int n, PKD_FIELD field, int iUnitSize, double dvFac, int bMarked) -> std::byte* {
    auto pBuff = static_cast<std::byte*>(vBuff);
    int iIndex = *piIndex;

    while (iIndex < n && iSize >= iUnitSize) {
        PARTICLE* p = pkd->Particle(iIndex++);
        if (bMarked && !p->bMarked) continue;
        if (field == PKD_FIELD::oPosition) {
            auto& d = *reinterpret_cast<TinyVector<double, 3>*>(pBuff);
            d = pkd->particles.position(p);
        }
        else if (field == PKD_FIELD::oVelocity) {
            auto& V = *reinterpret_cast<TinyVector<float, 3>*>(pBuff);
            V = pkd->particles.velocity(p) * dvFac;
        }
        else {
            const auto& src = pkd->particles.get<std::byte>(p, field);
            memcpy(pBuff, &src, iUnitSize);
        }
        pBuff += iUnitSize;
        iSize -= iUnitSize;
    }
    *piIndex = iIndex;
    return pBuff;
}

void pkdSendArray(PKD pkd, int iNode, PKD_FIELD field, int iUnitSize, double dvFac, int bMarked) {
#ifdef MPI_VERSION
    int iIndex = 0;
    pkd->mdl->Send(iNode, [&](std::span<std::byte> buffer) -> size_t {
        auto pStart = buffer.data();
        auto pEnd = pkdPackArray(pkd, buffer.size(), pStart, &iIndex, pkd->Local(), field, iUnitSize, dvFac, bMarked);
        return static_cast<size_t>(pEnd - pStart);
    });
#endif
}

/*****************************************************************************\
* Receive an array/vector from a specified node
\*****************************************************************************/

auto pkdRecvArray(PKD pkd, int iNode, void* pDest, int /*iUnitSize*/) -> void* {
#ifdef MPI_VERSION
    auto* pOut = static_cast<std::byte*>(pDest);
    pkd->mdl->Recv(iNode, [&](std::span<std::byte> buffer) -> size_t {
        if (buffer.empty()) return 0;
        std::memcpy(pOut, buffer.data(), buffer.size());
        pOut += buffer.size();
        return 1;
    });
    return pOut;
#else
    return pDest;
#endif
}

/*****************************************************************************\
*
\*****************************************************************************/

uint32_t pkdWriteFIO(PKD pkd, FIO fio, double dvFac, double dTuFac, Bound bnd, fioSpeciesList nWritten) {
    double dvFacGas = sqrt(dvFac);
    for (auto& p : pkd->particles) {
        writeParticle(pkd, fio, dvFac, dvFacGas, bnd, p);
        ++nWritten[p.species()];
    }
    return pkd->particles.size();
}

void pkdSetSoft(PKD pkd, double dSoft) {
    pkd->particles.SetSoft(dSoft);
}

void pkdSetSmooth(PKD pkd, double dSmooth) {
    for (auto& p : pkd->particles) {
        if ((p.is_gas() || p.is_star() || p.is_bh()) && p.ball() == 0.0) p.set_ball(dSmooth);
    }
}

void pkdPhysicalSoft(PKD pkd, double dSoftMax, double dFac, int bSoftMaxMul) {
    pkd->particles.PhysicalSoft(dSoftMax, dFac, bSoftMaxMul);
}

class SetMarkedHelper : public particleCACHEhelper {
    virtual void combine(void* dst, const void* src, const void* key) override final {
        auto p1 = static_cast<PARTICLE*>(dst);
        auto p2 = static_cast<const PARTICLE*>(src);
        if (p2->bMarked) p1->bMarked = 1;
#ifdef NN_FLAG_IN_PARTICLE
        if (p2->bNNflag) p1->bNNflag = 1;
#endif
    }

  public:
    explicit SetMarkedHelper(PKD pkd) : particleCACHEhelper(pkd, true) {}
};

void extensiveMarkerTest(PKD pkd, struct pkdTimestepParameters* ts, SPHOptions* SPHoptions) {
    std::stack<std::pair<int, int>> cellStack;

    // Add the toptree cell corresponding to the ROOT cell to the stack
    cellStack.push(std::make_pair(pkd->iTopTree[ROOT], pkd->Self()));
    int nParticles = 0;

    while (!cellStack.empty()) {
        auto [iCell, id] = cellStack.top();
        cellStack.pop();
        auto c = (id == pkd->Self()) ? pkd->tree[iCell] : pkd->tree[pkd->mdl->Fetch<KDN>(CID_CELL, iCell, id)];

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
        for (auto pj = c->lower(); pj <= c->upper(); ++pj) {
            auto p = (id == pkd->Self()) ? pkd->particles[pj] : pkd->particles[pkd->mdl->Fetch<PARTICLE>(CID_PARTICLE, pj, id)];
            ++nParticles;

            // Shortcut if necessary flag already set, as the test will pass in the end
            if ((SPHoptions->doSetDensityFlags && p.marked()) || (SPHoptions->doSetNNflags && p.NN_flag())) continue;

            // Get position and ball size of particle p
            auto pr = p.position();
            float fBallFactor = (SPHoptions->dofBallFactor) ? SPHoptions->fBallFactor : 1.0f;
            float pBall2 = std::min(SPHoptions->ballSizeLimit, p.ball() * fBallFactor);
            pBall2 *= pBall2;

            int needsCheck = 0;
            // Loop over all particles in my root
            auto rootc = pkd->tree[ROOT];
            for (auto qj = rootc->lower(); qj <= rootc->upper(); ++qj) {
                auto q = pkd->particles[qj];

                // Skip those that are not active
                if (!q.is_rung_range(ts->uRungLo, ts->uRungHi) && !(SPHoptions->useDensityFlags && q.marked()) && !(SPHoptions->useNNflags && q.NN_flag()))
                    continue;

                // Get position and ball size of particle q
                auto qr = q.position();
                float qBall2 = std::min(SPHoptions->ballSizeLimit, q.ball() * fBallFactor);
                qBall2 *= qBall2;

                // Calculate distance squared between particle p and q
                auto dist = pr - qr;
                float dist2 = dot(dist, dist);

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
        }
    }
    assert(nParticles == pkd->nGas);
}

void pkdGravAll(PKD pkd, struct pkdKickParameters* kick, struct pkdLightconeParameters* lc, struct pkdTimestepParameters* ts, double dTime, int nReps,
                int bPeriodic, int bGPU, int bEwald, int iRoot1, int iRoot2, double fEwCut, double fEwhCut, double dThetaMin, SPHOptions* SPHoptions,
                uint64_t* pnActive, double* pdPart, double* pdPartNumAccess, double* pdPartMissRatio, double* pdCell, double* pdCellNumAccess,
                double* pdCellMissRatio, double* pdFlop, uint64_t* pnRung) {

    double dActive;
    double dPartSum;
    double dCellSum;
    int i;

#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Gravity");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif

    /*
    ** Clear all the rung counters to be safe.
    */
    for (i = 0; i <= IRUNGMAX; ++i) pkd->nRung[i] = 0;

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    pkd->mdl->TimeReset();
#endif

    /*
    ** Set up Ewald tables and stuff.
    */
    if (bPeriodic && bEwald && SPHoptions->doGravity) {
        pkdEwaldInit(pkd, nReps, fEwCut, fEwhCut, bGPU); /* ignored in Flop count! */
    }
    pkdCopySPHOptionsToDevice(pkd, SPHoptions, bGPU);
    /*
    ** Start particle caching space (cell cache already active).
    */
    if (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) {
        pkd->mdl->CacheInitialize(CID_PARTICLE, std::make_shared<SetMarkedHelper>(pkd));
    }
    else {
        pkd->mdl->CacheInitialize(CID_PARTICLE, std::make_shared<particleCACHEhelper>(pkd));
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

    *pnActive = pkdGravWalk(pkd, kick, lc, ts, dTime, nReps, bPeriodic && bEwald, bGPU, iRoot1, iRoot2, 0, dThetaMin, pdFlop, &dPartSum, &dCellSum, SPHoptions);

    assert(pkd->nWpPending == 0);

    if (SPHoptions->doExtensiveILPTest && (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags)) {
        pkd->mdl->FlushCache(CID_PARTICLE);
        pkd->mdl->CacheBarrier(CID_PARTICLE);
        extensiveMarkerTest(pkd, ts, SPHoptions);
    }

    dActive = (double)(*pnActive);
    if (*pnActive) {
        *pdPart = dPartSum / dActive;
        *pdCell = dCellSum / dActive;
    }
    else {
        assert(dPartSum == 0 && dCellSum == 0);
        *pdPart = 0; /* for the statistics we don't count this processor, see pstGravity(). */
        *pdCell = 0;
    }
    /*
    ** Get caching statistics.
    */
    if (*pnActive) {
        *pdCellNumAccess = pkd->mdl->NumAccess(CID_CELL) / dActive;
        *pdPartNumAccess = pkd->mdl->NumAccess(CID_PARTICLE) / dActive;
    }
    else {
        *pdCellNumAccess = 0;
        *pdPartNumAccess = 0;
    }
    *pdCellMissRatio = 100.0 * pkd->mdl->MissRatio(CID_CELL);     /* as a percentage */
    *pdPartMissRatio = 100.0 * pkd->mdl->MissRatio(CID_PARTICLE); /* as a percentage */
    /*
    ** Output flops count in GFlops!
    */
    *pdFlop *= 1e-9;
    /*
    ** Stop particle caching space.
    */
    pkd->mdl->FinishCache(CID_PARTICLE);

    for (i = 0; i <= IRUNGMAX; ++i) pnRung[i] = pkd->nRung[i];

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
}

/* This became easier; we already calculated these values on the fly */
void pkdCalcEandL(PKD pkd, double& T, double& U, double& Eth, TinyVector<double, 3>& L, TinyVector<double, 3>& F, double& W) {
    T = pkd->dEnergyT;
    U = pkd->dEnergyU;
    W = pkd->dEnergyW;
    L = pkd->dEnergyL;
    F = pkd->dEnergyF;
    Eth = 0.0;
    if (pkd->particles.present(PKD_FIELD::oSph)) {
        for (auto& p : pkd->particles) {
            if (p.is_gas()) Eth += p.sph().E;
        }
    }
}

/*
** Drift particles whose Rung falls between uRungLo (large step) and uRungHi (small step) inclusive,
** and those whose destination activity flag is set.
**
** Note that the drift funtion no longer wraps the particles around the periodic "unit" cell. This is
** now done by Domain Decomposition only.
*/
void pkdDrift(PKD pkd, int iRoot, double dTime, double dDelta, double dDeltaVPred, double dDeltaTime, int bDoGas, int bGasEvolveDensity) {
    int i;
    TinyVector<double, 3> rfinal, r0, dr;
    int pLower, pUpper;

    if (iRoot >= 0) {
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
        if (pkd->decorators.present(DECORATOR_FIELD::oMFM)) {
            assert(pkd->particles.present(PKD_FIELD::oAcceleration));
            for (i = pLower; i <= pUpper; ++i) {
                auto p = pkd->particles[i];
                const auto& v = p.velocity();

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
                dMin = min(dMin, rfinal);
                dMax = max(dMax, rfinal);
            }
        }
        else {
            assert(pkd->particles.present(PKD_FIELD::oNewSph));
            for (i = pLower; i <= pUpper; ++i) {
                auto p = pkd->particles[i];
                const auto& v = p.velocity();
                if (bGasEvolveDensity && p.is_gas()) {
                    auto& NewSph = p.newsph();
                    // float dfBalldt = 1.0f / 3.0f * p.ball() * p.density() * NewSph.divv;
                    // p.set_ball(p.ball() + dDelta * dfBalldt);
                    p.set_density(p.density() - dDelta * p.density() * NewSph.divv);
                }

                r0 = p.position();
                p.set_position(rfinal = r0 + dDelta * v);
                assert(isfinite(rfinal[0]));
                assert(isfinite(rfinal[1]));
                assert(isfinite(rfinal[2]));
                dMin = min(dMin, rfinal);
                dMax = max(dMax, rfinal);
            }
        }
    }
    else {
        for (i = pLower; i <= pUpper; ++i) {
            auto p = pkd->particles[i];
            const auto& v = p.velocity();
            r0 = p.position();
            p.set_position(rfinal = r0 + dDelta * v);
            dMin = min(dMin, rfinal);
            dMax = max(dMax, rfinal);
        }
    }
    pkd->bnd = Bound(dMin, dMax);
    mdlDiag(pkd->mdl, "Out of pkdDrift\n");
}

void pkdReorderDecorator(PKD pkd) {
    auto ephem = pkd->Ephemeral<uint32_t>();

    // The ephemeral storage is populated as a map between decorators and their
    // particles. This way, there is no need to do complex searches later.
    for (uint32_t i = 0; i < pkd->Local(); ++i) {
        if (pkd->particles[i].have_decorator()) {
            ephem[pkd->particles[i].decorator_index()] = i;
        }
    }

    uint64_t iDecorator = 0;
    for (uint32_t i = 0; i < pkd->Local(); ++i) {
        if (pkd->particles[i].have_decorator()) {
            auto p = pkd->particles[i];
            auto prev_decor = p.decorator_field();
            auto new_decor = &pkd->decorators[iDecorator];
            if (prev_decor != new_decor) {
                // As we built the map before, we can simply get the particle
                // whose decorator is in the assigned position of p's decorator
                auto particle_to_swap = pkd->particles[ephem[iDecorator]];

                // Swap the decorators
                auto a = pkd->decorators[prev_decor];
                auto b = pkd->decorators[new_decor];
                swap(a, b);

                // Update the mapping
                ephem[p.decorator_index()] = ephem[iDecorator];

                // Update the particles' decorators pointers
                particle_to_swap.set_decorator(prev_decor);
                p.set_decorator(new_decor);
            }
            assert(p.decorator_index() == iDecorator);
            iDecorator++;
        }
    }
    // uint64_t iDecoratorCheck = 0;
    // for (uint32_t i=0; i<pkd->Local(); ++i) {
    //     if (pkd->particles[i].have_decorator()) {
    //         assert(pkd->particles[i].decorator_index()==iDecoratorCheck);
    //         iDecoratorCheck++;
    //     }
    // }
    return;
}

#ifdef OPTIM_REORDER_IN_NODES
/// @brief Order particles in buckets by type
/// @param pkd
/// After this function is complete, particles in buckets will be in the following order:
///     Gas, Star, Blackhole, Dark
void pkdReorderWithinNodes(PKD pkd) {
    for (int i = NRESERVED_NODES; i < pkd->Nodes(); i++) {
        auto node = pkd->tree[i];
        if (node->is_bucket()) { // We are in a bucket
            auto start = node->begin();
            auto i = std::partition(start, node->end(), [](auto& p) { return p.is_gas(); });
            node->Ngas() = i - start;
            start = i;
            // We perform another swap, just to have nGas->nStar->DM
            if (pkd->bStars) {
                i = std::partition(start, node->end(), [](auto& p) { return p.is_star(); });
                node->Nstar() = i - start;
                start = i;
            }
            if (pkd->bBlackholes) {
                i = std::partition(start, node->end(), [](auto& p) { return p.is_bh(); });
                node->Nbh() = i - start;
            }
        }
    }
    if (pkd->mdl->CacheStatus(CID_CELL)) pkd->mdl->FinishCache(CID_CELL);
    pkd->mdl->CacheInitialize(CID_CELL, std::make_shared<treeCACHEhelper>(pkd));
}
#endif

void pkdSetupInterpScale(PKD pkd, double dBoxSize, double mrMax) {
    const int nTable = 1000;
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    double dr, rt[nTable], at_inv[nTable];
    /*
    ** Setup lookup table.
    */
    dr = mrMax / (nTable - 1);
    for (int i = 0; i < nTable; ++i) {
        rt[i] = i * dr;
        at_inv[i] = 1.0 / csmComoveLookbackTime2Exp(pkd->csm, rt[i] / dLightSpeed);
    }
    pkd->interp_accel = gsl_interp_accel_alloc();
    pkd->interp_scale = gsl_spline_alloc(gsl_interp_cspline, nTable);
    gsl_spline_init(pkd->interp_scale, rt, at_inv, nTable);
}

/*
 * Stripped down versions of routines from master.c
 */
void pkdKickKDKOpen(PKD pkd, double dTime, double dDelta, uint8_t uRungLo, uint8_t uRungHi) {
    if (pkd->csm->val.bComove) {
        dDelta = csmComoveKickFac(pkd->csm, dTime, dDelta);
    }
    pkdKick(pkd, dTime, dDelta, 0, 0, 0, 0, uRungLo, uRungHi);
}

void pkdKickKDKClose(PKD pkd, double dTime, double dDelta, uint8_t uRungLo, uint8_t uRungHi) {
    if (pkd->csm->val.bComove) {
        dDelta = csmComoveKickFac(pkd->csm, dTime, dDelta);
    }
    pkdKick(pkd, dTime, dDelta, 0, 0, 0, 0, uRungLo, uRungHi);
}

void pkdKick(PKD pkd, double dTime, double dDelta, int bDoGas, double dDeltaVPred, double dDeltaU, double dDeltaUPred, uint8_t uRungLo, uint8_t uRungHi) {
    assert(pkd->particles.present(PKD_FIELD::oVelocity));
    assert(pkd->particles.present(PKD_FIELD::oAcceleration));

    if (bDoGas) { // pkd->param.bMeshlessHydro
        for (auto& p : pkd->particles) {
            if (!p.is_gas()) {
                if (p.is_rung_range(uRungLo, uRungHi)) {
                    auto& a = p.acceleration();
                    auto& v = p.velocity();
                    v += a * dDelta;
                }
            }
        }
    }
    else {
        for (auto& p : pkd->particles) {
            if (p.is_rung_range(uRungLo, uRungHi)) {
                auto& a = p.acceleration();
                auto& v = p.velocity();
                v += a * dDelta;
            }
        }
    }

    mdlDiag(pkd->mdl, "Done pkdkick\n");
}

/* Kick the tree at iRoot. */
void pkdKickTree(PKD pkd, double dTime, double dDelta, double dDeltaVPred, double dDeltaU, double dDeltaUPred, int iRoot) {
    /* Skip to local tree */
    auto c = pkd->tree[iRoot];
    while (c->is_remote()) c = pkd->tree[iRoot = c->lchild()];

    /* Now just kick all of the particles in the tree */
    for (auto& p : *c) {
        auto& v = p.velocity();
        auto& a = p.acceleration();
        v += a * dDelta;
        a = 0.0;
    }
}

void pkdInitCosmology(PKD pkd, const struct csmVariables* cosmo) {
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
void pkdInitLightcone(PKD pkd, int bBowtie, int bLightConeParticles, double dBoxSize, double dRedshiftLCP, double alphaLCP, blitz::TinyVector<double, 3> hLCP,
                      int mOctants) {
#ifdef __linux__
    uint64_t nPageSize = sysconf(_SC_PAGESIZE);
#else
    uint64_t nPageSize = 512;
#endif
    pkd->mLightConeOctants = mOctants;
    /*
    ** Initialize Lookup table for converting lightcone velocities to physical (sim units)
    */
    double dTimeLCP = csmExp2Time(pkd->csm, 1.0 / (1.0 + dRedshiftLCP));
    double mrLCP = dLightSpeedSim(dBoxSize) * csmComoveKickFac(pkd->csm, dTimeLCP, (csmExp2Time(pkd->csm, 1.0) - dTimeLCP));
    pkdSetupInterpScale(pkd, dBoxSize, mrLCP);
    /*
    ** Initialize light cone offsets.
    */
    if (pkd->Self() == 0) printf("alphaLCP = %f\n", alphaLCP);
    if (alphaLCP < 0) {
        initLightBallOffsets(pkd, mrLCP);
    }
    else {
        initLightConeOffsets(pkd, bBowtie, hLCP, alphaLCP, mrLCP);
    }
    /*
    ** allocate enough space for light cone particle output
    */
    uint64_t nLightConeBytes = (1024 * 1024 * 2);
    pkd->nLightConeMax = nLightConeBytes / sizeof(LIGHTCONEP);
    pkd->nLightCone = 0; // array of 8 octants
    auto mask = pkd->mLightConeOctants == 0 ? 1 : pkd->mLightConeOctants;
    for (auto octant = 0; octant < 8; ++octant) {
        if (bLightConeParticles && (mask & (1 << octant))) {
            void* v;
            if (posix_memalign(&v, nPageSize, nLightConeBytes)) pkd->pLightCone[octant] = NULL;
            else pkd->pLightCone[octant] = static_cast<LIGHTCONEP*>(v);
            mdlassert(pkd->mdl, pkd->pLightCone[octant] != NULL);
            io_init(&pkd->afiLightCone[octant], 8, 2 * 1024 * 1024, IO_AIO | IO_LIBAIO);
        }
        else {
            pkd->afiLightCone[octant].nBuffers = 0;
            pkd->pLightCone[octant] = NULL;
        }
        pkd->afiLightCone[octant].fd = -1;
    }
    pkd->pHealpixData = NULL;
}

void pkdZeroNewRung(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung) { /* JW: Ugly -- need to clean up */
    for (auto& p : pkd->particles) {
        if (p.is_active()) p.set_new_rung(0);
    }
}

void pkdChemCompInit(PKD pkd, struct inChemCompInit in) {

    for (auto& p : pkd->particles) {
        if (p.is_gas()) {
            auto& Sph = p.sph();
            float fMass = p.mass();

            if (Sph.ElemMass[ELEMENT_H] < 0.0f) {
                Sph.ElemMass[ELEMENT_H] = in.dInitialH * fMass;
                Sph.ElemMass[ELEMENT_He] = in.dInitialHe * fMass;
                Sph.ElemMass[ELEMENT_C] = in.dInitialC * fMass;
                Sph.ElemMass[ELEMENT_N] = in.dInitialN * fMass;
                Sph.ElemMass[ELEMENT_O] = in.dInitialO * fMass;
                Sph.ElemMass[ELEMENT_Ne] = in.dInitialNe * fMass;
                Sph.ElemMass[ELEMENT_Mg] = in.dInitialMg * fMass;
                Sph.ElemMass[ELEMENT_Si] = in.dInitialSi * fMass;
                Sph.ElemMass[ELEMENT_Fe] = in.dInitialFe * fMass;
            }
            if (Sph.fMetalMass < 0.0f) Sph.fMetalMass = in.dInitialMetallicity * fMass;
        }

#ifdef STELLAR_EVOLUTION
        else if (p.is_star()) {
            auto& Star = p.star();

            if (Star.ElemAbun[ELEMENT_H] < 0.0f) {
                Star.ElemAbun[ELEMENT_H] = in.dInitialH;
                Star.ElemAbun[ELEMENT_He] = in.dInitialHe;
                Star.ElemAbun[ELEMENT_C] = in.dInitialC;
                Star.ElemAbun[ELEMENT_N] = in.dInitialN;
                Star.ElemAbun[ELEMENT_O] = in.dInitialO;
                Star.ElemAbun[ELEMENT_Ne] = in.dInitialNe;
                Star.ElemAbun[ELEMENT_Mg] = in.dInitialMg;
                Star.ElemAbun[ELEMENT_Si] = in.dInitialSi;
                Star.ElemAbun[ELEMENT_Fe] = in.dInitialFe;
            }
            if (Star.fMetalAbun < 0.0f) Star.fMetalAbun = in.dInitialMetallicity;
        }
#endif // STELLAR_EVOLUTION
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
    T.d = fabs(dDelta) / dT;
    if (T.d <= 1.0) return 0;
    iRung = T.ieee.exponent - 1023; /* log2(d) */
    if (iRung > uMaxRung) return uMaxRung;
    else return iRung;
}

int pkdUpdateRung(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung, int iMaxRung, uint64_t* nRungCount) {
    int iTempRung;
    assert(!pkd->bNoParticleOrder);
    for (auto i = 0; i < iMaxRung; ++i) nRungCount[i] = 0;
    for (auto& p : pkd->particles) {
        if (p.is_active()) {
            if (p.new_rung() > iMaxRung) p.set_new_rung(iMaxRung);
            if (p.new_rung() >= uRung) p.set_rung(p.new_rung());
            else if (p.rung() > uRung) p.set_rung(uRung);
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

void pkdDeleteParticle(PKD pkd, particleStore::ParticleReference& p) {
    /* p.iOrder = -2 - p.iOrder; JW: Not needed -- just preserve iOrder */
    int pSpecies = p.species();
    p.set_class(0.0, 0.0, 0, 1.0, FIO_SPECIES_UNKNOWN); /* Special "DELETED" class == FIO_SPECIES_UNKNOWN */

    // IA: We copy the last particle into this position, the tree will no longer be valid!!!
    //
    // If we do something else with the tree before a reconstruction (such as finishing the particle loop), we need to be extra careful and:
    //   * check that the particle is marked
    //   * check that the particle type correspond to its placing inside the node
    //
    // Even with this, some weird bugs may appear, so seriously, be EXTRA careful!!
    // PARTICLE* lastp = pkd->Particle( pkd->Local() - 1);

    //
    // assert(!pkdIsDeleted(pkd,lastp));
    // IA: We can encounter a case where the last particle is the one being deleted; workaround:
    // while (pkdIsDeleted(pkd,lastp) || lastp==p){
    //   lastp--;
    //}

    // pkd->CopyParticle(p, lastp);
    // pkd->Local() -= 1;
    // pkd->TreeNode(ROOT)->pUpper -= 1;
    p.set_marked(false);
    p.set_rung(0);
    p.set_new_rung(0);
}

int pkdIsGas(PKD pkd, PARTICLE* p) {
    return pkdSpecies(pkd, p) == FIO_SPECIES_SPH;
}

int pkdIsDark(PKD pkd, PARTICLE* p) {
    return pkdSpecies(pkd, p) == FIO_SPECIES_DARK;
}

int pkdIsStar(PKD pkd, PARTICLE* p) {
    return pkdSpecies(pkd, p) == FIO_SPECIES_STAR;
}

int pkdIsBH(PKD pkd, PARTICLE* p) {
    return pkdSpecies(pkd, p) == FIO_SPECIES_BH;
}

void pkdOutPsGroup(PKD pkd, char* pszFileName, int iType) {
    FILE* fp;
    int i;

    if (iType == OUT_PSGROUP_STATS) {
        fp = fopen(pszFileName, "a+");
        assert(fp != NULL);
        struct psGroup* gd = pkd->psGroupTable.pGroup;

        for (i = 1; i < pkd->psGroupTable.nGroups; ++i) {
            if (gd[i].iPid != pkd->Self()) continue;
            fprintf(fp, "%d", gd[i].iGlobalId);
            fprintf(fp, " %10" PRIu64 "", gd[i].nTotal);
            fprintf(fp, " %12.8e", gd[i].fMass);
            fprintf(fp, " %12.8e", gd[i].fRMSRadius);
            fprintf(fp, " %12.8e", gd[i].r[0]);
            fprintf(fp, " %12.8e", gd[i].r[1]);
            fprintf(fp, " %12.8e", gd[i].r[2]);
            fprintf(fp, " %12.8e", gd[i].v[0]);
            fprintf(fp, " %12.8e", gd[i].v[1]);
            fprintf(fp, " %12.8e", gd[i].v[2]);
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
            fprintf(fp, "\n");
        }
        if (fclose(fp) == EOF) {
            perror("pkdOutPsGroup: could not close file");
            exit(1);
        }
    }
    else assert(0);
}

/*
** Initialize the EOS tables
*/
void pkdInitializeEOS(PKD pkd) {
    auto materials = pkd->particles.getMaterials(true);
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
                }
                else {
                    pkd->materials[iMat] = EOSinitMaterial(iMat, pkd->SPHoptions.dKpcUnit, pkd->SPHoptions.dMsolUnit, NULL);
                }
                if (pkd->SPHoptions.useIsentropic) {
                    EOSinitIsentropicLookup(pkd->materials[iMat], NULL);
                }
                if (pkd->SPHoptions.doShearStrengthModel && (EOSYieldModel(pkd->materials[iMat]) < 0)) {
                    if (EOSYieldModel(pkd->materials[iMat]) == -1) {
                        printf("Yield strength was requested but material %d has no or faulty parameters.", iMat);
                    }
                    else if (EOSYieldModel(pkd->materials[iMat]) == -2) {
                        printf("Yield strength was requested but for material %d this is not implemented yet.", iMat);
                    }
                    assert(0);
                }
            }
#else
            printf("Trying to initialize an EOSlib material, but EOSlib was not compiled in!\n");
            assert(0);
#endif
        }
    }
}

void pkdCopySPHOptionsToDevice(PKD pkd, SPHOptions* SPHoptions, int bGPU) {
    if (bGPU) {
#ifdef USE_CUDA
        // Only one thread needs to transfer the SPHoptions to the GPU
        if (pkd->mdl->Core() == 0) {
            SPHOptionsGPU SPHoptionsGPU;
            copySPHOptionsGPU(SPHoptions, &SPHoptionsGPU);
            pkd->cudaClient->setupSPHOptions(&SPHoptionsGPU);
        }
        pkd->mdl->ThreadBarrier();
#endif
    }
}
