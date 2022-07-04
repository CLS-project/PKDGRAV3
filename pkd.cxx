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
#ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
#else
    #define PRIu64 "llu"
#endif
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
extern "C" {
#include "core/healpix.h"
}
#ifdef COOLING
    #include "cooling/cooling.h"
#endif
#if ( defined(COOLING) || defined(GRACKLE) ) && defined(STAR_FORMATION)
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

//************************************************************************************************************************
//** The following routines are used to adjust the size of a PARTICLE and NODE based on the memory model
//************************************************************************************************************************

/* Add a NODE structure: assume double alignment */
int pkdContext::NodeAddStruct(int n) {
    int iOffset = iTreeNodeSize;
    assert(kdNodeListPRIVATE == nullptr);
    iTreeNodeSize += n;
    return iOffset;
}
/* Add n doubles to the node structure */
int pkdContext::NodeAddDouble(int n) {
    int iOffset = iTreeNodeSize;
    assert(kdNodeListPRIVATE == nullptr);
    assert((iOffset & (sizeof(double)-1)) == 0);
    iTreeNodeSize += sizeof(double) * n;
    return iOffset;
}
/* Add n floats to the node structure */
int pkdContext::NodeAddFloat(int n) {
    int iOffset = iTreeNodeSize;
    assert(kdNodeListPRIVATE == nullptr);
    assert((iOffset & (sizeof(float)-1)) == 0);
    iTreeNodeSize += sizeof(float) * n;
    return iOffset;
}
/* Add n 64-bit integers to the node structure */
int pkdContext::NodeAddInt64(int n) {
    int iOffset = iTreeNodeSize;
    assert(kdNodeListPRIVATE == nullptr);
    assert((iOffset & (sizeof(int64_t)-1)) == 0);
    iTreeNodeSize += sizeof(int64_t) * n;
    return iOffset;
}
/* Add n 32-bit integers to the node structure */
int pkdContext::NodeAddInt32(int n) {
    int iOffset = iTreeNodeSize;
    assert(kdNodeListPRIVATE == nullptr);
    assert((iOffset & (sizeof(int32_t)-1)) == 0);
    iTreeNodeSize += sizeof(int32_t) * n;
    return iOffset;
}

/* Add a structure: assume double alignment */
int pkdContext::ParticleAddStruct(int n) {
    int iOffset = iParticleSize;
    assert(ParticleBase() == nullptr);
    assert((iOffset & (sizeof(double)-1)) == 0);
    if (nParticleAlign < sizeof(double)) nParticleAlign = sizeof(double);
    iParticleSize += n;
    return iOffset;
}

/* Add n doubles to the particle structure */
int pkdContext::ParticleAddDouble(int n) {
    int iOffset = iParticleSize;
    assert(ParticleBase() == nullptr);
    if ( (iOffset & (sizeof(int64_t)-1)) != 0 ) {
        assert(iParticle32==0);
        iParticle32 = iParticleSize;
        iOffset = iParticleSize += sizeof(float);
    }
    assert((iOffset & (sizeof(double)-1)) == 0);
    if (nParticleAlign < sizeof(double)) nParticleAlign = sizeof(double);
    iParticleSize += sizeof(double) * n;
    return iOffset;
}

/* Add n floats to the particle structure */
int pkdContext::ParticleAddFloat(int n) {
    int iOffset = iParticleSize;
    assert(ParticleBase() == nullptr);
    if ( n==1 && iParticle32) {
        iOffset = iParticle32;
        iParticle32 = 0;
    }
    else {
        assert((iOffset & (sizeof(float)-1)) == 0);
        iParticleSize += sizeof(float) * n;
    }
    return iOffset;
}

/* Add n 64-bit integers to the particle structure */
int pkdContext::ParticleAddInt64(int n) {
    int iOffset = iParticleSize;
    assert(ParticleBase() == nullptr);
    if ( (iOffset & (sizeof(int64_t)-1)) != 0 ) {
        assert(iParticle32==0 );
        iParticle32 = iParticleSize;
        iOffset = iParticleSize += sizeof(float);
    }
    assert((iOffset & (sizeof(int64_t)-1)) == 0);
    if (nParticleAlign < sizeof(int64_t)) nParticleAlign = sizeof(int64_t);
    iParticleSize += sizeof(int64_t) * n;
    return iOffset;
}

/* Add n 32-bit integers to the particle structure */
int pkdContext::ParticleAddInt32(int n) {
    int iOffset = iParticleSize;
    assert(ParticleBase() == nullptr);
    if ( n==1 && iParticle32) {
        iOffset = iParticle32;
        iParticle32 = 0;
    }
    else {
        assert((iOffset & (sizeof(int32_t)-1)) == 0);
        iParticleSize += sizeof(int32_t) * n;
    }
    return iOffset;
}

//************************************************************************************************************************
//**
//************************************************************************************************************************

/* Extend the tree by adding more nodes */
void pkdContext::ExtendTree() {
    if ( nTreeTiles >= (1<<nTreeBitsHi) ) {
        fprintf(stderr, "ERROR: insufficent nodes available in tree build"
                "-- Increase nTreeBitsLo and/or nTreeBitsHi\n"
                "nTreeBitsLo=%d nTreeBitsHi=%d\n",
                nTreeBitsLo, nTreeBitsHi);
        assert( nTreeTiles < (1<<nTreeBitsHi) );
    }
    kdNodeListPRIVATE[nTreeTiles] = static_cast<KDN *>(mdlMalloc(mdl,(1<<nTreeBitsLo)*NodeSize()));
    mdlassert(mdl,kdNodeListPRIVATE[nTreeTiles] != NULL);
    ++nTreeTiles;
    nMaxNodes = (1<<nTreeBitsLo) * nTreeTiles;
}

static void firstTouch(uint64_t n,char *p) {
    if (n>4096) n-= 4096;
    while (n>=4096) {
        *p = 0;
        p += 4096;
        n -= 4096;
    }
}

/* Greatest common divisor */
static int gcd ( int a, int b ) {
    while ( a != 0 ) {
        int c = a;
        a = b%a;
        b = c;
    }
    return b;
}

static void initLightConeOffsets(PKD pkd) {
    BND bnd = {{0,0,0},{0.5,0.5,0.5}};
    double min2;
    int ix,iy,iz,nBox;

    /*
    ** Set up the light cone offsets such that they proceed from the inner 8
    ** unit boxes outward layer by layer so that we can skip checks of the
    ** outer layers if we want.
    */
    nBox = 0;
    for (ix=0; ix<=1; ++ix) {
        for (iy=0; iy<=1; ++iy) {
            for (iz=0; iz<=1; ++iz) {
                pkd->lcOffset0[nBox] = ix - 0.5;
                pkd->lcOffset1[nBox] = iy - 0.5;
                pkd->lcOffset2[nBox] = iz - 0.5;
                ++nBox;
            }
        }
    }
    assert(nBox == 8);
    for (ix=-1; ix<=2; ++ix) {
        for (iy=-1; iy<=2; ++iy) {
            for (iz=-1; iz<=2; ++iz) {
                if (ix>=0 && ix<=1 && iy>=0 && iy<=1 && iz>=0 && iz<=1)
                    continue;
                pkd->lcOffset0[nBox] = ix - 0.5;
                pkd->lcOffset1[nBox] = iy - 0.5;
                pkd->lcOffset2[nBox] = iz - 0.5;
                ++nBox;
            }
        }
    }
    assert(nBox == 64);
    for (ix=-2; ix<=3; ++ix) {
        for (iy=-2; iy<=3; ++iy) {
            for (iz=-2; iz<=3; ++iz) {
                if (ix>=-1 && ix<=2 && iy>=-1 && iy<=2 && iz>=-1 && iz<=2)
                    continue;
                double r[3] = {ix - 0.5, iy - 0.5, iz - 0.5};
                MINDIST(&bnd,r,min2);
                if (min2 < 9.0) {
                    pkd->lcOffset0[nBox] = r[0];
                    pkd->lcOffset1[nBox] = r[1];
                    pkd->lcOffset2[nBox] = r[2];
                    ++nBox;
                }
            }
        }
    }
    assert(nBox == 184);
}

pkdContext::pkdContext(mdl::mdlClass *mdl,
                       int nStore,uint64_t nMinTotalStore,uint64_t nMinEphemeral,uint32_t nEphemeralBytes,
                       int nTreeBitsLo, int nTreeBitsHi,
                       int iCacheSize,int iWorkQueueSize,int iCUDAQueueSize,double *fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,uint64_t nBH,
                       uint64_t mMemoryModel, int bLightCone, int bLightConeParticles) : mdl(mdl),
    csm(nullptr) {
    PARTICLE *p;
    uint32_t pi;
    int j,ism;

#ifdef __linux__
    uint64_t nPageSize = sysconf(_SC_PAGESIZE);
#else
    uint64_t nPageSize = 512;
#endif
    uint64_t nPageMask = nPageSize-1;

#define RANDOM_SEED 1
    srand(RANDOM_SEED);

    this->kdNodeListPRIVATE = NULL;
    this->nStore = 0; /* Set properly below */
    SetLocal(0);
    this->nDark = nDark;
    this->nGas = nGas;
    this->nBH = nBH;
    this->nStar = nStar;
    this->nRejects = 0;
    for (j=0; j<3; ++j) {
        this->fPeriod[j] = fPeriod[j];
    }

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

    if ( this->bNoParticleOrder )
        this->iParticleSize = this->nParticleAlign = sizeof(UPARTICLE);
    else
        this->iParticleSize = this->nParticleAlign = sizeof(PARTICLE);
    this->iParticle32 = 0;
    this->iTreeNodeSize = sizeof(KDN);

    if (!this->bIntegerPosition) this->oFieldOffset[oPosition] = ParticleAddDouble(3);
    if ( mMemoryModel & PKD_MODEL_PARTICLE_ID )
        this->oFieldOffset[oParticleID] = ParticleAddInt64();
    else
        this->oFieldOffset[oParticleID] = 0;
    /*
    ** Add a global group id. This is used when outputing an array of particles with
    ** one group assignment per particle.Usually only used in testing as it adds a
    ** 64 bit integer.
    */
    if ( mMemoryModel & PKD_MODEL_GLOBALGID )
        this->oFieldOffset[oGlobalGid] = ParticleAddInt64();
    else
        this->oFieldOffset[oGlobalGid] = 0;

    this->oFieldOffset[oVelocity] = 0;
    if ( mMemoryModel & PKD_MODEL_VELOCITY ) {
        if (sizeof(vel_t) == sizeof(double)) {
            this->oFieldOffset[oVelocity] = ParticleAddDouble(3);
        }
    }
    if (this->bIntegerPosition) this->oFieldOffset[oPosition] = ParticleAddInt32(3);
    if ( mMemoryModel & PKD_MODEL_RELAXATION )
        this->oFieldOffset[oRelaxation] = ParticleAddDouble();
    else
        this->oFieldOffset[oRelaxation] = 0;

    if ( mMemoryModel & PKD_MODEL_SPH )
#ifdef OPTIM_UNION_EXTRAFIELDS
        this->oFieldOffset[oSph] = ParticleAddStruct(sizeof(EXTRAFIELDS));
#else
        this->oFieldOffset[oSph] = ParticleAddStruct(sizeof(SPHFIELDS));
#endif
    else
        this->oFieldOffset[oSph] = 0;

    if ( mMemoryModel & PKD_MODEL_NEW_SPH )
        this->oFieldOffset[oNewSph] = ParticleAddStruct(sizeof(NEWSPHFIELDS));
    else
        this->oFieldOffset[oNewSph] = 0;

    if ( mMemoryModel & PKD_MODEL_STAR ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
        this->oFieldOffset[oStar] = 1;  // this value is of no relevance as long as it is >0
        if (!this->oFieldOffset[oSph])
            this->oFieldOffset[oSph] = ParticleAddStruct(sizeof(EXTRAFIELDS));
#else
        this->oFieldOffset[oStar] = ParticleAddStruct(sizeof(STARFIELDS));
#endif
    }
    else
        this->oFieldOffset[oStar] = 0;

#ifdef BLACKHOLES
    if ( mMemoryModel & PKD_MODEL_BH ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
        this->oFieldOffset[oBH] = 1;    // this value is of no relevance as long as it is >0
        if (!this->oFieldOffset[oSph])
            this->oFieldOffset[oSph] = ParticleAddStruct(sizeof(EXTRAFIELDS));
#else
        this->oFieldOffset[oBH] = ParticleAddStruct(sizeof(BHFIELDS));
#endif
    }
    else
        this->oFieldOffset[oBH] = 0;
#endif // BLACKHOLES

    if ( mMemoryModel & PKD_MODEL_VELSMOOTH )
        this->oFieldOffset[oVelSmooth] = ParticleAddStruct(sizeof(VELSMOOTH));
    else
        this->oFieldOffset[oVelSmooth] = 0;
    if ( mMemoryModel & PKD_MODEL_VELOCITY ) {
        if (sizeof(vel_t) == sizeof(float)) {
            this->oFieldOffset[oVelocity] = ParticleAddFloat(3);
        }
    }
    if ( mMemoryModel & PKD_MODEL_ACCELERATION )
        this->oFieldOffset[oAcceleration] = ParticleAddFloat(3);
    else
        this->oFieldOffset[oAcceleration] = 0;

    if ( mMemoryModel & PKD_MODEL_MASS )
        this->oFieldOffset[oMass] = ParticleAddFloat();
    else
        this->oFieldOffset[oMass] = 0;

    if ( mMemoryModel & PKD_MODEL_SOFTENING )
        this->oFieldOffset[oSoft] = ParticleAddFloat();
    else
        this->oFieldOffset[oSoft] = 0;

    if ( mMemoryModel & (PKD_MODEL_SPH|PKD_MODEL_NEW_SPH|PKD_MODEL_BALL) )
        this->oFieldOffset[oBall] = ParticleAddFloat();
    else this->oFieldOffset[oBall] = 0;
    if ( mMemoryModel & (PKD_MODEL_SPH|PKD_MODEL_NEW_SPH|PKD_MODEL_DENSITY) )
        this->oFieldOffset[oDensity] = ParticleAddFloat();
    else this->oFieldOffset[oDensity] = 0;

    this->oFieldOffset[oGroup] = 0;
    if ( (mMemoryModel & PKD_MODEL_GROUPS) && !this->bNoParticleOrder) {
        this->oFieldOffset[oGroup] = ParticleAddInt32();
    }
    else this->oFieldOffset[oGroup] = 0;

    if ( mMemoryModel & PKD_MODEL_POTENTIAL ) {
        this->oFieldOffset[oPotential] = ParticleAddFloat();
    }
    else this->oFieldOffset[oPotential] = 0;

    /*
    ** Tree node memory models
    */
    if (this->bIntegerPosition) this->oNodePosition = NodeAddInt32(3);
    else this->oNodePosition = NodeAddDouble(3);
    if ( mMemoryModel & PKD_MODEL_NODE_BND ) {
        if (this->bIntegerPosition) this->oNodeBnd  = NodeAddStruct(sizeof(IBND));
        else this->oNodeBnd  = NodeAddStruct(sizeof(BND));
    }
    else {
        this->oNodeBnd  = 0;
    }

    if ( mMemoryModel & PKD_MODEL_NODE_VBND ) {
        this->oNodeVBnd  = NodeAddStruct(sizeof(BND));
    }
    else {
        this->oNodeVBnd = 0;
    }

    this->oNodeVelocity = 0;
    if ( (mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(double))
        this->oNodeVelocity = NodeAddDouble(3);
    if ( mMemoryModel & (PKD_MODEL_SPH|PKD_MODEL_BH) ) {
#ifdef OPTIM_REORDER_IN_NODES
        this->oNodeNgas = NodeAddInt32();
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
        this->oNodeNstar = NodeAddInt32();
#endif
        this->oNodeNbh = NodeAddInt32();
#endif
    }
    /*
    ** Three extra bounds are required by the fast gas SPH code.
    */
    if ( mMemoryModel & PKD_MODEL_NODE_SPHBNDS ) {
        this->oNodeSphBounds = NodeAddStruct(sizeof(SPHBNDS));
    }
    else
        this->oNodeSphBounds = 0;

    if ( mMemoryModel & PKD_MODEL_NODE_MOMENT )
        this->oNodeMom = NodeAddStruct(sizeof(FMOMR));
    else
        this->oNodeMom = 0;

    /* The acceleration is required for the new time step criteria */
    if ( mMemoryModel & PKD_MODEL_NODE_ACCEL )
        this->oNodeAcceleration = NodeAddFloat(3);
    else
        this->oNodeAcceleration = 0;

    if ( (mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(float))
        this->oNodeVelocity = NodeAddFloat(3);
    /*
    ** N.B.: Update pkd->MaxNodeSize in pkd.h if you add fields.  We need to
    **       know the size of a node when setting up the pst.
    */
    assert(NodeSize() > 0);
    if (NodeSize() > pkdContext::MaxNodeSize()) {
        fprintf(stderr, "Node size is too large. Node size=%" PRIu64 ", max node size=%" PRIu64 "\n",
                (uint64_t)NodeSize(), (uint64_t)MaxNodeSize());
    }
    assert(NodeSize()<=MaxNodeSize());

    /* Align the particle size and the tree node, and store the tree node parameters */
    this->iParticleSize = (this->iParticleSize + this->nParticleAlign - 1 ) & ~(this->nParticleAlign-1);
    this->iTreeNodeSize = (this->iTreeNodeSize + sizeof(double) - 1 ) & ~(sizeof(double)-1);
    this->nTreeBitsLo = nTreeBitsLo;
    this->nTreeBitsHi = nTreeBitsHi;
    this->iTreeMask = (1<<this->nTreeBitsLo) - 1;

    /* Adjust nStore so that we use an integer number of pages */
    int n = nPageSize / gcd(nPageSize,this->iParticleSize);
    nStore += n; /* Not "n-1" here because we reserve one at the end */
    nStore -= nStore % n;
    assert( (uint64_t)nStore*this->iParticleSize % nPageSize == 0);
    --nStore;
    this->nStore = nStore;

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
    ** IMPORTANT: There is a whole lot of pointer math here. If you mess with this
    **            you better be sure you get it right or really bad things will happen.
    */
    this->nEphemeralBytes = nEphemeralBytes;
    uint64_t nBytesPerThread = ((nStore+1)*ParticleSize()+nPageMask) & ~nPageMask; // Constraint (d)
    uint64_t nBytesParticles = (uint64_t)mdlCores(this->mdl) * nBytesPerThread; // Constraint (a)
    uint64_t nBytesEphemeral = (uint64_t)mdlCores(this->mdl) * (nStore+1)*1ul*this->nEphemeralBytes; // Constraint (a)
    uint64_t nBytesTreeNodes = 0;
    if (nBytesEphemeral < nMinEphemeral) nBytesEphemeral = nMinEphemeral; // Constraint (b)
    if (nBytesParticles + nBytesEphemeral < nMinTotalStore) // Constraint (c)
        nBytesTreeNodes = nMinTotalStore - nBytesParticles - nBytesEphemeral;
    // Align to a even number of "tree tiles"
    uint64_t nTreeTileBytesPerNode = (1<<this->nTreeBitsLo)*this->iTreeNodeSize*mdlCores(this->mdl);
    uint64_t nTreeTiles = (uint64_t)ceil(1.0 * nBytesTreeNodes / nTreeTileBytesPerNode);
    nBytesTreeNodes = nTreeTiles * nTreeTileBytesPerNode;
    char *pParticles, *pEphemeral, *pTreeNodes;
    if (mdlCore(this->mdl)==0) {
        uint64_t nBytesTotal = nBytesParticles + nBytesEphemeral + nBytesTreeNodes;
        void *vParticles;
#ifdef _MSC_VER
        pParticles = _aligned_malloc(nBytesTotal,nPageSize);
#else
        if (posix_memalign(&vParticles,nPageSize,nBytesTotal)) pParticles = NULL;
        else pParticles = static_cast<char *>(vParticles);
#endif
        mdlassert(mdl,pParticles != NULL);
        pEphemeral = pParticles + nBytesParticles;
        pTreeNodes = pEphemeral + nBytesEphemeral;
    }
    else pParticles = pEphemeral = pTreeNodes = 0; // Ignore anyway in mdlSetArray() below
    pParticles = static_cast<char *>(mdlSetArray(this->mdl,1,nBytesPerThread,pParticles));
    pEphemeral = static_cast<char *>(mdlSetArray(this->mdl,1,nBytesEphemeral/mdlCores(this->mdl),pEphemeral));
    pTreeNodes = static_cast<char *>(mdlSetArray(this->mdl,1,nBytesTreeNodes/mdlCores(this->mdl),pTreeNodes));
    firstTouch(nBytesParticles/mdlCores(this->mdl),pParticles);
    firstTouch(nBytesEphemeral/mdlCores(this->mdl),pEphemeral);
    firstTouch(nBytesTreeNodes/mdlCores(this->mdl),pTreeNodes);
    this->pStorePRIVATE = (PARTICLE *)pParticles;
    this->pLite = pEphemeral;
    /*
    ** Now we setup the node storage for the tree.  This storage is no longer
    ** continguous as the MDL now supports non-contiguous arrays.  We allocate
    ** a single "tile" for the tree.  If this is not sufficient, then additional
    ** tiles are allocated dynamically.  The default parameters allow for 2^32
    ** nodes total which is the integer limit anyway. We may use the extra storage
    ** from above if constraint (c) could not otherwise be met.
    */
    this->kdNodeListPRIVATE = static_cast<KDN **>(mdlMalloc(this->mdl,(1<<this->nTreeBitsHi)*sizeof(KDN *)));
    mdlassert(mdl,this->kdNodeListPRIVATE != NULL);
    if (nTreeTiles) {
        this->nTreeTilesReserved = nTreeTiles;
        this->nTreeTiles = nTreeTiles;
        for (j=0; j<nTreeTiles; ++j) {
            this->kdNodeListPRIVATE[j] = reinterpret_cast<KDN *>(pTreeNodes);
            pTreeNodes += (1<<this->nTreeBitsLo)*this->iTreeNodeSize;
        }
    }
    else {
        this->nTreeTilesReserved = 0;
        this->kdNodeListPRIVATE[0] = static_cast<KDN *>(mdlMalloc(this->mdl,(1<<this->nTreeBitsLo)*this->iTreeNodeSize));
        mdlassert(mdl,this->kdNodeListPRIVATE[0] != NULL);
        this->nTreeTiles = 1;
    }
    this->nMaxNodes = (1<<this->nTreeBitsLo) * this->nTreeTiles;
    this->nNodes = 0;

    /*
    ** We also allocate a temporary particle used for swapping.  We need to do
    ** this now because the outside world can no longer know the size of a
    ** particle a priori.
    */
    this->pTempPRIVATE = static_cast<PARTICLE *>(malloc(ParticleSize()));
    mdlassert(mdl,this->pTempPRIVATE != NULL);
    /*
    ** Initialize light cone offsets.
    */
    initLightConeOffsets(this);
    /*
    ** allocate enough space for light cone particle output
    */
    uint64_t nLightConeBytes = (1024*1024*16);
    this->nLightConeMax = nLightConeBytes / sizeof(LIGHTCONEP);
    this->nLightCone = 0;
    if (bLightCone && bLightConeParticles) {
        void *v;
#ifdef _MSC_VER
        this->pLightCone = _aligned_malloc(nLightConeBytes, nPageSize);
#else
        if (posix_memalign(&v, nPageSize, nLightConeBytes)) this->pLightCone = NULL;
        else this->pLightCone = static_cast<LIGHTCONEP *>(v);
#endif
        mdlassert(mdl,this->pLightCone != NULL);
        io_init(&this->afiLightCone,8,2*1024*1024,IO_AIO|IO_LIBAIO);
    }
    else {
        this->afiLightCone.nBuffers = 0;
        this->pLightCone = NULL;
    }
    this->afiLightCone.fd = -1;
    this->pHealpixData = NULL;

#ifdef MDL_CACHE_SIZE
    if ( iCacheSize > 0 ) mdlSetCacheSize(this->mdl,iCacheSize);
#endif
    // This is cheeserific - chooses the largest specified

#if defined(USE_CUDA) || defined(USE_CL)
    this->cudaClient = new CudaClient(*this->mdl);
    mdlSetCudaBufferSize(this->mdl,PP_CUDA_MEMORY_LIMIT,PP_CUDA_MEMORY_LIMIT);
#endif
#if defined(USE_METAL)
    this->metalClient = new MetalClient(*this->mdl);
#endif
    mdlSetWorkQueueSize(this->mdl,iWorkQueueSize,iCUDAQueueSize);
    /*
    ** Initialize neighbor list pointer to NULL if present.
    */
    if (this->oFieldOffset[oSph]) {
        for (pi=0; pi<(this->nStore+1); ++pi) {
            p = Particle(pi);
            *pkd_pNeighborList(this,p) = NULL;
        }
    }

    /*
    ** Initialize global group id.
    */
    if (this->oFieldOffset[oGlobalGid]) {
        for (pi=0; pi<(this->nStore+1); ++pi) {
            p = Particle(pi);
            *pkdGlobalGid(this,p) = 0;
        }
    }

    /*
    ** We support up to 256 classes
    */
    ParticleClasses.reserve(PKD_MAX_CLASSES);

    this->fSoftFix = -1.0;
    this->fSoftFac = 1.0;
    this->fSoftMax = HUGE_VALF;
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

pkdContext::~pkdContext() {
    PARTICLE *p;
    char **ppCList;
    uint32_t pi;
    int ism;
    int i;

    if (kdNodeListPRIVATE) {
        /*
        ** Close caching space and free up nodes.
        */
        if (mdlCacheStatus(mdl,CID_CELL))      mdlFinishCache(mdl,CID_CELL);
        if (mdlCacheStatus(mdl,CID_CELL2))     mdlFinishCache(mdl,CID_CELL2);
        if (mdlCacheStatus(mdl,CID_PARTICLE2)) mdlFinishCache(mdl,CID_PARTICLE2);
        for ( i=nTreeTilesReserved; i<nTreeTiles; i++)
            mdlFree(mdl,kdNodeListPRIVATE[i]);
        mdlFree(mdl,kdNodeListPRIVATE);
    }
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
        SIMD_free(ewt.hx.f);
        SIMD_free(ewt.hy.f);
        SIMD_free(ewt.hz.f);
        SIMD_free(ewt.hCfac.f);
        SIMD_free(ewt.hSfac.f);
    }

    /*
    ** Free any neighbor lists that were left hanging around.
    */
    if (oFieldOffset[oSph]) {
        for (pi=0; pi<(FreeStore()+1); ++pi) {
            p = Particle(pi);
            if (pkdIsGas(this,p)) {
                ppCList = pkd_pNeighborList(this,p);
                if (*ppCList) {
                    free(*ppCList);
                    *ppCList = NULL;
                }
            }
        }
    }
    /* Only thread zero allocated this memory block  */
    mdlThreadBarrier(mdl);
    if (mdlCore(mdl) == 0) {
#ifdef _MSC_VER
        _aligned_free(ParticleBase());
#else
        mdlFree(mdl, ParticleBase());
#endif
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

void pkdSetClass( PKD pkd, float fMass, float fSoft, FIO_SPECIES eSpecies, PARTICLE *p ) {
    if ( pkd->oFieldOffset[oMass] ) {
        auto pMass = static_cast<float *>(pkdField(p,pkd->oFieldOffset[oMass]));
        *pMass = fMass;
        fMass = 0.0;
    }
    if ( pkd->oFieldOffset[oSoft] ) {
        auto pSoft = static_cast<float *>(pkdField(p,pkd->oFieldOffset[oSoft]));
        *pSoft = fSoft;
        fSoft = 0.0;
    }
    /* NOTE: The above can both be true, in which case a "zero" class is recorded */
    /* NOTE: Species is always part of the class table, so there will be at least one class per species */
    PARTCLASS newClass(fMass,fSoft,eSpecies);

    /* TODO: This is a linear search which is fine for a small number of classes */
    auto iclass = std::find(pkd->ParticleClasses.begin(),pkd->ParticleClasses.end(),newClass);
    if (iclass==pkd->ParticleClasses.end()) {
        assert( pkd->ParticleClasses.size() < PKD_MAX_CLASSES );
        p->iClass = pkd->ParticleClasses.size();
        pkd->ParticleClasses.emplace_back(newClass);
    }
    else p->iClass = std::distance(pkd->ParticleClasses.begin(),iclass);
    if (pkd->bNoParticleOrder) {
        assert(p->iClass==0);
    }
}

int pkdGetClasses( PKD pkd, int nMax, PARTCLASS *pClass ) {
    std::copy(pkd->ParticleClasses.begin(),pkd->ParticleClasses.end(),pClass);
    return pkd->ParticleClasses.size();
}

void pkdSetClasses( PKD pkd, int n, PARTCLASS *pClass, int bUpdate ) {
    uint8_t map[PKD_MAX_CLASSES];
    PARTICLE *p;

    if ( bUpdate && pkd->ParticleClasses.size() && !pkd->bNoParticleOrder) {
        /* Build a map from the old class to the new class */
        assert( n >= pkd->ParticleClasses.size() );
        for (auto i=0; i<pkd->ParticleClasses.size(); ++i) {
            auto jj = std::find(pClass,pClass+n,pkd->ParticleClasses[i]);
            map[i] = std::distance(pClass,jj);
        }

        /* Now update the class with the new value */
        for (auto i=0; i<pkd->Local(); ++i) {
            p = pkd->Particle(i);
            assert( p->iClass < pkd->ParticleClasses.size() );
            p->iClass = map[p->iClass];
        }
    }

    /* Finally, set the new class table */
    pkd->ParticleClasses.clear();
    pkd->ParticleClasses.insert(pkd->ParticleClasses.end(),pClass,pClass+n);
}

void pkdReadFIO(PKD pkd,FIO fio,uint64_t iFirst,int nLocal,double dvFac, double dTuFac) {
    int i,j;
    PARTICLE *p;
    BHFIELDS *pBH;
    float *pPot, dummypot;
    double r[3];
    double vel[3];
    float fMass, fSoft,fDensity,u,fMetals[ELEMENT_COUNT],fTimer;
    FIO_SPECIES eSpecies;
    uint64_t iParticleID;

    mdlassert(pkd->mdl,fio != NULL);

#ifdef USE_ITT
    __itt_domain *domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle *shMyTask = __itt_string_handle_create("Read");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
    pkd->ParticleClasses.clear();
    if (pkd->oFieldOffset[oStar]) {
        /* Make sure star class established -- how do all procs know of these classes? How do we ensure they agree on the class identifiers? */
        p = pkd->Particle(pkd->Local());
        pkdSetClass(pkd,0,0,FIO_SPECIES_STAR,p);
    }
#ifdef BLACKHOLES
    assert(pkd->oFieldOffset[oMass]);
    p = pkd->Particle( pkd->Local());
    pkdSetClass(pkd,0,0,FIO_SPECIES_BH, p);
#endif

    // Protect against uninitialized values
    fMass = 0.0f;
    fSoft = 0.0f;

    fioSeek(fio,iFirst,FIO_SPECIES_ALL);
    for (i=0; i<nLocal; ++i) {
        p = pkd->Particle(pkd->Local()+i);
        /*
        ** General initialization.
        */
        p->uRung =  0;
        if (!pkd->bNoParticleOrder) p->uNewRung = 0;
        p->bMarked = 1;
        pkdSetDensity(pkd,p,0.0);
        if (pkd->oFieldOffset[oBall]) pkdSetBall(pkd,p,0.0);
        /*
        ** Clear the accelerations so that the timestepping calculations do not
        ** get funny uninitialized values!
        */
        if ( pkd->oFieldOffset[oAcceleration] ) {
            float *a = pkdAccel(pkd,p);
            for (j=0; j<3; ++j) a[j] = 0.0;
        }
        if ( pkd->oFieldOffset[oPotential]) pPot = pkdPot(pkd,p);
        else pPot = &dummypot;
        pkdSetGroup(pkd,p,0);

        /* Initialize SPH fields if present */
        SPHFIELDS *pSph;
        if (pkd->oFieldOffset[oSph]) {
            pSph = static_cast<SPHFIELDS *>(pkdField(p,pkd->oFieldOffset[oSph]));
#ifndef OPTIM_REMOVE_UNUSED
            pSph->u = pSph->uPred = pSph->uDot = pSph->c = pSph->divv = pSph->BalsaraSwitch
                                                 = pSph->fMetals = pSph->diff = pSph->fMetalsPred = pSph->fMetalsDot = 0.0;
#endif
        }
        else pSph = NULL;

        /* Initialize New SPH fields if present */
        NEWSPHFIELDS *pNewSph;
        if (pkd->oFieldOffset[oNewSph]) {
            pNewSph = static_cast<NEWSPHFIELDS *>(pkdField(p,pkd->oFieldOffset[oNewSph]));
            pNewSph->u = pNewSph->uDot = pNewSph->divv = pNewSph->Omega = 0.0;
        }
        else pNewSph = NULL;

        /* Initialize Star fields if present */
        STARFIELDS *pStar;
        if (pkd->oFieldOffset[oStar]) {
            pStar = static_cast<STARFIELDS *>(pkdField(p,pkd->oFieldOffset[oStar]));
            pStar->fTimer = 0;
            /*      pStar->iGasOrder = IORDERMAX;*/
        }
        else pStar = NULL;

        eSpecies = fioSpecies(fio);
        switch (eSpecies) {
        case FIO_SPECIES_SPH:
            ;
            float afSphOtherData[2];
            fioReadSph(fio,&iParticleID,r,vel,&fMass,&fSoft,pPot,
                       &fDensity,&u,&fMetals[0],afSphOtherData);
            pkdSetClass(pkd,fMass,fSoft,eSpecies,p);
            pkdSetDensity(pkd,p,fDensity);
            if (pNewSph) {
                pNewSph->u = -u; /* Can't do conversion until density known */
            }
            else {
                assert(dTuFac>0.0);
                pkdSetBall(pkd,p,afSphOtherData[0]);
                if (pkd->oFieldOffset[oSph]) {
                    pSph = pkdSph(pkd, p);
#ifndef OPTIM_REMOVE_UNUSED
                    pSph->u = pSph->uPred = pSph->uDot = pSph->c = pSph->divv =
                            pSph->BalsaraSwitch = pSph->diff =
                                                      pSph->fMetals = pSph->fMetalsPred = pSph->fMetalsDot = 0.0;
#endif
                    for (j = 0; j < ELEMENT_COUNT; j++) pSph->afElemMass[j] = fMetals[j] * fMass;
#ifdef HAVE_METALLICITY
                    pSph->fMetalMass = afSphOtherData[1] * fMass;
#endif
                    // If the value is negative, means that it is a temperature
                    u = (u<0.0) ? -u*dTuFac : u;
#ifndef OPTIM_REMOVE_UNUSED
                    pSph->u = u * dTuFac;
#endif
                    /* IA: -unused- variables
                    pSph->fMetals = fMetals;
                      pSph->uPred = pSph->u;
                      pSph->fMetalsPred = pSph->fMetals;
                              */
                    pSph->vPred[0] = vel[0]*sqrt(dvFac);
                    pSph->vPred[1] = vel[1]*sqrt(dvFac);
                    pSph->vPred[2] = vel[2]*sqrt(dvFac);
                    pSph->Frho = 0.0;
                    pSph->Fmom[0] = 0.0;
                    pSph->Fmom[1] = 0.0;
                    pSph->Fmom[2] = 0.0;
                    pSph->Fene = 0.0;
                    pSph->E = u + 0.5*(pSph->vPred[0]*pSph->vPred[0] +
                                       pSph->vPred[1]*pSph->vPred[1] +
                                       pSph->vPred[2]*pSph->vPred[2]);
                    pSph->E *= fMass;
                    pSph->Uint = u*fMass;
                    assert(pSph->E>0);
                    pSph->mom[0] = fMass*vel[0]*sqrt(dvFac);
                    pSph->mom[1] = fMass*vel[1]*sqrt(dvFac);
                    pSph->mom[2] = fMass*vel[2]*sqrt(dvFac);
                    pSph->lastMom[0] = 0.; // vel[0];
                    pSph->lastMom[1] = 0.; //vel[1];
                    pSph->lastMom[2] = 0.; //vel[2];
                    pSph->lastE = pSph->E;
#ifdef ENTROPY_SWITCH
                    pSph->S = 0.0;
                    pSph->lastS = 0.0;
                    pSph->maxEkin = 0.0;
#endif
                    pSph->lastUint = pSph->Uint;
                    pSph->lastHubble = 0.0;
                    pSph->lastMass = fMass;
                    pSph->lastAcc[0] = 0.;
                    pSph->lastAcc[1] = 0.;
                    pSph->lastAcc[2] = 0.;
#ifndef USE_MFM
                    pSph->lastDrDotFrho[0] = 0.;
                    pSph->lastDrDotFrho[1] = 0.;
                    pSph->lastDrDotFrho[2] = 0.;
                    pSph->drDotFrho[0] = 0.;
                    pSph->drDotFrho[1] = 0.;
                    pSph->drDotFrho[2] = 0.;
#endif
                    //pSph->fLastBall = 0.0;
                    pSph->lastUpdateTime = -1.;
                    // pSph->nLastNeighs = 100;
#ifdef COOLING
                    pSph->lastCooling = 0.;
                    pSph->cooling_dudt = 0.;
#endif
#ifdef FEEDBACK
                    pSph->fAccFBEnergy = 0.;
#endif
                    pSph->uWake = 0;
                }
            }
            break;
        case FIO_SPECIES_DARK:
            fioReadDark(fio,&iParticleID,r,vel,&fMass,&fSoft,pPot,&fDensity);
            pkdSetClass(pkd,fMass,fSoft,eSpecies,p);
            pkdSetDensity(pkd,p,fDensity);
            break;
        case FIO_SPECIES_STAR:
            ;
            float afStarOtherData[4];
            fioReadStar(fio,&iParticleID,r,vel,&fMass,&fSoft,pPot,&fDensity,
                        fMetals,&fTimer,afStarOtherData);
            pkdSetClass(pkd,fMass,fSoft,eSpecies,p);
            pkdSetDensity(pkd,p,fDensity);
            if (pkd->oFieldOffset[oStar]) {
                pStar = pkdStar(pkd,p);
                pStar->fTimer = fTimer;
                // We avoid that star in the IC could explode
                pStar->hasExploded = 1;
#ifdef FEEDBACK
                pStar->fSNEfficiency = afStarOtherData[3];
#endif
#ifdef STELLAR_EVOLUTION
                for (j = 0; j < ELEMENT_COUNT; j++)
                    pStar->afElemAbun[j] = fMetals[j];
                pStar->fMetalAbun = afStarOtherData[0];
                pStar->fInitialMass = afStarOtherData[1];
                pStar->fLastEnrichTime = afStarOtherData[2];
#endif
            }
            break;
        case FIO_SPECIES_BH:
            pkdSetBall(pkd,p,pkdSoft(pkd,p));
            float otherData[3];
            fioReadBH(fio,&iParticleID,r,vel,&fMass,&fSoft,pPot,
                      &fDensity,otherData,&fTimer);
            pkdSetClass(pkd,fMass,fSoft,eSpecies,p);
            if (pkd->oFieldOffset[oBH]) {
                pBH = pkdBH(pkd,p);
                pBH->fTimer = fTimer;
                pBH->pLowPot = NULL;
                pBH->newPos[0] = -1;
                pBH->lastUpdateTime = -1.;
                pBH->dInternalMass = otherData[0];
                pBH->dAccretionRate = otherData[1];
                pBH->dAccEnergy = otherData[2];
                pBH->dFeedbackRate = 0.0;
            }
            break;
        default:
            fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
            assert(0);
        }

        for (j=0; j<3; ++j) {
            pkdSetPos(pkd,p,j,r[j]);
        }
        if (!pkd->bNoParticleOrder) p->iOrder = iFirst++;
        if (pkd->oFieldOffset[oParticleID]) *pkdParticleID(pkd,p) = iParticleID;

        if (pkd->oFieldOffset[oVelocity]) {
            if (!pkdIsGas(pkd,p)) {
                // IA: dvFac = a*a, and for the gas we already provide
                // the peculiar velocity in the IC
                for (j=0; j<3; ++j) pkdVel(pkd,p)[j] = vel[j]*dvFac;
            }
            else {
                for (j=0; j<3; ++j) pkdVel(pkd,p)[j] = vel[j]*sqrt(dvFac);
            }
        }

    }

    pkd->AddLocal(nLocal);
    pkd->nActive += nLocal;

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
}

void pkdCalcBound(PKD pkd,BND *pbnd) {
    double r[3],dMin[3],dMax[3];
    PARTICLE *p;
    int i = 0;
    int j;

    mdlassert(pkd->mdl,pkd->Local() > 0);
    p = pkd->Particle(i);
    for (j=0; j<3; ++j) {
        dMin[j] = dMax[j] = pkdPos(pkd,p,j);
    }
    for (++i; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        for (j=0; j<3; ++j) r[j] = pkdPos(pkd,p,j);
        pkdMinMax(r,dMin,dMax);
    }
    for (j=0; j<3; ++j) {
        pbnd->fCenter[j] = pkd->bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
        pbnd->fMax[j] = pkd->bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
    }
}

void pkdEnforcePeriodic(PKD pkd,BND *pbnd) {
    PARTICLE *p;
    double r;
    int i,j;
#if defined(USE_SIMD) && defined(__SSE2__)
    if (pkd->bIntegerPosition) {
        __m128i period = _mm_set1_epi32 (INTEGER_FACTOR);
        __m128i top = _mm_setr_epi32 ( (INTEGER_FACTOR/2)-1,(INTEGER_FACTOR/2)-1,(INTEGER_FACTOR/2)-1,0x7fffffff );
        __m128i bot = _mm_setr_epi32 (-(INTEGER_FACTOR/2),-(INTEGER_FACTOR/2),-(INTEGER_FACTOR/2),-0x80000000 );

        auto pPos = static_cast<char *>(pkdField(pkd->Particle(0),pkd->oFieldOffset[oPosition]));
        const int iSize = pkd->ParticleSize();
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
        for (i=0; i<pkd->Local(); ++i) {
            p = pkd->Particle(i);
            for (j=0; j<3; ++j) {
                r = pkdPos(pkd,p,j);
                if (r < pbnd->fCenter[j] - pbnd->fMax[j]) r += 2*pbnd->fMax[j];
                else if (r >= pbnd->fCenter[j] + pbnd->fMax[j]) r -= 2*pbnd->fMax[j];
                pkdSetPos(pkd,p,j,r);
                /*
                ** If it still doesn't lie in the "unit" cell then something has gone quite wrong with the
                ** simulation. Either we have a super fast particle or the initial condition is somehow not conforming
                ** to the specified periodic box in a gross way.
                */
                //      mdlassert(pkd->mdl,((r >= pbnd->fCenter[j] - pbnd->fMax[j])&&
                //      (r < pbnd->fCenter[j] + pbnd->fMax[j])));
            }
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
    PARTICLE *pi, *pj;
    pi = pkd->Particle(i);
    pj = pkd->Particle(j);
    PARTITION(pi<pj,pi<=pj,
              pi=pkd->Particle(++i),pj=pkd->Particle(--j),
              pkdSwapParticle(pkd,pi,pj),
              pkdPos(pkd,pi,d) >= fSplit,pkdPos(pkd,pj,d) < fSplit);
    return (i);
}


int pkdUpperPart(PKD pkd,int d,double fSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkd->Particle(i);
    pj = pkd->Particle(j);
    PARTITION(pi<pj,pi<=pj,
              pi=pkd->Particle(++i),pj=pkd->Particle(--j),
              pkdSwapParticle(pkd,pi,pj),
              pkdPos(pkd,pi,d) < fSplit,pkdPos(pkd,pj,d) >= fSplit);
    return (i);
}


int pkdLowerPartWrap(PKD pkd,int d,double fSplit1,double fSplit2,int i,int j) {
    PARTICLE *pi = pkd->Particle(i);
    PARTICLE *pj = pkd->Particle(j);

    if (fSplit1 > fSplit2) {
        PARTITION(pi<pj,pi<=pj,
                  pi=pkd->Particle(++i),pj=pkd->Particle(--j),
                  pkdSwapParticle(pkd,pi,pj),
                  (pkdPos(pkd,pi,d) < fSplit2 || pkdPos(pkd,pi,d) >= fSplit1),
                  (pkdPos(pkd,pj,d) >= fSplit2 && pkdPos(pkd,pj,d) < fSplit1));
    }
    else {
        PARTITION(pi<pj,pi<=pj,
                  pi=pkd->Particle(++i),pj=pkd->Particle(--j),
                  pkdSwapParticle(pkd,pi,pj),
                  (pkdPos(pkd,pi,d) < fSplit2 && pkdPos(pkd,pi,d) >= fSplit1),
                  (pkdPos(pkd,pj,d) >= fSplit2 || pkdPos(pkd,pj,d) < fSplit1));
    }
    return (i);
}


int pkdUpperPartWrap(PKD pkd,int d,double fSplit1,double fSplit2,int i,int j) {
    PARTICLE *pi = pkd->Particle(i);
    PARTICLE *pj = pkd->Particle(j);

    if (fSplit1 > fSplit2) {
        PARTITION(pi<pj,pi<=pj,
                  pi=pkd->Particle(++i),pj=pkd->Particle(--j),
                  pkdSwapParticle(pkd,pi,pj),
                  (pkdPos(pkd,pi,d) >= fSplit2 && pkdPos(pkd,pi,d) < fSplit1),
                  (pkdPos(pkd,pj,d) < fSplit2 || pkdPos(pkd,pj,d) >= fSplit1));
    }
    else {
        PARTITION(pi<pj,pi<=pj,
                  pi=pkd->Particle(++i),pj=pkd->Particle(--j),
                  pkdSwapParticle(pkd,pi,pj),
                  (pkdPos(pkd,pi,d) >= fSplit2 || pkdPos(pkd,pi,d) < fSplit1),
                  (pkdPos(pkd,pj,d) < fSplit2 && pkdPos(pkd,pj,d) >= fSplit1));
    }
    return (i);
}


int pkdLowerOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkd->Particle(i);
    pj = pkd->Particle(j);
    PARTITION(pi<pj,pi<=pj,
              pi=pkd->Particle(++i),pj=pkd->Particle(--j),
              pkdSwapParticle(pkd,pi,pj),
              pi->iOrder >= nOrdSplit,pj->iOrder < nOrdSplit);
    return (i);
}


int pkdUpperOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkd->Particle(i);
    pj = pkd->Particle(j);
    PARTITION(pi<pj,pi<=pj,
              pi=pkd->Particle(++i),pj=pkd->Particle(--j),
              pkdSwapParticle(pkd,pi,pj),
              pi->iOrder < nOrdSplit,pj->iOrder >= nOrdSplit);
    return (i);
}


int pkdActiveOrder(PKD pkd) {
    int i=0;
    int j=pkd->Local()-1;
    PARTICLE *pi, *pj;
    pi = pkd->Particle(i);
    pj = pkd->Particle(j);
    PARTITION(pi<pj,pi<=pj,
              pi=pkd->Particle(++i),pj=pkd->Particle(--j),
              pkdSwapParticle(pkd,pi,pj),
              pkdIsActive(pkd,pi),!pkdIsActive(pkd,pj));
    return (pkd->nActive = i);
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
            pkd->CopyParticle(pkd->Particle(iRejects+i),pkd->Particle(nSplit+i));
    }
    pkd->SetLocal(nSplit);
    return (pkd->nRejects);
}


int pkdSwapRejects(PKD pkd,int idSwap) {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;

    if (idSwap != -1) {
        nBuf = (pkdSwapSpace(pkd))*pkd->ParticleSize();
        nOutBytes = pkd->nRejects*pkd->ParticleSize();
        mdlassert(pkd->mdl,pkd->Local() + pkd->nRejects <= pkd->FreeStore());
        mdlSwap(pkd->mdl,idSwap,nBuf,pkd->Particle(pkd->Local()),
                nOutBytes,&nSndBytes,&nRcvBytes);
        pkd->AddLocal(nRcvBytes/pkd->ParticleSize());
        pkd->nRejects -= nSndBytes/pkd->ParticleSize();
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
    nBuf = pkd->FreeStore()*pkd->ParticleSize();
    nOutBytes = pkd->Local()*pkd->ParticleSize();
    mdlSwap(pkd->mdl,idSwap,nBuf,pkd->ParticleBase(), nOutBytes,
            &nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl,nSndBytes/pkd->ParticleSize() == pkd->Local());
    pkd->SetLocal(nRcvBytes/pkd->ParticleSize());
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
** cache element.  Normal code should call pkdTreeNode().
*/
void *pkdTreeNodeGetElement(void *vData,int i,int iDataSize) {
    auto pkd = static_cast<PKD>(vData);
    return pkd->TreeNode(i);
}

int pkdColOrdRejects(PKD pkd,uint64_t nOrdSplit,int iSplitSide) {
    int nSplit;
    if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkd->Local()-1);
    else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkd->Local()-1);
    return pkdColRejects(pkd,nSplit);
}

//static int cmpParticles(const void *pva,const void *pvb) {
//    PARTICLE *pa = (PARTICLE *)pva;
//    PARTICLE *pb = (PARTICLE *)pvb;
//    return(pa->iOrder - pb->iOrder);
//    }

void pkdLocalOrder(PKD pkd,uint64_t iMinOrder, uint64_t iMaxOrder) {
    int i;
    assert(pkd->Local() == iMaxOrder - iMinOrder + 1);
    for (i=0; i<pkd->Local(); ++i) {
        PARTICLE *p1 = pkd->Particle(i);
        assert(p1->iOrder >= iMinOrder && p1->iOrder <= iMaxOrder);
        while (p1->iOrder - iMinOrder !=  i) {
            PARTICLE *p2 = pkd->Particle(p1->iOrder-iMinOrder);
            pkdSwapParticle(pkd,p1,p2);
        }
    }
    /* Above replaces: qsort(pkd->ParticleBase(),pkd->Local(),pkd->ParticleSize(),cmpParticles); */
}

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
    nFileSize = pkd->ParticleSize() * pkd->Local();
    char *pBuffer = (char *)pkd->ParticleBase();
    while (nFileSize) {
        size_t count = nFileSize > MAX_IO_BUFFER_SIZE ? MAX_IO_BUFFER_SIZE : nFileSize;
        io_write(&info, pBuffer, count);
        pBuffer += count;
        nFileSize -= count;
    }
    io_close(&info);
}

void pkdRestore(PKD pkd,const char *fname) {
    asyncFileInfo info;
    size_t nFileSize;
    int fd;
    io_init(&info, IO_MAX_ASYNC_COUNT, 0, IO_AIO|IO_LIBAIO);
    fd = io_open(&info, fname);
    if (fd<0) {
        perror(fname);
        abort();
    }
    struct stat s;
    if ( fstat(fd,&s) != 0 ) {
        perror(fname);
        abort();
    }
    nFileSize = s.st_size;
    pkd->SetLocal(nFileSize / pkd->ParticleSize());
    char *pBuffer = (char *)pkd->ParticleBase();
    while (nFileSize) {
        size_t count = nFileSize > MAX_IO_BUFFER_SIZE ? MAX_IO_BUFFER_SIZE : nFileSize;
        io_read(&info, pBuffer, count);
        pBuffer += count;
        nFileSize -= count;
    }
    io_close(&info);
}

/*****************************************************************************\
* Write particles received from another node
\*****************************************************************************/

static void writeParticle(PKD pkd,FIO fio,double dvFac,double dvFacGas,BND *bnd,PARTICLE *p) {
    const BHFIELDS *pBH;
    float *pPot, dummypot;
    double v[3],r[3];
    float fMass, fSoft, fDensity,fMetals[ELEMENT_COUNT], fTimer, fBall;
    uint64_t iParticleID;
    int j;

    dummypot = 0.0;

    if ( pkd->oFieldOffset[oPotential]) pPot = pkdPot(pkd,p);
    else pPot = &dummypot;
    if (pkd->oFieldOffset[oVelocity]) {
        /* IA: the gas velocity in the code is v = a \dot x
         *  and the dm/star velocity v = a^2 \dot x
         *
         *  However, we save both as \dot x such that
         *  they can be be directly compared at the output
         */
        if (!pkdIsGas(pkd,p)) {
            vel_t *pV = pkdVel(pkd,p);
            v[0] = pV[0] * dvFac;
            v[1] = pV[1] * dvFac;
            v[2] = pV[2] * dvFac;
        }
        else {
            vel_t *pV = pkdVel(pkd,p);
            v[0] = pV[0] * dvFacGas;
            v[1] = pV[1] * dvFacGas;
            v[2] = pV[2] * dvFacGas;
        }
    }
    else v[0] = v[1] = v[2] = 0.0;

    /* Initialize SPH fields if present */
    const auto *pSph = static_cast<SPHFIELDS *>(pkd->oFieldOffset[oSph] ? pkdField(p,pkd->oFieldOffset[oSph]) : nullptr);
    const auto *pNewSph = static_cast<NEWSPHFIELDS *>(pkd->oFieldOffset[oNewSph] ? pkdField(p,pkd->oFieldOffset[oNewSph]) : nullptr);
    const auto *pStar = static_cast<STARFIELDS *>(pkd->oFieldOffset[oStar] ? pkdField(p,pkd->oFieldOffset[oStar]) : nullptr);
    fMass = pkdMass(pkd,p);
    fSoft = pkdSoft0(pkd,p);
    if (pkd->fSoftFix >= 0.0) fSoft = 0.0;
    if (pkd->oFieldOffset[oParticleID]) iParticleID = *pkdParticleID(pkd,p);
    else if (!pkd->bNoParticleOrder) iParticleID = p->iOrder;
    else iParticleID = 0;
    if (pkd->oFieldOffset[oDensity]) fDensity = pkdDensity(pkd,p);
    else fDensity = 0.0;

    r[0] = pkdPos(pkd,p,0);
    r[1] = pkdPos(pkd,p,1);
    r[2] = pkdPos(pkd,p,2);
    /* Enforce periodic boundaries */
    for (j=0; j<3; ++j) {
        if (r[j] < bnd->fCenter[j] - bnd->fMax[j]) r[j] += 2*bnd->fMax[j];
        else if (r[j] >= bnd->fCenter[j] + bnd->fMax[j]) r[j] -= 2*bnd->fMax[j];
        /*
        ** If it still doesn't lie in the "unit" cell then something has gone quite wrong with the
        ** simulation. Either we have a super fast particle or the initial condition is somehow not conforming
        ** to the specified periodic box in a gross way.
        */
        mdlassert(pkd->mdl,((r[j] >= bnd->fCenter[j] - bnd->fMax[j])&&
                            (r[j] < bnd->fCenter[j] + bnd->fMax[j])));

    }

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
    switch (pkdSpecies(pkd,p)) {
    case FIO_SPECIES_SPH:
        if (pkd->oFieldOffset[oNewSph]) {
            assert(pNewSph);
            assert(pkd->SPHoptions.TuFac > 0.0f);
            double T;
            float otherData[3];
            otherData[0] = otherData[1] = otherData[2] = 0.0f;
            T = EOSTofRhoU(fDensity, pNewSph->u, &pkd->SPHoptions);
            for (int k = 0; k < ELEMENT_COUNT; k++) fMetals[k] = 0.0f;
            fioWriteSph(fio,iParticleID,r,v,fMass,fSoft,*pPot,
                        fDensity,T,&fMetals[0],0.0f,T,otherData);
        }
        else {
            pSph = pkdSph(pkd,p);
            assert(pSph);
            {
#if defined(COOLING)
                const double dRedshift = dvFacGas - 1.;
                float temperature =  cooling_get_temperature(pkd, dRedshift, pkd->cooling, p, pSph);
#elif defined(GRACKLE)
                gr_float fDensity = pkdDensity(pkd,p);
                gr_float fMetalDensity = pSph->fMetalMass*pSph->omega;
                gr_float fSpecificUint = pSph->Uint/pkdMass(pkd,p);

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

                for (int k = 0; k < ELEMENT_COUNT; k++) fMetals[k] = pSph->afElemMass[k] / fMass;

#ifdef STAR_FORMATION
                float SFR = pSph->SFR;
#else
                float SFR=0.;
#endif

                fBall = pkdBall(pkd,p);
                float otherData[3];
                otherData[0] = SFR;
                // We may have problem if the number of groups increses more than 2^24, but should be enough
                otherData[1] = pkd->oFieldOffset[oGroup] ? (float)pkdGetGroup(pkd,p) : 0 ;
#ifdef HAVE_METALLICITY
                otherData[2] = pSph->fMetalMass / fMass;
#endif

                fioWriteSph(fio,iParticleID,r,v,fMass,fSoft,*pPot,
                            fDensity,pSph->Uint/fMass, &fMetals[0], fBall, temperature, &otherData[0]);
            }
        }
        break;
    case FIO_SPECIES_DARK: {
        float otherData[2];
        otherData[0] = pkd->oFieldOffset[oGroup] ? (float)pkdGetGroup(pkd,p) : 0 ;
        fioWriteDark(fio,iParticleID,r,v,fMass,fSoft,*pPot,fDensity, &otherData[0]);
    }
    break;
    case FIO_SPECIES_STAR:
        pStar = pkdStar(pkd,p);
#ifdef STELLAR_EVOLUTION
        for (int k = 0; k < ELEMENT_COUNT; k++) fMetals[k] = pStar->afElemAbun[k];
#endif
        float otherData[6];
        otherData[0] = pStar->fTimer;
        otherData[1] = pkd->oFieldOffset[oGroup] ? (float)pkdGetGroup(pkd,p) : 0 ;
#ifdef STELLAR_EVOLUTION
        otherData[2] = pStar->fMetalAbun;
        otherData[3] = pStar->fInitialMass;
        otherData[4] = pStar->fLastEnrichTime;
#endif
#ifdef FEEDBACK
        otherData[5] = pStar->fSNEfficiency;
#endif
        fioWriteStar(fio,iParticleID,r,v,fMass,fSoft,*pPot,fDensity,
                     &fMetals[0],&otherData[0]);
        break;
    case FIO_SPECIES_BH: {
        pBH = pkdBH(pkd,p);
        fTimer = pBH->fTimer;
        float otherData[6];
        otherData[0] = pBH->dInternalMass;
        otherData[1] = pBH->dAccretionRate;
        otherData[2] = pBH->dEddingtonRatio;
        otherData[3] = pBH->dFeedbackRate;
        otherData[4] = pBH->dAccEnergy;
        otherData[5] = pkd->oFieldOffset[oGroup] ? (float)pkdGetGroup(pkd,p) : 0 ;
        fioWriteBH(fio,iParticleID,r,v,fMass,fSoft,*pPot,fDensity,
                   otherData,fTimer);
    }
    break;
    case FIO_SPECIES_LAST:
        break;
    default:
        fprintf(stderr,"Unsupported particle type: %d\n",pkdSpecies(pkd,p));
        assert(0);
    }
}

struct packWriteCtx {
    PKD pkd;
    FIO fio;
    BND *bnd;
    double dvFac;
    double dTuFac;
    int iIndex;
};

static int unpackWrite(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packWriteCtx *ctx = (struct packWriteCtx *)vctx;
    PKD pkd = ctx->pkd;
    PARTICLE *p = (PARTICLE *)vBuff;
    int n = nSize / pkd->ParticleSize();
    int i;
    double dvFacGas = sqrt(ctx->dvFac);
    assert( n*pkd->ParticleSize() == nSize);
    for (i=0; i<n; ++i) {
        writeParticle(pkd,ctx->fio,ctx->dvFac,dvFacGas,ctx->bnd,pkd->ParticleGet(p,i));
    }
    return 1;
}

void pkdWriteFromNode(PKD pkd,int iNode, FIO fio,double dvFac,double dTuFac,BND *bnd) {
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
    int n = nSize / pkd->ParticleSize();
    if ( n > nLeft ) n = nLeft;
    nSize = n*pkd->ParticleSize();
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
    numPart_file[3] = nBH;
    numPart_file[4] = nStar;
    numPart_file[5] = 0;

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
    int iIndex;
    int field;
    int iUnitSize;
    int bMarked;
};

char *pkdPackArray(PKD pkd,int iSize,void *vBuff,int *piIndex,int n,int field,int iUnitSize,double dvFac,int bMarked) {
    auto pBuff = static_cast<char *>(vBuff);
    int oOffset = pkd->oFieldOffset[field];
    int iIndex = *piIndex;

    while (iIndex<n && iSize>=iUnitSize) {
        PARTICLE *p = pkd->Particle(iIndex++);
        if (bMarked && !p->bMarked) continue;
        if (field==oPosition) {
            double *d = (double *)pBuff;
            pkdGetPos1(pkd,p,d);
        }
        else if (field==oVelocity) {
            vel_t *v = pkdVel(pkd,p);
            float *V = (float *)pBuff;
            V[0] = v[0] * dvFac;
            V[1] = v[1] * dvFac;
            V[2] = v[2] * dvFac;
        }
        else {
            const char *src = (const char *)pkdField(p,oOffset);
            memcpy(pBuff,src,iUnitSize);
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
void pkdSendArray(PKD pkd, int iNode, int field, int iUnitSize,double dvFac,int bMarked) {
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

uint32_t pkdWriteFIO(PKD pkd,FIO fio,double dvFac,double dTuFac,BND *bnd) {
    PARTICLE *p;
    int i;
    uint32_t nCount;
    double dvFacGas = sqrt(dvFac);
    nCount = 0;
    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        writeParticle(pkd,fio,dvFac,dvFacGas,bnd,p);
        nCount++;
    }
    return nCount;
}

void pkdSetSoft(PKD pkd,double dSoft) {
    pkd->fSoftFix = dSoft;
}

void pkdSetSmooth(PKD pkd,double dSmooth) {
    PARTICLE *p;
    int i;
    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        if ((pkdIsGas(pkd,p)||pkdIsStar(pkd,p)||pkdIsBH(pkd,p)) && pkdBall(pkd,p)==0.0)
            pkdSetBall(pkd,p,dSmooth);
    }
}

void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul) {
    pkd->fSoftFac = dFac;
    pkd->fSoftMax = bSoftMaxMul ? HUGE_VALF : dSoftMax;
}

static void initSetMarked(void *vpkd, void *v) {}
static void combSetMarked(void *vpkd, void *v1, const void *v2) {
    PARTICLE *p1 = (PARTICLE *)v1;
    const PARTICLE *p2 = (const PARTICLE *)v2;
    if (p2->bMarked) p1->bMarked = 1;
}

void pkdGravAll(PKD pkd,
                struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
                double dTime,int nReps,int bPeriodic,
                int bEwald,int nGroup,int iRoot1, int iRoot2,
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
        pkdEwaldInit(pkd,nReps,fEwCut,fEwhCut); /* ignored in Flop count! */
    }
    /*
    ** Start particle caching space (cell cache already active).
    */
    if (SPHoptions->doSetDensityFlags) {
        mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,pkd->ParticleBase(),pkd->ParticleSize(),
                   pkd->Local(),NULL,initSetMarked,combSetMarked);
    }
    else {
        mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkd->ParticleBase(),pkd->ParticleSize(),
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

    *pnActive = pkdGravWalk(pkd,kick,lc,ts,
                            dTime,nReps,bPeriodic && bEwald,nGroup,
                            iRoot1,iRoot2,0,dThetaMin,pdFlop,&dPartSum,&dCellSum,SPHoptions);

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
void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double *L,double *F,double *W) {
    int i;

    *T = pkd->dEnergyT;
    *U = pkd->dEnergyU;
    *W = pkd->dEnergyW;
    for (i=0; i<3; ++i) {
        L[i] = pkd->dEnergyL[i];
        F[i] = pkd->dEnergyF[i];
    }
    *Eth = 0.0;
    if (pkd->oFieldOffset[oSph]) {
        int n = pkd->Local();
        for (i=0; i<n; ++i) {
            PARTICLE *p = pkd->Particle(i);
            if (pkdIsGas(pkd,p)) *Eth += pkdSph(pkd,p)->E;
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

void addToLightCone(PKD pkd,double *r,float fPot,PARTICLE *p,int bParticleOutput) {
    vel_t *v = pkdVel(pkd,p);
    if (pkd->afiLightCone.fd>0 && bParticleOutput) {
        LIGHTCONEP *pLC = pkd->pLightCone;
        pLC[pkd->nLightCone].pos[0] = r[0];
        pLC[pkd->nLightCone].pos[1] = r[1];
        pLC[pkd->nLightCone].pos[2] = r[2];
        pLC[pkd->nLightCone].vel[0] = v[0];
        pLC[pkd->nLightCone].vel[1] = v[1];
        pLC[pkd->nLightCone].vel[2] = v[2];
#ifdef POTENTIAL_IN_LIGHTCONE
        pLC[pkd->nLightCone].pot    = fPot;
#endif
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

#ifndef USE_SIMD_LC
#define NBOX 184
void pkdProcessLightCone(PKD pkd,PARTICLE *p,float fPot,double dLookbackFac,double dLookbackFacLCP,double dDriftDelta,double dKickDelta,double dBoxSize,int bLightConeParticles) {
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    const double mrLCP = dLightSpeed*dLookbackFacLCP;
    double vrx0[NBOX],vry0[NBOX],vrz0[NBOX];
    double vrx1[NBOX],vry1[NBOX],vrz1[NBOX];
    double mr0[NBOX],mr1[NBOX],x[NBOX];
    double r0[3],r1[3];
    double dlbt, dt, xStart, mr;
    vel_t *v;
    int j,k,iOct,nBox,bParticleOutput;
    struct {
        double dt;
        double fOffset;
        int jPlane;
    } isect[4], temp;

    /*
    ** Check all 184 by default.
    */
    nBox = NBOX;
    xStart = (dLookbackFac*dLightSpeed - 3.0)/(dKickDelta*dLightSpeed);
    if (xStart > 1) return;
    else if (xStart < 0) {
        xStart = (dLookbackFac*dLightSpeed - 2.0)/(dKickDelta*dLightSpeed);
        if (xStart >= 0) xStart = 0;
        else {
            /*
            ** Check only 64!
            */
            nBox = 64;
            xStart = (dLookbackFac*dLightSpeed - 1.0)/(dKickDelta*dLightSpeed);
            if (xStart >= 0) xStart = 0;
            else {
                /*
                ** Check only 8!
                */
                nBox = 8;
                xStart = 0;
            }
        }
    }

    v = pkdVel(pkd,p);
    pkdGetPos1(pkd,p,r0);
    for (j=0; j<3; ++j) {
        if (r0[j] < -0.5) r0[j] += 1.0;
        else if (r0[j] >= 0.5) r0[j] -= 1.0;
    }
    for (j=0; j<3; ++j) {
        isect[j].dt = (0.5 - r0[j])/v[j];
        if (isect[j].dt > 0.0) {
            /*
            ** Particle crosses the upper j-coordinate boundary of the unit cell at isect[j].dt.
            */
            isect[j].fOffset = -1.0;
        }
        else {
            /*
            ** Particle crosses the lower j-coordinate boundary of the unit cell at isect[j].dt.
            */
            isect[j].dt = (-0.5 - r0[j])/v[j];
            isect[j].fOffset = 1.0;
        }
        isect[j].jPlane = j;
    }
    isect[3].dt = dDriftDelta;
    isect[3].fOffset = 0.0;
    isect[3].jPlane = 3;
    /*
    ** Sort them!
    */
    if (isect[0].dt>isect[1].dt) {
        temp = isect[0];
        isect[0] = isect[1];
        isect[1] = temp;
    }
    if (isect[2].dt>isect[3].dt) {
        temp = isect[2];
        isect[2] = isect[3];
        isect[3] = temp;
    }
    temp = isect[1];
    isect[1] = isect[2];
    isect[2] = temp;
    if (isect[0].dt>isect[1].dt) {
        temp = isect[0];
        isect[0] = isect[1];
        isect[1] = temp;
    }
    if (isect[2].dt>isect[3].dt) {
        temp = isect[2];
        isect[2] = isect[3];
        isect[3] = temp;
    }
    if (isect[1].dt>isect[2].dt) {
        temp = isect[1];
        isect[1] = isect[2];
        isect[2] = temp;
    }

    for (k=0; k<4; ++k) {
        double dtApprox;
        if (k==0) {
            /*
            ** Check lightcone from 0 <= dt < isect[k].dt
            */
            dt = isect[k].dt;
            dtApprox = dt/dDriftDelta*dKickDelta;
            dlbt = dLookbackFac;
        }
        else {
            /*
            ** Check lightcone from isect[k-1].dt <= dt < isect[k].dt
            */
            dt = isect[k].dt - isect[k-1].dt;
            dtApprox = dt/dDriftDelta*dKickDelta;
            dlbt = dLookbackFac - dtApprox;
        }
        for (j=0; j<3; ++j) r1[j] = r0[j] + dt*v[j];
        for (iOct=0; iOct<nBox; ++iOct) {
            vrx0[iOct] = pkd->lcOffset0[iOct] + r0[0];
            vry0[iOct] = pkd->lcOffset1[iOct] + r0[1];
            vrz0[iOct] = pkd->lcOffset2[iOct] + r0[2];
            vrx1[iOct] = pkd->lcOffset0[iOct] + r1[0];
            vry1[iOct] = pkd->lcOffset1[iOct] + r1[1];
            vrz1[iOct] = pkd->lcOffset2[iOct] + r1[2];
            mr0[iOct] = sqrt(vrx0[iOct]*vrx0[iOct] + vry0[iOct]*vry0[iOct] + vrz0[iOct]*vrz0[iOct]);
            mr1[iOct] = sqrt(vrx1[iOct]*vrx1[iOct] + vry1[iOct]*vry1[iOct] + vrz1[iOct]*vrz1[iOct]);
            x[iOct] = (dLightSpeed*dlbt - mr0[iOct])/(dLightSpeed*dtApprox - mr0[iOct] + mr1[iOct]);
        }
        for (iOct=0; iOct<nBox; ++iOct) {
            if (x[iOct] >= xStart && x[iOct] < 1.0) {
                double r[3];
                /*
                ** Create a new light cone particle.
                */
                r[0] = (1-x[iOct])*vrx0[iOct] + x[iOct]*vrx1[iOct];
                r[1] = (1-x[iOct])*vry0[iOct] + x[iOct]*vry1[iOct];
                r[2] = (1-x[iOct])*vrz0[iOct] + x[iOct]*vrz1[iOct];
                mr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
                addToLightCone(pkd,r,fPot,p,bLightConeParticles && (mr <= mrLCP));
            }
        }
        if (isect[k].jPlane == 3) break;
        /*
        ** Now we need to reposition r0 to the new segment.
        */
        for (j=0; j<3; ++j) r0[j] = r1[j];
        r0[isect[k].jPlane] += isect[k].fOffset;
    }
}
#undef NBOX
#endif

/*
** Drift particles whose Rung falls between uRungLo (large step) and uRungHi (small step) inclusive,
** and those whose destination activity flag is set.
**
** Note that the drift funtion no longer wraps the particles around the periodic "unit" cell. This is
** now done by Domain Decomposition only.
*/
void pkdDrift(PKD pkd,int iRoot,double dTime,double dDelta,double dDeltaVPred,double dDeltaTime,int bDoGas) {
    PARTICLE *p;
    vel_t *v;
    int i,j;
    double rfinal[3],r0[3],dMin[3],dMax[3],dr[3];
    int pLower, pUpper;

    if (iRoot>=0) {
        KDN *pRoot = pkd->TreeNode(iRoot);
        pLower = pRoot->pLower;
        pUpper = pRoot->pUpper;
    }
    else {
        pLower = 0;
        pUpper = pkd->Local();
    }

    mdlDiag(pkd->mdl, "Into pkdDrift\n");
    assert(pkd->oFieldOffset[oVelocity]);

    for (j=0; j<3; ++j) {
        dMin[j] = pkd->bnd.fCenter[j] - pkd->bnd.fMax[j];
        dMax[j] = pkd->bnd.fCenter[j] + pkd->bnd.fMax[j];
    }
    /*
    ** Update particle positions
    */
    if (bDoGas) {
        if (pkd->oFieldOffset[oSph] /*pkd->param.bMeshlessHydro*/) {
            assert(pkd->oFieldOffset[oSph]);
            assert(pkd->oFieldOffset[oAcceleration]);
            for (i=pLower; i<=pUpper; ++i) {
                p = pkd->Particle(i);
                v = pkdVel(pkd,p);

                if (pkdIsGas(pkd,p)) {
                    for (j=0; j<3; ++j) {
                        // As for gas particles we use dx/dt = v/a, we must use the
                        // Kick factor provided by pkdgrav.
                        // See Stadel 2001, Appendix 3.8
                        dr[j] = v[j]*dDeltaVPred;

                    }
                }
                else {
                    for (j=0; j<3; ++j) {
                        dr[j] = v[j]*dDelta;
                    }
                }

#ifdef FORCE_1D
                dr[2] = 0.0;
                dr[1] = 0.0;
#endif
#ifdef FORCE_2D
                dr[2] = 0.0;
#endif



                pkdGetPos1(pkd,p,r0);
                for (j=0; j<3; ++j) {
                    pkdSetPos(pkd,p,j,rfinal[j] = r0[j] + dr[j]);
                }


                pkdMinMax(rfinal,dMin,dMax);
            }
        }
        else {
            assert(pkd->oFieldOffset[oNewSph]);;
            for (i=pLower; i<=pUpper; ++i) {
                p = pkd->Particle(i);
                v = pkdVel(pkd,p);
                // if (pkdIsGas(pkd,p)) {
                // NEWSPHFIELDS *NewSph = pkdNewSph(pkd,p);
                // float dfBalldt = 1.0f / 3.0f * pkdBall(pkd,p) * pkdDensity(pkd,p) * NewSph->divv;
                // pkdSetBall(pkd,p,pkdBall(pkd,p) + dDelta * dfBalldt);
                // }
                for (j=0; j<3; ++j) {
                    pkdSetPos(pkd,p,j,rfinal[j] = pkdPos(pkd,p,j) + dDelta*v[j]);
                }
                pkdMinMax(rfinal,dMin,dMax);
            }
        }
    }
    else {
        for (i=pLower; i<=pUpper; ++i) {
            p = pkd->Particle(i);
            v = pkdVel(pkd,p);
            pkdGetPos1(pkd,p,r0);
            for (j=0; j<3; ++j) {
                pkdSetPos(pkd,p,j,rfinal[j] = r0[j] + dDelta*v[j]);
                assert(isfinite(rfinal[j]));
            }
            pkdMinMax(rfinal,dMin,dMax);
        }
    }
    for (j=0; j<3; ++j) {
        pkd->bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
        pkd->bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
    }
    mdlDiag(pkd->mdl, "Out of pkdDrift\n");
}




#ifdef OPTIM_REORDER_IN_NODES
void pkdReorderWithinNodes(PKD pkd) {
    KDN *node;
    for (int i=NRESERVED_NODES; i<pkd->Nodes(); i++) {
        int nGas = 0;
        node = pkd->TreeNode(i);
        if (!node->iLower) { // We are in a bucket
            int start = node->pLower;
            for (int pj=node->pLower; pj<=node->pUpper; ++pj) {
                if (pkdIsGas(pkd, pkd->Particle(pj))) {
                    pkdSwapParticle(pkd, pkd->Particle(pj), pkd->Particle(start+nGas) );
                    nGas++;
                }

            }
            pkdNodeSetNgas(pkd, node, nGas);
            start += nGas;
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
            // We perform another swap, just to have nGas->nStar->DM
            int nStar = 0;
            for (int pj=start; pj<=node->pUpper; ++pj) {
                if (pkdIsStar(pkd, pkd->Particle(pj))) {
                    pkdSwapParticle(pkd, pkd->Particle(pj), pkd->Particle(start+nStar) );
                    nStar++;
                }

            }
            pkdNodeSetNstar(pkd, node, nStar);
            start += nStar;
#endif

            // We perform another swap, just to have nGas->nStar->BH->DM
            int nBH = 0;
            for (int pj=start; pj<=node->pUpper; ++pj) {
                if (pkdIsBH(pkd, pkd->Particle(pj))) {
                    pkdSwapParticle(pkd, pkd->Particle(pj), pkd->Particle(start+nBH) );
                    nBH++;
                }

            }
            pkdNodeSetNBH(pkd, node, nBH);

        }
    }
    /* // Check that this works
    for (int i=NRESERVED_NODES; i<pkd->nNodes; i++){
       node = pkd->TreeNode(i);
       if (!node->iLower){ // We are in a bucket
          if (pkdNodeNstar(pkd,node)>0){
          printf("Start node %d (%d, %d)\n",i, node->pLower, node->pUpper);
          //abort();
          for (int pj=node->pLower;pj<=node->pUpper;++pj) {
             if (pkdIsGas(pkd, pkd->Particle(pj)) ) printf("%d is Gas \n", pj);
             if (pkdIsStar(pkd, pkd->Particle(pj)) ) printf("%d is Star \n", pj);
             if (pkdIsDark(pkd, pkd->Particle(pj)) ) printf("%d is DM \n", pj);
          }
          }
       }
    }
    */
}
#endif


void pkdEndTimestepIntegration(PKD pkd, struct inEndTimestep in) {
    PARTICLE *p;
    SPHFIELDS *psph;
    int i;
    double pDelta, dScaleFactor, dHubble, pa[3];

    int bComove = pkd->csm->val.bComove;

    mdlDiag(pkd->mdl, "Into pkdComputePrimiteVars\n");
    assert(pkd->oFieldOffset[oVelocity]);
#ifndef USE_MFM
    assert(pkd->oFieldOffset[oMass]);
#endif

    if (bComove) {
        dScaleFactor = csmTime2Exp(pkd->csm,in.dTime);
        dHubble = csmTime2Hub(pkd->csm,in.dTime);
    }
    else {
        dScaleFactor = 1.0;
        dHubble = 0.0;
    }
#if defined(COOLING) || defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
    const double a_inv = 1./dScaleFactor;
    const double a_inv3 = a_inv*a_inv*a_inv;
#endif
#ifdef GRACKLE
    pkdGrackleUpdate(pkd, dScaleFactor, in.achCoolingTable, in.units);
#endif
    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        if (pkdIsGas(pkd,p) && pkdIsActive(pkd, p)  ) {
            psph = pkdSph(pkd, p);

            // ##### Add ejecta from stellar evolution
#ifdef STELLAR_EVOLUTION
            if (in.bChemEnrich && psph->fReceivedMass > 0.0f) {
                pkdAddStellarEjecta(pkd, p, psph, in.dConstGamma);
            }
#endif

#if defined(COOLING) || defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)
            float fMass = pkdMass(pkd, p);
            float fDens = pkdDensity(pkd, p);
            const float fDensPhys = fDens*a_inv3;
#endif
            if (in.dDelta > 0) {
                pDelta = in.dTime - psph->lastUpdateTime;
            }
            else {
                pDelta = 0.0;
            }

            // ##### Gravity
            hydroSourceGravity(pkd, p, psph,
                               pDelta, &pa[0], dScaleFactor, bComove);


            // ##### Expansion effects
            hydroSourceExpansion(pkd, p, psph,
                                 pDelta, dScaleFactor, dHubble, bComove, in.dConstGamma);



            // ##### Synchronize Uint, Etot (and possibly S)
            hydroSyncEnergies(pkd, p, psph, pa, in.dConstGamma);


            // ##### Cooling
#ifdef COOLING
            double dRedshift = 1./dScaleFactor - 1.;
            const float delta_redshift = -pDelta * dHubble * (dRedshift + 1.);
            cooling_cool_part(pkd, pkd->cooling, p, psph, pDelta, in.dTime, delta_redshift, dRedshift);
#endif
#ifdef GRACKLE
            pkdGrackleCooling(pkd, p, pDelta, in.dTuFac);
#endif

#ifdef FEEDBACK
            // ##### Apply feedback
            pkdAddFBEnergy(pkd, p, psph, in.dConstGamma);
#endif

            // ##### Effective Equation Of State
#ifdef COOLING
            double denMin = in.dCoolingFloorDen;
#ifdef STAR_FORMATION
            double minOverDens = in.dSFMinOverDensity;
#else
            double minOverDens = 57.7;
#endif
            if (pkd->csm->val.bComove) {
                double rhoCrit0 = 3. * pkd->csm->val.dHubble0 * pkd->csm->val.dHubble0 /
                                  (8. * M_PI);
                double denCosmoMin = rhoCrit0 * pkd->csm->val.dOmegab *
                                     minOverDens *
                                     a_inv3; // We do this in proper density

                denMin = ( denCosmoMin > denMin) ? denCosmoMin : denMin;
            }

            if ( (fDensPhys > denMin) &&
                    (psph->Uint < in.dCoolingFlooru*fMass ) ) {
                psph->Uint = in.dCoolingFlooru*fMass;
            }
#endif
#ifdef EEOS_POLYTROPE
            /* Second, the polytropic EoS */
            if (fDensPhys > in.dEOSPolyFloorDen) {

                const double minUint =  fMass * polytropicEnergyFloor(a_inv3, fDens,
                                        in.dEOSPolyFloorIndex, in.dEOSPolyFloorDen,  in.dEOSPolyFlooru);

                if (psph->Uint < minUint) psph->Uint = minUint;
            }
#endif

#ifdef  EEOS_JEANS
            const float f2Ball = 2.*pkdBall(pkd,p);
            const double Ujeans = fMass * jeansEnergyFloor(fDens, f2Ball, in.dConstGamma, in.dEOSNJeans);

            if (psph->Uint < Ujeans)
                psph->Uint = Ujeans;
#endif

            // Actually set the primitive variables
            hydroSetPrimitives(pkd, p, psph, in.dTuFac, in.dConstGamma);


            // Set 'last*' variables for next timestep
            hydroSetLastVars(pkd, p, psph, pa, dScaleFactor, in.dTime, in.dDelta, in.dConstGamma);


        }
        else if (pkdIsBH(pkd,p) && pkdIsActive(pkd,p)) {
#ifdef BLACKHOLES
            pkdBHIntegrate(pkd, p, in.dTime, in.dDelta, in.dBHRadiativeEff);
#endif
        }
    }

}




void pkdLightConeVel(PKD pkd,double dBoxSize) {
    const int nTable=1000;
    const double rMax=3.0;
    gsl_spline *scale;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    double dr,rt[nTable],at_inv[nTable];
    PARTICLE *p;
    vel_t *v;
    double dvFac,r2;
    const double dLightSpeed = dLightSpeedSim(dBoxSize);
    int i,j;

    assert(pkd->oFieldOffset[oVelocity]);
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
    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        v = pkdVel(pkd,p);

        r2 = 0.0;
        for (j=0; j<3; ++j) {
            r2 += pkdPos(pkd,p,j)*pkdPos(pkd,p,j);
        }
        /*
        ** Use r -> 1/a spline table.
        */
        dvFac = gsl_spline_eval(scale,sqrt(r2),acc);
        /* input velocities are momenta p = a^2*x_dot and we want v_pec = a*x_dot */
        for (j=0; j<3; ++j) v[j] *= dvFac;
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
    PARTICLE *p;
    vel_t *v;
    float *a;
    int i,j,n;

    assert(pkd->oFieldOffset[oVelocity]);
    assert(pkd->oFieldOffset[oAcceleration]);

    if (bDoGas) {
        if (1 /*pkd->param.bMeshlessHydro*/) {
            assert(pkd->oFieldOffset[oSph]);
            n = pkd->Local();
            for (i=0; i<n; ++i) {

                p = pkd->Particle(i);
                if (pkdIsRungRange(p,uRungLo,uRungHi)) {
                    v = pkdVel(pkd,p);
                    if (pkdIsGas(pkd,p)) {
                        //SPHFIELDS *sph = pkdSph(pkd,p);
                        //sph->vPred[0] = v[0];
                        //sph->vPred[1] = v[1];
                        //sph->vPred[2] = v[2];
                    }
                    else {
                        a = pkdAccel(pkd,p);
                        for (j=0; j<3; ++j) {
                            v[j] += a[j]*dDelta;
                        }
                    }
                }
            }
        }
        else {
            assert(pkd->oFieldOffset[oSph]);
            n = pkd->Local();
            for (i=0; i<n; ++i) {
                p = pkd->Particle(i);
                if (pkdIsRungRange(p,uRungLo,uRungHi)) {
                    a = pkdAccel(pkd,p);
                    v = pkdVel(pkd,p);
                    if (pkdIsGas(pkd,p)) {
#ifndef OPTIM_REMOVE_UNUSED
                        SPHFIELDS *sph = pkdSph(pkd,p);
                        for (j=0; j<3; ++j) { /* NB: Pred quantities must be done before std. */
                            sph->vPred[j] = v[j] + a[j]*dDeltaVPred;
                        }
                        sph->uPred = sph->u + sph->uDot*dDeltaUPred;
                        sph->u += sph->uDot*dDeltaU;
                        sph->fMetalsPred = sph->fMetals + sph->fMetalsDot*dDeltaUPred;
                        sph->fMetals += sph->fMetalsDot*dDeltaU;
#endif
                    }
                    for (j=0; j<3; ++j) {
                        v[j] += a[j]*dDelta;
                    }
                }
                else {
                    a = pkdAccel(pkd,p);
                    for (j=0; j<3; ++j) {
                        v[j] += a[j]*dDelta;
                    }
                }
            }
        }
    }
    else {
        n = pkd->Local();
        for (i=0; i<n; ++i) {
            p = pkd->Particle(i);
            if (pkdIsRungRange(p,uRungLo,uRungHi)) {
                a = pkdAccel(pkd,p);
                v = pkdVel(pkd,p);
                for (j=0; j<3; ++j) {
                    v[j] += a[j]*dDelta;
                }
            }
        }
    }

    mdlDiag(pkd->mdl, "Done pkdkick\n");
}

/* Kick the tree at iRoot. */
void pkdKickTree(PKD pkd,double dTime,double dDelta,double dDeltaVPred,double dDeltaU,double dDeltaUPred,int iRoot) {
    KDN *c;
    PARTICLE *p;
    vel_t *v;
    float *a;
    int i,j;

    /* Skip to local tree */
    c = pkd->TreeNode(iRoot);
    while (c->bRemote) c = pkd->TreeNode(iRoot = c->iLower);

    /* Now just kick all of the particles in the tree */
    for (i=c->pLower; i<=c->pUpper; ++i) {
        p = pkd->Particle(i);
        a = pkdAccel(pkd,p);
        v = pkdVel(pkd,p);
        for (j=0; j<3; ++j) {
            v[j] += a[j]*dDelta;
            a[j] = 0.0;
        }
    }
}

void pkdInitCosmology(PKD pkd, struct csmVariables *cosmo) {
    /*
    ** Need to be careful to correctly copy the cosmo
    ** parameters. This is very ugly!
    */
    if (pkd->csm) csmFinish(pkd->csm);
    csmInitialize(&pkd->csm);
    pkd->csm->val = *cosmo;
    if (pkd->csm->val.classData.bClass) {
        csmClassGslInitialize(pkd->csm);
    }
}

void pkdZeroNewRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung) {  /* JW: Ugly -- need to clean up */
    PARTICLE *p;
    int i;

    if (!pkd->bNoParticleOrder) for (i=0; i<pkd->Local(); ++i) {
            p = pkd->Particle(i);
            if ( !pkdIsActive(pkd,p) ) continue;
            p->uNewRung = 0;
        }
}

void pkdActiveRung(PKD pkd, int iRung, int bGreater) {
    pkd->uMinRungActive = iRung;
    pkd->uMaxRungActive = bGreater ? 255 : iRung;
}

void pkdCountRungs(PKD pkd,uint64_t *nRungs) {
    PARTICLE *p;
    int i;
    for (i=0; i<=MAX_RUNG; ++i) nRungs[i] = 0;

    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        ++nRungs[p->uRung];
    }
    for (i=0; i<=MAX_RUNG; ++i) pkd->nRung[i] = nRungs[i];
}

void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
                  double dDelta, int iMaxRung,
                  double dEta,double dVelFac,double dAccFac,
                  int bDoGravity,int bEpsAcc,double dhMinOverSoft) {
    PARTICLE *p;
    float *a;
    vel_t *v;
    int i,uNewRung;
    double vel;
    double acc;
    int j;
    double dT;
    double fSoft;

    assert(pkd->oFieldOffset[oVelocity]);
    assert(pkd->oFieldOffset[oAcceleration]);
    assert(!pkd->bNoParticleOrder);

    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        if (pkdIsActive(pkd,p)) {
            v = pkdVel(pkd,p);
            a = pkdAccel(pkd,p);
            fSoft = pkdSoft(pkd,p);
            vel = 0;
            acc = 0;
            for (j=0; j<3; j++) {
                vel += v[j]*v[j];
                acc += a[j]*a[j];
            }
            mdlassert(pkd->mdl,vel >= 0);
            vel = sqrt(vel)*dVelFac;
            mdlassert(pkd->mdl,acc >= 0);
            acc = sqrt(acc)*dAccFac;
            dT = FLOAT_MAXVAL;
            if (acc>0) {
                if (bEpsAcc) {
                    dT = dEta*sqrt(fSoft/acc);
                }
            }
            uNewRung = pkdDtToRung(dT,dDelta,iMaxRung);
            if (uNewRung > p->uNewRung) p->uNewRung = uNewRung;
        }
    }
}


void pkdSphStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
                double dDelta, int iMaxRung,double dEta, double dAccFac, double dEtaUDot) {
#ifndef OPTIM_REMOVE_UNUSED
    int i,j,uNewRung;
    double acc;
    double dtNew;
    int u1,u2,u3;

    assert(pkd->oFieldOffset[oAcceleration]);
    assert(pkd->oFieldOffset[oSph]);
    assert(!pkd->bNoParticleOrder);

    for (i=0; i<pkd->Local(); ++i) {
        PARICLE *p = pkd->Particle(i);
        if (pkdIsActive(pkd,p)) {
            if (pkdIsGas(pkd,p)) {
                u1 = p->uNewRung;
                FLOAT *a = pkdAccel(pkd,p);
                acc = 0;
                for (j=0; j<3; j++) {
                    acc += a[j]*a[j];
                }
                acc = sqrt(acc)*dAccFac;
                dtNew = FLOAT_MAXVAL;
                if (acc>0) dtNew = dEta*sqrt(pkdBall(pkd,p)/acc);
                u2 = pkdDtToRung(dtNew,dDelta,iMaxRung);
                float uDot = *pkd_uDot(pkd,p);
                u3=0;
                if (uDot < 0) {
                    double dtemp = dEtaUDot*(*pkd_u(pkd,p))/fabs(uDot);
                    if (dtemp < dtNew) dtNew = dtemp;
                    u3 = pkdDtToRung(dtemp,dDelta,iMaxRung);
                }
                uNewRung = pkdDtToRung(dtNew,dDelta,iMaxRung);
                if (uNewRung > p->uNewRung) p->uNewRung = uNewRung;
                if (!(p->iOrder%10000) || (p->uNewRung > 5 && !(p->iOrder%1000))) {
                    /*SPHFIELDS *sph = pkdSph(pkd,p);*/
                }
            }
        }
    }
#endif //OPTIM_REMOVE_UNUSED
}


void pkdChemCompInit(PKD pkd, struct inChemCompInit in) {

    for (int i = 0; i < pkd->Local(); ++i) {
        PARTICLE *p = pkd->Particle( i);

        if (pkdIsGas(pkd, p)) {
            SPHFIELDS *pSph = pkdSph(pkd,p);
            float fMass = pkdMass(pkd, p);

            if (pSph->afElemMass[ELEMENT_H] < 0.0f) {
                pSph->afElemMass[ELEMENT_H]  = in.dInitialH  * fMass;
#ifdef HAVE_HELIUM
                pSph->afElemMass[ELEMENT_He] = in.dInitialHe * fMass;
#endif
#ifdef HAVE_CARBON
                pSph->afElemMass[ELEMENT_C]  = in.dInitialC  * fMass;
#endif
#ifdef HAVE_NITROGEN
                pSph->afElemMass[ELEMENT_N]  = in.dInitialN  * fMass;
#endif
#ifdef HAVE_OXYGEN
                pSph->afElemMass[ELEMENT_O]  = in.dInitialO  * fMass;
#endif
#ifdef HAVE_NEON
                pSph->afElemMass[ELEMENT_Ne] = in.dInitialNe * fMass;
#endif
#ifdef HAVE_MAGNESIUM
                pSph->afElemMass[ELEMENT_Mg] = in.dInitialMg * fMass;
#endif
#ifdef HAVE_SILICON
                pSph->afElemMass[ELEMENT_Si] = in.dInitialSi * fMass;
#endif
#ifdef HAVE_IRON
                pSph->afElemMass[ELEMENT_Fe] = in.dInitialFe * fMass;
#endif
            }
#ifdef HAVE_METALLICITY
            if (pSph->fMetalMass < 0.0f)
                pSph->fMetalMass = in.dInitialMetallicity * fMass;
#endif
        }

#ifdef STELLAR_EVOLUTION
        else if (pkdIsStar(pkd, p)) {
            STARFIELDS *pStar = pkdStar(pkd, p);

            if (pStar->afElemAbun[ELEMENT_H] < 0.0f) {
                pStar->afElemAbun[ELEMENT_H]  = in.dInitialH;
#ifdef HAVE_HELIUM
                pStar->afElemAbun[ELEMENT_He] = in.dInitialHe;
#endif
#ifdef HAVE_CARBON
                pStar->afElemAbun[ELEMENT_C]  = in.dInitialC;
#endif
#ifdef HAVE_NITROGEN
                pStar->afElemAbun[ELEMENT_N]  = in.dInitialN;
#endif
#ifdef HAVE_OXYGEN
                pStar->afElemAbun[ELEMENT_O]  = in.dInitialO;
#endif
#ifdef HAVE_NEON
                pStar->afElemAbun[ELEMENT_Ne] = in.dInitialNe;
#endif
#ifdef HAVE_MAGNESIUM
                pStar->afElemAbun[ELEMENT_Mg] = in.dInitialMg;
#endif
#ifdef HAVE_SILICON
                pStar->afElemAbun[ELEMENT_Si] = in.dInitialSi;
#endif
#ifdef HAVE_IRON
                pStar->afElemAbun[ELEMENT_Fe] = in.dInitialFe;
#endif
            }
            if (pStar->fMetalAbun < 0.0f)
                pStar->fMetalAbun = in.dInitialMetallicity;
        }
#endif //STELLAR_EVOLUTION
    }
}

void pkdCorrectEnergy(PKD pkd, double dTuFac, double z, double dTime, int iDirection ) {
    /*PARTICLE *p;
    SPHFIELDS *sph;
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
    PARTICLE *p;
    int i;
    double dT;

    assert(!pkd->bNoParticleOrder);
    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        if (pkdIsActive(pkd,p)) {
            dT = dEta/sqrt(pkdDensity(pkd,p)*dRhoFac);
            p->uNewRung = pkdDtToRung(dT,dDelta,iMaxRung);
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
    PARTICLE *p;
    int i;
    int iTempRung;
    assert(!pkd->bNoParticleOrder);
    for (i=0; i<iMaxRung; ++i) nRungCount[i] = 0;
    for (i=0; i<pkd->Local(); ++i) {
        p = pkd->Particle(i);
        if ( pkdIsActive(pkd,p) ) {
            if ( p->uNewRung > iMaxRung ) p->uNewRung = iMaxRung;
            if ( p->uNewRung >= uRung ) p->uRung = p->uNewRung;
            else if ( p->uRung > uRung) p->uRung = uRung;
        }
        /*
        ** Now produce a count of particles in rungs.
        */
        nRungCount[p->uRung] += 1;
    }
    iTempRung = iMaxRung;
    while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    return iTempRung;
}

void pkdDeleteParticle(PKD pkd, PARTICLE *p) {
    /* p->iOrder = -2 - p->iOrder; JW: Not needed -- just preserve iOrder */
    int pSpecies = pkdSpecies(pkd,p);
    //pkdSetClass(pkd,pkdMass(pkd,p),pkdSoft(pkd,p),FIO_SPECIES_LAST,p); /* Special "DELETED" class == FIO_SPECIES_LAST */
    pkdSetClass(pkd,0.0,0.0,FIO_SPECIES_LAST,p); /* Special "DELETED" class == FIO_SPECIES_LAST */

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
    p->bMarked = 0;
}


/* IA: We replace the deleted particles with those at the end of the particle array (which are still valid), making the tree
 *  no longer usable, unless *extreme* care is taken.
 *
 *  We also update the number of particles in the master level, assuming that the values in PKD are correct (i.e., pkdDeleteParticle was called correctly)
 */
void pkdMoveDeletedParticles(PKD pkd, total_t *n, total_t *nGas, total_t *nDark, total_t *nStar, total_t *nBH) {
    int nLocal = pkd->Local();

    //printf("Deleting particles...\n");
    for (int i=nLocal-1; i>=0; i--) {
        PARTICLE *p = pkd->Particle(i);
        //printf("Checking %d \t %d \t",i, p);
        if (pkdIsDeleted(pkd,p)) {
            //printf("Deleting %d \n", i);
            PARTICLE *lastp = pkd->Particle( pkd->Local() - 1);

            // IA: We can encounter a case where the last particle is deleted; workaround:
            int back_position=0;
            while (pkdIsDeleted(pkd,lastp)) {
                //lastp--;  We can't do this because sizeof(PARTICLE) != pkd->ParticleSize() !!!!
                back_position++;
                lastp = pkd->Particle( pkd->Local()-1-back_position);
            }
            //printf("Had to go back to %d (back_position %d) %d \n", pkd->Local()-1-back_position,back_position, lastp);
            if (lastp>=p) { // Double check that this has the expected behaviour

                //printf("This is at the right of/at i\n");
                assert(back_position==0);
                // If we start from the end, back_position should always? be zero
                pkd->SetLocal(pkd->Local() - back_position+1);
                pkd->CopyParticle( p, lastp);

                assert(!pkdIsDeleted(pkd,p));
            }
            else { // All the particles between p (or more to the left) and the end of the array are being deleted, so no copying is needed
                //printf("This is at the left of i; ");
                pkd->SetLocal(pkd->Local() - back_position);
                i = pkd->Local();
                //printf("New nLocal %d \n", pkd->Local());
            }

            assert(lastp > pkd->Particle(0));


        }

    }

    pkd->TreeNode(ROOT)->pUpper = pkd->Local() - 1;

    *n  = pkd->Local();
    //printf("Final nLocal %d \n", *n);
    *nGas = pkd->nGas;
    *nDark = pkd->nDark;
    *nStar = pkd->nStar;
    *nBH = pkd->nBH;

}




void pkdNewParticle(PKD pkd, PARTICLE *p) {
    PARTICLE *newp;

    mdlassert(pkd->mdl,pkd->Local() < pkd->FreeStore());
    newp = pkd->Particle(pkd->Local());
    pkd->CopyParticle(newp,p);
    newp->iOrder = IORDERMAX;
    pkd->AddLocal(1);
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

void
pkdGetNParts(PKD pkd, struct outGetNParts *out ) {
    int pi;
    int n;
    int nGas;
    int nDark;
    int nStar;
    int nBH;
    total_t iMaxOrder;
    total_t iOrder;
    PARTICLE *p;

    n = 0;
    nGas = 0;
    nDark = 0;
    nStar = 0;
    nBH = 0;
    iMaxOrder = 0;
    for (pi = 0; pi < pkd->Local(); pi++) {
        p = pkd->Particle(pi);
        iOrder = p->iOrder;
        if (iOrder>iMaxOrder) iMaxOrder = iOrder;
        n++;
        if (pkdIsGas(pkd, p)) {
            ++nGas;
        }
        else if (pkdIsDark(pkd, p)) {
            ++nDark;
        }
        else if (pkdIsStar(pkd, p)) {
            ++nStar;
        }
        else if (pkdIsBH(pkd, p)) {
            ++nBH;
        }
    }

    out->n  = n;
    out->nGas = nGas;
    out->nDark = nDark;
    out->nStar = nStar;
    out->nBH = nBH;
    out->nMaxOrder = iMaxOrder;

    pkdSetNParts(pkd, nGas, nDark, nStar, nBH);
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

void pkdInitRelaxation(PKD pkd) {
    assert(pkd->oFieldOffset[oRelaxation]);
    for (int i=0; i<pkd->Local(); ++i) {
        PARTICLE *p = pkd->Particle(i);
        auto pRelax = static_cast<double *>(pkdField(p,pkd->oFieldOffset[oRelaxation]));
        *pRelax = 0.0;
    }
}

double pkdTotalMass(PKD pkd) {
    PARTICLE *p;
    double m;
    int i,n;

    m = 0.0;
    n = pkd->Local();
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        m += pkdMass(pkd,p);
    }
    return m;
}

uint8_t pkdGetMinDt(PKD pkd) {
    PARTICLE *p;
    uint8_t minDt = 0;
    int i, n;

    n = pkd->Local();
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        if (minDt < p->uNewRung) minDt = p->uNewRung;
    }
    return minDt;
}


void pkdSetGlobalDt(PKD pkd, uint8_t minDt) {
    PARTICLE *p;
    int i, n;

    n = pkd->Local();
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        p->uNewRung = minDt;
    }
}


/*
** This function checks the predicate and returns a new value based on the flags.
** setIfTrue:    >0 -> return true if the predicate is true
**               <0 -> "clear if true"
** clearIfFalse: >0 -> return false if the predicate is false
**               <0 -> "set if false"
** A value of zero for either results in no action for the "IfTrue" or "IfFalse" flags.
** Conflicting options (e.g., setIfTrue and setIfFalse) result in a toggle.
*/
static inline int isSelected( int predicate, int setIfTrue, int clearIfFalse, int value ) {
    int s = ((predicate!=0)&(setIfTrue>0)) | (!(predicate!=0)&(clearIfFalse==0));
    int c = ((predicate!=0)&(setIfTrue==0)) | (!(predicate!=0)&(clearIfFalse>0));
    return (~s&~c&value) | (s&~(c&value));
}

int pkdCountSelected(PKD pkd) {
    int i;
    int n=pkd->Local();
    int N=0;
    for ( i=0; i<n; i++ ) if (pkd->Particle(i)->bMarked) ++N;
    return N;
}

int pkdSelActive(PKD pkd, int setIfTrue, int clearIfFalse) {
    int i;
    int n=pkd->Local();
    PARTICLE *p;
    for ( i=0; i<n; i++ ) {
        p=pkd->Particle(i);
        p->bMarked = isSelected(pkdIsActive(pkd,p), setIfTrue, clearIfFalse, p->bMarked);
    }
    return n;
}
int pkdSelSpecies(PKD pkd,uint64_t mSpecies, int setIfTrue, int clearIfFalse) {
    int i;
    int n=pkd->Local();
    int N=0;
    if (mSpecies&(1<<FIO_SPECIES_ALL)) mSpecies = 0xffffffffu;
    for ( i=0; i<n; i++ ) {
        PARTICLE *p=pkd->Particle(i);
        p->bMarked = isSelected((1<<pkdSpecies(pkd,p)) & mSpecies,setIfTrue,clearIfFalse,p->bMarked);
        if (p->bMarked) ++N;
    }
    return N;
}
int pkdSelMass(PKD pkd,double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double m;
    int i,n,nSelected;

    n = pkd->Local();
    nSelected = 0;
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        m = pkdMass(pkd,p);
        p->bMarked = isSelected((m >= dMinMass && m <=dMaxMass),setIfTrue,clearIfFalse,p->bMarked);
        if ( p->bMarked ) nSelected++;
    }
    return nSelected;
}
int pkdSelById(PKD pkd,uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    int i,n,nSelected;

    n = pkd->Local();
    nSelected = 0;
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        p->bMarked = isSelected((p->iOrder >= idStart && p->iOrder <= idEnd),setIfTrue,clearIfFalse,p->bMarked);
        if ( p->bMarked ) nSelected++;
    }
    return nSelected;
}
int pkdSelPhaseDensity(PKD pkd,double dMinDensity, double dMaxDensity, int setIfTrue, int clearIfFalse ) {
    int n = pkd->Local();
    int nSelected = 0;
    for ( int i=0; i<n; i++ ) {
        auto p = pkd->Particle(i);
        auto pvel = static_cast<VELSMOOTH *>(pkdField(p,pkd->oFieldOffset[oVelSmooth]));
        float density = pkdDensity(pkd,p) * pow(pvel->veldisp2,-1.5);
        p->bMarked = isSelected((density >= dMinDensity && density <=dMaxDensity),setIfTrue,clearIfFalse,p->bMarked);
        if ( p->bMarked ) nSelected++;
    }
    return nSelected;
}
int pkdSelBox(PKD pkd,double *dCenter, double *dSize, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    int i,j,n,nSelected;
    int predicate;

    n = pkd->Local();
    nSelected = 0;
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        predicate = 1;
        for (j=0; j<3; j++ ) {
            double dx = dCenter[j] - pkdPos(pkd,p,j);
            predicate = predicate && dx < dSize[j] && dx >= -dSize[j];
        }
        p->bMarked = isSelected(predicate,setIfTrue,clearIfFalse,p->bMarked);
        if ( p->bMarked ) nSelected++;
    }
    return nSelected;
}
int pkdSelSphere(PKD pkd,double *r, double dRadius, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double d2,dx,dy,dz,dRadius2;
    int i,n,nSelected;

    n = pkd->Local();
    nSelected = 0;
    dRadius2 = dRadius*dRadius;
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        dx = r[0] - pkdPos(pkd,p,0);
        dy = r[1] - pkdPos(pkd,p,1);
        dz = r[2] - pkdPos(pkd,p,2);

        d2 = dx*dx + dy*dy + dz*dz;
        p->bMarked = isSelected((d2<=dRadius2),setIfTrue,clearIfFalse,p->bMarked);
        if ( p->bMarked ) nSelected++;
    }
    return nSelected;
}

/*
** Select all particles that fall inside a cylinder between two points P1 and P2
** with radius dRadius.
*/
int pkdSelCylinder(PKD pkd,double *dP1, double *dP2, double dRadius,
                   int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double dCyl[3], dPart[3];
    double dLength2, dRadius2, dL2;
    double pdotr;
    int i,j,n,nSelected;
    int predicate;

    dRadius2 = dRadius*dRadius;
    dLength2 = 0.0;
    for ( j=0; j<3; j++ ) {
        dCyl[j] = dP2[j] - dP1[j];
        dLength2 += dCyl[j] * dCyl[j];
    }
    n = pkd->Local();
    nSelected = 0;
    for ( i=0; i<n; i++ ) {
        p = pkd->Particle(i);
        pdotr = 0.0;
        for ( j=0; j<3; j++ ) {
            dPart[j] = pkdPos(pkd,p,j) - dP1[j];
            pdotr += dPart[j] * dCyl[j];
        }

        if ( pdotr < 0.0 || pdotr > dLength2 ) predicate = 0;
        else {
            dL2 = dPart[0]*dPart[0] + dPart[1]*dPart[1] + dPart[2]*dPart[2] - pdotr*pdotr/dLength2;
            predicate = (dL2 <= dRadius2);
        }
        p->bMarked = isSelected(predicate,setIfTrue,clearIfFalse,p->bMarked);
        if ( p->bMarked ) nSelected++;
    }
    return nSelected;
}
int pkdSelGroup(PKD pkd, int iGroup, int setIfTrue, int clearIfFalse) {
    int i;
    int n=pkd->Local();
    int N=0;
    for ( i=0; i<n; i++ ) {
        PARTICLE *p=pkd->Particle(i);
        p->bMarked = isSelected(pkdGetGroup(pkd,p)==iGroup,setIfTrue,clearIfFalse,p->bMarked);
        if (p->bMarked) ++N;
    }
    return N;
}
int pkdSelBlackholes(PKD pkd, int setIfTrue, int clearIfFalse) {
    int i;
    int n=pkd->Local();
    int N = 0;
    assert(pkd->oFieldOffset[oStar]);
    for ( i=0; i<n; i++ ) {
        PARTICLE *p=pkd->Particle(i);
        if (pkdIsStar(pkd, p)) {
            STARFIELDS *pStar = pkdStar(pkd,p);
            p->bMarked = isSelected(pStar->fTimer < 0,setIfTrue,clearIfFalse,p->bMarked);
        }
        else p->bMarked = 0;
        if (p->bMarked) ++N;
    }
    return N;
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
    int i,j,d,nOut;

    assert(!pkd->bNoParticleOrder); /* We need particle IDs */

    nOut = 0;
    for (i=0; i<pkd->Local(); ++i) {
        PARTICLE *p = pkd->Particle(i);
        float *pPot = pkdPot(pkd,p);
        for (j=0; j<nIn; ++j) {
            if (ID[j] == p->iOrder) {
                double r0[3];
                vel_t *v = pkdVel(pkd,p);
                pkdGetPos1(pkd,p,r0);
                out[nOut].id = p->iOrder;
                out[nOut].mass = pkdMass(pkd,p);
                out[nOut].phi = pPot ? *pPot : 0.0;
                for (d=0; d<3; ++d) {
                    out[nOut].r[d] = r0[d];
                    out[nOut].v[d] = v[d];
                }
                ++nOut;
                break;
            }
        }
    }
    return nOut;
}
