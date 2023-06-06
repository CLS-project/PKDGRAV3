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
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
    #include <malloc.h>
#endif
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"
#include "hydro/hydro.h"
#ifdef FEEDBACK
    #include "starformation/feedback.h"
#endif
#ifdef BLACKHOLES
    #include "blackhole/merger.h"
    #include "blackhole/evolve.h"
    #include "blackhole/step.h"
#endif
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#include <sys/stat.h>

using blitz::TinyVector;
using blitz::floor;
using blitz::dot;
using blitz::abs;
using blitz::all;
using blitz::any;

const int primes[1000] = {
    2,     3,     5,     7,    11,    13,    17,    19,    23,    29,
    31,    37,    41,    43,    47,    53,    59,    61,    67,    71,
    73,    79,    83,    89,    97,   101,   103,   107,   109,   113,
    127,   131,   137,   139,   149,   151,   157,   163,   167,   173,
    179,   181,   191,   193,   197,   199,   211,   223,   227,   229,
    233,   239,   241,   251,   257,   263,   269,   271,   277,   281,
    283,   293,   307,   311,   313,   317,   331,   337,   347,   349,
    353,   359,   367,   373,   379,   383,   389,   397,   401,   409,
    419,   421,   431,   433,   439,   443,   449,   457,   461,   463,
    467,   479,   487,   491,   499,   503,   509,   521,   523,   541,
    547,   557,   563,   569,   571,   577,   587,   593,   599,   601,
    607,   613,   617,   619,   631,   641,   643,   647,   653,   659,
    661,   673,   677,   683,   691,   701,   709,   719,   727,   733,
    739,   743,   751,   757,   761,   769,   773,   787,   797,   809,
    811,   821,   823,   827,   829,   839,   853,   857,   859,   863,
    877,   881,   883,   887,   907,   911,   919,   929,   937,   941,
    947,   953,   967,   971,   977,   983,   991,   997,  1009,  1013,
    1019,  1021,  1031,  1033,  1039,  1049,  1051,  1061,  1063,  1069,
    1087,  1091,  1093,  1097,  1103,  1109,  1117,  1123,  1129,  1151,
    1153,  1163,  1171,  1181,  1187,  1193,  1201,  1213,  1217,  1223,
    1229,  1231,  1237,  1249,  1259,  1277,  1279,  1283,  1289,  1291,
    1297,  1301,  1303,  1307,  1319,  1321,  1327,  1361,  1367,  1373,
    1381,  1399,  1409,  1423,  1427,  1429,  1433,  1439,  1447,  1451,
    1453,  1459,  1471,  1481,  1483,  1487,  1489,  1493,  1499,  1511,
    1523,  1531,  1543,  1549,  1553,  1559,  1567,  1571,  1579,  1583,
    1597,  1601,  1607,  1609,  1613,  1619,  1621,  1627,  1637,  1657,
    1663,  1667,  1669,  1693,  1697,  1699,  1709,  1721,  1723,  1733,
    1741,  1747,  1753,  1759,  1777,  1783,  1787,  1789,  1801,  1811,
    1823,  1831,  1847,  1861,  1867,  1871,  1873,  1877,  1879,  1889,
    1901,  1907,  1913,  1931,  1933,  1949,  1951,  1973,  1979,  1987,
    1993,  1997,  1999,  2003,  2011,  2017,  2027,  2029,  2039,  2053,
    2063,  2069,  2081,  2083,  2087,  2089,  2099,  2111,  2113,  2129,
    2131,  2137,  2141,  2143,  2153,  2161,  2179,  2203,  2207,  2213,
    2221,  2237,  2239,  2243,  2251,  2267,  2269,  2273,  2281,  2287,
    2293,  2297,  2309,  2311,  2333,  2339,  2341,  2347,  2351,  2357,
    2371,  2377,  2381,  2383,  2389,  2393,  2399,  2411,  2417,  2423,
    2437,  2441,  2447,  2459,  2467,  2473,  2477,  2503,  2521,  2531,
    2539,  2543,  2549,  2551,  2557,  2579,  2591,  2593,  2609,  2617,
    2621,  2633,  2647,  2657,  2659,  2663,  2671,  2677,  2683,  2687,
    2689,  2693,  2699,  2707,  2711,  2713,  2719,  2729,  2731,  2741,
    2749,  2753,  2767,  2777,  2789,  2791,  2797,  2801,  2803,  2819,
    2833,  2837,  2843,  2851,  2857,  2861,  2879,  2887,  2897,  2903,
    2909,  2917,  2927,  2939,  2953,  2957,  2963,  2969,  2971,  2999,
    3001,  3011,  3019,  3023,  3037,  3041,  3049,  3061,  3067,  3079,
    3083,  3089,  3109,  3119,  3121,  3137,  3163,  3167,  3169,  3181,
    3187,  3191,  3203,  3209,  3217,  3221,  3229,  3251,  3253,  3257,
    3259,  3271,  3299,  3301,  3307,  3313,  3319,  3323,  3329,  3331,
    3343,  3347,  3359,  3361,  3371,  3373,  3389,  3391,  3407,  3413,
    3433,  3449,  3457,  3461,  3463,  3467,  3469,  3491,  3499,  3511,
    3517,  3527,  3529,  3533,  3539,  3541,  3547,  3557,  3559,  3571,
    3581,  3583,  3593,  3607,  3613,  3617,  3623,  3631,  3637,  3643,
    3659,  3671,  3673,  3677,  3691,  3697,  3701,  3709,  3719,  3727,
    3733,  3739,  3761,  3767,  3769,  3779,  3793,  3797,  3803,  3821,
    3823,  3833,  3847,  3851,  3853,  3863,  3877,  3881,  3889,  3907,
    3911,  3917,  3919,  3923,  3929,  3931,  3943,  3947,  3967,  3989,
    4001,  4003,  4007,  4013,  4019,  4021,  4027,  4049,  4051,  4057,
    4073,  4079,  4091,  4093,  4099,  4111,  4127,  4129,  4133,  4139,
    4153,  4157,  4159,  4177,  4201,  4211,  4217,  4219,  4229,  4231,
    4241,  4243,  4253,  4259,  4261,  4271,  4273,  4283,  4289,  4297,
    4327,  4337,  4339,  4349,  4357,  4363,  4373,  4391,  4397,  4409,
    4421,  4423,  4441,  4447,  4451,  4457,  4463,  4481,  4483,  4493,
    4507,  4513,  4517,  4519,  4523,  4547,  4549,  4561,  4567,  4583,
    4591,  4597,  4603,  4621,  4637,  4639,  4643,  4649,  4651,  4657,
    4663,  4673,  4679,  4691,  4703,  4721,  4723,  4729,  4733,  4751,
    4759,  4783,  4787,  4789,  4793,  4799,  4801,  4813,  4817,  4831,
    4861,  4871,  4877,  4889,  4903,  4909,  4919,  4931,  4933,  4937,
    4943,  4951,  4957,  4967,  4969,  4973,  4987,  4993,  4999,  5003,
    5009,  5011,  5021,  5023,  5039,  5051,  5059,  5077,  5081,  5087,
    5099,  5101,  5107,  5113,  5119,  5147,  5153,  5167,  5171,  5179,
    5189,  5197,  5209,  5227,  5231,  5233,  5237,  5261,  5273,  5279,
    5281,  5297,  5303,  5309,  5323,  5333,  5347,  5351,  5381,  5387,
    5393,  5399,  5407,  5413,  5417,  5419,  5431,  5437,  5441,  5443,
    5449,  5471,  5477,  5479,  5483,  5501,  5503,  5507,  5519,  5521,
    5527,  5531,  5557,  5563,  5569,  5573,  5581,  5591,  5623,  5639,
    5641,  5647,  5651,  5653,  5657,  5659,  5669,  5683,  5689,  5693,
    5701,  5711,  5717,  5737,  5741,  5743,  5749,  5779,  5783,  5791,
    5801,  5807,  5813,  5821,  5827,  5839,  5843,  5849,  5851,  5857,
    5861,  5867,  5869,  5879,  5881,  5897,  5903,  5923,  5927,  5939,
    5953,  5981,  5987,  6007,  6011,  6029,  6037,  6043,  6047,  6053,
    6067,  6073,  6079,  6089,  6091,  6101,  6113,  6121,  6131,  6133,
    6143,  6151,  6163,  6173,  6197,  6199,  6203,  6211,  6217,  6221,
    6229,  6247,  6257,  6263,  6269,  6271,  6277,  6287,  6299,  6301,
    6311,  6317,  6323,  6329,  6337,  6343,  6353,  6359,  6361,  6367,
    6373,  6379,  6389,  6397,  6421,  6427,  6449,  6451,  6469,  6473,
    6481,  6491,  6521,  6529,  6547,  6551,  6553,  6563,  6569,  6571,
    6577,  6581,  6599,  6607,  6619,  6637,  6653,  6659,  6661,  6673,
    6679,  6689,  6691,  6701,  6703,  6709,  6719,  6733,  6737,  6761,
    6763,  6779,  6781,  6791,  6793,  6803,  6823,  6827,  6829,  6833,
    6841,  6857,  6863,  6869,  6871,  6883,  6899,  6907,  6911,  6917,
    6947,  6949,  6959,  6961,  6967,  6971,  6977,  6983,  6991,  6997,
    7001,  7013,  7019,  7027,  7039,  7043,  7057,  7069,  7079,  7103,
    7109,  7121,  7127,  7129,  7151,  7159,  7177,  7187,  7193,  7207,
    7211,  7213,  7219,  7229,  7237,  7243,  7247,  7253,  7283,  7297,
    7307,  7309,  7321,  7331,  7333,  7349,  7351,  7369,  7393,  7411,
    7417,  7433,  7451,  7457,  7459,  7477,  7481,  7487,  7489,  7499,
    7507,  7517,  7523,  7529,  7537,  7541,  7547,  7549,  7559,  7561,
    7573,  7577,  7583,  7589,  7591,  7603,  7607,  7621,  7639,  7643,
    7649,  7669,  7673,  7681,  7687,  7691,  7699,  7703,  7717,  7723,
    7727,  7741,  7753,  7757,  7759,  7789,  7793,  7817,  7823,  7829,
    7841,  7853,  7867,  7873,  7877,  7879,  7883,  7901,  7907,  7919
};


/*
** Assumes that p does not already occur in the hash table!!!
*/
void smHashAdd(SMX smx,void *p) {
    struct hashElement *t;
    uint32_t i = ((intptr_t)(p))%smx->nHash;
    if (!smx->pHash[i].p) {
        smx->pHash[i].p = p;
    }
    else {
        t = smx->pFreeHash;
        assert(t != NULL);
        smx->pFreeHash = t->coll;
        t->coll = smx->pHash[i].coll;
        smx->pHash[i].coll = t;
        t->p = p;
    }
}

/*
** Assumes that p is definitely in the hash table!!!
*/
void smHashDel(SMX smx,void *p) {
    struct hashElement *t,*tt;
    uint32_t i = ((intptr_t)(p))%smx->nHash;

    if (!smx->pHash[i].coll) {
        /*
        ** It has to be the first element.
        */
        smx->pHash[i].p = NULL;
    }
    else if (smx->pHash[i].p == p) {
        /*
        ** It is the first element, but there are others!
        */
        t = smx->pHash[i].coll;
        smx->pHash[i].coll = t->coll;
        smx->pHash[i].p = t->p;
        t->coll = smx->pFreeHash;
        smx->pFreeHash = t;
    }
    else {
        tt = &smx->pHash[i];
        while (tt->coll->p != p) tt = tt->coll;
        t = tt->coll;
        tt->coll = t->coll; /* unlink */
        t->coll = smx->pFreeHash;
        smx->pFreeHash = t;
    }
}


int smHashPresent(SMX smx,void *p) {
    struct hashElement *t;
    uint32_t i = ((intptr_t)(p))%smx->nHash;

    if (smx->pHash[i].p == p) return 1;
    t = smx->pHash[i].coll;
    while (t) {
        if (t->p == p) return 1;
        else t = t->coll;
    }
    return 0;
}


static int smInitializeBasic(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,int bSymmetric,int iSmoothType,int bMakeCache) {
    SMX smx;
    void (*initParticle)(void *,void *) = NULL;
    void (*pack)(void *,void *,const void *) = NULL;
    void (*unpack)(void *,void *,const void *) = NULL;
    void (*init)(void *,void *) = NULL;
    void (*flush)(void *,void *,const void *) = NULL;
    void (*comb)(void *,void *,const void *) = NULL;
    bool bPacked = false;
    uint32_t iPackSize = 0;
    uint32_t iFlushSize = 0;
    int i;

    smx = new struct smContext;
    assert(smx != NULL);
    smx->pSentinel = static_cast<PARTICLE *>(malloc(pkd->particles.ParticleSize()));
    assert(smx->pSentinel != NULL);
    smx->pkd = pkd;
    smx->fcnSmoothNode = NULL;
    smx->fcnSmoothGetNvars = NULL;
    smx->fcnSmoothFillBuffer = NULL;
    smx->fcnSmoothUpdate = NULL;
    if (smf != NULL) smf->pkd = pkd;
    smx->nSmooth = nSmooth;
    smx->bPeriodic = bPeriodic;
    smx->bSymmetric  = bSymmetric;
    smx->bSearchGasOnly = 0;
    smx->iSmoothType = iSmoothType;
    /*
    ** Initialize the context for compressed nearest neighbor lists.
    */
    smx->lcmp = lcodeInit(pkd->Threads(),pkd->Self(),pkd->Local(),nSmooth);

    switch (iSmoothType) {
    case SMX_NULL:
        smx->fcnSmooth = NullSmooth;
        initParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = NULL;
        break;
    case SMX_DENSITY:
        smx->fcnSmooth = bSymmetric?DensitySym:Density;
        initParticle = initDensity; /* Original Particle */
        init = initDensity; /* Cached copies */
        comb = combDensity;
        break;
    case SMX_DENSITY_F1:
        assert(!bSymmetric);
        smx->fcnSmooth = DensityF1;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        break;
    case SMX_DENSITY_M3:
        assert(!bSymmetric);
        assert(pkd->particles.present(PKD_FIELD::oGroup));
        smx->fcnSmooth = DensityM3;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        break;
    case SMX_GRADIENT_M3:
        assert(!bSymmetric);
        smx->fcnSmooth = LinkGradientM3;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        break;
    case SMX_HOP_LINK:
        assert(!bSymmetric);
        smx->fcnSmooth = LinkHopChains;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        break;
    case SMX_PRINTNN:
        smx->fcnSmooth = PrintNN;
        initParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = NULL;
        break;
    case SMX_HYDRO_DENSITY:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroDensity;
        initParticle = NULL; /* Original Particle */
        pack = packHydroDensity;
        unpack = unpackHydroDensity;
        init = NULL; /* Cached copies */
        comb = NULL;
        bPacked = true;
        iPackSize = sizeof(hydroDensityPack);
        break;
    case SMX_HYDRO_DENSITY_FINAL:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroDensityFinal;
        initParticle = NULL; /* Original Particle */
        pack = packHydroDensity;
        unpack = unpackHydroDensity;
        init = NULL; /* Cached copies */
        comb = NULL;
        bPacked = true;
        iPackSize = sizeof(hydroDensityPack);
        break;
    case SMX_HYDRO_GRADIENT:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroGradients;
        initParticle = NULL; /* Original Particle */
        pack = packHydroGradients;
        unpack = unpackHydroGradients;
        init = NULL; /* Cached copies */
        comb = NULL;
        bPacked = true;
        iPackSize = sizeof(hydroGradientsPack);
        break;
    case SMX_HYDRO_FLUX:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroRiemann;
        initParticle = initHydroFluxes; /* Original Particle */
        pack = packHydroFluxes;
        unpack = unpackHydroFluxes;
        init = initHydroFluxesCached; /* Cached copies */
        flush = flushHydroFluxes;
        comb = combHydroFluxes;
        bPacked = true;
        iPackSize = sizeof(hydroFluxesPack);
        iFlushSize = sizeof(hydroFluxesFlush);
        break;
#ifdef OPTIM_FLUX_VEC
    case SMX_HYDRO_FLUX_VEC:
        assert (pkd->particles.present(PKD_FIELD::oSph));
        smx->fcnSmoothNode = hydroRiemann_vec;
        smx->fcnSmoothGetNvars = hydroFluxGetNvars;
        smx->fcnSmoothFillBuffer = hydroFluxFillBuffer;
        smx->fcnSmoothUpdate = hydroFluxUpdateFromBuffer;
        initParticle = initHydroFluxes; /* Original Particle */
        pack = packHydroFluxes;
        unpack = unpackHydroFluxes;
        init = initHydroFluxesCached; /* Cached copies */
        flush = flushHydroFluxes;
        comb = combHydroFluxes;
        bPacked = true;
        iPackSize = sizeof(hydroFluxesPack);
        iFlushSize = sizeof(hydroFluxesFlush);
        break;
#endif
    case SMX_HYDRO_STEP:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroStep;
        initParticle = NULL; /* Original Particle */
        pack = packHydroStep;
        unpack = unpackHydroStep;
        init = initHydroStep; /* Cached copies */
        flush = flushHydroStep;
        comb = combHydroStep;
        bPacked = true;
        iPackSize = sizeof(hydroStepPack);
        iFlushSize = sizeof(hydroStepFlush);
        break;
#ifdef FEEDBACK
    case SMX_SN_FEEDBACK:
        smx->fcnSmooth = smSNFeedback;
        initParticle = NULL; /* Original Particle */
        pack = packSNFeedback;
        unpack = unpackSNFeedback;
        init = initSNFeedback; /* Cached copies */
        flush = flushSNFeedback;
        comb = combSNFeedback;
        smx->bSearchGasOnly = 1;
        bPacked = true;
        iPackSize = sizeof(snFeedbackPack);
        iFlushSize = sizeof(snFeedbackFlush);
        break;
#endif
#ifdef BLACKHOLES
    case SMX_BH_MERGER:
        smx->fcnSmooth = smBHmerger;
        initParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = combBHmerger;
        break;
    case SMX_BH_STEP:
        smx->fcnSmooth = smBHstep;
        initParticle = NULL; /* Original Particle */
        pack = packBHstep;
        unpack = unpackBHstep;
        init = NULL; /* Cached copies */
        comb = NULL;
        smx->bSearchGasOnly = 1;
        bPacked = true;
        iPackSize = sizeof(bhStepPack);
        break;
    case SMX_BH_DRIFT:
        smx->fcnSmooth = smBHevolve;
        initParticle = NULL; /* Original Particle */
        pack = packBHevolve;
        unpack = unpackBHevolve;
        init = initBHevolve; /* Cached copies */
        flush = flushBHevolve;
        comb = combBHevolve;
        smx->bSearchGasOnly = 1;
        bPacked = true;
        iPackSize = sizeof(bhEvolvePack);
        iFlushSize = sizeof(bhEvolveFlush);
        break;
#endif
#ifdef STELLAR_EVOLUTION
    case SMX_CHEM_ENRICHMENT:
        smx->fcnSmooth = smChemEnrich;
        initParticle = NULL;
        pack = packChemEnrich;
        unpack = unpackChemEnrich;
        init = initChemEnrich;
        flush = flushChemEnrich;
        comb = combChemEnrich;
        smx->bSearchGasOnly = 1;
        bPacked = true;
        iPackSize = sizeof(stevPack);
        iFlushSize = sizeof(stevFlush);
        break;
#endif
    case SMX_BALL:
        assert(pkd->particles.present(PKD_FIELD::oBall));
        smx->fcnSmooth = BallSmooth;
        initParticle = initBall; /* Original Particle */
        init = initBall; /* Cached copies */
        comb = NULL;
        break;

    default:
        assert(0);
    }

    if ( (smx->fcnSmoothNode != NULL) && ( (smx->fcnSmoothGetNvars == NULL) ||
                                           (smx->fcnSmoothFillBuffer == NULL) || (smx->fcnSmoothUpdate == NULL) ) ) {
        fprintf(stderr, "ERROR: Trying to use particle buffer in node smooth,"
                "but not all the required fuctions are set\n");
        abort();
    }
    /*
    ** Initialize the ACTIVE particles in the tree.
    ** There are other particles in the tree -- just not active.
    */
    auto pRoot = pkd->tree[ROOT];
    auto nTree = pRoot->count();
    if (initParticle != NULL) {
        for (auto &p : *pRoot) {
            if (p.is_active()) initParticle(pkd,&p);
        }
    }
    /*
    ** Start particle caching space (cell cache is already active).
    */
    if (bMakeCache) {
        smx->bOwnCache = 1;
        if (bPacked) {
            if (bSymmetric) {
                mdlPackedCacheCO(pkd->mdl,CID_PARTICLE,NULL,pkd->particles,nTree,
                                 pkd->particles.ParticleSize(),pkd,iPackSize,
                                 pack,unpack,iFlushSize,init,flush,comb);
            }
            else {
                mdlPackedCacheRO(pkd->mdl,CID_PARTICLE,NULL,pkd->particles,nTree,
                                 pkd->particles.ParticleSize(),pkd,iPackSize,
                                 pack,unpack);
            }
        }
        else if (bSymmetric) {
            mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,
                       pkd->particles,pkd->particles.ParticleSize(),
                       nTree,pkd,init,comb);
        }
        else {
            mdlROcache(pkd->mdl,CID_PARTICLE,NULL,
                       pkd->particles,pkd->particles.ParticleSize(),
                       nTree);
        }
    }
    else smx->bOwnCache = 0;
    /*
    ** Allocate Nearest-Neighbor List.
    */
    smx->nnListSize = 0;
    smx->nnListMax = NNLIST_INCREMENT;
    smx->nnList = static_cast<NN *>(malloc(smx->nnListMax*sizeof(NN)));
    assert(smx->nnList != NULL);

    /*
    ** Allocate priority queue.
    */
    smx->pq = static_cast<PQ *>(malloc(nSmooth*sizeof(PQ)));
    assert(smx->pq != NULL);
    PQ_INIT(smx->pq,nSmooth);
    /*
    ** Allocate hash table entries.
    ** The constant here just sets the hash table loading factor, for numbers larger than
    ** the 1000'th prime we end up using the result here as the hash table modulus.
    */
    smx->nHash = (int)floor(nSmooth*1.543765241931);
    for (i=0; i<1000; ++i) {
        if (primes[i] > smx->nHash) {
            smx->nHash = primes[i];
            break;
        }
    }
    smx->pHash = static_cast<struct hashElement *>(malloc((smx->nHash+nSmooth)*sizeof(struct hashElement)));
    assert(smx->pHash != NULL);
    for (i=0; i<smx->nHash; ++i) {
        smx->pHash[i].p = NULL;
        smx->pHash[i].coll = NULL;
    }
    /*
    ** set up the extra entries that may be needed for collision chains
    */
    smx->pFreeHash = &smx->pHash[i];
    for (; i<(smx->nHash+nSmooth-1); ++i) {
        smx->pHash[i].p = NULL;
        smx->pHash[i].coll = &smx->pHash[i+1];
    }
    smx->pHash[i].p = NULL;
    smx->pHash[i].coll = NULL;
    /*
    ** Allocate special stacks for searching within the tree.
    ** 1024 is more than enough.
    */
    smx->ST = new struct smContext::stStack[1024];
    assert(smx->ST != NULL);
    smx->S = new int[1024];
    assert(smx->S != NULL);
    /*
    ** Set up the sentinel particle with some very far away distance.
    ** This is used to initially load the priority queue and all references
    ** to this particle should end up being replaced in the priority queue
    ** as long as there are nSmooth particles.
    */
    auto sentinel = pkd->particles[smx->pSentinel];
    if (pkd->bIntegerPosition) sentinel.raw_position<int32_t>() = TinyVector<int32_t,3>(INT32_MAX,INT32_MAX,INT32_MAX);
    else sentinel.raw_position<double>() = TinyVector<double,3>(HUGE_VAL,HUGE_VAL,HUGE_VAL);
    /*
    ** Need to cast the pLite to an array of extra stuff.
    */
    assert(pkd->EphemeralBytes() >= sizeof(struct smExtraArray));
    smx->ea = (struct smExtraArray *)(pkd->pLite); /* Used only for SPH */
    *psmx = smx;
    return (1);
}

int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,int bSymmetric,int iSmoothType) {
    return smInitializeBasic(psmx,pkd,smf,nSmooth,bPeriodic,bSymmetric,iSmoothType,1);
}

int smInitializeRO(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,int iSmoothType) {
    return smInitializeBasic(psmx,pkd,smf,nSmooth,bPeriodic,0,iSmoothType,0);
}

void smFinish(SMX smx,SMF *smf) {
    char achOut[128];

    /*
     * Output statistics.
     */
    snprintf(achOut, sizeof(achOut), "Cell Accesses: %g\n",
             mdlNumAccess(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    snprintf(achOut, sizeof(achOut), "    Miss ratio: %g\n",
             mdlMissRatio(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    snprintf(achOut, sizeof(achOut), "Particle Accesses: %g\n",
             mdlNumAccess(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    snprintf(achOut, sizeof(achOut), "    Miss ratio: %g\n",
             mdlMissRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    /*
    ** Stop particle caching space.
    */
    if (smx->bOwnCache)
        mdlFinishCache(smx->pkd->mdl,CID_PARTICLE);
    /*
    ** Finish compressed lists.
    */
    lcodeFinish(smx->lcmp);
    /*
    ** Free up context storage.
    */
    delete [] smx->S;
    delete [] smx->ST;
    free(smx->pq);
    free(smx->nnList);
    free(smx->pHash);
    free(smx->pSentinel);
    delete smx;
}

static auto getCell(PKD pkd, int iCell, int id) {
    return (id==pkd->Self()) ? pkd->tree[iCell]
           : pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,CID_CELL,iCell,id))];
}

PQ *pqSearch(SMX smx,PQ *pq,TinyVector<double,3> r,int iRoot) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    int idSelf = smx->pkd->Self();
    struct smContext::stStack *S = smx->ST;
    int iCell,id;
    int sp = 0;
    // int pEnd, pj;
    double fDist2;
    TinyVector<double,3> dr;

    /* Start at the root node of the tree */
    auto kdn = getCell(pkd,pkd->iTopTree[iRoot],id = idSelf);
    while (1) {
        while (kdn->is_cell()) {
            auto [iLower,idLower,iUpper,idUpper] = kdn->get_child_cells(id);
            kdn = getCell(pkd,iLower,idLower);
            auto min1 = kdn->bound().mindist(r);
            kdn = getCell(pkd,iUpper,idUpper);
            auto min2 = kdn->bound().mindist(r);
            if (min1 < min2) {
                if (min1 > pq->fDist2) goto NoIntersect;
                S[sp].id = idUpper;
                S[sp].iCell = iUpper;
                S[sp].min = min2;
                ++sp;
                id = idLower;
                iCell = iLower;
                kdn = getCell(pkd,iCell,id);
            }
            else {
                if (min2 > pq->fDist2) goto NoIntersect;
                S[sp].id = idLower;
                S[sp].iCell = iLower;
                S[sp].min = min1;
                ++sp;
                id = idUpper;
                iCell = iUpper;
            }
        }
        /* Now at a bucket */
        if (id == idSelf ) {
#ifdef OPTIM_REORDER_IN_NODES
            auto pEnd = (smx->bSearchGasOnly) ? kdn->lower() + kdn->Ngas() : kdn->upper() + 1;
#else
            auto pEnd = kdn->upper() + 1;
#endif
            for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                auto p = pkd->particles[pj];
                if (!p.marked()) continue;
#ifndef OPTIM_REORDER_IN_NODES
                if (smx->bSearchGasOnly && !p.is_gas()) continue;
#endif
                dr = r - p.position();
                fDist2 = dot(dr,dr);
                if (fDist2 <= pq->fDist2) {
                    if (pq->iPid == idSelf) {
                        pkd->particles[pq->iIndex].set_marked(true);
                    }
                    else {
                        smHashDel(smx,pq->pPart);
                        mdlRelease(mdl,CID_PARTICLE,pq->pPart);
                        pq->iPid = idSelf;
                    }
                    pq->pPart = &p;
                    pq->fDist2 = fDist2;
                    pq->dr = dr;
                    pq->iIndex = pj;
                    p.set_marked(false); /* de-activate a particle that enters the queue */
                    PQ_REPLACE(pq);
                }
            }
        }
        else {
#ifdef OPTIM_REORDER_IN_NODES
            auto pEnd = (smx->bSearchGasOnly) ? kdn->lower() + kdn->Ngas() : kdn->upper() + 1;
#else
            auto pEnd = kdn->upper() + 1;
#endif
            for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                auto p = pkd->particles[static_cast<PARTICLE *>(mdlFetch(mdl,CID_PARTICLE,pj,id))];
#ifndef OPTIM_REORDER_IN_NODES
                if (smx->bSearchGasOnly && !p.is_gas()) continue;
#endif
                if (smHashPresent(smx,&p)) continue;
                dr = r - p.position();
                fDist2 = dot(dr,dr);
                if (fDist2 <= pq->fDist2) {
                    if (pq->iPid == idSelf) {
                        pkd->particles[pq->iIndex].set_marked(true);
                    }
                    else {
                        smHashDel(smx,pq->pPart);
                        mdlRelease(mdl,CID_PARTICLE,pq->pPart);
                    }
                    auto p = pkd->particles[static_cast<PARTICLE *>(mdlAcquire(mdl,CID_PARTICLE,pj,id))];
                    pq->pPart = &p;
                    pq->fDist2 = fDist2;
                    pq->dr = dr;
                    pq->iIndex = pj;
                    pq->iPid = id;
                    smHashAdd(smx,&p);
                    PQ_REPLACE(pq);
                }
            }
        }


NoIntersect:
        if (sp) {
            --sp;
            if (S[sp].min > pq->fDist2) goto NoIntersect;
            id = S[sp].id;
            iCell = S[sp].iCell;
            kdn = getCell(pkd,iCell,id);
        }
        else return pq;
    }
}

void smSmoothInitialize(SMX smx) {
    int i;
    /*
    ** Initialize the priority queue first.
    */
    auto sentinel = smx->pkd->particles[smx->pSentinel];
    for (i=0; i<smx->nSmooth; ++i) {
        smx->pq[i].pPart = smx->pSentinel;
        smx->pq[i].iIndex = smx->pkd->Local();
        smx->pq[i].iPid = smx->pkd->Self();
        smx->pq[i].dr = sentinel.position();
        smx->pq[i].fDist2 = dot(smx->pq[i].dr,smx->pq[i].dr);
    }
    for (i=0; i<3; ++i) smx->rLast[i] = 0.0;
}

void smSmoothFinish(SMX smx) {
    int i;
    /*
    ** Release acquired pointers and source-reactivate particles in prioq.
    */
    for (i=0; i<smx->nSmooth; ++i) {
        if (smx->pq[i].iPid == smx->pkd->Self()) {
            smx->pkd->particles[smx->pq[i].iIndex].set_marked(true);
        }
        else {
            smHashDel(smx,smx->pq[i].pPart);
            mdlRelease(smx->pkd->mdl,CID_PARTICLE,smx->pq[i].pPart);
        }
    }
}

float smSmoothSingle(SMX smx,SMF *smf,particleStore::ParticleReference &p,int iRoot1, int iRoot2) {
    PKD pkd = smx->pkd;
    int ix,iy,iz;
    TinyVector<double,3> r;
    double fBall;
    int j;
    PQ *pq;

    auto p_r = p.position();

    /*
    ** Correct distances and rebuild priority queue.
    */
    if (smx->bPeriodic) {
        for (j=0; j<3; ++j) {
            if (p_r[j] > smx->rLast[j] + 0.5 * pkd->fPeriod[j])
                smx->rLast[j] += pkd->fPeriod[j];
            else if (p_r[j] < smx->rLast[j] - 0.5 * pkd->fPeriod[j])
                smx->rLast[j] -= pkd->fPeriod[j];
        }
    }

    for (j=0; j<smx->nSmooth; ++j) {
        smx->pq[j].dr += p_r-smx->rLast;
        smx->pq[j].fDist2 = dot(smx->pq[j].dr,smx->pq[j].dr);
    }
    r = p_r;
    smx->rLast = p_r;

    PQ_BUILD(smx->pq,smx->nSmooth,pq);
    pq = pqSearch(smx,pq,r,iRoot1);
    if (iRoot2) pq = pqSearch(smx,pq,r,iRoot2);
    /*
    ** Search in replica boxes if it is required.
    */
//    printf("fPeriod %f %f %f \n",pkd->fPeriod[0],pkd->fPeriod[1],pkd->fPeriod[2] );
//    printf("x %f y %f z %f \n", p_r[0], p_r[1], p_r[2]);
    if (smx->bPeriodic) {
        fBall = sqrt(pq->fDist2);
        TinyVector<int,3> iStart = floor((p_r - fBall) / pkd->fPeriod + 0.5);
        TinyVector<int,3> iEnd   = floor((p_r + fBall) / pkd->fPeriod + 0.5);
        for (ix=iStart[0]; ix<=iEnd[0]; ++ix) {
            r[0] = p_r[0] - ix*pkd->fPeriod[0];
            for (iy=iStart[1]; iy<=iEnd[1]; ++iy) {
                r[1] = p_r[1] - iy*pkd->fPeriod[1];
                for (iz=iStart[2]; iz<=iEnd[2]; ++iz) {
                    r[2] = p_r[2] - iz*pkd->fPeriod[2];
                    if (ix || iy || iz) {
                        pq = pqSearch(smx,pq,r,iRoot1);
                        if (iRoot2) pq = pqSearch(smx,pq,r,iRoot2);
                    }
                }
            }
        }
    }
    fBall = sqrt(pq->fDist2);

    /* IA: I do not fully understand this fBall, so I will compute my own kernel length such that it encloses
     * all the particles in the neighbor list. This means that h > 0.5*max(dist). I have taken 0.501 as a safe
     * value, because 0.5 would exclude the furthest particle(s) */
//    int i;
//    fBall = 0.0;
//    for (i=0; i<smx->nSmooth; ++i){
//       if (fBall < smx->pq[i].fDist2) fBall = smx->pq[i].fDist2;
//    }
//    fBall = 0.50*sqrt(fBall);


    /*
    ** Apply smooth funtion to the neighbor list.
    */
    smx->fcnSmooth(&p,fBall,smx->nSmooth,smx->pq,smf);
    return fBall;
}

void smSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    float fBall;

    /*
    ** Initialize the bInactive flags for all local particles.
    */
    for (auto &p : pkd->particles) p.set_marked(true);
    smSmoothInitialize(smx);
    smf->pfDensity = NULL;
    switch (smx->iSmoothType) {
#ifdef BLACKHOLES
    case SMX_BH_DRIFT:
        /*
        ** For BH accretion, we use the ephemeral storage to keep track of accreting BHs
        */
        if (smf->bBHAccretion) {
            assert(pkd->EphemeralBytes() >= sizeof(remoteID));
            for (auto i = 0; i < pkd->Local(); ++i) {
                if (pkd->particles[i].is_gas()) {
                    auto &BHAccretor = static_cast<remoteID *>(pkd->pLite)[i];
                    BHAccretor.iPid = NOT_ACCRETED;
                }
            }
            pkd->mdl->CacheInitialize(CID_GROUP,NULL,pkd->pLite,pkd->Local(),
                                      std::make_shared<BHAccretionCache>());
        }

        for (auto &p : pkd->particles) {
            if (p.is_bh()) {
                smSmoothSingle(smx,smf,p,ROOT,0);
            }
        }

        if (smf->bBHAccretion) mdlFinishCache(pkd->mdl,CID_GROUP);
        break;
    case SMX_BH_STEP:
        for (auto &p : pkd->particles) {
            if (p.is_bh()) {
                smSmoothSingle(smx,smf,p,ROOT,0);
            }
        }
        break;
#endif
#ifdef FEEDBACK
    case SMX_SN_FEEDBACK:
        for (auto &p : pkd->particles) {
            if (p.is_star()) {
                auto &star = p.star();

                if (!star.bCCSNFBDone && ((smf->dTime - star.fTimer) > smf->dCCSNFBDelay)) {
                    smSmoothSingle(smx,smf,p, ROOT, 0);
                }

                if (!star.bSNIaFBDone && ((smf->dTime - star.fTimer) > smf->dSNIaFBDelay)) {
                    smSmoothSingle(smx,smf,p, ROOT, 0);
                }
            }
        }
        break;
#endif
#ifdef STELLAR_EVOLUTION
    case SMX_CHEM_ENRICHMENT:
        for (auto &p : pkd->particles) {
            if (p.is_star()) {
                auto &star = p.star();
                if (smf->dTime > star.fNextEnrichTime) {
                    smSmoothSingle(smx, smf, p, ROOT, 0);
                }
            }
        }
        break;
#endif
    default:
        for (auto &p : pkd->particles) {
            if (!smf->bMeshlessHydro ) {
                smSmoothSingle(smx,smf,p,ROOT,0);
                //p.set_ball(smSmoothSingle(smx,smf,p,ROOT,0));
            }
            else {
                if (p.is_active()) {
                    fBall = smSmoothSingle(smx,smf,p,ROOT,0);
                    if (smf->bUpdateBall) {
                        p.set_ball(fBall);
                    }
                }
            }
            /*
            ** Call mdlCacheCheck to make sure we are making progress!
            */
            mdlCacheCheck(pkd->mdl);
        }
    }
    smSmoothFinish(smx);
}

void smGather(SMX smx,double fBall2,TinyVector<double,3> r) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    int idSelf = pkd->Self();
    struct smContext::stStack *S = smx->ST;
    int iCell,id;
    int sp = 0;
    int pj, pEnd, nCnt;

    nCnt = smx->nnListSize;

    auto kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = idSelf);

    while (1) {
        auto min2 = kdn->bound().mindist(r);
        if (min2 > fBall2) {
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            int idUpper,iUpper;
            std::tie(iCell,id,iUpper,idUpper) = kdn->get_child_cells(id);
            kdn = getCell(pkd,iCell,id);
            S[sp].id = idUpper;
            S[sp].iCell = iUpper;
            S[sp].min = 0.0;
            ++sp;
            continue;
        }
        else {
            if (id == pkd->Self()) {
                pEnd = kdn->upper();
                for (pj=kdn->lower(); pj<=pEnd; ++pj) {
                    auto p = pkd->particles[pj];
                    if (!p.is_gas()) continue;
                    TinyVector<double,3> dr{r - p.position()};
                    double fDist2 = dot(dr,dr);
                    if (fDist2 <= fBall2) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            assert(smx->nnList != NULL);
                        }
                        smx->nnList[nCnt].fDist2 = fDist2;
                        smx->nnList[nCnt].dr = dr;
                        smx->nnList[nCnt].pPart = &p;
                        smx->nnList[nCnt].iIndex = pj;
                        smx->nnList[nCnt].iPid = idSelf;
                        ++nCnt;
                    }
                }
            }
            else {
                pEnd = kdn->upper();
                for (pj=kdn->lower(); pj<=pEnd; ++pj) {
                    auto p = pkd->particles[static_cast<PARTICLE *>(mdlFetch(mdl,CID_PARTICLE,pj,id))];
                    if (!p.is_gas()) continue;
                    TinyVector<double,3> dr{r - p.position()};
                    double fDist2 = dot(dr,dr);
                    if (fDist2 <= fBall2) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            assert(smx->nnList != NULL);
                        }
                        smx->nnList[nCnt].fDist2 = fDist2;
                        smx->nnList[nCnt].dr = dr;
                        smx->nnList[nCnt].pPart = static_cast<PARTICLE *>(mdlAcquire(mdl,CID_PARTICLE,pj,id));
                        smx->nnList[nCnt].iIndex = pj;
                        smx->nnList[nCnt].iPid = id;
                        ++nCnt;
                    }
                }
            }
        }
NoIntersect:
        if (sp) {
            --sp;
            id = S[sp].id;
            iCell = S[sp].iCell;
            kdn = getCell(pkd,iCell,id);
        }
        else break;
    }
    smx->nnListSize = nCnt;
}

void smDoGatherLocal(SMX smx,double fBall2,TinyVector<double,3> r,void (*Do)(SMX,PARTICLE *,double)) {
    PKD pkd = smx->pkd;
    int *S = smx->S;
    int sp = 0;
    int iCell;

    auto kdn = pkd->tree[iCell = ROOT];
    while (1) {
        auto min2 = kdn->bound().mindist(r);
        if (min2 > fBall2) {
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            S[sp++] = kdn->rchild();
            kdn = pkd->tree[iCell = kdn->lchild()];
            continue;
        }
        else {
            for (auto &p : *kdn) {
                TinyVector<double,3> dr{r - p.position()};
                double fDist2 = dot(dr,dr);
                if (fDist2 <= fBall2) {
                    Do(smx,&p,fDist2);
                }
            }
        }
NoIntersect:
        if (sp) kdn = pkd->tree[iCell = S[--sp]];
        else break;
    }
}


void smReSmoothSingle(SMX smx,SMF *smf,particleStore::ParticleReference &p,double fBall) {
    PKD pkd = smx->pkd;
    TinyVector<double,3> r;

    auto R = p.position();

    smx->nnListSize = 0;
    /*
    ** Note for implementing SLIDING PATCH, the offsets for particles are
    ** negative here, reflecting the relative +ve offset of the simulation
    ** volume.
    */
    if (smx->bPeriodic) {
        TinyVector<int,3> iStart = floor((R - fBall) / pkd->fPeriod + 0.5);
        TinyVector<int,3> iEnd   = floor((R + fBall) / pkd->fPeriod + 0.5);
        for (int ix=iStart[0]; ix<=iEnd[0]; ++ix) {
            r[0] = R[0] - ix*pkd->fPeriod[0];
            for (int iy=iStart[1]; iy<=iEnd[1]; ++iy) {
                r[1] = R[1] - iy*pkd->fPeriod[1];
                for (int iz=iStart[2]; iz<=iEnd[2]; ++iz) {
                    r[2] = R[2] - iz*pkd->fPeriod[2];
                    smGather(smx,fBall*fBall,r);
                }
            }
        }
    }
    else {
        smGather(smx,fBall*fBall,R);
    }
    /*
    ** Apply smooth funtion to the neighbor list.
    */
    smx->fcnSmooth(&p,fBall,smx->nnListSize,smx->nnList,smf);
    /*
    ** Release acquired pointers.
    */
    for (int i=0; i<smx->nnListSize; ++i) {
        if (smx->nnList[i].iPid != pkd->Self()) {
            mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[i].pPart);
        }
    }
}


int smReSmooth(SMX smx,SMF *smf, int iSmoothType) {
    PKD pkd = smx->pkd;
    int nSmoothed=0;

    smf->pfDensity = NULL;
    switch (iSmoothType) {
    case SMX_HYDRO_DENSITY:
        for (auto &p : pkd->particles) {
            if (p.is_active() && p.marked() && p.is_gas()) {
                smReSmoothSingle(smx,smf,p,p.ball());
                nSmoothed++;
            }
        }
        break;

#ifdef BLACKHOLES
    case SMX_BH_DRIFT:
        for (auto &p : pkd->particles) {
            if (p.is_bh()) {
                smReSmoothSingle(smx,smf,p,p.ball());
                nSmoothed++;
            }
        }
        break;

    case SMX_BH_STEP:
        for (auto &p : pkd->particles) {
            if (p.is_bh()) {
                smReSmoothSingle(smx,smf,p,2.*p.soft());
                nSmoothed++;
            }
        }
        break;
#endif

#ifdef FEEDBACK
    /* IA: If computing the hydrostep, we also do the smooth over the newly formed stars that has not yet exploded, such that
     *  they can increase the rung of the neighbouring gas particles before exploding
     */
    case SMX_HYDRO_STEP:
        for (auto &p : pkd->particles) {
            if (p.is_gas()) {
                if (p.is_active()) {
                    smReSmoothSingle(smx,smf,p,p.ball());
                    nSmoothed++;
                }
            }
            //else if (p.is_star()) {
            //    if (!(p.star().bCCSNFBDone && p.star().bSNIaFBDone)) {
            //         IA: In principle this does NOT improve the Isolated Galaxy case, as we wait until the end of the
            //         step to update the primitive variables
            //        smReSmoothSingle(smx,smf,p,p.ball());
            //        nSmoothed++;
            //    }
            //}
        }
        break;

    case SMX_SN_FEEDBACK:
        for (auto &p : pkd->particles) {
            if (p.is_star()) {
                const auto &star = p.star();
                if (!star.bCCSNFBDone && ((smf->dTime - star.fTimer) > smf->dCCSNFBDelay)) {
                    smSmoothSingle(smx,smf,p,ROOT,0);
                    nSmoothed++;
                }
                if (!star.bSNIaFBDone && ((smf->dTime - star.fTimer) > smf->dSNIaFBDelay)) {
                    smSmoothSingle(smx,smf,p,ROOT,0);
                    nSmoothed++;
                }
            }
        }
        break;
#endif
#ifdef STELLAR_EVOLUTION
    case SMX_CHEM_ENRICHMENT:
        for (auto &p : pkd->particles) {
            if (p.is_star()) {
                const auto &star = p.star();
                if (smf->dTime > star.fNextEnrichTime) {
                    smReSmoothSingle(smx,smf,p,p.ball());
                    nSmoothed++;
                }
            }
        }
        break;
#endif
    default:
        for (auto &p : pkd->particles) {
            if (p.is_active() && p.is_gas()) {
                smReSmoothSingle(smx,smf,p,p.ball());
                nSmoothed++;
            }
        }
    }

    return nSmoothed;
}

#ifdef OPTIM_SMOOTH_NODE

/* Allocate a buffer for N particles each one with nVar variables.
 * Also allocate an array of pointers to the beggining of each variable array
 *
 * If oldBuff is present, the new buffer is initialized with the contents of
 * oldBuff, with a size oldN
 */
void static inline allocNodeBuffer(const int N, const int nVar, my_real **p_buff,
                                   my_real ***p_ptrs, my_real *oldBuff, const int oldN) {

    assert(oldN < N);

    *p_buff = new (std::align_val_t(64)) my_real[N*nVar];
    assert(*p_buff!=NULL);
    if (oldBuff == NULL) {
        *p_ptrs = new (std::align_val_t(64)) my_real*[nVar];
        assert(*p_ptrs!=NULL);
    }

    myreal *buff = *p_buff;
    myreal **ptrs = *p_ptrs;

    // Fill pointer array
    for (int i=0; i<nVar; i++)
        ptrs[i] = &buff[i*N];

    // Fill if requested
    if (oldBuff != NULL)
        for (int var=0; var<nVar; var++)
            for (int i=0; i<oldN; i++)
                buff[var*N + i] =  oldBuff[var*oldN + i];

}

void static inline reallocNodeBuffer(const int N, const int nVar, my_real **p_buff,
                                     my_real ***p_ptrs, const int oldN) {

    my_real *tmp_buffer;

    allocNodeBuffer(N, nVar, &tmp_buffer, p_ptrs, *p_buff, oldN);

    delete [] *p_buff;

    *p_buff = tmp_buffer;
}

/* IA: In this version, we loop over the buckets, rather than over the particles.
 *
 * For each bucket, we look for all the surroiding buckets that may interact
 * with any particle in said bucket.
 *
 * Then, we put all those particles (including the own bucket)
 * in an interaction list.
 */
int smReSmoothNode(SMX smx,SMF *smf, int iSmoothType) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    int nCnt;
    int nSmoothed=0;

    smx->nnListSize = 0;
    int nnListMax_p = NNLIST_INCREMENT;

    NN *nnList_p;
    nnList_p = static_cast<NN *>(malloc(sizeof(NN)*nnListMax_p));

    // Here we store the pointers to the particles whose interaction
    // need to be computed
    std::vector<PARTICLE *> sinks;
    sinks.reserve(64); // At most, the size of the bucket

    /* For allowing vectorization, it is better to use an structure of
     *  arrays rather than an array of structures.
     *
     *  In our case, the structure is just an array of pointers to the locations
     *  in the buffer where a given array of variables starts.
     *
     *  Having everything in the same buffer can be advantageous as it should
     *  be in cache, but it was painful to code...
     */
    my_real *input_buffer = NULL;
    my_real **input_pointers = NULL;
    my_real *output_buffer = NULL;
    my_real **output_pointers = NULL;
    int inNvar, outNvar;
    if (smx->fcnSmoothGetNvars) {
        smx->fcnSmoothGetNvars(&inNvar, &outNvar);
        allocNodeBuffer(nnListMax_p, inNvar, &input_buffer, &input_pointers, NULL,0);
        allocNodeBuffer(nnListMax_p, outNvar, &output_buffer, &output_pointers, NULL,0);
    }


    double fBall_factor = 1.;
    if (iSmoothType==SMX_HYDRO_DENSITY)
        fBall_factor *= 1.2; // An small margin is kept in case fBall needs to increase


    for (auto i = NRESERVED_NODES; i < pkd->Nodes() - 1; ++i) {
        auto node = pkd->tree[i];
        if (node->is_bucket()) {

            // Prepare the interaction list

            auto bnd_node = node->bound();

            TinyVector<double,3> r{bnd_node.center()};

            // First, we add all the particles whose interactions need to be computed
            sinks.clear();
#ifdef OPTIM_REORDER_IN_NODES
            auto pEnd = node->lower() + node->Ngas();
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
            //if (iSmoothType==SMX_HYDRO_DENSITY) pEnd += node->Nstar();
#endif
#ifdef BLACKHOLES
            //if (iSmoothType==SMX_HYDRO_DENSITY) pEnd += node->Nbh();
#endif
#else // OPTIM_REORDER_IN_NODES
            auto pEnd = node->upper() + 1;
#endif

            TinyVector<double,3> fMax_shrink{0.};
            for (auto pj = node->lower(); pj < pEnd; ++pj) {
                auto p = pkd->particles[pj];

#ifdef OPTIM_AVOID_IS_ACTIVE
                int pIsActive = p.marked();
#else
                int pIsActive = p.is_active();
#endif


                if (pIsActive) {
                    if ( (iSmoothType==SMX_HYDRO_DENSITY) ||
                            (iSmoothType==SMX_HYDRO_DENSITY_FINAL) ) {

#ifndef OPTIM_AVOID_IS_ACTIVE
                        if (!p.marked())
                            continue;
#endif

#if defined(FEEDBACK) && !defined(STELLAR_EVOLUTION)
                        // If there is only feedback, once the particle explodes
                        // there is no need to updated its fBall.
                        // However, if we have stellar evolution, fBall needs to be
                        // always updated.
                        if (p.is_star() && p.star().bCCSNFBDone && p.star().bSNIaFBDone)
                            continue;
#endif

#ifndef OPTIM_REORDER_IN_NODES
                        // Explicit check of the types
                        if (!p.is_gas() && !p.is_star())
                            continue;
#endif
                    }
                    else {
#ifndef OPTIM_REORDER_IN_NODES
                        // Explicit check of the types
                        if (!p.is_gas()) continue;
#endif
                    } //SMX_HYDRO_DENSITY

                    fMax_shrink = blitz::max(fMax_shrink,abs(p.position() - bnd_node.center()) + p.ball()*fBall_factor);

                    sinks.push_back(&p);
                }
            }
            // There are no elligible particles in this bucket, go to the next
            if (sinks.empty()) continue;

            nCnt = 0;

            bnd_node.shrink(fMax_shrink);

            if (smx->bPeriodic) {
                TinyVector<int,3> iStart{floor((r - bnd_node.apothem()) / pkd->fPeriod + 0.5)};
                TinyVector<int,3> iEnd{floor((r + bnd_node.apothem()) / pkd->fPeriod + 0.5)};
                for (int ix=iStart[0]; ix<=iEnd[0]; ++ix) {
                    r[0] = bnd_node.center(0) - ix*pkd->fPeriod[0];
                    for (int iy=iStart[1]; iy<=iEnd[1]; ++iy) {
                        r[1] = bnd_node.center(1) - iy*pkd->fPeriod[1];
                        for (int iz=iStart[2]; iz<=iEnd[2]; ++iz) {
                            r[2] = bnd_node.center(2) - iz*pkd->fPeriod[2];
                            buildInteractionList(smx, smf, node, bnd_node, &nCnt, r, ix, iy, iz);
                        }
                    }
                }
            }
            else {
                buildInteractionList(smx, smf, node, bnd_node, &nCnt, r, 0, 0, 0);
            }


            // IA: Now we should have inside nnList all the particles in the
            //  bucket (sinks) and those of which can interact with them
            //  from other buckets (smx->nnList)
            //
            // We just have to proceed to compute the correct dx, dy, dz and
            // pass that nnList to the smoothfcn routine
            //
            // However, we have different options to do so:
            // 1) Naive: we pass the whole nnList
            //
            // 2) Sorting: we could follow Gonnet, 2007 (10.1002/jcc.20563) to
            // reduce the number of distance computations, but this is troublesome
            // in our case because:
            //      a) we need to compute the distance anyway for sorting
            //      b) we could sort relative to the cell, but this is suboptimal
            //      c) we are not computing cell-cell interactions, so there is
            //            no well-defined axis that could be used for projection
            //



            // For the smoothing length determination we can bypass the typical
            //  flow of calling fcnsmooth, as probably we have gathered more
            //  neighbours than needed and thus the iterative procedure should be
            //  faster
            if (iSmoothType==SMX_HYDRO_DENSITY) {
                hydroDensity_node(pkd, smf, bnd_node, sinks, smx->nnList, nCnt);
            }
            else {
                for (auto &P : sinks) {
                    auto partj = pkd->particles[P];
                    const float pBall2 = partj.ball()*partj.ball();
                    const TinyVector<double,3> dr_node{bnd_node.center() - partj.position()};

                    int nCnt_p = 0;
                    for (auto pk = 0; pk < nCnt; ++pk) {
                        if (P == smx->nnList[pk].pPart) continue;
                        const TinyVector<double,3> dr{smx->nnList[pk].dr - dr_node};
                        const double fDist2 = dot(dr,dr);
                        if (fDist2 < pBall2) {

                            // Reasons not to compute this interaction
                            if ( iSmoothType==SMX_HYDRO_FLUX ||
                                    iSmoothType==SMX_HYDRO_FLUX_VEC) {

                                const auto &qBall = smx->nnList[pk].fBall;
                                if (qBall*qBall < fDist2)
                                    continue;

#ifdef OPTIM_AVOID_IS_ACTIVE
                                int qIsActive = smx->nnList[pk].bMarked;
#else
                                abort(); // not supported
#endif

#ifdef OPTIM_NO_REDUNDANT_FLUXES
                                if (qIsActive) {
                                    if (dr[0] > 0) continue;
                                    else if (dr[0]==0) {
                                        if (dr[1] > 0) continue;
                                        else if (dr[1]==0) {
                                            if (dr[2] > 0) continue;
                                            else if (dr[2]==0) abort();
                                        }
                                    }
                                }
#endif
                            }

                            // Try pointer to pPart declared as restrict, to check if compiler does something better

                            if (nCnt_p >= nnListMax_p) {
                                nnListMax_p += NNLIST_INCREMENT;
                                nnList_p = static_cast<NN *>(realloc(nnList_p,nnListMax_p*sizeof(NN)));
                                assert(nnList_p != NULL);
                                if (smx->fcnSmoothGetNvars) {
                                    printf("WARNING: Increasing smoothNode buffer size to %d\n",
                                           nnListMax_p);
                                    int oldListMax = nnListMax_p - NNLIST_INCREMENT;
                                    reallocNodeBuffer(nnListMax_p, inNvar,
                                                      &input_buffer, &input_pointers, oldListMax);
                                    reallocNodeBuffer(nnListMax_p, outNvar,
                                                      &output_buffer, &output_pointers, oldListMax);
                                }
                            }

                            nnList_p[nCnt_p].pPart = (smx->nnList[pk].iPid == pkd->Self()) ?
                                                     smx->nnList[pk].pPart : static_cast<PARTICLE *>(mdlAcquire(mdl,CID_PARTICLE,smx->nnList[pk].iIndex,smx->nnList[pk].iPid));
                            nnList_p[nCnt_p].fDist2 = fDist2;
                            nnList_p[nCnt_p].dr = dr;
                            nnList_p[nCnt_p].iIndex = smx->nnList[pk].iIndex;
                            nnList_p[nCnt_p].iPid = smx->nnList[pk].iPid;

                            if (smx->fcnSmoothFillBuffer) {
                                PARTICLE *q = nnList_p[nCnt_p].pPart;

                                smx->fcnSmoothFillBuffer(input_pointers, q, nCnt_p,
                                                         fDist2, dr, smf);
                            }

                            ++nCnt_p;
                        }

                    }

                    //abort();
                    //printf("nCnt_p %d \n", nCnt_p);
                    //assert(nCnt_p<200);

                    if (smx->fcnSmoothNode) {
                        smx->fcnSmoothNode(&partj,partj.ball(),nCnt_p,
                                           input_pointers, output_pointers, smf);
                        for (auto pk = 0; pk < nCnt_p; ++pk) {
                            smx->fcnSmoothUpdate(output_pointers,input_pointers,
                                                 &partj, nnList_p[pk].pPart, pk, smf);
                        }
                    }
                    else {
                        smx->fcnSmooth(&partj,partj.ball(),nCnt_p,nnList_p,smf);
                    }

                    for (auto pk = 0; pk < nCnt_p; ++pk) {
                        if (nnList_p[pk].iPid != pkd->Self()) {
                            mdlRelease(pkd->mdl,CID_PARTICLE,nnList_p[pk].pPart);
                        }
                    }
                }
            }

            nSmoothed += sinks.size();

        }
    }
    if (smx->fcnSmoothGetNvars) {
        delete [] input_buffer;
        delete [] input_pointers;
        delete [] output_buffer;
        delete [] output_pointers;
    }
    free(nnList_p);
    //printf("nSmoothed %d \n", nSmoothed);
    return nSmoothed;
}



void buildInteractionList(SMX smx, SMF *smf, KDN *node, Bound bnd_node, int *nCnt_tot, TinyVector<double,3> r, int ix, int iy, int iz) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    int id, sp, iCell;
    struct smContext::stStack *S = smx->ST;
    int nCnt = *nCnt_tot;

    // We look for the biggest node that encloses the needed domain
    // We can only take advantage of this if we are are in the original cell
    auto kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = pkd->Self());

    //  Now we start the walk as usual
    sp = 0;
    while (1) {
        auto bnd = kdn->bound();
        if (any(abs(bnd.center()-r) - bnd.apothem() - bnd_node.apothem() > 0.)) goto NoIntersect;

        /*
        ** We have an intersection to test.
        */
        if (kdn->is_cell()) {
            int idUpper,iUpper;
            std::tie(iCell,id,iUpper,idUpper) = kdn->get_child_cells(id);
            kdn = getCell(pkd,iCell,id);
            S[sp].id = idUpper;
            S[sp].iCell = iUpper;
            S[sp].min = 0.0;
            ++sp;
            continue;
        }
        else {
            if (id == pkd->Self()) {
#ifdef OPTIM_REORDER_IN_NODES
                auto pEnd = kdn->lower() + kdn->Ngas();
#else
                auto pEnd = kdn->upper() + 1;
#endif
                //printf("pEnd %d \n", pEnd);
                for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                    auto p = pkd->particles[pj];
#ifndef OPTIM_REORDER_IN_NODES
                    if (!p.is_gas()) continue;
#endif
                    const TinyVector<double,3> dr {r - p.position()};
                    if (all(abs(dr) <= bnd_node.apothem())) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            //printf("realloc \n");
                            assert(smx->nnList != NULL);
                        }
                        smx->nnList[nCnt].fDist2 = dot(dr,dr);
                        smx->nnList[nCnt].dr = dr;
                        smx->nnList[nCnt].pPart = &p;
                        smx->nnList[nCnt].fBall = p.ball();
                        smx->nnList[nCnt].bMarked = p.marked();
                        smx->nnList[nCnt].iIndex = pj;
                        smx->nnList[nCnt].iPid = pkd->Self();
                        ++nCnt;
                    }
                }
            }
            else {
#ifdef OPTIM_REORDER_IN_NODES
                auto pEnd = kdn->lower() + kdn->Ngas();
#else
                auto pEnd = kdn->upper() + 1;
#endif
                for (auto pj = kdn->lower(); pj < pEnd; ++pj) {
                    auto p = pkd->particles[static_cast<PARTICLE *>(mdlFetch(mdl,CID_PARTICLE,pj,id))];
#ifndef OPTIM_REORDER_IN_NODES
                    if (!p.is_gas()) continue;
#endif
                    const TinyVector<double,3> dr {r - p.position()};
                    if (all(abs(dr) <= bnd_node.apothem())) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            //printf("realloc \n");
                            assert(smx->nnList != NULL);
                        }

                        smx->nnList[nCnt].fDist2 = dot(dr,dr);
                        smx->nnList[nCnt].dr = dr;
                        smx->nnList[nCnt].fBall = p.ball();
                        smx->nnList[nCnt].bMarked = p.marked();
                        // The actual particle information fill be acquired later
                        smx->nnList[nCnt].pPart = NULL;
                        smx->nnList[nCnt].iIndex = pj;
                        smx->nnList[nCnt].iPid = id;
                        ++nCnt;
                    }
                }
            }
        }
NoIntersect:
        if (sp) {
            --sp;
            id = S[sp].id;
            iCell = S[sp].iCell;
            kdn = getCell(pkd,iCell,id);
        }
        else break;
    }

    *nCnt_tot = nCnt;

}
#endif //OPTIM_SMOOTH_NODE
