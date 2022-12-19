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
#endif
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#include <sys/stat.h>

#define ASSERT_CONCAT_(a, b) a##b
#define ASSERT_CONCAT(a, b) ASSERT_CONCAT_(a, b)
/* These can't be used after statements in c89. */
#ifdef __COUNTER__
#define STATIC_ASSERT(e,m)                      \
    { enum { ASSERT_CONCAT(static_assert_, __COUNTER__) = 1/(!!(e)) }; }
#else
/* This can't be used twice on the same line so ensure if using in headers
 * that the headers are not included twice (by wrapping in #ifndef...#endif)
 * Note it doesn't cause an issue when used on same line of separate modules
 * compiled with gcc -combine -fwhole-program.  */
#define STATIC_ASSERT(e,m) \
    { enum { ASSERT_CONCAT(assert_line_, __LINE__) = 1/(!!(e)) }; }
#endif


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
    void (*init)(void *,void *) = NULL;
    void (*comb)(void *,void *,const void *) = NULL;
    int i,pi,j;
    int nTree;

    smx = new struct smContext;
    assert(smx != NULL);
    smx->pSentinel = static_cast<PARTICLE *>(malloc(pkd->ParticleSize()));
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
        smx->fcnPost = NULL;
        break;
    case SMX_DENSITY:
        smx->fcnSmooth = bSymmetric?DensitySym:Density;
        initParticle = initDensity; /* Original Particle */
        init = initDensity; /* Cached copies */
        comb = combDensity;
        smx->fcnPost = NULL;
        break;
    case SMX_DENSITY_F1:
        assert(!bSymmetric);
        smx->fcnSmooth = DensityF1;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_DENSITY_M3:
        assert(!bSymmetric);
        assert(pkd->particles.present(PKD_FIELD::oGroup));
        smx->fcnSmooth = DensityM3;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_GRADIENT_M3:
        assert(!bSymmetric);
        smx->fcnSmooth = LinkGradientM3;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_HOP_LINK:
        assert(!bSymmetric);
        smx->fcnSmooth = LinkHopChains;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_PRINTNN:
        smx->fcnSmooth = PrintNN;
        initParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_HYDRO_DENSITY:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroDensity;
        initParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_HYDRO_GRADIENT:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroGradients;
        initParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_HYDRO_FLUX:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroRiemann;
        initParticle = initHydroFluxes; /* Original Particle */
        init = initHydroFluxesCached; /* Cached copies */
        comb = combThirdHydroLoop;
        smx->fcnPost = NULL;
        break;
    case SMX_HYDRO_FLUX_VEC:
        assert (pkd->particles.present(PKD_FIELD::oSph));
        smx->fcnSmoothNode = hydroRiemann_vec;
        smx->fcnSmoothGetNvars = hydroFluxGetNvars;
        smx->fcnSmoothFillBuffer = hydroFluxFillBuffer;
        smx->fcnSmoothUpdate = hydroFluxUpdateFromBuffer;
        initParticle = initHydroFluxes; /* Original Particle */
        init = initHydroFluxesCached; /* Cached copies */
        comb = combThirdHydroLoop;
        smx->fcnPost = NULL;
        break;
    case SMX_HYDRO_STEP:
        assert( pkd->particles.present(PKD_FIELD::oSph) ); /* Validate memory model */
        smx->fcnSmooth = hydroStep;
        initParticle = initHydroStep; /* Original Particle */
        init = initHydroStep; /* Cached copies */
        comb = combHydroStep;
        smx->fcnPost = NULL;
        break;
#ifdef FEEDBACK
    case SMX_SN_FEEDBACK:
        smx->fcnSmooth = smSNFeedback;
        initParticle = NULL; /* Original Particle */
        init = initSNFeedback; /* Cached copies */
        comb = combSNFeedback;
        smx->fcnPost = NULL;
        break;
#endif
#ifdef BLACKHOLES
    case SMX_BH_MERGER:
        smx->fcnSmooth = smBHmerger;
        initParticle = NULL; /* Original Particle */
        init = NULL; /* Cached copies */
        comb = combBHmerger;
        smx->fcnPost = NULL;
        break;
    case SMX_BH_DRIFT:
        smx->fcnSmooth = smBHevolve;
        initParticle = NULL; /* Original Particle */
        init = initBHevolve; /* Cached copies */
        comb = combBHevolve;
        smx->fcnPost = NULL;
        break;
#endif
#ifdef STELLAR_EVOLUTION
    case SMX_CHEM_ENRICHMENT:
        smx->fcnSmooth = smChemEnrich;
        initParticle = NULL;
        init = initChemEnrich;
        comb = combChemEnrich;
        smx->fcnPost = NULL;
        break;
#endif
#ifndef OPTIM_REMOVE_UNUSED
    case SMX_DENDVDX:
        assert( pkd->oFieldOffset[oSph] ); /* Validate memory model */
        smx->fcnSmooth = DenDVDX;
        initParticle = NULL; /* Original Particle */
        init = initDenDVDX; /* Cached copies */
        comb = combDenDVDX;
        smx->fcnPost = NULL;
        break;
    case SMX_SPHFORCES:
        assert( pkd->oFieldOffset[oSph] ); /* Validate memory model */
        assert( pkd->oFieldOffset[oAcceleration] ); /* Validate memory model */
        smx->fcnSmooth = SphForces;
        initParticle = initSphForcesParticle; /* Original Particle */
        init = initSphForces; /* Cached copies */
        comb = combSphForces;
        smx->fcnPost = NULL;
        break;
    case SMX_DIST_DELETED_GAS:
        assert(bSymmetric != 0);
        smx->fcnSmooth = DistDeletedGas;
        initParticle = NULL;
        init = initDistDeletedGas;
        comb = combDistDeletedGas;
        smx->fcnPost = NULL;
        break;
    case SMX_DIST_SN_ENERGY:
        assert(bSymmetric != 0);
        smx->fcnSmooth = DistSNEnergy;
        initParticle = NULL;
        init = initDistSNEnergy;
        comb = combDistSNEnergy;
        smx->fcnPost = NULL;
        break;
    case SMX_MEANVEL:
        assert( pkd->oFieldOffset[oVelSmooth]); /* Validate memory model */
        smx->fcnSmooth = bSymmetric?MeanVelSym:MeanVel;
        initParticle = initMeanVel; /* Original Particle */
        init = initMeanVel; /* Cached copies */
        comb = combMeanVel;
        smx->fcnPost = NULL;
        break;
    case SMX_DIVV:
        assert( pkd->oFieldOffset[oVelSmooth]); /* Validate memory model */
        smx->fcnSmooth = bSymmetric?DivvSym:Divv;
        initParticle = initDivv; /* Original Particle */
        init = initDivv; /* Cached copies */
        comb = combDivv;
        smx->fcnPost = NULL;
        break;
    case SMX_VELDISP2:
        assert( pkd->oFieldOffset[oVelSmooth]); /* Validate memory model */
        smx->fcnSmooth = bSymmetric?VelDisp2Sym:VelDisp2;
        initParticle = initVelDisp2; /* Original Particle */
        init = initVelDisp2; /* Cached copies */
        comb = combVelDisp2;
        smx->fcnPost = NULL;
        break;
    case SMX_FOF:
        assert(bSymmetric == 0);
        smx->fcnSmooth = NULL;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
    case SMX_RELAXATION:
        assert( pkd->oFieldOffset[oRelaxation]); /* Validate memory model */
        assert(bSymmetric == 0);
        smx->fcnSmooth = AddRelaxation;
        initParticle = NULL;
        init = NULL;
        comb = NULL;
        smx->fcnPost = NULL;
        break;
#endif //OPTIM_REMOVE_UNUSED
    case SMX_BALL:
        assert(pkd->particles.present(PKD_FIELD::oBall));
        smx->fcnSmooth = BallSmooth;
        initParticle = initBall; /* Original Particle */
        init = initBall; /* Cached copies */
        comb = NULL;
        smx->fcnPost = NULL;
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
    nTree = pkd->TreeNode(ROOT)->pUpper + 1;
    if (initParticle != NULL) {
        for (pi=0; pi<nTree; ++pi) {
            PARTICLE *p = pkd->Particle(pi);
            /*if (TYPETest(p,smx->eParticleTypes))*/
            if (pkdIsActive(pkd,p)) initParticle(pkd,p);
        }
    }
    /*
    ** Start particle caching space (cell cache is already active).
    */
    if (bMakeCache) {
        smx->bOwnCache = 1;
        if (bSymmetric) {
            mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,
                       pkd->ParticleBase(),pkd->ParticleSize(),
                       nTree,pkd,init,comb);
        }
        else {
            mdlROcache(pkd->mdl,CID_PARTICLE,NULL,
                       pkd->ParticleBase(),pkd->ParticleSize(),
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
    for (j=0; j<3; ++j) {
        if (pkd->bIntegerPosition) pkdSetPosRaw(pkd,smx->pSentinel,j,INT32_MAX);
        else pkdSetPos(pkd,smx->pSentinel,j,HUGE_VAL);
    }
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
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;
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
    ** Now do any post calculations, these ususlly involve some sort of
    ** normalizations of the smoothed quantities, usually division by
    ** the local density! Do NOT put kernel normalizations in here as
    ** these do not depend purely on local properties in the case of
    ** "Gather-Scatter" kernel.
    */
    if (smx->fcnPost != NULL) {
        for (pi=0; pi<pkd->Local(); ++pi) {
            p = pkd->Particle(pi);
            smx->fcnPost(pkd,p,smf);
        }
    }
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

static KDN *getCell(PKD pkd, int iCell, int id) {
    if (id==pkd->Self()) return pkd->TreeNode(iCell);
    return static_cast<KDN *>(mdlFetch(pkd->mdl,CID_CELL,iCell,id));
}

PQ *pqSearch(SMX smx,PQ *pq,double r[3],int iRoot) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    KDN *kdn;
    int idSelf = smx->pkd->Self();
    struct smContext::stStack *S = smx->ST;
    double min1,min2;
    double p_r[3];
    int iCell,id;
    int sp = 0;
    BND bnd;
    int pEnd, pj;
    PARTICLE *p;
    double dx,dy,dz,fDist2;

    /* Start at the root node of the tree */
    kdn = getCell(pkd,pkd->iTopTree[iRoot],id = idSelf);
    while (1) {
        while (kdn->iLower) {
            int idLower,iLower,idUpper,iUpper;
            pkdGetChildCells(kdn,id,idLower,iLower,idUpper,iUpper);
            kdn = getCell(pkd,iLower,idLower);
            bnd = pkdNodeGetBnd(pkd, kdn);
            MINDIST(&bnd,r,min1);
            kdn = getCell(pkd,iUpper,idUpper);
            bnd = pkdNodeGetBnd(pkd, kdn);
            MINDIST(&bnd,r,min2);
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
            pEnd = kdn->pUpper;
            for (pj=kdn->pLower; pj<=pEnd; ++pj) {
                p = pkd->Particle(pj);
                if (!p->bMarked) continue;
                pkdGetPos1(pkd,p,p_r);
                dx = r[0] - p_r[0];
                dy = r[1] - p_r[1];
                dz = r[2] - p_r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= pq->fDist2) {
                    if (pq->iPid == idSelf) {
                        pkd->Particle(pq->iIndex)->bMarked = 1;
                    }
                    else {
                        smHashDel(smx,pq->pPart);
                        mdlRelease(mdl,CID_PARTICLE,pq->pPart);
                        pq->iPid = idSelf;
                    }
                    pq->pPart = p;
                    pq->fDist2 = fDist2;
                    pq->dx = dx;
                    pq->dy = dy;
                    pq->dz = dz;
                    pq->iIndex = pj;
                    p->bMarked = 0; /* de-activate a particle that enters the queue */
                    PQ_REPLACE(pq);
                }
            }
        }
        else {
            pEnd = kdn->pUpper;
            for (pj=kdn->pLower; pj<=pEnd; ++pj) {
                p = static_cast<PARTICLE *>(mdlFetch(mdl,CID_PARTICLE,pj,id));
                if (smHashPresent(smx,p)) continue;
                pkdGetPos1(pkd,p,p_r);
                dx = r[0] - p_r[0];
                dy = r[1] - p_r[1];
                dz = r[2] - p_r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= pq->fDist2) {
                    if (pq->iPid == idSelf) {
                        pkd->Particle(pq->iIndex)->bMarked = 1;
                    }
                    else {
                        smHashDel(smx,pq->pPart);
                        mdlRelease(mdl,CID_PARTICLE,pq->pPart);
                    }
                    p = static_cast<PARTICLE *>(mdlAcquire(mdl,CID_PARTICLE,pj,id));
                    pq->pPart = p;
                    pq->fDist2 = fDist2;
                    pq->dx = dx;
                    pq->dy = dy;
                    pq->dz = dz;
                    pq->iIndex = pj;
                    pq->iPid = id;
                    smHashAdd(smx,p);
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
    for (i=0; i<smx->nSmooth; ++i) {
        smx->pq[i].pPart = smx->pSentinel;
        smx->pq[i].iIndex = smx->pkd->Local();
        smx->pq[i].iPid = smx->pkd->Self();
        smx->pq[i].dx = pkdPos(smx->pkd,smx->pSentinel,0);
        smx->pq[i].dy = pkdPos(smx->pkd,smx->pSentinel,1);
        smx->pq[i].dz = pkdPos(smx->pkd,smx->pSentinel,2);
        smx->pq[i].fDist2 = pow(smx->pq[i].dx,2) + pow(smx->pq[i].dy,2) +
                            pow(smx->pq[i].dz,2);
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
            smx->pkd->Particle(smx->pq[i].iIndex)->bMarked = 1;
        }
        else {
            smHashDel(smx,smx->pq[i].pPart);
            mdlRelease(smx->pkd->mdl,CID_PARTICLE,smx->pq[i].pPart);
        }
    }
}

float smSmoothSingle(SMX smx,SMF *smf,PARTICLE *p,int iRoot1, int iRoot2) {
    PKD pkd = smx->pkd;
    int ix,iy,iz;
    double p_r[3];
    double r[3],fBall;
    int iStart[3],iEnd[3];
    int j;
    PQ *pq;

    pkdGetPos1(pkd,p,p_r);

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
        smx->pq[j].dx += p_r[0]-smx->rLast[0];
        smx->pq[j].dy += p_r[1]-smx->rLast[1];
        smx->pq[j].dz += p_r[2]-smx->rLast[2];
        smx->pq[j].fDist2 = pow(smx->pq[j].dx,2) + pow(smx->pq[j].dy,2) +
                            pow(smx->pq[j].dz,2);
    }
    for (j=0; j<3; ++j) smx->rLast[j] = r[j] = p_r[j];

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
        for (j=0; j<3; ++j) {
            iStart[j] = d2i(floor((p_r[j] - fBall)/pkd->fPeriod[j] + 0.5));
            iEnd[j] = d2i(floor((p_r[j] + fBall)/pkd->fPeriod[j] + 0.5));
        }
        for (ix=iStart[0]; ix<=iEnd[0]; ++ix) {
            r[0] = p_r[0] - ix*pkd->fPeriod[0];
            for (iy=iStart[1]; iy<=iEnd[1]; ++iy) {
                r[1] = p_r[1] - iy*pkd->fPeriod[1];
                for (iz=iStart[2]; iz<=iEnd[2]; ++iz) {
                    r[2] = p_r[2] - iz*pkd->fPeriod[2];
//    printf("\t ix %d iy %d iz %d \n", ix, iy, iz);
//    printf("\t x %f y %f z %f \n", r[0], r[1], r[2]);
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
    smx->fcnSmooth(p,fBall,smx->nSmooth,smx->pq,smf);
    return fBall;
}

void smSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;
    float fBall;

    /*
    ** Initialize the bInactive flags for all local particles.
    */
    for (pi=0; pi<pkd->Local(); ++pi) {
        p = pkd->Particle(pi);
        p->bMarked = 1;
    }
    smSmoothInitialize(smx);
    smf->pfDensity = NULL;
    for (pi=0; pi<pkd->Local(); ++pi) {
        p = pkd->Particle(pi);
        if (!smf->bMeshlessHydro ) {
            smSmoothSingle(smx,smf,p,ROOT,0);
            //pkdSetBall(pkd,p,smSmoothSingle(smx,smf,p,ROOT,0));
        }
        else {
            if (pkdIsActive(pkd,p)) {
                fBall = smSmoothSingle(smx,smf,p,ROOT,0);
                if (smf->bUpdateBall) {
                    pkdSetBall(pkd,p,fBall);
                }
                /*
                    smSmoothFinish(smx);
                    for (int pj=0;pj<pkd->Local();++pj) {
                    PARTICLE *p2 = pkd->Particle(pj);
                    p2->bMarked = 1;
                    }
                    smSmoothInitialize(smx);
                    smf->pfDensity = NULL;
                */
            }
        }
        /*
        ** Call mdlCacheCheck to make sure we are making progress!
        */
        mdlCacheCheck(pkd->mdl);
    }
    smSmoothFinish(smx);
}

void smGather(SMX smx,double fBall2,double r[3], PARTICLE *pp) {
    PARTICLE *p;
    KDN *kdn;
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    int idSelf = pkd->Self();
    struct smContext::stStack *S = smx->ST;
    double min2;
    int iCell,id;
    int sp = 0;
    BND bnd;
    double p_r[3];
    double dx, dy, dz, fDist2;
    int pj, pEnd, nCnt;

    nCnt = smx->nnListSize;

    kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = idSelf);

    while (1) {
        bnd = pkdNodeGetBnd(pkd, kdn);
        MINDIST(&bnd,r,min2);
        if (min2 > fBall2) {
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->iLower) {
            int idUpper,iUpper;
            pkdGetChildCells(kdn,id,id,iCell,idUpper,iUpper);
            kdn = getCell(pkd,iCell,id);
            S[sp].id = idUpper;
            S[sp].iCell = iUpper;
            S[sp].min = 0.0;
            ++sp;
            continue;
        }
        else {
            if (id == pkd->Self()) {
                pEnd = kdn->pUpper;
                for (pj=kdn->pLower; pj<=pEnd; ++pj) {
                    p = pkd->Particle(pj);
                    if (!pkdIsGas(pkd,p)) continue;
                    pkdGetPos1(pkd,p,p_r);
                    dx = r[0] - p_r[0];
                    dy = r[1] - p_r[1];
                    dz = r[2] - p_r[2];
                    fDist2 = dx*dx + dy*dy + dz*dz;
                    if (fDist2 <= fBall2) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            assert(smx->nnList != NULL);
                        }
                        smx->nnList[nCnt].fDist2 = fDist2;
                        smx->nnList[nCnt].dx = dx;
                        smx->nnList[nCnt].dy = dy;
                        smx->nnList[nCnt].dz = dz;
                        smx->nnList[nCnt].pPart = p;
                        smx->nnList[nCnt].iIndex = pj;
                        smx->nnList[nCnt].iPid = idSelf;
                        ++nCnt;
                    }
                }
            }
            else {
                pEnd = kdn->pUpper;
                for (pj=kdn->pLower; pj<=pEnd; ++pj) {
                    p = static_cast<PARTICLE *>(mdlFetch(mdl,CID_PARTICLE,pj,id));
                    if (!pkdIsGas(pkd,p)) continue;
                    pkdGetPos1(pkd,p,p_r);
                    dx = r[0] - p_r[0];
                    dy = r[1] - p_r[1];
                    dz = r[2] - p_r[2];
                    fDist2 = dx*dx + dy*dy + dz*dz;
                    if (fDist2 <= fBall2) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            assert(smx->nnList != NULL);
                        }
                        smx->nnList[nCnt].fDist2 = fDist2;
                        smx->nnList[nCnt].dx = dx;
                        smx->nnList[nCnt].dy = dy;
                        smx->nnList[nCnt].dz = dz;
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

void smDoGatherLocal(SMX smx,double fBall2,double r[3],void (*Do)(SMX,PARTICLE *,double)) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    double p_r[3];
    double min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,pEnd;
    BND bnd;

    kdn = pkd->TreeNode(iCell = ROOT);
    while (1) {
        bnd = pkdNodeGetBnd(pkd, kdn);
        MINDIST(&bnd,r,min2);
        if (min2 > fBall2) {
            goto NoIntersect;
        }
        /*
        ** We have an intersection to test.
        */
        if (kdn->iLower) {
            kdn = pkd->TreeNode(iCell = kdn->iLower);
            S[sp++] = iCell+1;
            continue;
        }
        else {
            pEnd = kdn->pUpper;
            for (pj=kdn->pLower; pj<=pEnd; ++pj) {
                p = pkd->Particle(pj);
                pkdGetPos1(pkd,p,p_r);
                dx = r[0] - p_r[0];
                dy = r[1] - p_r[1];
                dz = r[2] - p_r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2) {
                    Do(smx,p,fDist2);
                }
            }
        }
NoIntersect:
        if (sp) kdn = pkd->TreeNode(iCell = S[--sp]);
        else break;
    }
}


void smReSmoothSingle(SMX smx,SMF *smf,PARTICLE *p,double fBall) {
    PKD pkd = smx->pkd;
    double R[3], r[3];
    int iStart[3],iEnd[3];
    int i,j;
    int ix,iy,iz;

    pkdGetPos1(pkd,p,R);

    smx->nnListSize = 0;
    /*
    ** Note for implementing SLIDING PATCH, the offsets for particles are
    ** negative here, reflecting the relative +ve offset of the simulation
    ** volume.
    */
    if (smx->bPeriodic) {
        for (j=0; j<3; ++j) {
            iStart[j] = d2i(floor((R[j] - fBall)/pkd->fPeriod[j] + 0.5));
            iEnd[j] = d2i(floor((R[j] + fBall)/pkd->fPeriod[j] + 0.5));
        }
        for (ix=iStart[0]; ix<=iEnd[0]; ++ix) {
            r[0] = R[0] - ix*pkd->fPeriod[0];
            for (iy=iStart[1]; iy<=iEnd[1]; ++iy) {
                r[1] = R[1] - iy*pkd->fPeriod[1];
                for (iz=iStart[2]; iz<=iEnd[2]; ++iz) {
                    r[2] = R[2] - iz*pkd->fPeriod[2];
                    smGather(smx,fBall*fBall,r, p);
                }
            }
        }
    }
    else {
        smGather(smx,fBall*fBall,R, p);
    }
    /*
    ** Apply smooth funtion to the neighbor list.
    */
    smx->fcnSmooth(p,0.5*fBall,smx->nnListSize,smx->nnList,smf);
    /*
    ** Release acquired pointers.
    */
    for (i=0; i<smx->nnListSize; ++i) {
        if (smx->nnList[i].iPid != pkd->Self()) {
            mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[i].pPart);
        }
    }
}


int  smReSmooth(SMX smx,SMF *smf, int iSmoothType) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi, nSmoothed=0;

    smf->pfDensity = NULL;
    switch (iSmoothType) {
    case SMX_HYDRO_DENSITY:
        for (pi=0; pi<pkd->Local(); ++pi) {
            p = pkd->Particle(pi);
#ifdef FEEDBACK
            // We follow the density of stars that has not yet exploded to have a proper fBall
            if (pkdIsActive(pkd,p) && p->bMarked && (pkdIsGas(pkd,p) || pkdIsStar(pkd,p))) {
                if (pkdIsStar(pkd,p) && (pkdStar(pkd,p)->hasExploded==1)) continue;
#else
            if (pkdIsActive(pkd,p) && p->bMarked && pkdIsGas(pkd,p)) {
#endif


                //if (pkdIsStar(pkd,p)) printf("%d \n",pkdStar(pkd,p)->hasExploded);

                smReSmoothSingle(smx,smf,p,2.*pkdBall(pkd,p));
                nSmoothed++;
            }
        }
        break;
    case SMX_BH_DRIFT:
        for (pi=0; pi<pkd->Local(); ++pi) {
            p = pkd->Particle(pi);
            if (pkdIsBH(pkd,p)) {
                smReSmoothSingle(smx,smf,p,2.*pkdBall(pkd,p));
                nSmoothed++;
            }
        }
        break;

#ifdef FEEDBACK
    /* IA: If computing the hydrostep, we also do the smooth over the newly formed stars that has not yet exploded, such that
     *  they can increase the rung of the neighbouring gas particles before exploding
     */
    case SMX_HYDRO_STEP:
        for (pi=0; pi<pkd->Local(); ++pi) {
            p = pkd->Particle(pi);
            if (pkdIsGas(pkd,p)) {
                if (pkdIsActive(pkd,p)) {
                    smReSmoothSingle(smx,smf,p, 2.*pkdBall(pkd,p));
                    nSmoothed++;
                }
            }
            if (pkdIsStar(pkd,p)) {
                if (pkdStar(pkd,p)->hasExploded==0) {
                    // IA: In principle this does NOT improve the Isolated Galaxy case, as we wait until the end of the
                    // step to update the primitive variables
                    //if ( (smf->dTime/*+pkd->param.dDelta/(1<<p->uRung)*/-pkdStar(pkd,p)->fTimer) < 0.95*pkd->param.dFeedbackDelay)
                    //if (pkdIsStar(pkd,p)) printf("SN dt\n");
                    //smReSmoothSingle(smx,smf,p, 2.*pkdBall(pkd,p));
                    //nSmoothed++;
                    //}
                }
            }
        }
        break;

    case SMX_SN_FEEDBACK:

        for (pi=0; pi<pkd->Local(); ++pi) {
            p = pkd->Particle(pi);
            if (pkdIsStar(pkd,p)) {
                if ( (pkdStar(pkd,p)->hasExploded == 0) &&
                        ((smf->dTime-pkdStar(pkd,p)->fTimer) > smf->dSNFBDelay) ) {
                    smReSmoothSingle(smx,smf,p, 2.*pkdBall(pkd,p));
                    pkdStar(pkd,p)->hasExploded = 1;

                    nSmoothed++;
                }
            }
        }
        break;
#endif
#ifdef STELLAR_EVOLUTION
    case SMX_CHEM_ENRICHMENT:
        for (pi = 0; pi < pkd->Local(); ++pi) {
            p = pkd->Particle( pi);
            if (pkdIsStar(pkd, p)) {
                STARFIELDS *pStar = pkdStar(pkd, p);
                if ((float)smf->dTime > pStar->fNextEnrichTime) {
                    smReSmoothSingle(smx, smf, p, 2.0 * pkdBall(pkd, p));
                    nSmoothed++;
                }
            }
        }
        break;
#endif
    default:
        for (pi=0; pi<pkd->Local(); ++pi) {
            p = pkd->Particle(pi);
            if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
                smReSmoothSingle(smx,smf,p, 2.*pkdBall(pkd,p));
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

    *p_buff = (my_real *) _mm_malloc(N*nVar*sizeof(my_real), 64);
    assert(*p_buff!=NULL);
    if (oldBuff == NULL) {
        *p_ptrs = (my_real **)_mm_malloc(nVar*sizeof(my_real *), 64);
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

    _mm_free(*p_buff);

    *p_buff = tmp_buffer;
}

/* IA: In this version, we loop over the buckets, rather than over the particles.
 *
 * For each bucket, we look for all the surroiding buckets that may interact
 * with any particle in said bucket.
 *
 * Then, we put all those particles (including the own bucket)
 * in a interaction list.
 */
int  smReSmoothNode(SMX smx,SMF *smf, int iSmoothType) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    int pj, pk, nCnt;
    double dx, dy, dz;
    double fDist2;
    PARTICLE *p;
    int nSmoothed=0;

    smx->nnListSize = 0;
    int nnListMax_p = NNLIST_INCREMENT;
    KDN *node;
    BND bnd_node;

    NN *nnList_p;
    nnList_p = static_cast<NN *>(malloc(sizeof(NN)*nnListMax_p));

    // Here we store the pointers to the particle whose interaction
    // need to be computed
    PARTICLE **sinks;
    sinks = static_cast<PARTICLE **>(malloc(64*sizeof(PARTICLE *))); // At most, the size of the bucket

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




    for (int i=NRESERVED_NODES; i<pkd->Nodes()-1; i++) {
        node = pkd->TreeNode(i);
        if (!node->iLower) { // We are in a bucket

            // Prepare the interaction list

            bnd_node = pkdNodeGetBnd(pkd,node);


            //printf("fBall %e nodeBall %e \n", fBall, pkdNodeBall(pkd,node));
            // Size of the ball that contains all possible particles
            // interacting with this bucket

            double r[3];
            r[0] = bnd_node.fCenter[0];
            r[1] = bnd_node.fCenter[1];
            r[2] = bnd_node.fCenter[2];

            // First, we add all the particles whose interactions need to be computed
            int nActive = 0;
            float nodeBall = 0.;
#ifdef OPTIM_REORDER_IN_NODES
            int pEnd = node->pLower + pkdNodeNgas(pkd,node);
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
            if (iSmoothType==SMX_HYDRO_DENSITY) pEnd += pkdNodeNstar(pkd,node);
#endif
#ifdef BLACKHOLES
            if (iSmoothType==SMX_HYDRO_DENSITY) pEnd += pkdNodeNbh(pkd,node);
#endif
#else // OPTIM_REORDER_IN_NODES
            int pEnd = node->pUpper+1;
#endif

            double fMax_shrink[3] = {0.,0.,0.};
            for (pj=node->pLower; pj<pEnd; ++pj) {
                p = pkd->Particle(pj);

#ifdef OPTIM_AVOID_IS_ACTIVE
                int pIsActive = p->bMarked;
#else
                int pIsActive = pkdIsActive(pkd,p);
#endif


                if (pIsActive) {
                    if (iSmoothType==SMX_HYDRO_DENSITY) {

#ifndef OPTIM_AVOID_IS_ACTIVE
                        if (!p->bMarked)
                            continue;
#endif

#if defined(FEEDBACK) && !defined(STELLAR_EVOLUTION)
                        // If there is only feedback, once the particle explodes
                        // there is no need to updated its fBall.
                        // However, if we have stellar evolution, fBall needs to be
                        // always updated.
                        if (pkdIsStar(pkd,p) && (pkdStar(pkd,p)->hasExploded==1))
                            continue;
#endif

#ifndef OPTIM_REORDER_IN_NODES
                        // Explicit check of the types
                        if (!pkdIsGas(pkd,p) && !pkdIsStar(pkd,p))
                            continue;
#endif
                    }
                    else {
#ifndef OPTIM_REORDER_IN_NODES
                        // Explicit check of the types
                        if (!pkdIsGas(pkd,p)) continue;
#endif
                    } //SMX_HYDRO_DENSITY

                    for (int j=0; j<3; j++) {
                        const double disp = fabs(pkdPos(pkd,p,j) - bnd_node.fCenter[j]) + pkdBall(pkd,p)*2.;
                        fMax_shrink[j] = (disp > fMax_shrink[j]) ? disp : fMax_shrink[j];
                    }

                    if (nodeBall<pkdBall(pkd,p)) nodeBall=pkdBall(pkd,p);
                    sinks[nActive] = p;
                    nActive++;
                }
            }
            // There are no elligibles particle in this bucket, go to the next
            if (nActive==0) continue;

            //printf("nodeBall %e nActive %d \n", nodeBall, nActive);
            nCnt = 0;
            //printf("%e %e \n", 2.*nodeBall, pkdNodeBall(pkd,node));
            int nCnt_own = nActive;
            //printf("start node %d %d \n", pkd->Self(), i);

            // Remember! pkdBall gives HALF the radius of the enclosing sphere!
            nodeBall *= 2.;


            for (int j=0; j<3; j++)
                bnd_node.fMax[j] = fMax_shrink[j];

            if (smx->bPeriodic) {
                double iStart[3], iEnd[3];
                for (int j=0; j<3; ++j) {
                    iStart[j] = d2i(floor((r[j] - bnd_node.fMax[j])/pkd->fPeriod[j] + 0.5));
                    iEnd[j] = d2i(floor((r[j] + bnd_node.fMax[j])/pkd->fPeriod[j] + 0.5));
                }
                for (int ix=iStart[0]; ix<=iEnd[0]; ++ix) {
                    r[0] = bnd_node.fCenter[0] - ix*pkd->fPeriod[0];
                    for (int iy=iStart[1]; iy<=iEnd[1]; ++iy) {
                        r[1] = bnd_node.fCenter[1] - iy*pkd->fPeriod[1];
                        for (int iz=iStart[2]; iz<=iEnd[2]; ++iz) {
                            r[2] = bnd_node.fCenter[2] - iz*pkd->fPeriod[2];
                            buildInteractionList(smx, smf, node, bnd_node, &nCnt, r, ix, iy, iz);
                        }
                    }
                }
            }
            else {
                buildInteractionList(smx, smf, node, bnd_node, &nCnt, r, 0, 0, 0);
            }



            //printf("interaction list completed nCnt %d nCnt_own %d nActive  %d \n", nCnt, nCnt_own, nActive);

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
                hydroDensity_node(pkd, smf, bnd_node, sinks, smx->nnList,
                                  nCnt_own, nCnt);
            }
            else {
                for (pj=0; pj<nCnt_own; pj++) {
                    PARTICLE *partj = sinks[pj];
                    float fBall2_p = 4.*pkdBall(pkd,partj)*pkdBall(pkd,partj);
                    double dx_node = -pkdPos(pkd,partj,0)+bnd_node.fCenter[0];
                    double dy_node = -pkdPos(pkd,partj,1)+bnd_node.fCenter[1];
                    double dz_node = -pkdPos(pkd,partj,2)+bnd_node.fCenter[2];

                    int nCnt_p = 0;
                    for (pk=0; pk<nCnt; pk++) {
                        dx = -dx_node + smx->nnList[pk].dx;
                        dy = -dy_node + smx->nnList[pk].dy;
                        dz = -dz_node + smx->nnList[pk].dz;

                        fDist2 = dx*dx + dy*dy + dz*dz;
                        if (fDist2 < fBall2_p) {
                            PARTICLE *q = smx->nnList[pk].pPart;
                            float qh = pkdBall(pkd,q);
                            if (fDist2==0.)
                                continue;
                            if (*pkdParticleID(pkd,partj) == *pkdParticleID(pkd,q))
                                continue;

                            // Reasons not to compute this interaction
                            if ( iSmoothType==SMX_HYDRO_FLUX ||
                                    iSmoothType==SMX_HYDRO_FLUX_VEC) {

                                if (4.*qh*qh < fDist2)
                                    continue;

#ifdef OPTIM_AVOID_IS_ACTIVE
                                int qIsActive = q->bMarked;
#else
                                int qIsActive = pkdIsActive(pkd,q);
#endif

#ifdef OPTIM_NO_REDUNDANT_FLUXES
                                if (qIsActive) {
                                    if (dx > 0) continue;
                                    else if (dx==0) {
                                        if (dy > 0) continue;
                                        else if (dy==0) {
                                            if (dz > 0) continue;
                                            else if (dz==0) abort();
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

                            nnList_p[nCnt_p].fDist2 = fDist2;
                            nnList_p[nCnt_p].dx = dx;
                            nnList_p[nCnt_p].dy = dy;
                            nnList_p[nCnt_p].dz = dz;
                            nnList_p[nCnt_p].pPart = smx->nnList[pk].pPart;
                            nnList_p[nCnt_p].iIndex = smx->nnList[pk].iIndex;
                            nnList_p[nCnt_p].iPid = smx->nnList[pk].iPid;

                            if (smx->fcnSmoothFillBuffer) {
                                PARTICLE *q = smx->nnList[pk].pPart;

                                smx->fcnSmoothFillBuffer(input_pointers, q, nCnt_p,
                                                         fDist2, dx, dy, dz, smf);
                            }


                            nCnt_p++;
                        }

                    }

                    //abort();
                    //printf("nCnt_p %d \n", nCnt_p);
                    //assert(nCnt_p<200);

                    if (smx->fcnSmoothNode) {
                        smx->fcnSmoothNode(partj,pkdBall(pkd,partj),nCnt_p,
                                           input_pointers, output_pointers, smf);
                        for (pk=0; pk<nCnt_p; pk++) {
                            smx->fcnSmoothUpdate(output_pointers,input_pointers,
                                                 partj, nnList_p[pk].pPart, pk, smf);
                        }
                    }
                    else {
                        smx->fcnSmooth(partj,pkdBall(pkd,partj),nCnt_p,nnList_p,smf);
                    }
                }
            }



            nSmoothed += nCnt_own;

            for (pk=0; pk<nCnt; ++pk) {
                if (smx->nnList[pk].iPid != pkd->Self()) {
                    mdlRelease(mdl,CID_PARTICLE,smx->nnList[pk].pPart);
                }
            }

            //printf("end node %d %d \n", pkd->Self(), i);
        }
    }
    if (smx->fcnSmoothGetNvars) {
        _mm_free(input_buffer);
        _mm_free(input_pointers);
        _mm_free(output_buffer);
        _mm_free(output_pointers);
    }
    free(nnList_p);
    free(sinks);
    //printf("nSmoothed %d \n", nSmoothed);
    return nSmoothed;
}



void buildInteractionList(SMX smx, SMF *smf, KDN *node, BND bnd_node, int *nCnt_tot, double r[3], int ix, int iy, int iz) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    int id, sp, iCell, pEnd, pj;
    double dx, dy, dz, p_r[3], fDist2;
    KDN *kdn;
    BND bnd;
    struct smContext::stStack *S = smx->ST;
    int nCnt = *nCnt_tot;

    // We look for the biggest node that encloses the needed domain
    id = pkd->Self();

// We can only take advantage of this if we are are in the original cell
    kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = pkd->Self());

    //  Now we start the walk as usual
    sp = 0;
    while (1) {
        bnd = pkdNodeGetBnd(pkd, kdn);
        for (int bnd_j=0; bnd_j<3; bnd_j++) {
            if (fabs(bnd.fCenter[bnd_j]-r[bnd_j]) - bnd.fMax[bnd_j] - bnd_node.fMax[bnd_j] > 0. ) goto NoIntersect;
        }


        /*
        ** We have an intersection to test.
        */
        if (kdn->iLower) {
            int idUpper,iUpper;
            pkdGetChildCells(kdn,id,id,iCell,idUpper,iUpper);
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
                pEnd = kdn->pLower+pkdNodeNgas(pkd,kdn);
#else
                pEnd = kdn->pUpper+1;
#endif
                //printf("pEnd %d \n", pEnd);
                for (pj=kdn->pLower; pj<pEnd; ++pj) {
                    p = pkd->Particle(pj);
#ifndef OPTIM_REORDER_IN_NODES
                    if (!pkdIsGas(pkd,p)) continue;
#endif
                    pkdGetPos1(pkd,p,p_r);
                    dx = r[0] - p_r[0];
                    dy = r[1] - p_r[1];
                    dz = r[2] - p_r[2];
                    if (fabs(dx) <= bnd_node.fMax[0] &&
                            fabs(dy) <= bnd_node.fMax[1] &&
                            fabs(dz) <= bnd_node.fMax[2] ) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            //printf("realloc \n");
                            assert(smx->nnList != NULL);
                        }
                        fDist2 = dx*dx + dy*dy + dz*dz;
                        smx->nnList[nCnt].fDist2 = fDist2;
                        smx->nnList[nCnt].dx = dx;
                        smx->nnList[nCnt].dy = dy;
                        smx->nnList[nCnt].dz = dz;
                        smx->nnList[nCnt].pPart = p;
                        smx->nnList[nCnt].iIndex = pj;
                        smx->nnList[nCnt].iPid = pkd->Self();
                        ++nCnt;
                    }
                }
            }
            else {
#ifdef OPTIM_REORDER_IN_NODES
                pEnd = kdn->pLower+pkdNodeNgas(pkd,kdn);
#else
                pEnd = kdn->pUpper+1;
#endif
                for (pj=kdn->pLower; pj<pEnd; ++pj) {
                    p = static_cast<PARTICLE *>(mdlFetch(mdl,CID_PARTICLE,pj,id));
#ifndef OPTIM_REORDER_IN_NODES
                    if (!pkdIsGas(pkd,p)) continue;
#endif
                    pkdGetPos1(pkd,p,p_r);
                    dx = r[0] - p_r[0];
                    dy = r[1] - p_r[1];
                    dz = r[2] - p_r[2];
                    if (fabs(dx) <= bnd_node.fMax[0] &&
                            fabs(dy) <= bnd_node.fMax[1] &&
                            fabs(dz) <= bnd_node.fMax[2] ) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            //printf("realloc \n");
                            assert(smx->nnList != NULL);
                        }

                        fDist2 = dx*dx + dy*dy + dz*dz;
                        smx->nnList[nCnt].fDist2 = fDist2;
                        smx->nnList[nCnt].dx = dx;
                        smx->nnList[nCnt].dy = dy;
                        smx->nnList[nCnt].dz = dz;
                        smx->nnList[nCnt].pPart = static_cast<PARTICLE *>(mdlAcquire(mdl,CID_PARTICLE,pj,id));

                        // This should be faster regarding caching and memory transfer, but the call to pkdIsActive can be a bottleneck here!
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
