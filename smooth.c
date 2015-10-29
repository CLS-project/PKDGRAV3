#ifdef HAVE_CONFIG_H
#include "config.h"
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
#include "rbtree.h"
#include <sys/stat.h>

#define ASSERT_CONCAT_(a, b) a##b
#define ASSERT_CONCAT(a, b) ASSERT_CONCAT_(a, b)
/* These can't be used after statements in c89. */
#ifdef __COUNTER__
#define STATIC_ASSERT(e,m)						\
    { enum { ASSERT_CONCAT(static_assert_, __COUNTER__) = 1/(!!(e)) }; }
#else
/* This can't be used twice on the same line so ensure if using in headers
 * that the headers are not included twice (by wrapping in #ifndef...#endif)
 * Note it doesn't cause an issue when used on same line of separate modules
 * compiled with gcc -combine -fwhole-program.  */
#define STATIC_ASSERT(e,m) \
    { enum { ASSERT_CONCAT(assert_line_, __LINE__) = 1/(!!(e)) }; }
#endif

#if defined(USE_SIMD) && defined(__SSE2__)
static inline FLOAT mindist(const BND *bnd,const FLOAT *pos) {
#ifdef __AVX__
    typedef union {
	uint64_t i[4];
	__m256d p;
	} vint64;
    static const vint64 isignmask = {{0x8000000000000000,0x8000000000000000,0x8000000000000000,0x8000000000000000}};
    static const vint64 zero = {{0,0,0,0}};
    __m256i justthree = {-1,-1,-1,0};
    __m256d fCenter = _mm256_maskload_pd(bnd->fCenter,justthree);
    __m256d fMax = _mm256_maskload_pd(bnd->fMax,justthree);
    __m256d ppos = _mm256_maskload_pd(pos,justthree);
    __m128d d;
    __m256d m = _mm256_max_pd(zero.p,_mm256_sub_pd(_mm256_andnot_pd(isignmask.p,_mm256_sub_pd(fCenter,ppos)),fMax));
    m = _mm256_mul_pd(m,m);
    d = _mm_hadd_pd(_mm256_extractf128_pd(m,1),_mm256_castpd256_pd128(m));
    return _mm_cvtsd_f64(_mm_hadd_pd(d,d));
#else
    typedef union {
	uint64_t i[2];
	__m128d p;
	} vint64;
    static const vint64 isignmask = {{0x8000000000000000,0x8000000000000000}};
    static const vint64 zero = {{0,0}};
    __m128d m2,m1;
    m1 = _mm_max_sd(zero.p,_mm_sub_sd(_mm_andnot_pd(isignmask.p,_mm_sub_sd(_mm_load_sd(bnd->fCenter+0),_mm_load_sd(pos+0))),_mm_load_sd(bnd->fMax+0)));
    m2 = _mm_mul_sd(m1,m1);
    m1 = _mm_max_sd(zero.p,_mm_sub_sd(_mm_andnot_pd(isignmask.p,_mm_sub_sd(_mm_load_sd(bnd->fCenter+1),_mm_load_sd(pos+1))),_mm_load_sd(bnd->fMax+1)));
    m2 = _mm_add_sd(m2,_mm_mul_sd(m1,m1));
    m1 = _mm_max_sd(zero.p,_mm_sub_sd(_mm_andnot_pd(isignmask.p,_mm_sub_sd(_mm_load_sd(bnd->fCenter+2),_mm_load_sd(pos+2))),_mm_load_sd(bnd->fMax+2)));
    m2 = _mm_add_sd(m2,_mm_mul_sd(m1,m1));
    return _mm_cvtsd_f64(m2);
#endif
    }
#define MINDIST(bnd,pos,min2) ((min2) = mindist(bnd,pos))
#else
#define MINDIST(bnd,pos,min2) {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;					\
    for (BND_j=0;BND_j<3;++BND_j) {\
	BND_dMin = fabs((bnd)->fCenter[BND_j] - (pos)[BND_j]) - (bnd)->fMax[BND_j]; \
	if (BND_dMin > 0) (min2) += BND_dMin*BND_dMin;			\
	}\
    }
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
    void (*comb)(void *,void *,void *) = NULL;
    int i,pi,j;
    int nTree;
    int iTopDepth;

    smx = malloc(sizeof(struct smContext));
    assert(smx != NULL);
    smx->pSentinel = malloc(pkdParticleSize(pkd));
    assert(smx->pSentinel != NULL);
    smx->pkd = pkd;
    if (smf != NULL) smf->pkd = pkd;
    smx->nSmooth = nSmooth;
    smx->bPeriodic = bPeriodic;
    /*
    ** Initialize the context for compressed nearest neighbor lists.
    */
    smx->lcmp = lcodeInit(pkd->nThreads,pkd->idSelf,pkd->nLocal,nSmooth);

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
	smx->fcnSmooth = DensityF1;
	initParticle = initDensity; /* Original Particle */
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
    case SMX_DENSITY_M3:
	assert(pkd->oGroup);
	smx->fcnSmooth = DensityM3;
	initParticle = initDensity; /* Original Particle */
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
    case SMX_GRADIENT_M3:
	smx->fcnSmooth = LinkGradientM3;
	initParticle = NULL;
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
    case SMX_HOP_LINK:
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
    case SMX_DENDVDX:
	assert( pkd->oSph ); /* Validate memory model */
	smx->fcnSmooth = DenDVDX;
	initParticle = NULL; /* Original Particle */
	init = initDenDVDX; /* Cached copies */
	comb = combDenDVDX;
	smx->fcnPost = NULL;
	break;
    case SMX_SPHFORCES:
	assert( pkd->oSph ); /* Validate memory model */
	assert( pkd->oAcceleration ); /* Validate memory model */
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
	assert( pkd->oVelSmooth); /* Validate memory model */
	smx->fcnSmooth = bSymmetric?MeanVelSym:MeanVel;
	initParticle = initMeanVel; /* Original Particle */
	init = initMeanVel; /* Cached copies */
	comb = combMeanVel;
	smx->fcnPost = NULL;
	break;
    case SMX_DIVV:
	assert( pkd->oVelSmooth); /* Validate memory model */
	smx->fcnSmooth = bSymmetric?DivvSym:Divv;
	initParticle = initDivv; /* Original Particle */
	init = initDivv; /* Cached copies */
	comb = combDivv;
	smx->fcnPost = NULL;
	break;
    case SMX_VELDISP2:
	assert( pkd->oVelSmooth); /* Validate memory model */
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
	assert( pkd->oRelaxation); /* Validate memory model */
	assert(bSymmetric == 0);
	smx->fcnSmooth = AddRelaxation;
	initParticle = NULL;
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
#ifdef SYMBA
    case SMX_SYMBA:
	assert(bSymmetric == 0);
	smx->fcnSmooth = DrmininDrift;
	initParticle = NULL;
	init = NULL;
	comb = NULL;
	smx->fcnPost = NULL;
	break;
#endif /* SYMBA */

    default:
	assert(0);
    }
    /*
    ** Initialize the ACTIVE particles in the tree.
    ** There are other particles in the tree -- just not active.
    */
    nTree = pkdTreeNode(pkd,ROOT)->pUpper + 1;
    if (initParticle != NULL) {
	for (pi=0;pi<nTree;++pi) {
	    PARTICLE *p = pkdParticle(pkd,pi);
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
		       pkdParticleBase(pkd),pkdParticleSize(pkd),
		       nTree,pkd,init,comb);
	    }
	else {
	    mdlROcache(pkd->mdl,CID_PARTICLE,NULL,
		       pkdParticleBase(pkd),pkdParticleSize(pkd),
		       nTree);
	    }
	}
    else smx->bOwnCache = 0;
    /*
    ** Allocate Nearest-Neighbor List.
    */
    smx->nnListSize = 0;
    smx->nnListMax = NNLIST_INCREMENT;
    smx->nnList = malloc(smx->nnListMax*sizeof(NN));
    assert(smx->nnList != NULL);

    /*
    ** Allocate priority queue.
    */
    smx->pq = malloc(nSmooth*sizeof(PQ));
    assert(smx->pq != NULL);
    PQ_INIT(smx->pq,nSmooth);
    /*
    ** Allocate hash table entries.
    ** The constant here just sets the hash table loading factor, for numbers larger than
    ** the 1000'th prime we end up using the result here as the hash table modulus.
    */
    smx->nHash = (int)floor(nSmooth*1.543765241931);
    for (i=0;i<1000;++i) {
	if (primes[i] > smx->nHash) {
	    smx->nHash = primes[i];
	    break;
	}
    }
    smx->pHash = malloc((smx->nHash+nSmooth)*sizeof(struct hashElement));
    assert(smx->pHash != NULL);
    for (i=0;i<smx->nHash;++i) {
	smx->pHash[i].p = NULL;
	smx->pHash[i].coll = NULL;
    }
    /*
    ** set up the extra entries that may be needed for collision chains
    */
    smx->pFreeHash = &smx->pHash[i];
    for (;i<(smx->nHash+nSmooth-1);++i) {
	smx->pHash[i].p = NULL;
	smx->pHash[i].coll = &smx->pHash[i+1];
    }
    smx->pHash[i].p = NULL;
    smx->pHash[i].coll = NULL;
    /*
    ** Allocate special stacks for searching within the tree.
    ** 1024 is more than enough.
    */
    smx->ST = malloc(1024*sizeof(struct stStack));
    assert(smx->ST != NULL);
    /*
    ** Set up the sentinel particle with some very far away distance.
    ** This is used to initially load the priority queue and all references
    ** to this particle should end up being replaced in the priority queue
    ** as long as there are nSmooth particles set bSrcActive=1.
    */
    for (j=0;j<3;++j) {
#ifdef INTEGER_POSITION
	smx->pSentinel->rPRIVATE[j] = 0x7fffffff;
#else
	pkdSetPos(smx->pSentinel->r,j,HUGE_VAL);
#endif
    }
    smx->pSentinel->bSrcActive = 1;
    smx->pSentinel->bDstActive = 0;
    /*
    ** Need to cast the pLite to an array of extra stuff.
    */
    STATIC_ASSERT(EPHEMERAL_BYTES >= sizeof(struct smExtraArray),
	"The pLite structure is not large enough");
    STATIC_ASSERT(EPHEMERAL_BYTES >= sizeof(struct smGroupArray)/4 + sizeof(uint32_t),
	"The pLite structure is not large enough");
    smx->ea = (struct smExtraArray *)(pkd->pLite); /* Used only for SPH */
    smx->pl = (uint32_t *)(pkd->pLite); /* Used for Hop: particle links */
    smx->ga = (struct smGroupArray *)(smx->pl + pkd->nLocal + 1); /* Hop: Skeleton groups nGroups << nLocal/4 */
    *psmx = smx;
    return(1);
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
    sprintf(achOut, "Cell Accesses: %g\n",
	    mdlNumAccess(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Miss ratio: %g\n",
	    mdlMissRatio(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(smx->pkd->mdl,CID_CELL));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "Particle Accesses: %g\n",
	    mdlNumAccess(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Miss ratio: %g\n",
	    mdlMissRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(smx->pkd->mdl,CID_PARTICLE));
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
	for (pi=0;pi<pkd->nLocal;++pi) {
	    p = pkdParticle(pkd,pi);
	    if ( pkdIsSrcActive(p,0,MAX_RUNG) && pkdIsDstActive(p,0,MAX_RUNG) )
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
    free(smx->ST);
    free(smx->pq);
    free(smx->nnList);
    free(smx->pHash);
    free(smx->pSentinel);
    free(smx);
}

static KDN *getCell(PKD pkd, int iCell, int id) {
    if (id==pkd->idSelf) return pkdTreeNode(pkd,iCell);
    return mdlFetch(pkd->mdl,CID_CELL,iCell,id);
    }

PQ *pqSearch(SMX smx,PQ *pq,FLOAT r[3]) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    KDN *kdn;
    int idSelf = smx->pkd->idSelf;
    struct stStack *S = smx->ST;
    FLOAT dMin,min1,min2;
    double p_r[3];
    int j,iCell,id;
    int sp = 0;
    const BND *bnd;
    int pEnd, pj;
    PARTICLE *p;
    FLOAT dx,dy,dz,fDist2;

    /* Start at the root node of the tree */
    kdn = getCell(pkd,ROOT,id = 0);
    while (1) {
	while (kdn->iLower) {
	    int idLower,iLower,idUpper,iUpper;
	    pkdGetChildCells(kdn,id,idLower,iLower,idUpper,iUpper);
	    kdn = getCell(pkd,iLower,idLower);
            bnd = pkdNodeBnd(pkd, kdn);
	    MINDIST(bnd,r,min1);
	    kdn = getCell(pkd,iUpper,idUpper);
            bnd = pkdNodeBnd(pkd, kdn);
	    MINDIST(bnd,r,min2);
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
	    for (pj=kdn->pLower;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		if (!p->bMarked) continue;
		pkdGetPos1(p->r,p_r);
		dx = r[0] - p_r[0];
		dy = r[1] - p_r[1];
		dz = r[2] - p_r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= pq->fDist2) {
		    if (pq->iPid == idSelf) {
			pkdParticle(pkd,pq->iIndex)->bMarked = 1;
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
	    for (pj=kdn->pLower;pj<=pEnd;++pj) {
		p = mdlFetch(mdl,CID_PARTICLE,pj,id);
		if (!p->bSrcActive || smHashPresent(smx,p)) continue;
		pkdGetPos1(p->r,p_r);
		dx = r[0] - p_r[0];
		dy = r[1] - p_r[1];
		dz = r[2] - p_r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= pq->fDist2) {
		    if (pq->iPid == idSelf) {
			pkdParticle(pkd,pq->iIndex)->bMarked = 1;
			}
		    else {
			smHashDel(smx,pq->pPart);
			mdlRelease(mdl,CID_PARTICLE,pq->pPart);
			}
		    p = mdlAcquire(mdl,CID_PARTICLE,pj,id);
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
    for (i=0;i<smx->nSmooth;++i) {
	smx->pq[i].pPart = smx->pSentinel;
	smx->pq[i].iIndex = smx->pkd->nLocal;
	smx->pq[i].iPid = smx->pkd->idSelf;
	smx->pq[i].dx = pkdPos(smx->pSentinel->r,0);
	smx->pq[i].dy = pkdPos(smx->pSentinel->r,1);
	smx->pq[i].dz = pkdPos(smx->pSentinel->r,2);
	smx->pq[i].fDist2 = pow(smx->pq[i].dx,2) + pow(smx->pq[i].dy,2) + 
	    pow(smx->pq[i].dz,2);
	}
    for (i=0;i<3;++i) smx->rLast[i] = 0.0;
    }

void smSmoothFinish(SMX smx) {
    int i;
    /*
    ** Release acquired pointers and source-reactivate particles in prioq.
    */
    for (i=0;i<smx->nSmooth;++i) {
	if (smx->pq[i].iPid == smx->pkd->idSelf) {
	    pkdParticle(smx->pkd,smx->pq[i].iIndex)->bMarked = 1;
	    }
	else {
	    smHashDel(smx,smx->pq[i].pPart);
	    mdlRelease(smx->pkd->mdl,CID_PARTICLE,smx->pq[i].pPart);
	    }
	}
    }

float smSmoothSingle(SMX smx,SMF *smf,PARTICLE *p) {
    PKD pkd = smx->pkd;
    int ix,iy,iz;
    double p_r[3];
    FLOAT r[3],fBall;
    int iStart[3],iEnd[3];
    int j;
    PQ *pq;

    pkdGetPos1(p->r,p_r);

    /*
    ** Correct distances and rebuild priority queue.
    */
    if (smx->bPeriodic) {
	for (j=0;j<3;++j) {
	    if (p_r[j] > smx->rLast[j] + 0.5 * pkd->fPeriod[j])
		smx->rLast[j] += pkd->fPeriod[j];
	    else if (p_r[j] < smx->rLast[j] - 0.5 * pkd->fPeriod[j])
		smx->rLast[j] -= pkd->fPeriod[j];
	    }
	}

    for (j=0;j<smx->nSmooth;++j) {
	smx->pq[j].dx += p_r[0]-smx->rLast[0];
	smx->pq[j].dy += p_r[1]-smx->rLast[1];
	smx->pq[j].dz += p_r[2]-smx->rLast[2];
	smx->pq[j].fDist2 = pow(smx->pq[j].dx,2) + pow(smx->pq[j].dy,2) + 
	    pow(smx->pq[j].dz,2);
	}
    for (j=0;j<3;++j) smx->rLast[j] = r[j] = p_r[j];

    PQ_BUILD(smx->pq,smx->nSmooth,pq);
    pq = pqSearch(smx,pq,r);
    /*
    ** Search in replica boxes if it is required.
    */
    if (smx->bPeriodic) {
	fBall = sqrt(pq->fDist2);
	for (j=0;j<3;++j) {
	    iStart[j] = d2i(floor((p_r[j] - fBall)/pkd->fPeriod[j] + 0.5));
	    iEnd[j] = d2i(floor((p_r[j] + fBall)/pkd->fPeriod[j] + 0.5));
	    }
	for (ix=iStart[0];ix<=iEnd[0];++ix) {
	    r[0] = p_r[0] - ix*pkd->fPeriod[0];
	    for (iy=iStart[1];iy<=iEnd[1];++iy) {
		r[1] = p_r[1] - iy*pkd->fPeriod[1];
		for (iz=iStart[2];iz<=iEnd[2];++iz) {
		    r[2] = p_r[2] - iz*pkd->fPeriod[2];
		    if (ix || iy || iz) {
			pq = pqSearch(smx,pq,r);
			}
		    }
		}
	    }
	}
    fBall = sqrt(pq->fDist2);
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

    /*
    ** Initialize the bInactive flags for all local particles.
    */
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	p->bMarked = p->bSrcActive;
    }
    smSmoothInitialize(smx);
    smf->pfDensity = NULL;
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;
	pkdSetBall(pkd,p,smSmoothSingle(smx,smf,p));

	/*
	** Call mdlCacheCheck to make sure we are making progress!
	*/
	mdlCacheCheck(pkd->mdl);
    }
    smSmoothFinish(smx);
}


void UpdateSphBounds(SMX smx) {
#if 0
    PKD pkd = smx->pkd;
    PARTICLE *p;
    KDN *pkdn,*p1,*p2;
    SPHBNDS *bn,*b1,*b2;
    double Bmin[3];
    double Bmax[3];
    double BImin[3];
    double BImax[3];
    uint32_t iNode;
    int nDepth,d,pj;

    assert(pkd->oNodeSphBounds);
    nDepth = 1;
    iNode = ROOT;
    while (1) {
	while (pkdTreeNode(pkd,iNode)->iLower) {
	    iNode = pkdTreeNode(pkd,iNode)->iLower;
	    ++nDepth;
	}
	/*
	** Now calculate all bucket quantities!
	*/
	pkdn = pkdTreeNode(pkd,iNode);
	/*
	** Update bucket fast gas bounds.
	** Default bounds always makes the cell look infinitely far away, regardless from where.
	*/
	for (d=0;d<3;++d) Bmin[d] = HUGE_VAL;
	for (d=0;d<3;++d) Bmax[d] = -HUGE_VAL;
	for (d=0;d<3;++d) BImin[d] = HUGE_VAL;
	for (d=0;d<3;++d) BImax[d] = -HUGE_VAL;
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
	    p = pkdParticle(pkd,pj);
	    if (pkdIsGas(pkd,p)) {
		if (smx->ea[pj].bDone) {
		    /*
		    ** Use the actual currently calculated fBall.
		    */
		    for (d=0;d<3;++d) Bmin[d] = fmin(Bmin[d],p->r[d] - p->fBall);
		    for (d=0;d<3;++d) Bmax[d] = fmax(Bmax[d],p->r[d] + p->fBall);
		    /*
		    ** If the particle is inactive and we have calculated it already, then we
		    ** actually don't want it included in our BI ball bound. No test needed here.
		    ** if (!pkdActive(pkd,p)) Ignore this particle in the BI ball bound;
		    */
		}
		else {
		    /*
		    ** Use the old fBall, but increased by the factor given by the maximum growth of hsmooth.
		    */
		    for (d=0;d<3;++d) Bmin[d] = fmin(Bmin[d],p->r[d] - (1+pkd->param.ddHonHLimit)*p->fBall);
		    for (d=0;d<3;++d) Bmax[d] = fmax(Bmax[d],p->r[d] + (1+pkd->param.ddHonHLimit)*p->fBall);
		    if (!pkdIsActive(pkd,p)) {
			for (d=0;d<3;++d) BImin[d] = fmin(BImin[d],p->r[d] - (1+pkd->param.ddHonHLimit)*p->fBall);
			for (d=0;d<3;++d) BImax[d] = fmax(BImax[d],p->r[d] + (1+pkd->param.ddHonHLimit)*p->fBall);
		    }
		}
	    }
	    /*
	    ** Note that minimums can always safely be increased and maximums safely decreased in parallel, even on 
	    ** a shared memory machine, without needing locking since these bounds should always simply be seen as being
	    ** a conservative bound on the particles in the algorithms. This is true AS LONG AS a double precision store
	    ** operation is atomic (i.e., that the individual bytes are not written one-by-one). We take advantage of 
	    ** this fact in the fast gas algorithm where we locally reduce the bounds to exclude particles which have 
	    ** already been completed in the direct neighbor search phase.
	    */
	    bn = pkdNodeSphBounds(pkd,pkdn);
	    for (d=0;d<3;++d) bn->B.min[d] = Bmin[d];
	    for (d=0;d<3;++d) bn->B.max[d] = Bmax[d];
	    for (d=0;d<3;++d) bn->BI.min[d] = BImin[d];
	    for (d=0;d<3;++d) bn->BI.max[d] = BImax[d];
	}
	/*
	** Finished with the bucket, move onto the next one,
	** or to the parent.
	*/
	while (iNode & 1) {
	    iNode = pkdTreeNode(pkd,iNode)->iParent;
	    --nDepth;
	    if (!iNode) {
		assert(nDepth == 0);
		return;	/* exit point!!! */
	    }
	    pkdn = pkdTreeNode(pkd,iNode);
	    p1 = pkdTreeNode(pkd,pkdn->iLower);
	    p2 = pkdTreeNode(pkd,pkdn->iLower + 1);
	    b1 = pkdNodeSphBounds(pkd,p1);
	    b2 = pkdNodeSphBounds(pkd,p2);
	    bn = pkdNodeSphBounds(pkd,pkdn);
	    for (d=0;d<3;++d) bn->B.min[d] = fmin(b1->B.min[d],b2->B.min[d]);
	    for (d=0;d<3;++d) bn->B.max[d] = fmax(b1->B.max[d],b2->B.max[d]);
	    for (d=0;d<3;++d) bn->BI.min[d] = fmin(b1->BI.min[d],b2->BI.min[d]);
	    for (d=0;d<3;++d) bn->BI.max[d] = fmax(b1->BI.max[d],b2->BI.max[d]);
	}
	++iNode;
    }
#endif
}

#if 0
static inline int iOpenInactive(PKD pkd,KDN *k,CELT *check,KDN **pc,PARTICLE **pp) {
    PARTICLE *p;
    KDN *c;
    SPHBNDS *bk,*bc;
    double sk,sc;  /* size proxies for the cells */
    int iCell,iPart,j;
        
    bk = pkdNodeSphBounds(pkd,k);
    if (check->iCell < 0) {
	iPart = -check->iCell;
	if (check->id == pkd->idSelf) {
	    p = pkdParticle(pkd,iPart);
	}
	else {
	    p = mdlAcquire(pkd->mdl,CID_PARTICLE,iPart,check->id);
	}
	*pc = NULL;
	*pp = p;
	for (j=0;j<3;++j) {
	    /*
	    ** If any of the dimensions show no overlap, then the 2 bounds do not overlap.
	    */
	    if (bk->BI.max[j] < p->r[j]+check->rOffset[j] || 
		bk->BI.min[j] > p->r[j]+check->rOffset[j]) return(10); /* ignore particle */
	}
	if (k->iLower) return(0); /* particle stays on the checklist (open k) */
	else return(1);  /* test particle by particle */
    }
    else {
	iCell = check->iCell;
	if (check->id == pkd->idSelf) {
	    c = pkdTreeNode(pkd,iCell);
	}
	else if (check->id < 0) {
	    c = pkdTopNode(pkd,iCell);
	    assert(c->iLower != 0);
	}
	else {
	    c = mdlAcquire(pkd->mdl,CID_CELL,iCell,check->id);
	}
	*pc = c;
	*pp = NULL;
	bc = pkdNodeSphBounds(pkd,c);
	for (j=0;j<3;++j) {
	    /*
	    ** If any of the dimensions show no overlap, then the 2 bounds do not overlap.
	    */
	    if (bk->BI.max[j] < bc->A.min[j]+check->rOffset[j] || 
		bk->BI.min[j] > bc->A.max[j]+check->rOffset[j]) return(10); /* ignore cell */
	}
	if (k->iLower) {
	    /*
	    ** Open the larger of the 2 cells. We use a Manhatten metric to define the size.
	    */
	    sk = 0; sc = 0;
	    for (j=0;j<3;++j) {
		sk += bk->BI.max[j] - bk->BI.min[j];  
		sc += bc->A.max[j] - bc->A.min[j];
	    }
	    if (sk > sc) return(0); /* cell stays on checklist (open cell k) */
	    else if (c->iLower) return(3); /* open cell c */
	    else return(2); /* open the bucket c */
	}
	else if (c->iLower) return(3); /* open cell c */
	else return(2); /* open the bucket c */
    }
}
#endif

/*
** Returns the number of elements added to the do queue.
*/
uint32_t BoundWalkInactive(SMX smx) {
#if 0
    PKD pkd = smx->pkd;
    PARTICLE *p,*pp;
    KDN *k,*c,*kSib;
    SPHBNDS *bn;
    double BImin[3];
    double BImax[3];
    FLOAT rOffset[3];
    FLOAT d2;
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxInitCheck,nCheck;
    int iCell,iSib,iCheckCell;
    int i,ii,j,pj,d;
    int iOpen;
    uint32_t nDo;

    assert(pkd->oNodeSphBounds);
    nDo = 0;
    /*
    ** Allocate Checklist.
    */
    nMaxInitCheck = 3;
    nMaxInitCheck = nMaxInitCheck*nMaxInitCheck*nMaxInitCheck;	/* all replicas */
    iCell = pkd->iTopRoot;
    while ((iCell = pkdTopNode(pkd,iCell)->iParent)) ++nMaxInitCheck; /* all top tree siblings */
    assert(nMaxInitCheck < pkd->nMaxCheck);  /* we should definitely have enough to cover us here! */
    nCheck = 0;
    iStack = -1;
    /*
    ** First we add any replicas of the entire box
    ** to the Checklist.
    */
    if (smx->bPeriodic) {
	for (ix=-1;ix<=1;++ix) {
	    rOffset[0] = ix*pkd->fPeriod[0];
	    for (iy=-1;iy<=1;++iy) {
		rOffset[1] = iy*pkd->fPeriod[1];
		for (iz=-1;iz<=1;++iz) {
		    rOffset[2] = iz*pkd->fPeriod[2];
		    bRep = ix || iy || iz;
		    if (bRep) {
			pkd->Check[nCheck].iCell = ROOT;
			/* If leaf of top tree, use root of
			   local tree.
			*/
			if (pkdTopNode(pkd,ROOT)->iLower) {
			    pkd->Check[nCheck].id = -1;
			}
			else {
			    pkd->Check[nCheck].id = pkdTopNode(pkd,ROOT)->pLower;
			}
			for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = rOffset[j];
			++nCheck;
		    }
		}
	    }
	}
    }
    /*
    ** This adds all siblings of a chain leading from the local tree leaf in the top
    ** tree up to the ROOT of the top tree.
    */
    iCell = pkd->iTopRoot;
    iSib = SIBLING(iCell);
    while (iSib) {
	if (pkdTopNode(pkd,iSib)->iLower) {
	    pkd->Check[nCheck].iCell = iSib;
	    pkd->Check[nCheck].id = -1;
	}
	else {
	    /* If leaf of top tree, use root of local tree */
	    pkd->Check[nCheck].iCell = ROOT;
	    pkd->Check[nCheck].id = pkdTopNode(pkd,iSib)->pLower;
	}
	for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
	++nCheck;
	iCell = pkdTopNode(pkd,iCell)->iParent;
	iSib = SIBLING(iCell);
    }
    /*
    ** We are now going to work on the local tree.
    ** Make iCell point to the root of the tree again.
    */
    k = pkdTreeNode(pkd,iCell = ROOT);
    /*
    ** If iCell is now a bucket we add it to its own checklist.
    */
    if (k->iLower == 0) {
	pkd->Check[nCheck].iCell = iCell;
	pkd->Check[nCheck].id = pkd->idSelf;
	for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
	++nCheck;
    }
    while (1) {
	while (1) {
	    /*
	    ** Process the Checklist.
	    */
/*
	    printf("\nCELL:%d ",iCell);
*/
	    ii = 0;
	    for (i=0;i<nCheck;++i) {
		iOpen = iOpenInactive(pkd,k,&pkd->Check[i],&c,&pp);
/*
		printf("%1d",iOpen);
*/
		switch (iOpen) {
		case 0:
		    /*
		    ** This checkcell stays on the checklist.
		    */
		    pkd->Check[ii++] = pkd->Check[i];
		    break;
		case 1:
		    /*
		    ** We check individual particles against each other here.
		    */
		    for (d=0;d<3;++d) BImin[d] = HUGE_VAL;
		    for (d=0;d<3;++d) BImax[d] = -HUGE_VAL;		    
		    for (pj=k->pLower;pj<=k->pUpper;++pj) {
			p = pkdParticle(pkd,pj);
			/*
			** Have we potentially missed this particle?
			*/
			if (pkdIsGas(pkd,p) && !pkdIsActive(pkd,p) && !smx->ea[pj].bDone) {
			    d2 = 0;
			    for (j=0;j<3;++j) {
				d2 += pow(p->r[j] - (pp->r[j] + pkd->Check[i].rOffset[j]),2);
			    }
			    if (d2 < pow((1+pkd->param.ddHonHLimit)*p->fBall,2)) {
				/*
				** Needs an updated density, so add it to the end of the 
				** do queue.
				*/
				assert(nDo < pkd->nLocal);
				smx->ea[nDo++].iIndex = pj;
/*
				printf("Added inactive:%d due to active:%d\n",pj,-pkd->Check[i].iCell);
*/
				/*
				** Set the particle as done here, so that we don't try to add it again!
				*/
				smx->ea[pj].bDone = 1;
			    }
			    else {
				/*
				** The particle remains not done, and we can use it to calculate the inactive
				** ball bounds.
				*/
				for (d=0;d<3;++d) BImin[d] = fmin(BImin[d],p->r[d] - (1+pkd->param.ddHonHLimit)*p->fBall);
				for (d=0;d<3;++d) BImax[d] = fmax(BImax[d],p->r[d] + (1+pkd->param.ddHonHLimit)*p->fBall);
			    }
			}
		    }
		    /*
		    ** If any particles were added to the Do queue in this part we can shrink
		    ** the ball bounds for inactives for this bucket at this stage! This will tend to make
		    ** the further checks trivial once all inactives in a bucket are done.
		    */
		    bn = pkdNodeSphBounds(pkd,k);
		    for (d=0;d<3;++d) bn->BI.min[d] = BImin[d];
		    for (d=0;d<3;++d) bn->BI.max[d] = BImax[d];
		    break;
		case 2:
		    /*
		    ** Now I am trying to open a bucket, which means I place each of its particles
		    ** on the checklist. These are marked by a negative cell id.
		    */
		    if (nCheck + (c->pUpper - c->pLower + 1) > pkd->nMaxCheck) {
			pkd->nMaxCheck += (c->pUpper - c->pLower + 1) + 1000;
			pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			assert(pkd->Check != NULL);
			for (ism=0;ism<pkd->nMaxStack;++ism) {
			    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->S[ism].Check != NULL);
			}
			printf("Case 2: CPU:%d increased checklist size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		    }
		    for (pj=c->pLower;pj<=c->pUpper;++pj) {
			if (pkd->Check[i].id == pkd->idSelf) p = pkdParticle(pkd,pj);
			else p = mdlAcquire(pkd->mdl,CID_PARTICLE,pj,pkd->Check[i].id);
			/*
			** Only add those particle which we really need to check here!
			*/
			if (pkdIsGas(pkd,p) && pkdIsActive(pkd,p)) {
			    pkd->Check[nCheck] = pkd->Check[i];
			    pkd->Check[nCheck].iCell = -pj;
			    ++nCheck;
			}
			if (pkd->Check[i].id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
		    }
		    break;
		case 3:
		    /*
		    ** Open the cell.
		    ** Here we ASSUME that the children of
		    ** c are all in sequential memory order!
		    ** (the new tree build assures this)
		    ** (also true for the top tree)
		    ** We could do a prefetch here for non-local
		    ** cells.
		    */
		    if (nCheck + 2 > pkd->nMaxCheck) {
			pkd->nMaxCheck += 1000;
			pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			assert(pkd->Check != NULL);
			for (ism=0;ism<pkd->nMaxStack;++ism) {
			    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->S[ism].Check != NULL);
			}
			printf("Case 3: CPU:%d increased checklist size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		    }
		    iCheckCell = c->iLower;
		    pkd->Check[nCheck] = pkd->Check[i];
		    pkd->Check[nCheck+1] = pkd->Check[i];
		    /*
		    ** If we are opening a leaf of the top tree
		    ** we need to correctly set the processor id.
		    ** (this is a bit tricky)
		    */
		    if (pkd->Check[i].id < 0) {
			if (pkdTopNode(pkd,iCheckCell)->pLower >= 0) {
			    pkd->Check[nCheck].iCell = ROOT;
			    pkd->Check[nCheck].id = pkdTopNode(pkd,iCheckCell)->pLower;
			}
			else {
			    pkd->Check[nCheck].iCell = iCheckCell;
			    assert(pkd->Check[nCheck].id == -1);
			}
			if (pkdTopNode(pkd,iCheckCell+1)->pLower >= 0) {
			    pkd->Check[nCheck+1].iCell = ROOT;
			    pkd->Check[nCheck+1].id = pkdTopNode(pkd,iCheckCell+1)->pLower;
			}
			else {
			    pkd->Check[nCheck+1].iCell = iCheckCell+1;
			    assert(pkd->Check[nCheck+1].id == -1);
			}
		    }
		    else {
			pkd->Check[nCheck].iCell = iCheckCell;
			pkd->Check[nCheck+1].iCell = iCheckCell+1;
			assert(pkd->Check[nCheck].id == pkd->Check[i].id);
			assert(pkd->Check[nCheck+1].id == pkd->Check[i].id);
		    }
		    nCheck += 2;
		    break;
		case 10:
		    /*
		    ** This checkcell is removed from the checklist since it has no overlap with the current cell.
		    */
		    break;		
		}
		if (pkd->Check[i].id >= 0 && pkd->Check[i].id != pkd->idSelf) {
		    if (c) mdlRelease(pkd->mdl,CID_CELL,c);
		    if (pp) mdlRelease(pkd->mdl,CID_PARTICLE,pp);
		}
	    }
	    nCheck = ii;
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
	    if (!k->iLower) {
		break;
	    }
	    k = pkdTreeNode(pkd,iCell = k->iLower);
	    /*
	    ** Make sure all the check lists are long enough to handle 2 more cells.
	    */
	    if (nCheck+1 >= pkd->nMaxCheck) {
		pkd->nMaxCheck += 1000;
		pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
		assert(pkd->Check != NULL);
		for (ism=0;ism<pkd->nMaxStack;++ism) {
		    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
		    assert(pkd->S[ism].Check != NULL);
		    }
		printf("F CPU:%d increased check list size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		}
	    /*
	    ** Check iCell has any inactive.
	    */
	    if (k->nActive < (k->pUpper-k->pLower+1)) {
		/*
		** iCell has inactives, continue processing it.
		** Put the sibling onto the checklist.
		*/
		pkd->Check[nCheck].iCell = iCell + 1;
		pkd->Check[nCheck].id = pkd->idSelf;
		for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
		++nCheck;
		/*
		** If iCell is now a bucket we add it to its own checklist.
		*/
		if (k->iLower == 0) {
		    pkd->Check[nCheck].iCell = iCell;
		    pkd->Check[nCheck].id = pkd->idSelf;
		    for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
		    ++nCheck;
		}
		/*
		** Test whether the sibling has inactives as well.
		** If not we don't push it onto the stack, but we
		** have to be careful to not pop the stack when we
		** hit the sibling. See the goto "InactiveAscend" below
		** for how this is done.
		*/
		kSib = pkdTreeNode(pkd,iCell+1);
		if (kSib->nActive < (kSib->pUpper-kSib->pLower+1)) {
		    /*
		    ** Sibling has inactives as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.
		    */
		    ++iStack;
		    assert(iStack < pkd->nMaxStack);
		    pkd->S[iStack].nCheck = nCheck;
		    /*
		    ** Maybe use a memcpy here!
		    ** for (i=0;i<nCheck;++i) pkd->S[iStack].Check[i] = pkd->Check[i];
		    */
		    memcpy(pkd->S[iStack].Check,pkd->Check,nCheck*sizeof(CELT));
		    pkd->S[iStack].Check[nCheck-1].iCell = iCell;
		    }
		}
	    else {
		/*
		** Skip iCell, but add it to the Checklist.
		** No need to push anything onto the stack here.
		*/
		pkd->Check[nCheck].iCell = iCell;
		pkd->Check[nCheck].id = pkd->idSelf;
		for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
		++nCheck;
		/*
		** Move onto processing the sibling.
		*/
		k = pkdTreeNode(pkd,++iCell);
		}
	    }
	/*
	** Now the interaction list should be complete and the
	** Checklist should be empty!
	*/
	assert(nCheck == 0);

	while (iCell & 1) {
	InactiveAscend:
	    k = pkdTreeNode(pkd,iCell = k->iParent);
	    if (!iCell) {
		/*
		** Make sure stack is empty.
		*/
		assert(iStack == -1);
		return(nDo);
		}
	    }
	k = pkdTreeNode(pkd,++iCell);
	if (k->nActive == (k->pUpper-k->pLower+1)) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	nCheck = pkd->S[iStack].nCheck;
	/*
	** Use a memcpy here. This is where we would win with lists since we could simply take the
	** pointer to this list. We would have to link the old checklist into the freelist.
	** for (i=0;i<nCheck;++i) Check[i] = S[iStack].Check[i];
	*/
	memcpy(pkd->Check,pkd->S[iStack].Check,nCheck*sizeof(CELT));
	--iStack;
	}
#else
    return 0;
#endif
    }

#if 0
static inline int iOpenActive(PKD pkd,KDN *k,CELT *check,KDN **pc,PARTICLE **pp) {
    PARTICLE *p;
    KDN *c;
    SPHBNDS *bk,*bc;
    double sk,sc;  /* size proxies for the cells */
    int iCell,iPart,j;
    double dx0,dy0,dz0,dx1,dy1,dz1,xc,yc,zc,mink2;
        
    bk = pkdNodeSphBounds(pkd,k);
    if (check->iCell < 0) {
	iPart = -check->iCell;
	if (check->id == pkd->idSelf) {
	    p = pkdParticle(pkd,iPart);
	}
	else {
	    p = mdlAcquire(pkd->mdl,CID_PARTICLE,iPart,check->id);
	}
	*pc = NULL;
	*pp = p;
	xc = p->r[0] + check->rOffset[0];
	yc = p->r[1] + check->rOffset[1];
	zc = p->r[2] + check->rOffset[2];
	dx0 = xc - bk->A.max[0];
	dx1 = xc - bk->A.min[0];
	dy0 = yc - bk->A.max[1];
	dy1 = yc - bk->A.min[1];
	dz0 = zc - bk->A.max[2];
	dz1 = zc - bk->A.min[2];
	mink2 = ((dx0>0)?dx0*dx0:0) + ((dx1<0)?dx1*dx1:0) +
	    ((dy0>0)?dy0*dy0:0) + ((dy1<0)?dy1*dy1:0) + 
	    ((dz0>0)?dz0*dz0:0) + ((dz1<0)?dz1*dz1:0);
	/*
	** We have to be extra careful here when comparing to the fBalls of remote particles, for which
	** we don't exactly know if they have an updated softening or not. Safe is to multiply all by 
	** the maximum hsmooth growth factor regardless of their state.
	*/
	if (mink2 > pow((1+pkd->param.ddHonHLimit)*p->fBall,2)) return(10);
	if (k->iLower) return(0); /* particle stays on the checklist (open k) */
	else return(1);  /* test particle by particle */
    }
    else {
	iCell = check->iCell;
	if (check->id == pkd->idSelf) {
	    c = pkdTreeNode(pkd,iCell);
	}
	else if (check->id < 0) {
	    c = pkdTopNode(pkd,iCell);
	    assert(c->iLower != 0);
	}
	else {
	    c = mdlAcquire(pkd->mdl,CID_CELL,iCell,check->id);
	}
	*pc = c;
	*pp = NULL;
	bc = pkdNodeSphBounds(pkd,c);
	for (j=0;j<3;++j) {
	    /*
	    ** If any of the dimensions show no overlap, then the 2 bounds do not overlap.
	    */
	    if (bk->A.max[j] < bc->B.min[j]+check->rOffset[j] || 
		bk->A.min[j] > bc->B.max[j]+check->rOffset[j]) return(10); /* ignore cell */
	}
	if (k->iLower) {
	    /*
	    ** Open the larger of the 2 cells. We use a Manhatten metric to define the size.
	    */
	    sk = 0; sc = 0;
	    for (j=0;j<3;++j) {
		sk += bk->A.max[j] - bk->A.min[j];  
		sc += bc->B.max[j] - bc->B.min[j];
	    }
	    if (sk > sc) return(0); /* cell stays on checklist (open cell k) */
	    else if (c->iLower) return(3); /* open cell c */
	    else return(2); /* open the bucket c */
	}
	else if (c->iLower) return(3); /* open cell c */
	else return(2); /* open the bucket c */
    }
}
#endif

void BoundWalkActive(SMX smx,LIST **ppList,int *pnMaxpList) {
#if 0
    PKD pkd = smx->pkd;
    PARTICLE *p,*pp;
    KDN *k,*c,*kSib;
    FLOAT rOffset[3];
    FLOAT d2;
    int iStack,ism;
    int ix,iy,iz,bRep;
    int nMaxInitCheck,nCheck;
    int iCell,iSib,iCheckCell;
    int i,ii,j,pj;
    int iOpen;
    int nList;
    char **ppCList;
    uint32_t iIndex,iPid;

    assert(pkd->oNodeSphBounds);
    /*
    ** Allocate Checklist.
    */
    nMaxInitCheck = 3;
    nMaxInitCheck = nMaxInitCheck*nMaxInitCheck*nMaxInitCheck;	/* all replicas */
    iCell = pkd->iTopRoot;

    while ((iCell = pkdTopNode(pkd,iCell)->iParent)) ++nMaxInitCheck; /* all top tree siblings */
    assert(nMaxInitCheck < pkd->nMaxCheck);  /* we should definitely have enough to cover us here! */

    nCheck = 0;
    iStack = -1;
    /*
    ** First we add the replicas of the entire box
    ** to the Checklist.
    */
    if (smx->bPeriodic) {
	for (ix=-1;ix<=1;++ix) {
	    rOffset[0] = ix*pkd->fPeriod[0];
	    for (iy=-1;iy<=1;++iy) {
		rOffset[1] = iy*pkd->fPeriod[1];
		for (iz=-1;iz<=1;++iz) {
		    rOffset[2] = iz*pkd->fPeriod[2];
		    bRep = ix || iy || iz;
		    if (bRep) {
			pkd->Check[nCheck].iCell = ROOT;
			/* If leaf of top tree, use root of
			   local tree.
			*/
			if (pkdTopNode(pkd,ROOT)->iLower) {
			    pkd->Check[nCheck].id = -1;
			}
			else {
			    pkd->Check[nCheck].id = pkdTopNode(pkd,ROOT)->pLower;
			}
			for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = rOffset[j];
			++nCheck;
		    }
		}
	    }
	}
    }
    /*
    ** This adds all siblings of a chain leading from the local tree leaf in the top
    ** tree up to the ROOT of the top tree.
    */
    iCell = pkd->iTopRoot;
    iSib = SIBLING(iCell);
    while (iSib) {
	if (pkdTopNode(pkd,iSib)->iLower) {
	    pkd->Check[nCheck].iCell = iSib;
	    pkd->Check[nCheck].id = -1;
	}
	else {
	    /* If leaf of top tree, use root of local tree */
	    pkd->Check[nCheck].iCell = ROOT;
	    pkd->Check[nCheck].id = pkdTopNode(pkd,iSib)->pLower;
	}
	for (j=0;j<3;++j) pkd->Check[nCheck].rOffset[j] = 0.0;
	++nCheck;
	iCell = pkdTopNode(pkd,iCell)->iParent;
	iSib = SIBLING(iCell);
    }
    /*
    ** We are now going to work on the local tree.
    ** Make iCell point to the root of the tree again.
    */
    k = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
	while (1) {
	    /*
	    ** Process the Checklist.
	    */
/*
	    printf("\nCELL:%d ",iCell);
*/
	    ii = 0;
	    for (i=0;i<nCheck;++i) {
		iOpen = iOpenActive(pkd,k,&pkd->Check[i],&c,&pp);
/*
		printf("%1d",iOpen);
*/
		switch (iOpen) {
		case 0:
		    /*
		    ** This checkcell stays on the checklist.
		    */
		    pkd->Check[ii++] = pkd->Check[i];
		    break;
		case 1:
		    iIndex = -pkd->Check[i].iCell;
		    iPid = pkd->Check[i].id;
		    /*
		    ** We check individual particles against each other here.
		    */
		    for (pj=k->pLower;pj<=k->pUpper;++pj) {
			p = pkdParticle(pkd,pj);
			/*
			** Is it possible that this particle is not on our list?
			*/
			if (pkdIsGas(pkd,p) && pkdIsActive(pkd,p)) {
			    d2 = 0;
			    for (j=0;j<3;++j) {
				d2 += pow(p->r[j] - (pp->r[j] + pkd->Check[i].rOffset[j]),2);
			    }
			    /*
			    ** Here the remote processor may have updated pp->fBall, but we have no easy way of 
			    ** knowing this. We assume that it hasn't been updated, but if it has been then at 
			    ** worst we are being overly conservative by an extra factor of dHonHlimit.
			    ** We can cull the extra remote particles from the active lists in the next fast gas
			    ** phase, although this is only really needed for iterative SPH calculations, such as
			    ** radiative transfer.
			    */
			    if (d2 < pow((1+pkd->param.ddHonHLimit)*pp->fBall,2)) {
				/*
				** Check if particle pp (a remote particle) is on the list of particle p (local and active).
				** If not, then add pp to p's list.
				*/
				ppCList = pkd_pNeighborList(pkd,p);
				assert(*ppCList != NULL);
				if (!bInList(smx->lcmp,*ppCList,iIndex,iPid)) {
				    /*
				    ** Add particle pp to the neighbor list of p!
				    ** Unfortunately this means throwing away the old compressed list and creating a new one.
				    ** Adding to an already compressed list would be a useful function, but we still would
				    ** usually have to allocate/reallocate new storage for the compressed list of particle pp.
				    */
				    lcodeDecode(smx->lcmp,*ppCList,ppList,pnMaxpList,&nList);
				    assert(*ppCList != NULL);
				    free(*ppCList);
				    if (nList == *pnMaxpList) {
					*pnMaxpList *= 2;
					*ppList = realloc(*ppList,(*pnMaxpList)*sizeof(LIST));
					assert(*ppList != NULL);
				    }
				    (*ppList)[nList].iIndex = iIndex;
				    (*ppList)[nList].iPid = iPid;
				    ++nList;
				    /*
				    ** Compress the final list here as much as possible.
				    */
				    qsort(*ppList,nList,sizeof(LIST),lcodeCmpList);
				    lcodeEncode(smx->lcmp,*ppList,nList,ppCList);
				}
			    }
			}
		    }
		    break;
		case 2:
		    /*
		    ** Now I am trying to open a bucket, which means I place each of its particles
		    ** on the checklist. These are marked by a negative cell id.
		    */
		    if (nCheck + (c->pUpper - c->pLower + 1) > pkd->nMaxCheck) {
			pkd->nMaxCheck += (c->pUpper - c->pLower + 1) + 1000;
			pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			assert(pkd->Check != NULL);
			for (ism=0;ism<pkd->nMaxStack;++ism) {
			    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->S[ism].Check != NULL);
			}
			printf("Case 2: CPU:%d increased checklist size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		    }
		    for (pj=c->pLower;pj<=c->pUpper;++pj) {
			if (pkd->Check[i].id == pkd->idSelf) p = pkdParticle(pkd,pj);
			else p = mdlAcquire(pkd->mdl,CID_PARTICLE,pj,pkd->Check[i].id);
			/*
			** Only add those particle which we really need to check here!
			*/
			if (pkdIsGas(pkd,p) && pkdIsActive(pkd,p)) {
			    pkd->Check[nCheck] = pkd->Check[i];
			    pkd->Check[nCheck].iCell = -pj;
			    ++nCheck;
			}
			if (pkd->Check[i].id != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,p);
		    }
		    break;
		case 3:
		    /*
		    ** Open the cell.
		    ** Here we ASSUME that the children of
		    ** c are all in sequential memory order!
		    ** (the new tree build assures this)
		    ** (also true for the top tree)
		    ** We could do a prefetch here for non-local
		    ** cells.
		    */
		    if (nCheck + 2 > pkd->nMaxCheck) {
			pkd->nMaxCheck += 1000;
			pkd->Check = realloc(pkd->Check,pkd->nMaxCheck*sizeof(CELT));
			assert(pkd->Check != NULL);
			for (ism=0;ism<pkd->nMaxStack;++ism) {
			    pkd->S[ism].Check = realloc(pkd->S[ism].Check,pkd->nMaxCheck*sizeof(CELT));
			    assert(pkd->S[ism].Check != NULL);
			}
			printf("Case 3: CPU:%d increased checklist size to %d\n",mdlSelf(pkd->mdl),pkd->nMaxCheck);
		    }
		    iCheckCell = c->iLower;
		    pkd->Check[nCheck] = pkd->Check[i];
		    pkd->Check[nCheck+1] = pkd->Check[i];
		    /*
		    ** If we are opening a leaf of the top tree
		    ** we need to correctly set the processor id.
		    ** (this is a bit tricky)
		    */
		    if (pkd->Check[i].id < 0) {
			if (pkdTopNode(pkd,iCheckCell)->pLower >= 0) {
			    /* no need to search local trees again */
			    if ((iPid = pkdTopNode(pkd,iCheckCell)->pLower) != pkd->idSelf) {
				pkd->Check[nCheck].iCell = ROOT;
				pkd->Check[nCheck].id = iPid;
				++nCheck;
			    }
			}
			else {
			    pkd->Check[nCheck].iCell = iCheckCell;
			    assert(pkd->Check[nCheck].id == -1);
			    ++nCheck;
			}
			if (pkdTopNode(pkd,iCheckCell+1)->pLower >= 0) {
			    /* no need to search local trees again */
			    if ((iPid = pkdTopNode(pkd,iCheckCell+1)->pLower) != pkd->idSelf) {
				pkd->Check[nCheck].iCell = ROOT;
				pkd->Check[nCheck].id = iPid;
				++nCheck;
			    }
			}
			else {
			    pkd->Check[nCheck].iCell = iCheckCell+1;
			    assert(pkd->Check[nCheck].id == -1);
			    ++nCheck;
			}
		    }
		    else {
			pkd->Check[nCheck].iCell = iCheckCell;
			pkd->Check[nCheck+1].iCell = iCheckCell+1;
			assert(pkd->Check[nCheck].id == pkd->Check[i].id);
			assert(pkd->Check[nCheck+1].id == pkd->Check[i].id);
			nCheck += 2;
		    }
		    break;
		case 10:
		    /*
		    ** This checkcell is removed from the checklist since it has no overlap with the current cell.
		    */
		    break;		
		}
		if (pkd->Check[i].id >= 0 && pkd->Check[i].id != pkd->idSelf) mdlRelease(pkd->mdl,CID_CELL,c);
	    }
	    nCheck = ii;
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
	    if (!k->iLower) break;
	    k = pkdTreeNode(pkd,iCell = k->iLower);
	    /*
	    ** Check iCell has any active. It is actually active gas that we want to check here, but this
	    ** will have to do. All bounds checking is done against the active gas particles, so there are
	    ** only a few extra tests due to opening cells which may have actives, but no active gas particles.
	    */
	    if (k->nActive) {
		/*
		** iCell has actives, continue processing it.
		** Do NOT put the sibling onto the checklist, since we only want to check remote cells!
		*/
		/*
		** Test whether the sibling has inactives as well.
		** If not we don't push it onto the stack, but we
		** have to be careful to not pop the stack when we
		** hit the sibling. See the goto "InactiveAscend" below
		** for how this is done.
		*/
		kSib = pkdTreeNode(pkd,iCell+1);
		if (kSib->nActive) {
		    /*
		    ** Sibling has actives as well.
		    ** Push Checklist for the sibling onto the stack
		    ** before proceeding deeper in the tree.
		    */
		    ++iStack;
		    assert(iStack < pkd->nMaxStack);
		    pkd->S[iStack].nCheck = nCheck;
		    /*
		    ** Maybe use a memcpy here!
		    ** for (i=0;i<nCheck;++i) pkd->S[iStack].Check[i] = pkd->Check[i];
		    */
		    memcpy(pkd->S[iStack].Check,pkd->Check,nCheck*sizeof(CELT));
		    }
		}
	    else {
		/*
		** Skip iCell and do NOT add it to the Checklist.
		** No need to push anything onto the stack here.
		*/
		/*
		** Move onto processing the sibling.
		*/
		k = pkdTreeNode(pkd,++iCell);
		}
	    }
	/*
	** Now the interaction list should be complete and the
	** Checklist should be empty!
	*/
	assert(nCheck == 0);

	while (iCell & 1) {
	InactiveAscend:
	    k = pkdTreeNode(pkd,iCell = k->iParent);
	    if (!iCell) {
		/*
		** Make sure stack is empty.
		*/
		assert(iStack == -1);
		return;
		}
	    }
	k = pkdTreeNode(pkd,++iCell);
	if (!k->nActive) goto InactiveAscend;
	/*
	** Pop the Checklist from the top of the stack,
	** also getting the state of the interaction list.
	*/
	nCheck = pkd->S[iStack].nCheck;
	/*
	** Use a memcpy here. This is where we would win with lists since we could simply take the
	** pointer to this list. We would have to link the old checklist into the freelist.
	** for (i=0;i<nCheck;++i) Check[i] = S[iStack].Check[i];
	*/
	memcpy(pkd->Check,pkd->S[iStack].Check,nCheck*sizeof(CELT));
	--iStack;
	}
#endif
    }


void DoLocalSearch(SMX smx,SMF *smf,PARTICLE *p,double *rLast) {
    PKD pkd = smx->pkd;
    PQ *pq;
    double p_r[3], r[3];
    float fBall;
    int i,j,bDone,ix,iy,iz;

    pkdGetPos1(p->r,p_r);

    /*
    ** Correct distances and rebuild priority queue.
    */
    for (i=0;i<smx->nSmooth;++i) {
	smx->pq[i].dx += p_r[0]-rLast[0];
	smx->pq[i].dy += p_r[1]-rLast[1];
	smx->pq[i].dz += p_r[2]-rLast[2];
	smx->pq[i].fDist2 = pow(smx->pq[i].dx,2) + pow(smx->pq[i].dy,2) + 
	    pow(smx->pq[i].dz,2);
    }
    for (j=0;j<3;++j) rLast[j] = p_r[j];
    PQ_BUILD(smx->pq,smx->nSmooth,pq);
    pq = pqSearch(smx,pq,p_r);
    /*
    ** Search in replica boxes if it is required.
    */
    if (smx->bPeriodic) {
	for (ix=-1;ix<=1;++ix) {
	    r[0] = p_r[0] - ix*pkd->fPeriod[0];
	    for (iy=-1;iy<=1;++iy) {
		r[1] = p_r[1] - iy*pkd->fPeriod[1];
		for (iz=-1;iz<=1;++iz) {
		    r[2] = p_r[2] - iz*pkd->fPeriod[2];
		    if (ix || iy || iz) {
			pq = pqSearch(smx,pq,r);
		    }
		}
	    }
	}
    }
    fBall = sqrtf(pq->fDist2);
    pkdSetBall(pkd,p,fBall);
    /*
    ** Apply smooth funtion to the neighbor list.
    */
    smx->fcnSmooth(p,fBall,smx->nSmooth,smx->pq,smf);
    /*
    ** Call mdlCacheCheck to make sure we are making progress!
    */
    mdlCacheCheck(pkd->mdl);
}


/*
** New fast gas method uses only read-only caching and an efficient method to determine nearest
** neighbor lists for the the subsequent fast gas phase 2 code which calculates the SPH forces
** in a momentum conserving way. The new method is also extremely conducive to iterative SPH 
** calculations. All particles which interact with actives are found on the nearest neighbor list 
** of the active particles ONLY. If particles are mutual active nearest neighbors then only one
** of them is found on the interaction list of the other, EXCEPT in the case of a remote neighbor
** where both directions of the interaction must be calculated (with the advantage that we never 
** need the combiner cache of course). - Joachim Stadel
*/
void smFastGasPhase1(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p,*pp;
    FLOAT rLast[3];
    int pi,i,j,ii;
    uint32_t uHead,uTail;
    LIST *pList;
    int nMaxpList;
    int nList;
    int nInitList;
    char **ppCList;

    assert(pkd->oSph); /* Validate memory model */
    assert(pkd->oNodeSphBounds); /* Validate memory model */
    /*
    ** Initialize a default sized list. We will make this 2*nSmooth to start with.
    ** tList is a temporary list used to remove elements from a neighbor's own list.
    */
    nMaxpList = 2*pkd->param.nSmooth;
    pList = malloc(nMaxpList*sizeof(LIST));
    assert(pList != NULL);
    /*
    ** Initialize the bInactive and bDone flags for all local particles.
    ** Set DstActive to initially be all rung-active gas particles.
    */
    uHead = 0;
    uTail = 0;
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	p->bMarked = p->bSrcActive;
	smx->ea[pi].bDone = 0;
	if (pkdIsGas(pkd,p) && pkdIsActive(pkd,p)) {
	    /*
	    ** Clear the neighbor list of this active gas particle!
	    */
	    ppCList = pkd_pNeighborList(pkd,p);
	    if (*ppCList) {
		free(*ppCList);
		*ppCList = NULL;
	    }
	    /*
	    ** Place it on the do queue.
	    */
	    smx->ea[uTail++].iIndex = pi;
	}
    }
    if (uTail == 0) return;  /* no active particles??? */
    else if (uTail == pkd->nLocal) uTail = 0;  /* wrap uTail around in this case */
    smx->ea[pkd->nLocal].bDone = 1;  /* initialize for Sentinel, but this is not really needed */
    /*
    ** Initialize the priority queue first.
    */
    for (i=0;i<smx->nSmooth;++i) {
	smx->pq[i].pPart = smx->pSentinel;
	smx->pq[i].iIndex = pkd->nLocal;
	smx->pq[i].iPid = pkd->idSelf;
	smx->pq[i].dx = pkdPos(smx->pSentinel->r,0);
	smx->pq[i].dy = pkdPos(smx->pSentinel->r,1);
	smx->pq[i].dz = pkdPos(smx->pSentinel->r,2);
	smx->pq[i].fDist2 = pow(smx->pq[i].dx,2) + pow(smx->pq[i].dy,2) + 
	    pow(smx->pq[i].dz,2);
    }
    for (j=0;j<3;++j) rLast[j] = 0.0;

    do {
	/*
	** Remove leading element from the do queue.
	*/
	assert(uHead < pkd->nLocal);
	pi = smx->ea[uHead++].iIndex;  
	if (uHead == pkd->nLocal) uHead = 0;
	p = pkdParticle(pkd,pi);
	DoLocalSearch(smx,smf,p,rLast);
	/*
	** Mark this particle as done.
	*/
	smx->ea[pi].bDone = 1;
	if (pkdIsActive(pkd,p)) {
	    nList = 0;
	    /*
	    ** If this active particle already had inactive neighbors added to
	    ** its list then decode the list and start adding to it.
	    */
	    ppCList = pkd_pNeighborList(pkd,p);
	    if (*ppCList) {
		lcodeDecode(smx->lcmp,*ppCList,&pList,&nMaxpList,&nList);
		free(*ppCList);
		*ppCList = NULL;
	    }
	    nInitList = nList;
	    /*
	    ** Loop through the neighbors adding particles to the head of the 
	    ** do queue which are both rung INACTIVE, local, not done, and also DstActive.
	    */
	    for (i=0;i<smx->nSmooth;++i) {
		if (smx->pq[i].iPid == pkd->idSelf) {
		    /*
		    ** We skip the self interaction here, and do not add 
		    ** it to the compressed list.
		    */
		    if (smx->pq[i].iIndex == pi) continue;
		    pp = pkdParticle(pkd,smx->pq[i].iIndex);
		    if (pkdIsGas(pkd,pp)){
			if (!pkdIsActive(pkd,pp)) {
			    /* add this inactive particle to the list of p, but it may
			    ** already have been added IF p started with a list! */
			    for (ii=0;ii<nInitList;++ii) {
				if (pList[ii].iIndex == smx->pq[i].iIndex && 
				    pList[ii].iPid == smx->pq[i].iPid) goto DontAddToList;
			    }
			    if (nList == nMaxpList) {
				nMaxpList *= 2;
				pList = realloc(pList,nMaxpList*sizeof(LIST));
				assert(pList != NULL);
			    }
			    pList[nList].iIndex = smx->pq[i].iIndex;
			    pList[nList].iPid = smx->pq[i].iPid;
			    ++nList;
			DontAddToList:
			    if (!smx->ea[smx->pq[i].iIndex].bDone) {
				/*
				** Needs an updated density, so add it to the head of the 
				** do queue.
				*/
				if (uHead == 0) uHead = pkd->nLocal;
				smx->ea[--uHead].iIndex = smx->pq[i].iIndex;
			    }
			}
			else {
			    /*
			    ** pp is an active particle!
			    */
			    if (smx->ea[smx->pq[i].iIndex].bDone) {
				/*
				** If this (pi) particle was already added to the list of particle pj, then this
				** was not correct, and we should remove pi from pj's list if we care about the
				** ordering of the hsmooths among local particles (which we don't) as long as only
				** one of the 2 particles (both active) has the other on its list then all is good. 
				** If it is not in the list of pj, then we can simply add pj to the list of pi.
				*/
				ppCList = pkd_pNeighborList(pkd,pp);
/*
				printf("%p---------------check local list----------------- %d\n",*ppCList,pi);
				lcodeDecode(smx->lcmp,*ppCList,&tList,&nMaxtList,&ntList);
				printf("%d--",smx->pq[i].iIndex);
				lcodePrintList(tList,ntList);
*/		
				if (*ppCList) {
				    /*
				    ** pp certainly does have a list.
				    */
				    if (!bInListLocal(smx->lcmp,*ppCList,pi)) {
					if (nList == nMaxpList) {
					    nMaxpList *= 2;
					    pList = realloc(pList,nMaxpList*sizeof(LIST));
					    assert(pList != NULL);
					}
					pList[nList].iIndex = smx->pq[i].iIndex;
					pList[nList].iPid = smx->pq[i].iPid;
					++nList;
				    }
				}
			    }
			    else {
				/*
				** We can't be sure about this particle since it has not had its up-to-date
				** fBall calculated yet. We add it to p's particle list, but check later if this 
				** needs to be corrected. If pp finds p as a neighbor, then the test must be performed.
				*/
				if (nList == nMaxpList) {
				    nMaxpList *= 2;
				    pList = realloc(pList,nMaxpList*sizeof(LIST));
				    assert(pList != NULL);
				}
				pList[nList].iIndex = smx->pq[i].iIndex;
				pList[nList].iPid = smx->pq[i].iPid;
				++nList;
			    }
			}
		    } /* end of if (pkdIsGas(pkd,pp)) */
		    else {
			/*
			** We should never get here if the bSrcActive is set to 1 for all gas particles!
			*/
			assert(pkdIsGas(pkd,pp));
		    } /* end of !pkdIsGas(pkd,pp) */
		}
		else {
		    /*
		    ** It is a remote particle and needs to be done at any case.
		    */
		    if (nList == nMaxpList) {
			nMaxpList *= 2;
			pList = realloc(pList,nMaxpList*sizeof(LIST));
			assert(pList != NULL);
		    }
		    pList[nList].iIndex = smx->pq[i].iIndex;
		    pList[nList].iPid = smx->pq[i].iPid;
		    ++nList;
		}
	    }
	    /*
	    ** Compress the final list here as much as possible.
	    */
	    qsort(pList,nList,sizeof(LIST),lcodeCmpList);
	    ppCList = pkd_pNeighborList(pkd,p);
	    /*
	    ** Free the old list if it exists!
	    */
	    if (*ppCList) {
		free(*ppCList);
		*ppCList = NULL;
	    }
	    if (nList) {
		lcodeEncode(smx->lcmp,pList,nList,ppCList);
	    }
/*
	    printf("%p-%d--",*ppCList,pi);
	    lcodePrintList(pList,nList);
*/
	}
	else {
	    /*
	    ** We have to check if we need to add this inactive to any of the local active's
	    ** lists. The reason for this is that we mark this inactive as done, and this 
	    ** means it will no longer be found in the subsequent search over inactives.
	    */
	    for (i=0;i<smx->nSmooth;++i) {
		if (smx->pq[i].iPid == pkd->idSelf) {
		    /*
		    ** We skip the self interaction here, and do not add 
		    ** it to the compressed list.
		    */
		    if (smx->pq[i].iIndex == pi) continue;
		    pp = pkdParticle(pkd,smx->pq[i].iIndex);
		    if (pkdIsGas(pkd,pp) && pkdIsActive(pkd,pp)) {
			/*
			** pp is an active particle!
			*/
			ppCList = pkd_pNeighborList(pkd,pp);
			if (*ppCList == NULL) {
			    /*
			    ** Start a new compressed list for particle pp.
			    */
			    pList[0].iIndex = pi;
			    pList[0].iPid = pkd->idSelf;
			    /*
			    ** Compress the final list here as much as possible.
			    */
			    lcodeEncode(smx->lcmp,pList,1,ppCList);
			}
			else if (!bInListLocal(smx->lcmp,*ppCList,pi)) {
			    /*
			    ** Add this particle to the neighbor list of pp!
			    ** Unfortunately this means throwing away the old compressed list and creating a new one.
			    ** Adding to an already compressed list would be a useful function, but we still would
			    ** usually have to allocate/reallocate new storage for the compressed list of particle pp.
			    */
			    lcodeDecode(smx->lcmp,*ppCList,&pList,&nMaxpList,&nList);
			    free(*ppCList);
			    if (nList == nMaxpList) {
				nMaxpList *= 2;
				pList = realloc(pList,nMaxpList*sizeof(LIST));
				assert(pList != NULL);
			    }
			    pList[nList].iIndex = pi;
			    pList[nList].iPid = pkd->idSelf;
			    ++nList;
			    /*
			    ** Compress the final list here as much as possible.
			    */
			    qsort(pList,nList,sizeof(LIST),lcodeCmpList);
			    lcodeEncode(smx->lcmp,pList,nList,ppCList);
			}
		    } /* end of if (pkdIsGas(pkd,pp) && pkdIsActive(pkd,pp)) */
		}
	    }
	}	
    } while (uHead != uTail);
    /*
    ** Update all local sph bounds, making sure we only increase current minimums and decrease current
    ** maximums (only that is safe).
    */
    UpdateSphBounds(smx);
    /*
    ** Start of ball bound searching with the newly updated local bounds.
    ** First we make all local inactive which haven't been done search for active particles within their
    ** dHonHlimit extended ball bounds. We need to calculate new densities for these inactives and add them
    ** to the lists of the local actives (can't add them to the lists of remote actives - these actives must 
    ** find them on their own in the last pass if need be).
    */
    uHead = 0;
    uTail = BoundWalkInactive(smx);
    assert(uTail < pkd->nLocal);
/*
    printf("After BoundWalkInactive added %d new inactive particles. They get new densities now.\n",uTail);
*/
    /*
    ** New inactive particles have been added to the list.
    */
    while (uHead != uTail) {
	/*
	** Remove leading element from the do queue.
	*/
	assert(uHead < pkd->nLocal);
	pi = smx->ea[uHead++].iIndex;  
	if (uHead == pkd->nLocal) uHead = 0;
	p = pkdParticle(pkd,pi);
	DoLocalSearch(smx,smf,p,rLast);
	/*
	** Note that this particle is already marked as done by BoundWalkInactive!
	*/
	for (i=0;i<smx->nSmooth;++i) {
	    if (smx->pq[i].iPid == pkd->idSelf) {
		pp = pkdParticle(pkd,smx->pq[i].iIndex);
		if (pkdIsGas(pkd,pp) && pkdIsActive(pkd,pp)) {
		    /*
		    ** Add this particle to the neighbor list of pp!
		    ** Unfortunately this means throwing away the old compressed list and creating a new one.
		    ** Adding to an already compressed list would be a useful function, but we still would
		    ** usually have to allocate/reallocate new storage for the compressed list of particle pp.
		    */
		    ppCList = pkd_pNeighborList(pkd,pp);
		    assert(smx->ea[smx->pq[i].iIndex].bDone);
		    assert(*ppCList != NULL);
		    lcodeDecode(smx->lcmp,*ppCList,&pList,&nMaxpList,&nList);
		    free(*ppCList);
		    if (nList == nMaxpList) {
			nMaxpList *= 2;
			pList = realloc(pList,nMaxpList*sizeof(LIST));
			assert(pList != NULL);
		    }
		    pList[nList].iIndex = pi;
		    pList[nList].iPid = pkd->idSelf;
		    ++nList;
		    /*
		    ** Compress the final list here as much as possible.
		    */
		    qsort(pList,nList,sizeof(LIST),lcodeCmpList);
/*
		    lcodePrintList(pList,nList);
*/
		    lcodeEncode(smx->lcmp,pList,nList,ppCList);
		}
	    }
	}
    }
    /*
    ** Update all local sph bounds, making sure we only increase current minimums and decrease current
    ** maximums (only that is safe). We are really only interested in the ball bound for all particles (B) from here
    ** on now, although the completed inactives will again cause a shrink in BI ball bounds.
    */
    UpdateSphBounds(smx);
    /*
    ** In this last pass we make all local active search for REMOTE actives AND inactives which have not already
    ** been added to our lists. If they haven't, then add them, their densities will have been updated already.
    ** This is the only part where extra particles, other than true nearest neighbors at this step are added to 
    ** local lists. So only remote neighbors are a superset of the actual remote neighbors.
    */
    BoundWalkActive(smx,&pList,&nMaxpList);
    /*
    ** Release acquired pointers and source-reactivate particles in prioq.
    */
    for (i=0;i<smx->nSmooth;++i) {
	if (smx->pq[i].iPid == pkd->idSelf) {
	    pkdParticle(pkd,smx->pq[i].iIndex)->bMarked = 1;
	}
	else {
	    smHashDel(smx,smx->pq[i].pPart);
	    mdlRelease(pkd->mdl,CID_PARTICLE,smx->pq[i].pPart);
	}
    }    
    free(pList);
}


void smFastGasPhase2(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p,*pp;
    float pBall, ppBall;
    LIST *pList;
    double dx,dy,dz,fDist2;
    int nMaxpList,nList,i,nCnt,pi;
    char **ppCList;

    assert(pkd->oSph); /* Validate memory model */
    /*
    ** Initialize a default sized list. We will make this 2*nSmooth to start with.
    */
    nMaxpList = 2*pkd->param.nSmooth;
    pList = malloc(nMaxpList*sizeof(LIST));
    assert(pList != NULL);

    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if (pkdIsGas(pkd,p) && pkdIsActive(pkd,p)) {
	    ppCList = pkd_pNeighborList(pkd,p);
	    if (*ppCList) {
		lcodeDecode(smx->lcmp,*ppCList,&pList,&nMaxpList,&nList);
		if (nList > smx->nnListMax) {
		    smx->nnListMax = nList + NNLIST_INCREMENT;
		    smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
		    assert(smx->nnList != NULL);
		}
		nCnt = 0;
		for (i=0;i<nList;++i) {
		    if (pList[i].iPid == pkd->idSelf) {
			pp = pkdParticle(pkd,pList[i].iIndex);
		    }
		    else {
			pp = mdlAcquire(pkd->mdl,CID_PARTICLE,pList[i].iIndex,pList[i].iPid);
		    }
		    dx = pkdPos(p->r,0) - pkdPos(pp->r,0);
		    dy = pkdPos(p->r,1) - pkdPos(pp->r,1);
		    dz = pkdPos(p->r,2) - pkdPos(pp->r,2);
		    if (smx->bPeriodic) {
			/*
			** Correct for periodic boundaries.
			*/
			if (dx > 0.5*pkd->fPeriod[0]) dx -= pkd->fPeriod[0];
			else if (dx < -0.5*pkd->fPeriod[0]) dx += pkd->fPeriod[0];
			if (dy > 0.5*pkd->fPeriod[1]) dy -= pkd->fPeriod[1];
			else if (dy < -0.5*pkd->fPeriod[1]) dy += pkd->fPeriod[1];
			if (dz > 0.5*pkd->fPeriod[2]) dz -= pkd->fPeriod[2];
			else if (dz < -0.5*pkd->fPeriod[2]) dz += pkd->fPeriod[2];
		    }
		    fDist2 = dx*dx + dy*dy + dz*dz;
		    pBall= pkdBall(pkd,p);
		    ppBall= pkdBall(pkd,pp);
		    if (fDist2 < pBall*pBall || fDist2 < ppBall*ppBall) {
			smx->nnList[nCnt].fDist2 = fDist2;
			smx->nnList[nCnt].dx = dx;
			smx->nnList[nCnt].dy = dy;
			smx->nnList[nCnt].dz = dz;
			smx->nnList[nCnt].pPart = pp;
			smx->nnList[nCnt].iIndex = pList[i].iIndex;
			smx->nnList[nCnt].iPid = pList[i].iPid;
			++nCnt;
		    }
		    else {
			/*
			** This particle is not a neighbor in any sense now and can be
			** removed from the list. But for now we don't bother!
			*/
			if (pList[i].iPid != pkd->idSelf) mdlRelease(pkd->mdl,CID_PARTICLE,pp);
		    }
		} /* end of for (i=0;i<nList;++i) */
		/*
		** Apply smooth funtion to the neighbor list.
		*/
		smx->fcnSmooth(p,pkdBall(pkd,p),nCnt,smx->nnList,smf);
		/*
		** Release acquired pointers.
		*/
		for (i=0;i<nCnt;++i) {
		    if (smx->nnList[i].iPid != pkd->idSelf) {
			mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[i].pPart);
		    }
		}
		/*
		** Call mdlCacheCheck to make sure we are making progress!
		*/
		mdlCacheCheck(pkd->mdl);
	    } /* end of if (*ppCList) */
	} /* end of p is active... */
    }
}


void pkdFastGasCleanup(PKD pkd) {
    PARTICLE *p;
    uint32_t pi;
    char **ppCList;

    assert(pkd->oSph); /* Validate memory model */
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	ppCList = pkd_pNeighborList(pkd,p);
	if (*ppCList) {
	    free(*ppCList);
	    *ppCList = NULL;
	}
    }
}

void smGather(SMX smx,FLOAT fBall2,FLOAT r[3]) {
    KDN *kdn;
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    int idSelf = pkd->idSelf;
    struct stStack *S = smx->ST;
    FLOAT min2;
    int iCell,id;
    int sp = 0;
    const BND *bnd;
    PARTICLE *p;
    double p_r[3];
    FLOAT dx, dy, dz, fDist2;
    int pj, pEnd, nCnt;

    kdn = getCell(pkd,iCell=ROOT,id = 0);
    while (1) {
        bnd = pkdNodeBnd(pkd, kdn);
	MINDIST(bnd,r,min2);
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
	    if (id == pkd->idSelf) {
		pEnd = kdn->pUpper;
		for (pj=kdn->pLower;pj<=pEnd;++pj) {
		    p = pkdParticle(pkd,pj);
		    if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
		    pkdGetPos1(p->r,p_r);
		    dx = r[0] - p_r[0];
		    dy = r[1] - p_r[1];
		    dz = r[2] - p_r[2];
		    fDist2 = dx*dx + dy*dy + dz*dz;
		    if (fDist2 <= fBall2) {
			if (nCnt >= smx->nnListMax) {
			    smx->nnListMax += NNLIST_INCREMENT;
			    smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
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
		for (pj=kdn->pLower;pj<=pEnd;++pj) {
		    p = mdlFetch(mdl,CID_PARTICLE,pj,id);
		    if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
		    pkdGetPos1(p->r,p_r);
		    dx = r[0] - p_r[0];
		    dy = r[1] - p_r[1];
		    dz = r[2] - p_r[2];
		    fDist2 = dx*dx + dy*dy + dz*dz;
		    if (fDist2 <= fBall2) {
			if (nCnt >= smx->nnListMax) {
			    smx->nnListMax += NNLIST_INCREMENT;
			    smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
			    assert(smx->nnList != NULL);
			    }
			smx->nnList[nCnt].fDist2 = fDist2;
			smx->nnList[nCnt].dx = dx;
			smx->nnList[nCnt].dy = dy;
			smx->nnList[nCnt].dz = dz;
			smx->nnList[nCnt].pPart = mdlAcquire(mdl,CID_PARTICLE,pj,id);
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
}

void smDoGatherLocal(SMX smx,FLOAT fBall2,FLOAT r[3],void (*Do)(SMX,PARTICLE *,FLOAT)) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    double p_r[3];
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,pEnd;
    const BND *bnd;

    kdn = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
        bnd = pkdNodeBnd(pkd, kdn);
	MINDIST(bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	}
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iLower);
	    S[sp++] = iCell+1;
	    continue;
	}
	else {
	    pEnd = kdn->pUpper;
	    for (pj=kdn->pLower;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
		pkdGetPos1(p->r,p_r);
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
	if (sp) kdn = pkdTreeNode(pkd,iCell = S[--sp]);
	else break;
    }
}


void smReSmoothSingle(SMX smx,SMF *smf,PARTICLE *p,FLOAT fBall) {
    PKD pkd = smx->pkd;
    double R[3], r[3];
    int iStart[3],iEnd[3];
    int i,j;
    int ix,iy,iz;

    pkdGetPos1(p->r,R);

    smx->nnListSize = 0;
    /*
    ** Note for implementing SLIDING PATCH, the offsets for particles are
    ** negative here, reflecting the relative +ve offset of the simulation
    ** volume.
    */
    if (smx->bPeriodic) {
	for (j=0;j<3;++j) {
	    iStart[j] = d2i(floor((R[j] - fBall)/pkd->fPeriod[j] + 0.5));
	    iEnd[j] = d2i(floor((R[j] + fBall)/pkd->fPeriod[j] + 0.5));
	}
	for (ix=iStart[0];ix<=iEnd[0];++ix) {
	    r[0] = R[0] - ix*pkd->fPeriod[0];
	    for (iy=iStart[1];iy<=iEnd[1];++iy) {
		r[1] = R[1] - iy*pkd->fPeriod[1];
		for (iz=iStart[2];iz<=iEnd[2];++iz) {
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
    smx->fcnSmooth(p,fBall,smx->nnListSize,smx->nnList,smf);
    /*
    ** Release acquired pointers.
    */
    for (i=0;i<smx->nnListSize;++i) {
	if (smx->nnList[i].iPid != pkd->idSelf) {
	    mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[i].pPart);
	}
    }
}

void smReSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;

    smf->pfDensity = NULL;
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if ( pkdIsDstActive(p,0,MAX_RUNG) && pkdIsSrcActive(p,0,MAX_RUNG) )
	    smReSmoothSingle(smx,smf,p,pkdBall(pkd,p));
    }
}

FLOAT phase_dist(PKD pkd,double dvTau2,PARTICLE *pa,PARTICLE *pb,double H) {
    int j;
    FLOAT dx,dv,dx2,dv2;
    vel_t *va, *vb;

    assert(pkd->oGroup); /* Validate memory model */
    assert(pkd->oVelocity); /* Validate memory model */

    va = pkdVel(pkd,pa);
    vb = pkdVel(pkd,pb);

    dx2=0.0;
    for (j=0;j<3;++j) {
	dx = pkdPos(pa->r,j) - pkdPos(pb->r,j);
	dx2 += dx*dx;
    }
    dx2 /= pkdBall(pkd,pa);   /* this is actually fBall2! */
    dv2 = 0.0;
    if (dvTau2 > 0) {
	for (j=0;j<3;++j) {
	    dv = (va[j] - vb[j]) + H*(pkdPos(pa->r,j) - pkdPos(pb->r,j));
	    dv2 += dv*dv;
	}
	dv2 /= dvTau2;
    }
    return(dx2 + dv2);
}


FLOAT corrPos(FLOAT c,FLOAT r,FLOAT l) {

    FLOAT d;

    d = r-c;
    if (d >= 0.5*l) return r - l;
    else if (d < -0.5*l) return r + l;
    else return r;
}


FLOAT PutInBox(FLOAT r,FLOAT l) {
    if (r < -0.5*l) return r + l;
    else if (r >= 0.5*l) return r - l;
    else return r;
}

typedef struct {
    RB_NODE node;
    FOFRM   data;
} RM_NODE;

typedef struct protoGroup {
    int iId;
    int nMembers;
    RB_TREE treeRemoteMembers;
} FOFPG;

/*
**  Copy the given tree (recursively) to an array
*/
int copy_rm( FOFRM *rm,RB_NODE *node) {
    int iCount = 0;
    if ( node != NULL ) {
	RM_NODE *rmnode = (RM_NODE *)(node);
	iCount = copy_rm( rm, node->link[0] );
	rm += iCount;
	*rm++ = rmnode->data;
	iCount++;
	iCount += copy_rm( rm, node->link[1] );
    }
    return iCount;
}

int CmpRMs(void *ctx,const void *v1,const void *v2) {
    FOFRM *rm1 = (FOFRM *)v1;
    FOFRM *rm2 = (FOFRM *)v2;
    if (rm1->iPid != rm2->iPid) return (rm1->iPid - rm2->iPid);
    else return (rm1->iIndex - rm2->iIndex);
}

void smFof(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    PARTICLE *p;
    float *pPot;
    int32_t *pGroup;
    int32_t *pPartGroup;
    vel_t *v;
    RB_TYPE rm_type;
    FOFRM   rm_data;
    FOFPG* protoGroup;

    int pi,pn,pnn,nCnt,i,j,k;
    int nRmListSize,nRmCnt,iRmIndex;
    int nFifo, iHead, iTail, iMaxGroups, iGroup;
    int *Fifo;
    int iStart[3],iEnd[3];
    int ix,iy,iz;
    int idSelf = pkd->idSelf;

    FLOAT R[3],r[3],l[3],lx,ly,lz,fBall,fBall2Max,fvBall2;
    FLOAT fMass;
    int nTree,tmp;

    assert(pkd->oGroup); /* Validate memory model */
    assert(pkd->oVelocity); /* Validate memory model */
    assert(pkd->oPotential); /* Validate memory model */

    if (smx->bPeriodic) {
	lx = pkd->fPeriod[0];
	ly = pkd->fPeriod[1];
	lz = pkd->fPeriod[2];
    }
    else {
	lx = FLOAT_MAXVAL;
	ly = FLOAT_MAXVAL;
	lz = FLOAT_MAXVAL;
    }
    l[0] = lx;
    l[1] = ly;
    l[2] = lz;

    iHead = 0;
    iTail = 0;
    tmp = pkd->nDark+pkd->nGas+pkd->nStar;

    nTree = pkdTreeNode(pkd,ROOT)->pUpper + 1;
    iMaxGroups = nTree+1;
    nFifo = nTree;
    Fifo = (int *)malloc(nFifo*sizeof(int));
    assert(Fifo != NULL);

    /*used something smaller than FOFGD here to reduce memory usage*/
    protoGroup = (FOFPG *)malloc(iMaxGroups*sizeof(FOFPG));
    assert(protoGroup != NULL);

    /* This is the "empty" group */
    iGroup = 0;
    protoGroup[iGroup].nMembers = 0;
    protoGroup[iGroup].iId = iGroup;
    protoGroup[iGroup].treeRemoteMembers = NULL;

    nRmListSize = 0;
    rb_type_create(&rm_type,sizeof(RM_NODE),0,CmpRMs,0,0);

    pkd->nGroups = 0;
    pkd->nMaxRm = 0;
    fBall2Max = 0.0;

    /* set spatial linking length for each particle (mass dependend when bTauAbs=0) */
    for (pn=0;pn<nTree;pn++) {
	p = pkdParticle(pkd,pn);
	fMass = pkdMass(pkd,p);
	pGroup = pkdInt32(p,pkd->oGroup);
	*pGroup = 0; 
	if (smf->bTauAbs) {
	    pkdSetBall(pkd,p,smf->dTau2);
	    /* enforce a real space linking length smaller than the mean particle separation at all times :*/
	    if (smf->dTau2 > 0.2f*powf(fMass,2.0/3.0) )
		pkdSetBall(pkd,p,0.2f*powf(fMass,2.0/3.0));
	}
	else {
	    pkdSetBall(pkd,p,smf->dTau2*powf(fMass,2.0/3.0));
	}
	if (pkdBall(pkd,p) > fBall2Max) fBall2Max = pkdBall(pkd,p);
    }

    /* the velocity space linking length is the same for all particles. Zero means: do regular 3D FOF */
    fvBall2 = smf->dVTau2;
 
    /* Have to restart particle chache, since we will need the updated p->fBall now */
    mdlFinishCache(mdl,CID_PARTICLE);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));

    /* Starting FOF search now... */
    for (pn=0;pn<nTree;pn++) {
	p = pkdParticle(pkd,pn);
	if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;
	pGroup = pkdInt32(p,pkd->oGroup);
	if (*pGroup ) continue;
	iGroup++;
	assert(iGroup < iMaxGroups);
	protoGroup[iGroup].nMembers = 0;
	protoGroup[iGroup].iId = iGroup;
	protoGroup[iGroup].treeRemoteMembers = NULL;
	nRmCnt = 0;
	/*
	** Mark particle and add it to the do-fifo
	*/
	*pGroup = iGroup;
	Fifo[iTail] = pn; iTail++;
	if (iTail == nFifo) iTail = 0;
	while (iHead != iTail) {
	    pi = Fifo[iHead];iHead++;
	    p = pkdParticle(pkd,pi);
	    if (iHead == nFifo) iHead=0;
	    /*
	    ** Do a Ball Gather at the radius p->fBall
	    */
	    pkdGetPos1(p->r,R);
	    smx->nnListSize = 0;
	    fBall = sqrt(pkdBall(pkd,p));
	    if (smx->bPeriodic) {
		for (j=0;j<3;++j) {
		    iStart[j] = d2i(floor((R[j] - fBall)/pkd->fPeriod[j] + 0.5));
		    iEnd[j] = d2i(floor((R[j] + fBall)/pkd->fPeriod[j] + 0.5));
		}
		for (ix=iStart[0];ix<=iEnd[0];++ix) {
		    r[0] = R[0] - ix*pkd->fPeriod[0];
		    for (iy=iStart[1];iy<=iEnd[1];++iy) {
			r[1] = R[1] - iy*pkd->fPeriod[1];
			for (iz=iStart[2];iz<=iEnd[2];++iz) {
			    r[2] = R[2] - iz*pkd->fPeriod[2];
			    smGather(smx,pkdBall(pkd,p),r);
			}
		    }
		}
	    }
	    else {
		smGather(smx,pkdBall(pkd,p),R);
	    }
	    nCnt = smx->nnListSize;
	    for (pnn=0;pnn<nCnt;++pnn ) {
		if (smx->nnList[pnn].iPid == idSelf) { /* Local neighbors: */

		    /* Do not add particles which already are in a group*/
		    pPartGroup = pkdInt32(smx->nnList[pnn].pPart,pkd->oGroup);
		    if (*pPartGroup) continue;
		
		    if (fvBall2 > 0.0) {
			/* Check if this particle is close enough in velocity space  */
			if (phase_dist(pkd,fvBall2,p,smx->nnList[pnn].pPart,smf->H) > 1.0) continue;
		    }
		    /*
		    **  Mark particle and add it to the do-fifo
		    */
		    *pPartGroup = iGroup;
		    Fifo[iTail] = smx->nnList[pnn].iIndex;iTail++;
		    if (iTail == nFifo) iTail = 0;
	      
		} else {	 /* Nonlocal neighbors: */
		
		    /* Make remote member linking symmetric by using smaller linking length if different: */
		    if (fvBall2 > 0.0) { /* Check velocity space distance */
			if (phase_dist(pkd,fvBall2,p,smx->nnList[pnn].pPart,smf->H) > 1.0 ||
			    phase_dist(pkd,fvBall2,smx->nnList[pnn].pPart,p,smf->H) > 1.0) continue;
		    }
		    else { /* real space distance */
			if (smx->nnList[pnn].fDist2 > pkdBall(pkd,smx->nnList[pnn].pPart)) continue;
		    }
		
		    /* Add to remote member (RM) list if new */
		    rm_data.iIndex = smx->nnList[pnn].iIndex;
		    rm_data.iPid = smx->nnList[pnn].iPid;

		    if ( rb_insert(&rm_type,&protoGroup[iGroup].treeRemoteMembers,&rm_data) ) {
			nRmCnt++;
			nRmListSize++;
		    }
		}
	    }
	}
	if ( nRmCnt > pkd->nMaxRm ) pkd->nMaxRm = nRmCnt;
    }
    free(Fifo);
    
    /*
    ** Now we can already reject groups which are too small, if they are entirely local
    */
    for (pn=0; pn<nTree ; pn++) {
	p = pkdParticle(pkd,pn);
	pGroup = pkdInt32(p,pkd->oGroup);
	if (*pGroup >=0 && *pGroup < iMaxGroups)
	    ++(protoGroup[*pGroup].nMembers);
	else
	    printf("ERROR: idSelf=%i , p->pGroup=%i too large. iMaxGroups=%i \n",pkd->idSelf,*pGroup,iMaxGroups);
    }
    /*
    ** Create a remapping and give unique local Ids !
    */
    iMaxGroups = iGroup;
    iGroup= 1 + idSelf;
    pkd->nGroups = 0;
    protoGroup[0].iId = tmp;
    for (i=1;i<=iMaxGroups;i++) {
	protoGroup[i].iId = iGroup;
	if (protoGroup[i].nMembers < smf->nMinMembers && protoGroup[i].treeRemoteMembers == NULL) {
	    protoGroup[i].iId = tmp;
	}
	else {
	    iGroup += pkd->nThreads;
	    ++pkd->nGroups;
	}
    }
    /*
    ** Update the particle groups ids.
    */
    for (pi=0;pi<nTree;pi++) {
	p = pkdParticle(pkd,pi);
	pGroup = pkdInt32(p,pkd->oGroup);
	*pGroup = protoGroup[*pGroup].iId;
    }
    /*
    ** Allocate the remote members array
    */
    pkd->nRm = nRmListSize;
    pkd->remoteMember = mdlMalloc(mdl,(pkd->nRm+1)*sizeof(FOFRM));
    iRmIndex = 0;
    /*
    ** Allocate memory for group data
    */
    pkd->groupData = (FOFGD *) malloc((1+pkd->nGroups)*sizeof(FOFGD));
    assert(pkd->groupData != NULL);
    k=1;
    for (i=0;i<pkd->nGroups;i++) {
	while (protoGroup[k].iId == tmp) k++;
	pkd->groupData[i].iGlobalId = protoGroup[k].iId;
	pkd->groupData[i].iLocalId = protoGroup[k].iId;
	pkd->groupData[i].nLocal = protoGroup[k].nMembers;
	pkd->groupData[i].iFirstRm = iRmIndex;
	pkd->groupData[i].nRemoteMembers = copy_rm(pkd->remoteMember+iRmIndex,protoGroup[k].treeRemoteMembers);
	iRmIndex += pkd->groupData[i].nRemoteMembers;
	rb_free(&rm_type, &protoGroup[k].treeRemoteMembers);
	k++;
	pkd->groupData[i].bMyGroup = 1;
	pkd->groupData[i].fMass = 0.0;
	for (j=0;j<3;j++) {
	    pkd->groupData[i].rcom[j] = 0.0;
	    pkd->groupData[i].r[j] = 0.0;
	    pkd->groupData[i].v[j] = 0.0;
	}
	pkd->groupData[i].fRMSRadius = 0.0;
	pkd->groupData[i].potordenmax = -1.0;
    }
    free(protoGroup);
    /* Sanity check: the list size should match the number of elements copied */
    assert( iRmIndex == nRmListSize );
    /*
    ** Calculate local group properties
    */
    for (pi=0;pi<nTree;++pi) {
	p = pkdParticle(pkd,pi);
	pGroup = pkdInt32(p,pkd->oGroup);
	fMass = pkdMass(pkd,p);
	v = pkdVel(pkd,p);
	if (*pGroup != tmp) {
	    i = (*pGroup - 1 - pkd->idSelf)/pkd->nThreads;
	    for (j=0;j<3;j++) {
		if (pkd->groupData[i].fMass > 0.0)
		    r[j] = corrPos(pkd->groupData[i].rcom[j]/pkd->groupData[i].fMass,pkdPos(p->r,j),l[j]);
		else  r[j] = pkdPos(p->r,j);
		pkd->groupData[i].rcom[j] += r[j]*fMass;
		pkd->groupData[i].fRMSRadius +=  r[j]*r[j]*fMass;
		pkd->groupData[i].v[j] += v[j]*fMass;
	    }
	    pkd->groupData[i].fMass += fMass;
	    if(smf->iCenterType == 1){ /* maximum of fabs(potential) is stored in case 1 (centered on potential)*/
		pPot = pkdPot(pkd,p);
		if ( fabs(*pPot) > pkd->groupData[i].potordenmax) {
		    pkd->groupData[i].potordenmax = fabs(*pPot);
		    for (j=0;j<3;j++) pkd->groupData[i].r[j] = r[j];
		}
	    } else {
		if (pkdDensity(pkd,p) > pkd->groupData[i].potordenmax) {
		    pkd->groupData[i].potordenmax = pkdDensity(pkd,p);
		    for (j=0;j<3;j++) pkd->groupData[i].r[j] = r[j];
		}
	    }
	}
    }
}


int CmpGroups(const void *v1,const void *v2) {
    FOFGD *g1 = (FOFGD *)v1;
    FOFGD *g2 = (FOFGD *)v2;
    return g1->nTotal - g2->nTotal;
}

/*
 * This will expand the FIFO by doubling nFifo, reallocating the memory and
 * moving the head-to-end part to the end (and updating iHead)
 *                     tail     nFifo
 * Before: |XXXXXXXXXXX|XXXXXXXX|
 *         0           head
 *
 *                     tail                           nFifo
 * After:  |XXXXXXXXXXX|                     |XXXXXXXX|
 *                                           head
 */
static void expandFIFO( int *iHead, int *iTail, int *nFifo, FOFRM **rmFifo ) {
    int nMove;
    assert(*iHead == *iTail);
    nMove = *nFifo - *iHead;
    *nFifo *= 2;
    *rmFifo = (FOFRM *)realloc(*rmFifo,*nFifo*sizeof(FOFRM));
    assert(*rmFifo != NULL);
    memmove(*rmFifo + *nFifo - nMove, *rmFifo + *iHead, nMove*sizeof(FOFRM));
    *iHead = *nFifo - nMove;
    }


int smGroupMerge(SMF *smf,int bPeriodic) {
    PKD pkd = smf->pkd;
    MDL mdl = smf->pkd->mdl;
    PARTICLE *p;
    PARTICLE *pPart;
    int32_t *pGroup;
    int32_t *pPartGroup, iPartGroup;
    FLOAT l[3], rtmp,rcomsubgroup;
    int pi,id,i,j,k,index,listSize, sgListSize, lsgListSize;
    int nLSubGroups,nSubGroups,nMyGroups;
    int iHead, iTail, nFifo,tmp, nTree;

    FOFGD *sG;
    FOFRM *rmFifo;
    FOFRM *remoteRM;
    FOFRM rm;
    FOFGD **subGroup; /* Array of group data pointers */
    FOFGD **lSubGroup; /* Array of group data pointers */

    if (bPeriodic) {
	for (j=0;j<3;j++) l[j] = pkd->fPeriod[j];
    }
    else {
	for (j=0;j<3;j++) l[j] = FLOAT_MAXVAL;
    }

    tmp = pkd->nDark+pkd->nGas+pkd->nStar;
    nTree = pkdTreeNode(pkd,ROOT)->pUpper + 1;
    nFifo = 30*pkd->nMaxRm + 1;
    sgListSize = 10*pkd->nThreads;
    lsgListSize = 10*pkd->nThreads;

    subGroup = (FOFGD **)malloc(sgListSize*sizeof(FOFGD *));
    assert(subGroup != NULL);

    lSubGroup = (FOFGD **)malloc(lsgListSize*sizeof(FOFGD *));
    assert(lSubGroup != NULL);
    rmFifo = (FOFRM *)malloc(nFifo*sizeof(FOFRM));
    assert(rmFifo != NULL);

    /*
    ** Start RO particle cache.
    */
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
    /*
    ** Start CO group data cache.
    */
    /*printf( "Processor %d cache has %d entries\n", mdlSelf(mdl), pkd->nGroups );*/
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->groupData,sizeof(FOFGD), pkd->nGroups,pkd,initGroupMerge,combGroupMerge);
    /*
    ** Start RO remote member cache.
    */
    mdlROcache(mdl,CID_RM,NULL,pkd->remoteMember,sizeof(FOFRM),pkd->nRm);

    for (i=0; i < pkd->nGroups ;i++) {
	nSubGroups = 0;
	nLSubGroups = 0;
	iHead = 0;
	iTail = 0;
	pkd->groupData[i].nTotal = pkd->groupData[i].nLocal;
	if (pkd->groupData[i].bMyGroup == 0)	goto NextGroup;
	if (pkd->groupData[i].nRemoteMembers && pkd->groupData[i].bMyGroup) {

	    /* Add all remote members to the Fifo: */
	  for (j=pkd->groupData[i].iFirstRm; j < pkd->groupData[i].iFirstRm + pkd->groupData[i].nRemoteMembers ;j++){
	      rmFifo[iTail] = pkd->remoteMember[j];		
	      if (++iTail == nFifo) iTail = 0;
	      if (iTail == iHead) expandFIFO(&iHead,&iTail,&nFifo,&rmFifo);
	    }
	    while (iHead != iTail) {
		rm = rmFifo[iHead];
		if (++iHead == nFifo) iHead = 0;
		if (rm.iPid == pkd->idSelf) {
		    pPart = pkdParticle(pkd,rm.iIndex);
		    pPartGroup = pkdInt32(pPart,pkd->oGroup);
		    /* Local: Have I got this group already? If RM not in a group, ignore it */
		    if (*pPartGroup == pkd->groupData[i].iLocalId || *pPartGroup == tmp
			|| *pPartGroup == 0 ) goto NextRemoteMember;
		    for (k=0; k < nLSubGroups ;k++) {
			if (lSubGroup[k]->iLocalId == *pPartGroup) {
			    goto NextRemoteMember;
			}
		    }
		    /* Local: New subgroup found, add to list: */
		    sG = pkd->groupData + (*pPartGroup - 1 - pkd->idSelf)/pkd->nThreads;
		    if (nLSubGroups >= lsgListSize) {
			lsgListSize = d2i(lsgListSize * 1.5);
			lSubGroup = (FOFGD **)realloc(lSubGroup, lsgListSize*sizeof(FOFGD *));
			assert(lSubGroup != NULL);
		    }
		    lSubGroup[nLSubGroups++] = sG;
		    if (sG->potordenmax > pkd->groupData[i].potordenmax) {
			pkd->groupData[i].bMyGroup = 0;
			goto NextGroup;
		    }
		    if (sG->iLocalId < pkd->groupData[i].iGlobalId)
			pkd->groupData[i].iGlobalId = sG->iLocalId;
		    pkd->groupData[i].nTotal += sG->nLocal;
		    /* Add all its remote members to the Fifo: */
		    for (j=sG->iFirstRm; j< sG->iFirstRm + sG->nRemoteMembers ;j++) {
		      rmFifo[iTail] = pkd->remoteMember[j];			
		      if (++iTail == nFifo) iTail = 0;
		      if (iTail == iHead) expandFIFO(&iHead,&iTail,&nFifo,&rmFifo);
		    }
		}
		else {
		    pPart = mdlAcquire(mdl,CID_PARTICLE,rm.iIndex,rm.iPid);
		    iPartGroup = * pkdInt32(pPart,pkd->oGroup);
		    mdlRelease(mdl,CID_PARTICLE,pPart);
		    
		    /* Remote: ignore if not in a group */
		    if (iPartGroup == tmp) {
			goto NextRemoteMember;
		    }
		    /* Remote: Have I got this group already? */
		    for (k=0; k < nSubGroups ;++k) {
			if (iPartGroup == subGroup[k]->iLocalId) {
			    goto NextRemoteMember;
			}
		    }
		    /* Remote: New subgroup found, add to list: */
		    index = (iPartGroup - 1 - rm.iPid)/pkd->nThreads ;
		    sG = mdlAcquire(mdl,CID_GROUP,index,rm.iPid);
		    
		    if (nSubGroups >= sgListSize) {
			sgListSize = d2i(sgListSize * 1.5);
			subGroup = (FOFGD **)realloc(subGroup, sgListSize*sizeof(FOFGD *));
			assert(subGroup != NULL);
		    }
		    subGroup[nSubGroups++] = sG;
		    if (sG->potordenmax > pkd->groupData[i].potordenmax) {
			pkd->groupData[i].bMyGroup = 0;
			goto NextGroup;
		    }
		    if (sG->iLocalId < pkd->groupData[i].iGlobalId) 
			pkd->groupData[i].iGlobalId = sG->iLocalId;
		    pkd->groupData[i].nTotal += sG->nLocal;
		    /* Add all its remote members to the Fifo: */
		    for (j=sG->iFirstRm; j < sG->iFirstRm + sG->nRemoteMembers ;j++) {
			remoteRM = mdlAcquire(mdl,CID_RM,j,rm.iPid);
			rmFifo[iTail] = *remoteRM;
			if (++iTail == nFifo) iTail = 0;
			if (iTail == iHead) expandFIFO(&iHead,&iTail,&nFifo,&rmFifo);
			mdlRelease(mdl,CID_RM,remoteRM);
		    }
		}
	    NextRemoteMember:
		;
	    }
	    if (pkd->groupData[i].nTotal < smf->nMinMembers) {
		/*
		** Nonlocal group too small:
		*/
		pkd->groupData[i].iGlobalId = 0;
		pkd->groupData[i].bMyGroup = 0;
		for (k=0;k < nSubGroups;++k) {
		    subGroup[k]->iGlobalId = 0;
		    subGroup[k]->bMyGroup = 0;
		}
		for (k=0;k < nLSubGroups;++k) {
		    lSubGroup[k]->iGlobalId = 0;
		    lSubGroup[k]->bMyGroup = 0;
		}
	    }
	    else {
		/*
		** Nonlocal group big enough: calculate properties
		*/
		for (k=0;k<nSubGroups + nLSubGroups;++k) {
		    if(k < nSubGroups) sG = subGroup[k];
		    else sG = lSubGroup[k-nSubGroups];
		    sG->iGlobalId = pkd->groupData[i].iGlobalId;
		    sG->bMyGroup = 0;		    
		    for (j=0;j<3;j++) { /** including check for groups split across period bounds **/
		      rtmp = sG->rcom[j]/sG->fMass;      /* might be across period */
		      rcomsubgroup = corrPos(pkd->groupData[i].rcom[j]/pkd->groupData[i].fMass,rtmp,l[j]);      		
		      pkd->groupData[i].rcom[j] += rcomsubgroup*sG->fMass;
		      if (rcomsubgroup != rtmp){  /** shift across periodic boundary **/
			sG->fRMSRadius  += sG->fMass * (pow(rcomsubgroup,2)-pow(rtmp,2));
		      } 
		      pkd->groupData[i].v[j] += sG->v[j];
		    }
		    pkd->groupData[i].fRMSRadius += sG->fRMSRadius;  
		    pkd->groupData[i].fMass += sG->fMass;
		}
	    }
	NextGroup:
	    /*
	    ** Release non-local pointers.
	    */
	    for (k=0;k < nSubGroups;k++) {
		mdlRelease(mdl,CID_GROUP,subGroup[k]);
	    }
	}
    }
    mdlFinishCache(mdl,CID_PARTICLE);
    mdlFinishCache(mdl,CID_GROUP);
    mdlFinishCache(mdl,CID_RM);

    free(subGroup);
    free(lSubGroup);
    free(rmFifo);
    mdlFree(mdl,pkd->remoteMember);
    pkd->nRm = 0;

    for (pi=0;pi<nTree ;pi++){
	p = pkdParticle(pkd,pi);
	pGroup = pkdInt32(p,pkd->oGroup);
	index = (*pGroup - 1 - pkd->idSelf)/pkd->nThreads ;
	if(index >= 0 && index < pkd->nGroups )
	    *pGroup = pkd->groupData[index].iGlobalId;
	else
	    *pGroup = 0;
    }

    /*
    ** Here the embarasssingly parallel output of particle group ids,
    ** local density and group links were written. This was removed now.
    **
    ** Output of particle arrays with local density and group ids have to be added 
    ** elsewhere. Group link output is not needed, links can be generated from the
    ** group id arrays in postprocessing quite easily.
    **
    ** JD -- Feb 24, 2009
    */

    /* Move real groups to low memory and normalize their properties. */
    nMyGroups=0;
    for (i=0; i< pkd->nGroups;i++) {
	if (pkd->groupData[i].bMyGroup && pkd->groupData[i].iGlobalId != 0) {
	    for (j=0;j<3;j++)
		pkd->groupData[i].rcom[j] /= pkd->groupData[i].fMass;
	    /* 
	    ** Do not calculate fDeltaR2 with the corrected positions!
	    */
	    pkd->groupData[i].fRMSRadius /= pkd->groupData[i].fMass;
	    pkd->groupData[i].fRMSRadius -= pkd->groupData[i].rcom[0]*pkd->groupData[i].rcom[0]
		+ pkd->groupData[i].rcom[1]*pkd->groupData[i].rcom[1]
		+ pkd->groupData[i].rcom[2]*pkd->groupData[i].rcom[2];
	    pkd->groupData[i].fRMSRadius = sqrt(pkd->groupData[i].fRMSRadius);
	    /* 
	    ** Now put all the positions back into the box and normalise the rest
	    */
	    for (j=0;j<3;j++) {
		pkd->groupData[i].rcom[j] = PutInBox(pkd->groupData[i].rcom[j],l[j]);
		pkd->groupData[i].r[j] = PutInBox(pkd->groupData[i].r[j],l[j]);
		pkd->groupData[i].v[j] /= pkd->groupData[i].fMass;
	    }
	    pkd->groupData[nMyGroups] = pkd->groupData[i];
	    nMyGroups++;
	}
    }
    if (nMyGroups == pkd->nGroups) {
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,(nMyGroups+1)*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
    }
    pkd->groupData[nMyGroups].bMyGroup = 0;

    /* Below the master node collects all group information, to write it out later.
    **
    ** Ideally this would be replaced by a scheme where every node writes out the data
    ** of his own groups, resp. sends them to his I/O node to do so for him.
    **
    ** But since the amount of data is rather small (e.g. it gave ascii files of around 50 MB for VL-2)
    ** the current serial output of group data is not expected to become a bottleneck anytime soon.
    */

    /* Start RO group data cache and master reads and saves all the group data. */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->groupData,sizeof(FOFGD), nMyGroups + 1);
    if (pkd->idSelf == 0) {
	listSize = pkd->nThreads*(pkd->nGroups+1);
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,listSize*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
	for (id=1; id < pkd->nThreads; id++) {
	    index = 0;
	    while (1) {
		sG = mdlAcquire(mdl,CID_GROUP,index,id);
		mdlRelease(mdl,CID_GROUP,sG);
		if (sG->bMyGroup != 0) {
		    if (nMyGroups >= listSize-1) {
			listSize *= 2;
			pkd->groupData = (FOFGD *) realloc(pkd->groupData, listSize*sizeof(FOFGD));
			assert(pkd->groupData != NULL);
		    }
		    pkd->groupData[nMyGroups] = *sG;
		    index++;
		    nMyGroups++;
		}
		else {
		    break;
		}
	    }
	}
	pkd->groupData[nMyGroups].bMyGroup = 0;
    }
    mdlFinishCache(mdl,CID_GROUP);
    if (pkd->idSelf != 0) {
	nMyGroups = 0;
    }
    else {
	/*
	** Master orders groupData
	*/
	qsort(pkd->groupData,nMyGroups,sizeof(FOFGD), CmpGroups);
    }
    pkd->nGroups = nMyGroups;
    return nMyGroups;
}


#if 0
void DoBins(SMX smx,PARTICLE *p,FLOAT fDist2) {
    FOFGD *pCurGrp = smx->pCurrentGroup;
    FOFBIN *pBin;
    int iStart,nBins,pid;

    iStart = pCurGrp->first.FirstBin;
    nBins = pCurGrp->num.nBin;
    pid = pCurGrp->pid;

    /* JS: calc this in an arithmetic way. */
    k = 0;
    while (fDist2 > pkd->groupBin[iBin+k].fRadius*pkd->groupBin[iBin+k].fRadius) {
	k++;
	if ( k == smf->nBins) goto nextParticle;
    }

    

    fMass = pkdMass(pkd,p);
    if (pid != smx->pkd->idSelf) {
	pBin = mdlAcquire(smx->pkd->mdl,CID_BIN,iBin,pid);
    } else {
	pBin = pkd->groupBin[iBin];
    }
    pBin->nMembers++;
    pBin->fMassInBin += fMass;
    if (pid != smx->pkd->idSelf) mdlRelease(smx->pkd->mdl,CID_BIN,pBin);
}
#endif

int smGroupProfiles(SMX smx, SMF *smf, int nTotalGroups) {
    assert(0);
    return 0;
#if 0
    PKD pkd = smf->pkd;
    MDL mdl = smf->pkd->mdl;
    PARTICLE *p;
    int32_t *pBin, *pPartBin;
    double *v;
    double dx2;
    FLOAT l[3],L[3],r[3],relvel[3],com[3];
    FLOAT rvir,Mvir,fBall,lastbin,fMass,fAvgDens;
    int pn,i,j,k,iBin,nBins,nTree,index,nCnt,pnn;
    int iStart[3],iEnd[3];
    int ix,iy,iz;
    FOFGD *gdp;
    FOFBIN *bin;

    if (nTotalGroups==0) return 0;

    assert(pkd->oGroup); /* Validate memory model */
    assert(pkd->oVelocity); /* Validate memory model */
    if (smx->bPeriodic) {
	for (j=0;j<3;j++) l[j] = pkd->fPeriod[j];
    }
    else {
	for (j=0;j<3;j++) l[j] = FLOAT_MAXVAL;
    }
    nTree = pkdTreeNode(pkd,ROOT)->pUpper + 1;
    /*
    ** Start RO group data cache and read all if you are not master
    **
    ** Every node needs to know all group positions and properties to set up the bins. 
    ** Only konwing the properties of ones own groups would not be enough,
    ** because often local particles contribute mass to a bin around
    ** a group which belonged to another node.
    **
    ** This is another reason to keep having the master collect everything, so it can be
    ** distributed to everybody here. 
    */
    if (pkd->idSelf != 0) {
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,nTotalGroups*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
    }
    mdlROcache(mdl,CID_GROUP,NULL,pkd->groupData,sizeof(FOFGD),pkd->nGroups);
    if (pkd->idSelf != 0) {
	for (i=0; i< nTotalGroups; i++) {
	    gdp = mdlAcquire(mdl,CID_GROUP,i,0);
	    mdlRelease(mdl,CID_GROUP,gdp);
	    pkd->groupData[i] = *gdp;
	}
    }
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Allocate memory for the bins.
    */
    nBins = nTotalGroups*smf->nBins;
    if ( pkd->groupBin != NULL ) free(pkd->groupBin);
    pkd->groupBin = (FOFBIN *) malloc( (nBins)*sizeof(FOFBIN) );
    assert(pkd->groupBin != NULL);
    /*
    ** Initalize bin array
    */
    iBin = 0;
    for (i=0; i< nTotalGroups; i++) {
	if (smf->bLogBins == 2) { /* logarithmic bins with same fixed non-comoving size for all groups */
	    lastbin = smf->binFactor/smf->a; /* NB: lastbin is actually the first bin in this case */
	}
	else {
	    /* estimate virial radius, assuming isothermal shperes */
	    fAvgDens = 0.5*pkd->groupData[i].fMass
		*0.238732414/pow(pkd->groupData[i].fRMSRadius,3.0); /*using half mass radius and assuming spherical halos*/
	    lastbin = pow(fAvgDens,0.5)*pkd->groupData[i].fRMSRadius*smf->binFactor;
	}
	for (j=0; j < smf->nBins; j++) {
	    pkd->groupBin[iBin].nMembers = 0;
	    if (smf->bLogBins == 1) {/* logarithmic bins */
		dx2 = smf->nBins-(j+1);
		pkd->groupBin[iBin].fRadius = smf->fMinRadius*pow(lastbin/smf->fMinRadius,((float) j)/( (float)smf->nBins-1.0));
	    }
	    else if (smf->bLogBins == 2) {/* logarithmic bins with same fixed non-comoving size for all groups */
		pkd->groupBin[iBin].fRadius = lastbin*pow(10.0, 0.1*j);
	    }
	    else { /* linear bins */
		dx2 = j+1;
		pkd->groupBin[iBin].fRadius = lastbin*dx2/smf->nBins;
	    }
	    pkd->groupBin[iBin].fMassInBin = 0.0;
	    iBin++;
	}
    }
    /*
    ** Add local particles to their corresponding bins
    */
    for (index=0; index< nTotalGroups; index++) {
	k = index*smf->nBins + smf->nBins-1;
	fBall = pkd->groupBin[k].fRadius;
      
	for (j = 0; j < 3; j++) {
	    if(smf->iCenterType == 0)
		com[j] = pkd->groupData[index].rcom[j];
	    else
		com[j] = pkd->groupData[index].r[j];
	}
	smx->nnListSize = 0;
	if (smx->bPeriodic) {
	    for (j=0;j<3;++j) {
		iStart[j] = floor((com[j] - fBall)/pkd->fPeriod[j] + 0.5);
		iEnd[j] = floor((com[j] + fBall)/pkd->fPeriod[j] + 0.5);
	    }
	    for (ix=iStart[0];ix<=iEnd[0];++ix) {
		r[0] = com[0] - ix*pkd->fPeriod[0];
		for (iy=iStart[1];iy<=iEnd[1];++iy) {
		    r[1] = com[1] - iy*pkd->fPeriod[1];
		    for (iz=iStart[2];iz<=iEnd[2];++iz) {
			r[2] = com[2] - iz*pkd->fPeriod[2];
			smGatherLocal(smx,fBall*fBall,r,DoBins);
		    }
		}
	    }
	}
	else {
	    smGatherLocal(smx,fBall*fBall,com,DoBins);
	}
    }
    /*
    ** Start CO group profiles cache.
    */
    mdlCOcache(mdl,CID_BIN,NULL,pkd->groupBin,sizeof(FOFBIN),nBins,pkd,initGroupBins,combGroupBins);
    if (pkd->idSelf != 0) {
	for (i=0; i< nBins; i++) {
	    if (pkd->groupBin[i].fMassInBin > 0.0) {
		bin = mdlAcquire(mdl,CID_BIN,i,0);
		*bin = pkd->groupBin[i];
		mdlRelease(mdl,CID_BIN,bin);
	    }
	}
    }
    mdlFinishCache(mdl,CID_BIN);
    if (pkd->idSelf != 0) {
	free(pkd->groupData);
	free(pkd->groupBin);
    }
    pkd->nBins =  nBins;
    return nBins;
#endif
}


