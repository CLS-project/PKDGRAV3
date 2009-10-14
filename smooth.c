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
#include "smoothfcn.h"
#include "rbtree.h"
#include <sys/stat.h>

const char *smooth_module_id = "$Id$";
const char *smooth_h_module_id = SMOOTH_H_MODULE_ID;


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


int smInitialize(SMX *psmx,PKD pkd,SMF *smf,int nSmooth,int bPeriodic,int bSymmetric,int iSmoothType) {
    SMX smx;
    void (*initParticle)(void *,void *) = NULL;
    void (*init)(void *,void *) = NULL;
    void (*comb)(void *,void *,void *) = NULL;
    int pi,j;
    int nTree;
    int iTopDepth;

    smx = malloc(sizeof(struct smContext));
    assert(smx != NULL);
    smx->pkd = pkd;
    if (smf != NULL) smf->pkd = pkd;
#ifdef BADSMOOTH
    nSmooth=nSmooth*3;
#endif
    smx->nSmooth = nSmooth;
    smx->bPeriodic = bPeriodic;

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
    for (i=0;i<1000,++i) {
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
	smx->pHash[i].coll = smx->pHash[i+1].coll;
    }
    smx->pHash[i].p = NULL;
    smx->pHash[i].coll = NULL;
    /*
    ** Allocate special stacks for searching.
    ** There is a mistake here, since I use these stacks for the remote trees as well.
    ** This can be easily fixed, but a hack for now.
    */
    smx->S = malloc(1024*sizeof(int));
    assert(smx->S != NULL);
    smx->Smin = malloc(1024*sizeof(FLOAT));
    assert(smx->Smin != NULL);
    /*
    ** Allocate special stacks for searching within the top tree.
    ** Calculate the number of levels in the top tree.
    */
    iTopDepth = 1+(int)ceil(log((double)smx->pkd->nThreads)/log(2.0));
    smx->ST = malloc(iTopDepth*sizeof(int));
    assert(smx->ST != NULL);
    smx->SminT = malloc(iTopDepth*sizeof(FLOAT));
    assert(smx->SminT != NULL);
    /*
    ** Set up the sentinel particle with some very far away distance.
    ** This is used to initially load the priority queue and all references
    ** to this particle should end up being replaced in the priority queue
    ** as long as there are nSmooth particles set bSrcActive=1.
    */
    for (j=0;j<3;++j) {
	smx->pSentinel.r[j] = HUGE_VAL;
    }
    smx->pSentinel.bSrcActive = 1;
    smx->pSentinel.bDstActive = 0;
    /*
    ** Need to cast the pLite to an array for char for local flags.
    */
    smx->bInactive = UNION_CAST(pkd->pLite,PLITE *,char *);
   
    *psmx = smx;
    return(1);
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
    sprintf(achOut, "    Min ratio: %g\n",
	    mdlMinRatio(smx->pkd->mdl,CID_CELL));
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
    sprintf(achOut, "    Min ratio: %g\n",
	    mdlMinRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    sprintf(achOut, "    Coll ratio: %g\n",
	    mdlCollRatio(smx->pkd->mdl,CID_PARTICLE));
    mdlDiag(smx->pkd->mdl, achOut);
    /*
    ** Stop particle caching space.
    */
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
    ** Free up context storage.
    */
    free(smx->S);
    free(smx->Smin);
    free(smx->ST);
    free(smx->SminT);
    free(smx->pq);
    free(smx->nnList);
    free(smx->pHash);
    free(smx);
}


/*
** This function performs a local nearest neighbor search.
*/
PQ *pqSearchLocal(SMX smx,PQ *pq,FLOAT r[3],int *pbDone) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    FLOAT dx,dy,dz,dMin,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int i,j,pj,pEnd,iCell,iParent;
    int sp = 0;
    int sm = 0;
    int idSelf = pkd->idSelf;

    *pbDone = 1;	/* assume that we will complete the search */
    /*
    ** We don't perform containment tests except at the
    ** root, so that the pbDone flag can be correctly
    ** set.
    */
    kdn = pkdTreeNode(pkd,iCell = ROOT);
    S[sp] = iCell;
    /*
    ** Start of PRIOQ searching loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iLower);
	    MINDIST(kdn->bnd,r,min1);
	    kdn = pkdTreeNode(pkd,++iCell);
	    MINDIST(kdn->bnd,r,min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		kdn = pkdTreeNode(pkd,--iCell);
		if (min1 >= pq->fDist2) goto NotContained;
	    }
	    else {
		Smin[sm++] = min1;
		if (min2 >= pq->fDist2) goto NotContained;
	    }
	}
	pEnd = kdn->pUpper;
	for (pj=kdn->pLower;pj<=pEnd;++pj) {
	    if (smx->bInactive[pj]) continue;
	    p = pkdParticle(pkd,pj);
	    dx = r[0] - p->r[0];
	    dy = r[1] - p->r[1];
	    dz = r[2] - p->r[2];
	    fDist2 = dx*dx + dy*dy + dz*dz;
	    if (fDist2 < pq->fDist2) {
		if (pq->iPid == idSelf) {
		    smx->bInactive[pq->iIndex] = 0;
		} 
		else {
		    smHashDel(smx,pq->pPart);
		    mdlRelease(pkd->mdl,CID_PARTICLE,pq->pPart);
		    pq->iPid = idSelf;
		}
		pq->pPart = p;
		pq->fDist2 = fDist2;
		pq->dx = dx;
		pq->dy = dy;
		pq->dz = dz;
		pq->iIndex = pj;
		smx->nInactive[pj] = 1; /* de-activate a particle that enters the queue */
		PQ_REPLACE(pq);
	    }
	}
    NoIntersect:
	while (iCell == S[sp]) {
	    if (sp) {
		--sp;
		kdn = pkdTreeNode(pkd,iCell = kdn->iParent);
	    }
	    else {
		/*
		** Containment Test!
		*/
		for (j=0;j<3;++j) {
		    dMin = kdn->bnd.fMax[j] -
			fabs(kdn->bnd.fCenter[j] - r[j]);
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = kdn->iParent;
			if (!iParent) {
			    *pbDone = 0;		/* EXIT, not contained! */
			    break;
			}
			S[sp] = iParent;
			goto NotContained;
		    }
		}
		return pq;
	    }
	}
    NotContained:
	kdn = pkdTreeNode(pkd,iCell ^= 1);
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
	    MINDIST(kdn->bnd,r,min2);
	}
	if (min2 >= pq->fDist2) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iParent);
	    goto NoIntersect;
	}
	S[++sp] = iCell;
    }
}



PQ *pqSearchRemote(SMX smx,PQ *pq,int id,FLOAT r[3]) {
    PKD pkd = smx->pkd;
    MDL mdl = smx->pkd->mdl;
    KDN *kdn;
    PARTICLE *p;
    KDN *pkdn,*pkdu;
    FLOAT dx,dy,dz,min1,min2,fDist2;
    FLOAT *Smin = smx->Smin;
    int *S = smx->S;
    int pj,pEnd,iCell;
    int sp = 0;
    int sm = 0;
    int idSelf = pkd->idSelf;

    assert(id != idSelf);
    kdn = pkdTreeNode(pkd,iCell = ROOT);
    S[sp] = iCell;
    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    /*
    ** Start of PRIOQ searching loop.
    */
    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (pkdn->iLower) {
	    kdn = pkdTreeNode(pkd,iCell = pkdn->iLower);
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    MINDIST(pkdn->bnd,r,min1);
	    kdn = pkdTreeNode(pkd,++iCell);
	    pkdu = mdlAquire(mdl,CID_CELL,iCell,id);
	    MINDIST(pkdu->bnd,r,min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		kdn = pkdTreeNode(pkd,--iCell);
		mdlRelease(mdl,CID_CELL,pkdu);
		if (min1 >= pq->fDist2) goto NotContained;
	    }
	    else {
		Smin[sm++] = min1;
		mdlRelease(mdl,CID_CELL,pkdn);
		pkdn = pkdu;
		if (min2 >= pq->fDist2) goto NotContained;
	    }
	}
	pEnd = pkdn->pUpper;
	for (pj=pkdn->pLower;pj<=pEnd;++pj) {
	    p = mdlAquire(mdl,CID_PARTICLE,pj,id);
	    if (smHashPresent(p)) continue;
	    if (!p->bSrcActive) {
		mdlRelease(mdl,CID_PARTICLE,p);
		continue;
	    }
	    dx = r[0] - p->r[0];
	    dy = r[1] - p->r[1];
	    dz = r[2] - p->r[2];
	    fDist2 = dx*dx + dy*dy + dz*dz;
	    if (fDist2 < pq->fDist2) {
		if (pq->iPid == idSelf) {
		    smx->bInactive[pq->iIndex] = 0;
		}
		else {
		    smHashDel(smx,pq->pPart);
		    mdlRelease(mdl,CID_PARTICLE,pq->pPart);
		}
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
	    else mdlRelease(mdl,CID_PARTICLE,p);
	}
    NoIntersect:
	while (iCell == S[sp]) {
	    if (!sp) {
		mdlRelease(mdl,CID_CELL,pkdn);
		return pq;
	    }
	    --sp;
	    kdn = pkdTreeNode(pkd,iCell = pkdn->iParent);
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	}
    NotContained:
	kdn = pkdTreeNode(pkd,iCell ^= 1);
	mdlRelease(mdl,CID_CELL,pkdn);
	pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
	    MINDIST(pkdn->bnd,r,min2);
	}
	if (min2 >= pq->fDist2) {
	    kdn = pkdTreeNode(pkd,iCell = pkdn->iParent);
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    goto NoIntersect;
	}
	S[++sp] = iCell;
    }
}


PQ *pqSearch(SMX smx,PQ *pq,FLOAT r[3],int bReplica,int *pbDone) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    int idSelf = smx->pkd->idSelf;
    FLOAT *Smin = smx->SminT;
    int *S = smx->ST;
    FLOAT dMin,min1,min2;
    int j,iCell,id,iParent;
    int sp = 0;
    int sm = 0;

    *pbDone = 0;
    if (bReplica) kdn = pkdTopNode(pkd,iCell = ROOT);
    else {
	kdn = pkdTopNode(pkd,iCell = pkd->iTopRoot);
	assert(kdn->pLower == idSelf);
    }
    if (iCell != ROOT) S[sp] = kdn->iParent;
    else S[sp] = iCell;

    while (1) {
	/*
	** Descend to bucket via the closest cell at each level.
	*/
	while (kdn->iLower) {
	    kdn = pkdTopNode(pkd,iCell = kdn->iLower);
	    MINDIST(kdn->bnd,r,min1);
	    kdn = pkdTopNode(pkd,++iCell);
	    MINDIST(kdn->bnd,r,min2);
	    if (min1 < min2) {
		Smin[sm++] = min2;
		kdn = pkdTopNode(pkd,--iCell);
		if (min1 >= pq->fDist2) goto NotContained;
	    }
	    else {
		Smin[sm++] = min1;
		if (min2 >= pq->fDist2) goto NotContained;
	    }
	}
	id = kdn->pLower;	/* this is the thread id in LTT */
	if (id == pkd->idSelf) {
	    pq = pqSearchLocal(smx,pq,r,pbDone);
	    if (*pbDone) return(pq);
	}
	else {
	    pq = pqSearchRemote(smx,pq,id,r);
	}
    NoIntersect:
	while (iCell == S[sp]) {
	    if (sp) {
		--sp;
		kdn = pkdTopNode(pkd,iCell = kdn->iParent);
	    }
	    else if (!bReplica) {
		/*
		** Containment Test!
		*/
		for (j=0;j<3;++j) {
		    dMin = kdn->bnd.fMax[j] -
			fabs(kdn->bnd.fCenter[j] - r[j]);
		    if (dMin*dMin < pq->fDist2 || dMin < 0) {
			iParent = kdn->iParent;
			if (!iParent) {
			    *pbDone = 0;
			    return pq;
			}
			S[sp] = iParent;
			goto NotContained;
		    }
		}
		*pbDone = 1;
		return pq;
	    }
	    else return pq;
	}
    NotContained:
	kdn = pkdTopNode(pkd,iCell ^= 1);
	/*
	** Intersection Test. (ball-test)
	*/
	if (sm) min2 = Smin[--sm];
	else {
	    MINDIST(kdn->bnd,r,min2);
	}
	if (min2 >= pq->fDist2) {
	    kdn = pkdTopNode(pkd,iCell = kdn->iParent);
	    goto NoIntersect;
	}
	S[++sp] = iCell;
    }
}


void smSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    PQ *pq = smx->pq;
    FLOAT r[3],fBall;
    FLOAT rLast[3];
    int iStart[3],iEnd[3];
    int pi,i,j,bDone;
    int ix,iy,iz;

    /*
    ** Initialize the bInactive flags for all local particles.
    */
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	smx->bInactive[pi] = (p->bSrcActive)?0:1;
    }
    smx->nInactive[pkd->nLocal] = 0;  /* initialize for Sentinel, but this is not really needed */
    /*
    ** Initialize the priority queue first.
    */
    for (i=0;i<smx->nSmooth;++i) {
	smx->pq[i].pPart = &smx->pSentinel;
	smx->pq[i].iIndex = pkd->nLocal;
	smx->pq[i].iPid = pkd->idSelf;
	smx->pq[i].dx = smx->pSentinel.r[0];
	smx->pq[i].dy = smx->pSentinel.r[1];
	smx->pq[i].dz = smx->pSentinel.r[2];
	smx->pq[i].fDist2 = pow(smx->pq[i].dx,2) + pow(smx->pq[i].dy,2) + 
	    pow(smx->pq[i].dz,2);
    }
    for (j=0;j<3;++j) rLast[j] = 0.0;

    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;
	/*
	** Correct distances and rebuild priority queue.
	*/
	for (i=0;i<smx->nSmooth;++i) {
	    smx->pq[i].dx += p->r[0]-rLast[0];
	    smx->pq[i].dy += p->r[1]-rLast[1];
	    smx->pq[i].dz += p->r[2]-rLast[2];
	    smx->pq[i].fDist2 = pow(smx->pq[i].dx,2) + pow(smx->pq[i].dy,2) + 
		pow(smx->pq[i].dz,2);
	}
	for (j=0;j<3;++j) rLast[j] = p->r[j];
	PQ_BUILD(smx->pq,smx->nSmooth,pq);

	pq = pqSearch(smx,pq,p->r,0,&bDone);
	/*
	** Search in replica boxes if it is required.
	*/
	if (!bDone && smx->bPeriodic) {
	    fBall = sqrt(pq->fDist2);
	    for (j=0;j<3;++j) {
		iStart[j] = floor((p->r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		iEnd[j] = floor((p->r[j] + fBall)/pkd->fPeriod[j] + 0.5);
	    }
	    for (ix=iStart[0];ix<=iEnd[0];++ix) {
		r[0] = p->r[0] - ix*pkd->fPeriod[0];
		for (iy=iStart[1];iy<=iEnd[1];++iy) {
		    r[1] = p->r[1] - iy*pkd->fPeriod[1];
		    for (iz=iStart[2];iz<=iEnd[2];++iz) {
			r[2] = p->r[2] - iz*pkd->fPeriod[2];
			if (ix || iy || iz) {
			    pq = pqSearch(smx,pq,r,1,&bDone);
			}
		    }
		}
	    }
	}
	p->fBall = sqrt(pq->fDist2);
	/*
	** Apply smooth funtion to the neighbor list.
	*/
	smx->fcnSmooth(p,smx->nSmooth,smx->pq,smf);
	/*
	** Call mdlCacheCheck to make sure we are making progress!
	*/
	mdlCacheCheck(pkd->mdl);
    }
    /*
    ** Release aquired pointers and source-reactivate particles in prioq.
    */
    for (i=0;i<smx->nSmooth;++i) {
	if (smx->pq[i].iPid == pkd->idSelf) {
	    smx->bInactive[smx->pq[i].iIndex] = 0;
	}
	else {
	    PQ_HASHDEL(smx->pq[i].pPart);
	    mdlRelease(pkd->mdl,CID_PARTICLE,smx->pq[i].pPart);
	}
    }
}


void smGatherLocal(SMX smx,FLOAT fBall2,FLOAT r[3]) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,nCnt,pEnd;
    int idSelf = pkd->idSelf;

    nCnt = smx->nnListSize;
    kdn = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
	MINDIST(kdn->bnd,r,min2);
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
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
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
		    smx->nnList[nCnt].pj = pj;
		    smx->nnList[nCnt].pid = idSelf;
		    ++nCnt;
		}
	    }
	}
    NoIntersect:
	if (sp) kdn = pkdTreeNode(pkd,iCell = S[--sp]);
	else break;
    }
    smx->nnListSize = nCnt;
}


void smGatherRemote(SMX smx,FLOAT fBall2,FLOAT r[3],int id) {
    MDL mdl = smx->pkd->mdl;
    KDN *pkdn;
    PARTICLE *pp;
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int pj,nCnt,pEnd;
    int iCell;

    assert(id != smx->pkd->idSelf);
    nCnt = smx->nnListSize;
    iCell = ROOT;
    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
    while (1) {
	MINDIST(pkdn->bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	}
	/*
	** We have an intersection to test.
	*/
	if (pkdn->iLower) {
	    iCell = pkdn->iLower;
	    S[sp++] = iCell+1;
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	    continue;
	}
	else {
	    pEnd = pkdn->pUpper;
	    for (pj=pkdn->pLower;pj<=pEnd;++pj) {
		pp = mdlAquire(mdl,CID_PARTICLE,pj,id);
		if ( !pkdIsSrcActive(pp,0,MAX_RUNG) ) {
		    mdlRelease(mdl,CID_PARTICLE,pp);
		    continue;
		}
		dx = r[0] - pp->r[0];
		dy = r[1] - pp->r[1];
		dz = r[2] - pp->r[2];
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
		    smx->nnList[nCnt].pPart = pp;
		    smx->nnList[nCnt].pj = pj;
		    smx->nnList[nCnt].pid = id;
		    ++nCnt;
		}
		else mdlRelease(mdl,CID_PARTICLE,pp);
	    }
	}
    NoIntersect:
	if (sp) {
	    iCell = S[--sp];
	    mdlRelease(mdl,CID_CELL,pkdn);
	    pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
	}
	else break;
    }
    mdlRelease(mdl,CID_CELL,pkdn);
    smx->nnListSize = nCnt;
}


void smGather(SMX smx,FLOAT fBall2,FLOAT r[3]) {
    KDN *kdn;
    PKD pkd = smx->pkd;
    int idSelf = pkd->idSelf;
    int *S = smx->ST;
    FLOAT min2;
    int iCell,id;
    int sp = 0;

    kdn = pkdTopNode(pkd,iCell = ROOT);
    while (1) {
	MINDIST(kdn->bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	}
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    kdn = pkdTopNode(pkd,iCell = kdn->iLower);
	    S[sp++] = iCell+1;
	    continue;
	}
	else {
	    id = kdn->pLower; /* this is the thread id in LTT */
	    if (id != idSelf) {
		smGatherRemote(smx,fBall2,r,id);
	    }
	    else {
		smGatherLocal(smx,fBall2,r);
	    }
	}
    NoIntersect:
	if (sp) kdn = pkdTopNode(pkd,iCell = S[--sp]);
	else break;
    }
}

void smDoGatherLocal(SMX smx,FLOAT fBall2,FLOAT r[3],void (*Do)(SMX,PARTICLE *,FLOAT)) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,nCnt,pEnd;
    int idSelf = pkd->idSelf;

    nCnt = smx->nnListSize;
    kdn = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
	MINDIST(kdn->bnd,r,min2);
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
		dx = r[0] - p->r[0];
		dy = r[1] - p->r[1];
		dz = r[2] - p->r[2];
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


void smReSmoothOne(SMX smx,SMF *smf,void *p,FLOAT *R,FLOAT fBall) {
    PKD pkd = smx->pkd;
    FLOAT r[3];
    int iStart[3],iEnd[3];
    int i,j;
    int ix,iy,iz;

    smx->nnListSize = 0;
    /*
    ** Note for implementing SLIDING PATCH, the offsets for particles are
    ** negative here, reflecting the relative +ve offset of the simulation
    ** volume.
    */
    if (smx->bPeriodic) {
	for (j=0;j<3;++j) {
	    iStart[j] = floor((R[j] - fBall)/pkd->fPeriod[j] + 0.5);
	    iEnd[j] = floor((R[j] + fBall)/pkd->fPeriod[j] + 0.5);
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
    smx->fcnSmooth(p,smx->nnListSize,smx->nnList,smf);
    /*
    ** Release aquired pointers.
    */
    for (i=0;i<smx->nnListSize;++i) {
	if (smx->nnList[i].bRemote) {
	    mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[i].pPart);
	}
    }
}

void smReSmooth(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;

    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if ( pkdIsDstActive(p,0,MAX_RUNG) && pkdIsSrcActive(p,0,MAX_RUNG) )
	    smReSmoothOne(smx,smf,p,p->r,p->fBall);
    }
}


FLOAT phase_dist(PKD pkd,double dvTau2,PARTICLE *pa,PARTICLE *pb,double H) {
    int j;
    FLOAT dx,dv,dx2,dv2;
    double *va, *vb;

    assert(pkd->oGroup); /* Validate memory model */
    assert(pkd->oVelocity); /* Validate memory model */

    va = pkdVel(pkd,pa);
    vb = pkdVel(pkd,pb);

    dx2=0.0;
    for (j=0;j<3;++j) {
	dx = pa->r[j] - pb->r[j];
	dx2 += dx*dx;
    }
    dx2 /= pa->fBall;   /* this is actually fBall2! */
    dv2 = 0.0;
    if (dvTau2 > 0) {
	for (j=0;j<3;++j) {
	    dv = (va[j] - vb[j]) + H*(pa->r[j] - pb->r[j]);
	    dv2 += dv*dv;
	}
	dv2 /= dvTau2;
    }
    return(dx2 + dv2);
}


FLOAT corrPos(FLOAT com,FLOAT r,FLOAT l) {
    if (com > 0.25*l && r < -0.25*l) return r + l;
    else if (com < -0.25*l && r > 0.25*l) return r - l;
    else return r;
}


FLOAT PutInBox(FLOAT r,FLOAT l) {
    if (r < -0.5*l) return r + l;
    else if (r > 0.5*l) return r - l;
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
    int32_t *pBin;
    int32_t *pGroup;
    int32_t *pPartGroup;
    double *v;
    RB_TYPE rm_type;
    FOFRM   rm_data;
    FOFPG* protoGroup;
    FOFBIN *bin;

    int pi,pn,pnn,nCnt,i,j,k;
    int nRmListSize,nRmCnt,iRmIndex;
    int nFifo, iHead, iTail, iMaxGroups, iGroup;
    int *Fifo;
    int iStart[3],iEnd[3];
    int ix,iy,iz;
    int idSelf = pkd->idSelf;

    FLOAT r[3],l[3],relpos[3],lx,ly,lz,fBall,fBall2Max,rho,fvBall2;
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

    /* set spatial linking lenght for each particle (mass dependend when bTauAbs=0) */
    for (pn=0;pn<nTree;pn++) {
	p = pkdParticle(pkd,pn);
	fMass = pkdMass(pkd,p);
	pGroup = pkdInt32(p,pkd->oGroup);
	*pGroup = 0; 
	if (smf->bTauAbs) {
	    p->fBall = smf->dTau2;
	    /* enforce a real space linking length smaller than the mean particle separation at all times :*/
	    if (smf->dTau2 > 0.2*pow(fMass,0.6666) )
		p->fBall = 0.2*pow(fMass,0.6666);
	}
	else {
	    p->fBall = smf->dTau2*pow(fMass,0.6666);
	}
	if (p->fBall > fBall2Max) fBall2Max = p->fBall;
    }

    /* the velocity space linking lenght is the same for all particles. Zero means: do regular 3D FOF */
    fvBall2 = smf->dVTau2;
 
    /* Have to restart particle chache, since we will need the updated p->fBall now */
    mdlFinishCache(mdl,CID_PARTICLE);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),nTree);

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
	    smx->nnListSize = 0;
	    fBall = sqrt(p->fBall);
	    if (smx->bPeriodic) {
		for (j=0;j<3;++j) {
		    iStart[j] = floor((p->r[j] - fBall)/pkd->fPeriod[j] + 0.5);
		    iEnd[j] = floor((p->r[j] + fBall)/pkd->fPeriod[j] + 0.5);
		}
		for (ix=iStart[0];ix<=iEnd[0];++ix) {
		    r[0] = p->r[0] - ix*pkd->fPeriod[0];
		    for (iy=iStart[1];iy<=iEnd[1];++iy) {
			r[1] = p->r[1] - iy*pkd->fPeriod[1];
			for (iz=iStart[2];iz<=iEnd[2];++iz) {
			    r[2] = p->r[2] - iz*pkd->fPeriod[2];
			    smGather(smx,p->fBall,r);
			}
		    }
		}
	    }
	    else {
		smGather(smx,p->fBall,p->r);
	    }
	    nCnt = smx->nnListSize;
	    for (pnn=0;pnn<nCnt;++pnn ) {
		if (smx->nnList[pnn].pid == idSelf) { /* Local neighbors: */

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
		    Fifo[iTail] = smx->nnList[pnn].pj;iTail++;
		    if (iTail == nFifo) iTail = 0;
	      
		} else {	 /* Nonlocal neighbors: */
		
		    /* Make remote member linking symmetric by using smaller linking length if different: */
		    if (fvBall2 > 0.0) { /* Check velocity space distance */
			if (phase_dist(pkd,fvBall2,p,smx->nnList[pnn].pPart,smf->H) > 1.0 ||
			    phase_dist(pkd,fvBall2,smx->nnList[pnn].pPart,p,smf->H) > 1.0) continue;
		    }
		    else { /* real space distance */
			if (smx->nnList[pnn].fDist2 > smx->nnList[pnn].pPart->fBall) continue;
		    }
		
		    /* Add to remote member (RM) list if new */
		    rm_data.iIndex = smx->nnList[pnn].pj;
		    rm_data.iPid = smx->nnList[pnn].pid;

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
		    r[j] = corrPos(pkd->groupData[i].rcom[j]/pkd->groupData[i].fMass,p->r[j],l[j]);
		else  r[j] = p->r[j];
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
		if (p->fDensity > pkd->groupData[i].potordenmax) {
		    pkd->groupData[i].potordenmax = p->fDensity;
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

static void mktmpdir( const char *dirname ) {
    struct stat s;
    if ( stat(dirname,&s) == 0 ) {
	if ( S_ISDIR(s.st_mode) )
	    return;
    }
    mkdir( dirname, 0700 );
}

int smGroupMerge(SMF *smf,int bPeriodic) {
    PKD pkd = smf->pkd;
    MDL mdl = smf->pkd->mdl;
    PARTICLE *p;
    PARTICLE *pPart;
    int32_t *pBin, *pGroup;
    int32_t *pPartGroup, iPartGroup;
    FLOAT l[3], r,min,max,corr;
    int pi,id,i,j,k,index,listSize, sgListSize, lsgListSize;
    int nLSubGroups,nSubGroups,nMyGroups;
    int iHead, iTail, nFifo,tmp, nTree;
    FILE * pFile; /* files for parallel output of group ids, links and densities*/
    FILE * lFile;
    FILE * dFile;
    char filename [30];

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
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd), nTree);
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
		iTail++;
	    }
	    while (iHead != iTail) {
		rm = rmFifo[iHead];
		iHead++;
		if (iHead == nFifo) iHead = 0;
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
			lsgListSize *= 1.5;
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
			rmFifo[iTail++] = pkd->remoteMember[j];
			if (iTail == nFifo) iTail = 0;
		    }
		}
		else {
		    pPart = mdlAquire(mdl,CID_PARTICLE,rm.iIndex,rm.iPid);
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
		    sG = mdlAquire(mdl,CID_GROUP,index,rm.iPid);
		    
		    if (nSubGroups >= sgListSize) {
			sgListSize *= 1.5;
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
			remoteRM = mdlAquire(mdl,CID_RM,j,rm.iPid);
			rmFifo[iTail++] = *remoteRM;
			if (iTail == nFifo) iTail = 0;
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
		    for (j=0;j<3;j++) {
			pkd->groupData[i].rcom[j] += sG->rcom[j];
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
    /*

    /* Start RO group data cache and master reads and saves all the group data. */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->groupData,sizeof(FOFGD), nMyGroups + 1);
    if (pkd->idSelf == 0) {
	listSize = pkd->nThreads*(pkd->nGroups+1);
	pkd->groupData = (FOFGD *) realloc(pkd->groupData,listSize*sizeof(FOFGD));
	assert(pkd->groupData != NULL);
	for (id=1; id < pkd->nThreads; id++) {
	    index = 0;
	    while (1) {
		sG = mdlAquire(mdl,CID_GROUP,index,id);
		mdlRelease(mdl,CID_GROUP,sG);
		if (sG->bMyGroup != 0) {
		    if (nMyGroups >= listSize-1) {
			listSize *= 2.0;
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
	pBin = mdlAquire(smx->pkd->mdl,CID_BIN,iBin,pid);
    } else {
	pBin = pkd->groupBin[iBin];
    }
    pBin->nMembers++;
    pBin->fMassInBin += fMass;
    if (pid != smx->pkd->idSelf) mdlRelease(smx->pkd->mdl,CID_BIN,pBin);
}


int smGroupProfiles(SMX smx, SMF *smf, int nTotalGroups) {
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
	    gdp = mdlAquire(mdl,CID_GROUP,i,0);
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
		bin = mdlAquire(mdl,CID_BIN,i,0);
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
}
