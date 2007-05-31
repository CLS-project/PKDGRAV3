#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include <assert.h>
#include "mdl.h"
#include "pst.h"
#include "pkd.h"	
#include "outtype.h"
#include "smooth.h"


void pstAddServices(PST pst,MDL mdl)
    {
    int nThreads,nCell;

    nThreads = mdlThreads(mdl);
    mdlAddService(mdl,PST_SETADD,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetAdd,
		  sizeof(struct inSetAdd),0);
    mdlAddService(mdl,PST_LEVELIZE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstLevelize,
		  sizeof(struct inLevelize),0);
    mdlAddService(mdl,PST_READTIPSY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstReadTipsy,
		  sizeof(struct inReadTipsy),0);
    mdlAddService(mdl,PST_DOMAINDECOMP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDomainDecomp,
		  sizeof(struct inDomainDecomp),0);
    mdlAddService(mdl,PST_CALCBOUND,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcBound,
		  0,sizeof(struct outCalcBound));
    mdlAddService(mdl,PST_RUNGDDWEIGHT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstRungDDWeight,
		  sizeof(struct inRungDDWeight),0);
    mdlAddService(mdl,PST_WEIGHT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWeight,
		  sizeof(struct inWeight),sizeof(struct outWeight));
    mdlAddService(mdl,PST_COUNTVA,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCountVA,
		  sizeof(struct inCountVA),sizeof(struct outCountVA));
    mdlAddService(mdl,PST_WEIGHTWRAP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWeightWrap,
		  sizeof(struct inWeightWrap),sizeof(struct outWeightWrap));
    mdlAddService(mdl,PST_ORDWEIGHT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOrdWeight,
		  sizeof(struct inOrdWeight),sizeof(struct outOrdWeight));
    mdlAddService(mdl,PST_FREESTORE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFreeStore,
		  0,sizeof(struct outFreeStore));
    mdlAddService(mdl,PST_COLREJECTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstColRejects,
		  0,nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_SWAPREJECTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSwapRejects,
		  nThreads*sizeof(int),nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_COLORDREJECTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstColOrdRejects,
		  sizeof(struct inColOrdRejects),nThreads*sizeof(OREJ));
    mdlAddService(mdl,PST_DOMAINORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDomainOrder,
		  sizeof(struct inDomainOrder),0);
    mdlAddService(mdl,PST_LOCALORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstLocalOrder,
		  0,0);
    mdlAddService(mdl,PST_OUTARRAY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOutArray,
		  sizeof(struct inOutArray),0);
    mdlAddService(mdl,PST_OUTVECTOR,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOutVector,
		  sizeof(struct inOutVector),0);
    mdlAddService(mdl,PST_WRITETIPSY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWriteTipsy,
		  sizeof(struct inWriteTipsy),0);
    /*
    ** Calculate the number of levels in the top tree and use it to 
    ** define the size of the messages.
    */
    nCell = 1<<(1+(int)ceil(log((double)nThreads)/log(2.0)));
    mdlAddService(mdl,PST_BUILDTREE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstBuildTree,
		  sizeof(struct inBuildTree),nCell*sizeof(KDN));
    mdlAddService(mdl,PST_DISTRIBCELLS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDistribCells,
		  nCell*sizeof(KDN),0);
#ifdef GASOLINE
    mdlAddService(mdl,PST_CALCBOUNDBALL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcBoundBall,
		  sizeof(struct inCalcBoundBall),nCell*sizeof(BND));
    mdlAddService(mdl,PST_DISTRIBBOUNDBALL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDistribBoundBall,
		  nCell*sizeof(BND),0);
#endif
    mdlAddService(mdl,PST_CALCROOT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcRoot,
		  0,sizeof(struct ioCalcRoot));
    mdlAddService(mdl,PST_DISTRIBROOT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDistribRoot,
		  sizeof(struct ioCalcRoot),0);
    mdlAddService(mdl,PST_SMOOTH,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSmooth,
		  sizeof(struct inSmooth),0);
    mdlAddService(mdl,PST_GRAVITY,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGravity,
		  sizeof(struct inGravity),sizeof(struct outGravity));
    mdlAddService(mdl,PST_CALCEANDL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCalcEandL,
		  0,sizeof(struct outCalcEandL));
    mdlAddService(mdl,PST_DRIFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDrift,
		  sizeof(struct inDrift),0);
    mdlAddService(mdl,PST_DRIFTINACTIVE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDriftInactive,
		  sizeof(struct inDrift),0);
    mdlAddService(mdl,PST_CACHEBARRIER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCacheBarrier,
		  0,0);
    mdlAddService(mdl,PST_STEPVERYACTIVE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstStepVeryActiveKDK,
		  sizeof(struct inStepVeryActive),
		  sizeof(struct outStepVeryActive));
#ifdef HERMITE
    mdlAddService(mdl,PST_STEPVERYACTIVEH,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstStepVeryActiveHermite,
		  sizeof(struct inStepVeryActiveH),
		  sizeof(struct outStepVeryActiveH));
    mdlAddService(mdl,PST_COPY0,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCopy0,
		  sizeof(struct inCopy0),0);
    mdlAddService(mdl,PST_PREDICTOR,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPredictor,
		  sizeof(struct inPredictor),0);
    mdlAddService(mdl,PST_CORRECTOR,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCorrector,
		  sizeof(struct inCorrector),0);
    mdlAddService(mdl,PST_SUNCORRECTOR,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSunCorrector,
		  sizeof(struct inSunCorrector),0);
    mdlAddService(mdl,PST_PREDICTORINACTIVE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPredictorInactive,
		  sizeof(struct inPredictorInactive),0);
    mdlAddService(mdl,PST_AARSETHSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstAarsethStep,
		  sizeof(struct inAarsethStep), 0);
    mdlAddService(mdl,PST_FIRSTDT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFirstDt,0,0);
#endif
    mdlAddService(mdl,PST_ROPARTICLECACHE,pst,
		  (void (*)(void *,void *,int,void *,int *))pstROParticleCache,
		  0,0);
    mdlAddService(mdl,PST_PARTICLECACHEFINISH,pst,
		  (void (*)(void *,void *,int,void *,int *))pstParticleCacheFinish,
		  0,0);
    mdlAddService(mdl,PST_KICK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstKick,
		  sizeof(struct inKick),sizeof(struct outKick));
    mdlAddService(mdl,PST_SETSOFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetSoft,
		  sizeof(struct inSetSoft),0);
#ifdef CHANGESOFT
    mdlAddService(mdl,PST_PHYSICALSOFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPhysicalSoft,
		  sizeof(struct inPhysicalSoft),0);
    mdlAddService(mdl,PST_PREVARIABLESOFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPreVariableSoft,
		  0,0);
    mdlAddService(mdl,PST_POSTVARIABLESOFT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstPostVariableSoft,
		  sizeof(struct inPostVariableSoft),0);
#endif
    mdlAddService(mdl,PST_SETTOTAL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetTotal,
		  0,sizeof(struct outSetTotal));
    mdlAddService(mdl,PST_SETWRITESTART,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetWriteStart,
		  sizeof(struct inSetWriteStart),0);
    mdlAddService(mdl,PST_ONENODEREADINIT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstOneNodeReadInit,
		  sizeof(struct inReadTipsy), nThreads*sizeof(int));
    mdlAddService(mdl,PST_SWAPALL,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSwapAll,
		  sizeof(int),0);
    mdlAddService(mdl,PST_ACTIVETYPEORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveTypeOrder,
		  sizeof(struct inActiveTypeOrder),sizeof(int));
    mdlAddService(mdl,PST_ACTIVEORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveOrder,
		  0,sizeof(int));
    mdlAddService(mdl,PST_INITSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitStep,
		  sizeof(struct inInitStep),0);
    mdlAddService(mdl,PST_SETRUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetRung,
		  sizeof(struct inSetRung),0);
    mdlAddService(mdl,PST_ACTIVERUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveRung,
		  sizeof(struct inActiveRung),sizeof(int));
    mdlAddService(mdl,PST_CURRRUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCurrRung,
		  sizeof(struct inCurrRung),sizeof(struct outCurrRung));
    mdlAddService(mdl,PST_DENSITYSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDensityStep,
		  sizeof(struct inDensityStep),0);
    mdlAddService(mdl,PST_GETMAP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetMap,
		  sizeof(struct inGetMap),nThreads*sizeof(int));
    mdlAddService(mdl,PST_GRAVSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGravStep,
		  sizeof(struct inGravStep), 0);
    mdlAddService(mdl,PST_ACCELSTEP,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstAccelStep,
		  sizeof(struct inAccelStep), 0);
    mdlAddService(mdl,PST_ACTIVEEXACTTYPE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveExactType,
		  sizeof(struct inActiveType),sizeof(int));
    mdlAddService(mdl,PST_ACTIVETYPE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveType,
		  sizeof(struct inActiveType),sizeof(int));
    mdlAddService(mdl,PST_SETTYPE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetType,
		  sizeof(struct inActiveType),sizeof(int));
    mdlAddService(mdl,PST_RESETTYPE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstResetType,
		  sizeof(struct inActiveType),sizeof(int));
    mdlAddService(mdl,PST_COUNTTYPE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstCountType,
		  sizeof(struct inActiveType),sizeof(int));
    mdlAddService(mdl,PST_ACTIVEMASKRUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveMaskRung,
		  sizeof(struct inActiveType),sizeof(int));
    mdlAddService(mdl,PST_ACTIVETYPERUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstActiveTypeRung,
		  sizeof(struct inActiveType),sizeof(int));
    mdlAddService(mdl,PST_SETPARTICLETYPES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetParticleTypes,
		  sizeof(struct inSetParticleTypes),0);
    mdlAddService(mdl,PST_RESMOOTH,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstReSmooth,
		  sizeof(struct inReSmooth),0);
#ifdef GASOLINE
    mdlAddService(mdl,PST_RESMOOTHWALK,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstReSmoothWalk,
		  sizeof(struct inReSmooth),0);
#endif
    mdlAddService(mdl,PST_DTTORUNG,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDtToRung,
		  sizeof(struct inDtToRung),sizeof(struct outDtToRung));
    mdlAddService(mdl,PST_INITDT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitDt,
		  sizeof(struct inInitDt),0);
    mdlAddService(mdl,PST_COLNPARTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstColNParts,
		  0,nThreads*sizeof(struct outColNParts));
    mdlAddService(mdl,PST_NEWORDER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstNewOrder,
		  nThreads*sizeof(int),0);
    mdlAddService(mdl,PST_SETNPARTS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSetNParts,
		  sizeof(struct inSetNParts),0);
    mdlAddService(mdl,PST_CLEARTIMER,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstClearTimer,
		  sizeof(struct inClearTimer),0);
    mdlAddService(mdl,PST_FOF,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFof, 
		  sizeof(struct inFof),0);
    mdlAddService(mdl,PST_GROUPMERGE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGroupMerge, 
		  sizeof(struct inGroupMerge),sizeof(int));
    mdlAddService(mdl,PST_GROUPPROFILES,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGroupProfiles,
		  sizeof(struct inGroupProfiles),sizeof(int));
#ifdef RELAXATION
    mdlAddService(mdl,PST_INITRELAXATION,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstInitRelaxation,0,0);
#endif
#ifdef USE_MDL_IO
    mdlAddService(mdl,PST_FINDIOS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstFindIOS,
		  sizeof(struct inFindIOS),sizeof(struct outFindIOS));
    mdlAddService(mdl,PST_STARTIO,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstStartIO,
		  sizeof(struct inStartIO),0);
#endif

#ifdef PLANETS 
    mdlAddService(mdl,PST_READSS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstReadSS,
				  sizeof(struct inReadSS),0);
    mdlAddService(mdl,PST_WRITESS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstWriteSS,
				  sizeof(struct inWriteSS),0);
    mdlAddService(mdl,PST_SUNINDIRECT,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstSunIndirect,
		  sizeof(struct inSunIndirect),sizeof(struct outSunIndirect));
    mdlAddService(mdl,PST_GRAVSUN,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGravSun,
		  sizeof(struct inGravSun),0);
    mdlAddService(mdl,PST_NEXTCOLLISION,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstNextCollision,
				  0,sizeof(struct outNextCollision));
    mdlAddService(mdl,PST_GETCOLLIDERINFO,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetColliderInfo,
	    sizeof(struct inGetColliderInfo), sizeof(struct outGetColliderInfo));
    mdlAddService(mdl,PST_DOCOLLISION,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstDoCollision,
		  sizeof(struct inDoCollision),sizeof(struct outDoCollision));
    mdlAddService(mdl,PST_RESETCOLLIDERS,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstResetColliders,
		  sizeof(struct inResetColliders),0);
    mdlAddService(mdl,PST_GETVARIABLEVERYACTIVE,pst,
		  (void (*)(void *,void *,int,void *,int *)) pstGetVariableVeryActive,
				  0,sizeof(struct outGetVariableVeryActive));
#endif
    }

void pstInitialize(PST *ppst,MDL mdl,LCL *plcl)
    {
    PST pst;

    pst = (PST)malloc(sizeof(struct pstContext));
    mdlassert(mdl,pst != NULL);
    *ppst = pst;
    pst->plcl = plcl;
    pst->mdl = mdl;
    pst->idSelf = mdlSelf(mdl);
    pst->pstLower = NULL;
    pst->idUpper = -1;	/* invalidate upper 'id' */
    pst->nLeaves = 1;
    pst->nLower = 0;
    pst->nUpper = 0;
    pst->iSplitDim = -1;
    pst->iVASplitSide = 0;
    pst->nLowTot = 0;
    pst->nHighTot = 0;
    }


void pstFinish(PST pst)
    {
    PST pstKill;

    while (pst) {
	pstKill = pst;
	if(pst->nLeaves == 1 && pst->plcl->pkd)
	    pkdFinish(pst->plcl->pkd);
	pst = pst->pstLower;
	free(pstKill);
	}
    }


void pstSetAdd(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    PST pstNew;
    struct inSetAdd *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSetAdd));
    if (pst->nLeaves > 1) {
	++pst->nLeaves;
	if (pst->nLower > pst->nUpper) {
	    ++pst->nUpper;
	    mdlReqService(pst->mdl,pst->idUpper,PST_SETADD,in,nIn);
	    mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	    }
	else {
	    ++pst->nLower;
	    pstSetAdd(pst->pstLower,in,nIn,NULL,NULL);
	    }
	}
    else {
	++pst->nLeaves;
	++pst->nLower;
	++pst->nUpper;
	pstInitialize(&pstNew,pst->mdl,pst->plcl);
	pst->pstLower = pstNew;
	pst->idUpper = in->id;
	}
    if (pnOut) *pnOut = 0;
    }


void pstLevelize(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inLevelize *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inLevelize));
    pst->iLvl = in->iLvl;
    if (pst->nLeaves > 1) {
	++in->iLvl;
	mdlReqService(pst->mdl,pst->idUpper,PST_LEVELIZE,in,nIn);
	pstLevelize(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	plcl->nPstLvl = pst->iLvl+1;
	}
    if (pnOut) *pnOut = 0;
    }

void pstGetMap(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inGetMap *in = vin;
    int *out = vout;
    int *tmp,i;

    mdlassert(pst->mdl,nIn == sizeof(struct inGetMap));
    if (pst->nLeaves > 1) {
#ifdef HAVE_ALLOCA
	tmp = alloca(mdlThreads(pst->mdl)*sizeof(int));
#else
	tmp = malloc(mdlThreads(pst->mdl)*sizeof(int));
#endif
	mdlassert(pst->mdl,tmp != NULL);
	pstGetMap(pst->pstLower,in,nIn,vout,pnOut);
	in->nStart += pst->nLower;
	mdlReqService(pst->mdl,pst->idUpper,PST_GETMAP,in,nIn);
	mdlGetReply(pst->mdl,pst->idUpper,tmp,pnOut);
	for (i=0;i<pst->nUpper;++i) {
	    out[in->nStart+i] = tmp[in->nStart+i];
	    }
	in->nStart -= pst->nLower;
#ifndef HAVE_ALLOCA
	free(tmp);
#endif
	}
    else {
	out[in->nStart] = pst->idSelf;
	}
    if (pnOut) *pnOut = mdlThreads(pst->mdl)*sizeof(int);
    }

void pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inReadTipsy *in = vin;
    int *pout = vout;
    int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
    int *ptmp;
    int nThreads;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(struct inReadTipsy));
    nThreads = mdlThreads(pst->mdl);
    nFileStart = in->nFileStart;
    nFileEnd = in->nFileEnd;
    nFileTotal = nFileEnd - nFileStart + 1;
    if (pst->nLeaves > 1) {
	nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
	in->nFileStart = nFileSplit;
	mdlReqService(pst->mdl,pst->idUpper,PST_ONENODEREADINIT,in,nIn);
	in->nFileStart = nFileStart;
	in->nFileEnd = nFileSplit - 1;
	pstOneNodeReadInit(pst->pstLower,in,nIn,vout,pnOut);
	in->nFileEnd = nFileEnd;
	ptmp = malloc(nThreads*sizeof(*ptmp));
	mdlassert(pst->mdl,ptmp != NULL);
	mdlGetReply(pst->mdl,pst->idUpper,ptmp,pnOut);
	for(i = 0; i < nThreads; i++) {
	    if(ptmp[i] != -1)
		pout[i] = ptmp[i];
	    }
	free(ptmp);
	}
    else {
	for(i = 0; i < nThreads; i++)
	    pout[i] = -1;
	/*
	** Determine the size of the local particle store.
	*/
	nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
		
	pkdInitialize(&plcl->pkd,pst->mdl,nStore,in->fPeriod,
		      in->nDark,in->nGas,in->nStar);
	pout[pst->idSelf] = nFileTotal;
	}
    if (pnOut) *pnOut = nThreads*sizeof(*pout);
    }

void pstReadTipsy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inReadTipsy *in = vin;
    int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
    char achInFile[PST_FILENAME_SIZE];
    char achOutName[PST_FILENAME_SIZE];

    mdlassert(pst->mdl,nIn == sizeof(struct inReadTipsy));
    nFileStart = in->nFileStart;
    nFileEnd = in->nFileEnd;
    nFileTotal = nFileEnd - nFileStart + 1;
    if (pst->nLeaves > 1) {
	nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
	in->nFileStart = nFileSplit;
	mdlReqService(pst->mdl,pst->idUpper,PST_READTIPSY,in,nIn);
	in->nFileStart = nFileStart;
	in->nFileEnd = nFileSplit - 1;
	pstReadTipsy(pst->pstLower,in,nIn,NULL,NULL);
	in->nFileEnd = nFileEnd;
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	/*
	** Add the local Data Path to the provided filename.
	*/
	achInFile[0] = 0;
	achOutName[0] = 0;
	if (plcl->pszDataPath) {
	    strcat(achInFile,plcl->pszDataPath);
	    strcat(achInFile,"/");
	    strcat(achOutName,plcl->pszDataPath);
	    strcat(achOutName,"/");
	    }
	strcat(achInFile,in->achInFile);
	strcat(achOutName,in->achOutName);
	/*
	** Determine the size of the local particle store.
	*/
	nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
	pkdInitialize(&plcl->pkd,pst->mdl,nStore,in->fPeriod,
		      in->nDark,in->nGas,in->nStar);
	pkdReadTipsy(plcl->pkd,achInFile,achOutName,nFileStart,nFileTotal,in->bStandard,
		     in->dvFac,in->dTuFac,in->bDoublePos);
	}
    if (pnOut) *pnOut = 0;
    }


int _pstRejMatch(PST pst,int n1,OREJ *p1,int n2,OREJ *p2,int *pidSwap)
    {
    int id,i,i1=-1,i2=-1,nLarge,id1,id2;
    int s1,s2,r1,r2;

    /*
    ** Check to see if there is enough space...
    */
    s1 = 0;
    r1 = 0;
    for (i=0;i<n1;++i) {
	s1 += p1[i].nSpace;
	r1 += p1[i].nRejects;
	}
    s2 = 0;
    r2 = 0;
    for (i=0;i<n2;++i) {
	s2 += p2[i].nSpace;
	r2 += p2[i].nRejects;
	}
    mdlassert(pst->mdl,r1 <= s2);
    mdlassert(pst->mdl,r2 <= s1);
    /*
    ** First invalidate the pidSwap array.
    */
    for (id=0;id<mdlThreads(pst->mdl);++id) pidSwap[id] = -1;
    /*
    ** Now map largest nReject of p1 to largest nSpace of p2.
    */
    while (1) {
	nLarge = 0;
	for (i=0;i<n1;++i) {
	    if (p1[i].nRejects > nLarge) {
		nLarge = p1[i].nRejects;
		i1 = i;
		}
	    }
	if (nLarge == 0) break;
	nLarge = 0;
	for (i=0;i<n2;++i) {
	    if (p2[i].nSpace > nLarge) {
		nLarge = p2[i].nSpace;
		i2 = i;
		}
	    }
	if (nLarge == 0) break;
	p1[i1].nRejects = 0;
	p1[i1].nSpace = 0;
	p2[i2].nRejects = 0;
	p2[i2].nSpace = 0;
	id1 = p1[i1].id;
	id2 = p2[i2].id;
	pidSwap[id1] = id2;
	pidSwap[id2] = id1;
	}
    /*
    ** Now map largest nReject of p2 to largest nSpace of p1.
    ** However, already mapped stuff is ignored, by the above!
    */
    while (1) {
	nLarge = 0;
	for (i=0;i<n2;++i) {
	    if (p2[i].nRejects > nLarge) {
		nLarge = p2[i].nRejects;
		i2 = i;
		}
	    }
	if (nLarge == 0) break;
	nLarge = 0;
	for (i=0;i<n1;++i) {
	    if (p1[i].nSpace > nLarge) {
		nLarge = p1[i].nSpace;
		i1 = i;
		}
	    }
	if (nLarge == 0) break;
	p1[i1].nRejects = 0;
	p1[i1].nSpace = 0;
	p2[i2].nRejects = 0;
	p2[i2].nSpace = 0;
	id1 = p1[i1].id;
	id2 = p2[i2].id;
	pidSwap[id1] = id2;
	pidSwap[id2] = id1;
	}
    for (i=0;i<mdlThreads(pst->mdl);++i)
	if (pidSwap[i] != -1) return(1);
    return(0);
    }


#define MAX_ITTR	64
#define EPS_BOUND	0.01
#define MASS_EPS	1e-11
#define PARANOID_CHECK

void _pstRootSplit(PST pst,int iSplitDim,int bDoRootFind,int bDoSplitDimFind,
		   int bSplitVA) {
    int NUM_SAFETY = 4;			/* slop space when filling up memory */
    int nSafeTot;				/* total slop space we have to play with */
    int margin;					/* more slop */
    int d,ittr,nOut;
    int nLow=-1,nHigh=-1,nLowerStore,nUpperStore;
    FLOAT fLow,fHigh;
    FLOAT fl,fu,fm=-1,fmm;
    struct outFreeStore outFree;
    struct inWeight inWt;
    struct inWeightWrap inWtWrap;
    struct outWeight outWtLow;
    struct outWeight outWtHigh;
    struct outWeightWrap outWtWrLow;
    struct outWeightWrap outWtWrHigh;
    struct inCountVA inCtVA;
    struct outCountVA outCtVA;
    OREJ *pLowerRej,*pUpperRej;
    int *pidSwap,iRet;
    int nLowTot,nHighTot;
    int nLast;					/* number of particles at the last split
						   iteration */
    int nDiff=0;				/* Difference between one iteration and the next */
    char ach[256];				/* Debug */

    mdlTimer t;

    int pFlag;					/* 0 => we are splitting all particles at once. 
						   1 => we first split active, and then inactive. */
    int nTotalActive;
    int dBnd;

    mdlZeroTimer(pst->mdl,&t);
    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
    pstFreeStore(pst->pstLower,NULL,0,&outFree,NULL);
    nLowerStore = outFree.nFreeStore;
    mdlGetReply(pst->mdl,pst->idUpper,&outFree,NULL);
    nUpperStore = outFree.nFreeStore;

    mdlprintf(pst->mdl,"_pstRootSplit: id %d Level %d\n",pst->idSelf,pst->iLvl);
    mdlPrintTimer(pst->mdl,"TIME START _pstRootSplit ",&t);
    mdlprintf(pst->mdl,"_pstRootSplit: fA0 %f fIA0 %f RS? %d DC? %d\n",
	      pst->fSplit,pst->fSplitInactive,bDoRootFind, bDoSplitDimFind);
    /* Debug */
    /*
      sprintf(ach,"id: %d _pstRootSplit\n", pst->idSelf );
      mdlDiag(pst->mdl,ach);
    */

    if(bDoSplitDimFind || pst->iSplitDim == -1) {
	pst->iSplitDim = iSplitDim;
	}

    d = dBnd = pst->iSplitDim;

    mdlassert(pst->mdl,d < 3);

    fl = pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd];
    fu = pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd];
    fm = pst->fSplit;
    ittr = -1;
	
    if (bDoRootFind || fm<fl || fm>fu || (bSplitVA && (pst->iVASplitSide == 0))) {
	/*
	** First order the particles into active/inactive order...
	*/
	int nActiveOrder;
	pstActiveOrder(pst, NULL, 0, &nActiveOrder, NULL); /* SOON NO MORE ACTIVE ORDER */

	fmm = (fl + fu)/2;
	/*
	 * First find total number of active particles.
	 */
	inWt.iSplitDim = d;
	inWt.fSplit = fmm;
	inWt.ittr = 0;
	inWt.iSplitSide = 1;
	inWt.pFlag = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	nTotalActive = outWtLow.nLow + outWtHigh.nLow
	    + outWtLow.nHigh + outWtHigh.nHigh;
	mdlassert(pst->mdl,nActiveOrder == nTotalActive);
	pFlag = 1;
	if (nTotalActive <=1) {
	    pFlag = 0;			/* Divide them all */
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	    nTotalActive = outWtLow.nLow + outWtHigh.nLow
		+ outWtLow.nHigh + outWtHigh.nHigh;
	    }

	/*
	** Now start the ROOT finder based on balancing active weight ALONE!
	** (unless pFlag == 0)
	*/
	ittr = 0;
	while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
	    fm = fmm;
	    inWt.iSplitDim = d;
	    inWt.fSplit = fm;
	    inWt.ittr = ittr;
	    inWt.iSplitSide = 1;
	    inWt.pFlag = pFlag;
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	    /*
	    ** Add lower and Upper subsets weights and numbers
	    */
	    nLow = outWtLow.nLow + outWtHigh.nLow;
	    nHigh = outWtLow.nHigh + outWtHigh.nHigh;
	    fLow = outWtLow.fLow + outWtHigh.fLow;
	    fHigh = outWtLow.fHigh + outWtHigh.fHigh;
	    /*
	      printf("ittr:%d l:%d u:%d lw:%f uw:%f\n",ittr,nLow,nHigh,fLow,fHigh);
	    */
	    if(nLow == 1 && nHigh == 1) /* break on trivial case */
		break;
	    if(pFlag) {			/* split on work */
		if (fLow/pst->nLower > fHigh/pst->nUpper) fu = fm;
		else if (fLow/pst->nLower < fHigh/pst->nUpper) fl = fm;
		else break;
		}
	    else {				/* split on number */
		if (nLow/(double)pst->nLower >
		    nHigh/(double)pst->nUpper) fu = fm;
		else if (nLow/(double)pst->nLower <
			 nHigh/(double)pst->nUpper) fl = fm;
		else break;
		}
	    fmm = (fl + fu)/2;
	    ++ittr;
	    }
	/*
	** We have now found a new split. We need to decide which side of this split
	** we want to put the very active particles.
	*/
	inCtVA.iSplitDim = d;
	inCtVA.fSplit = fm;
	pstCountVA(pst,&inCtVA,sizeof(inCtVA),&outCtVA,NULL);
	if (outCtVA.nLow == 0 && outCtVA.nHigh == 0) {
	    pst->iVASplitSide = 0;
	    }
	else if (outCtVA.nHigh > outCtVA.nLow) {
	    pst->iVASplitSide = -1;
	    }
	else {
	    pst->iVASplitSide = 1;
	    }
	mdlPrintTimer(pst->mdl,"TIME active split _pstRootSplit ",&t);
	}

    pst->fSplit = fm;
	
    mdlprintf(pst->mdl, "id: %d (%d) Chose split: %f (%f,%f) %d %d\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd], pst->nLower, pst->nUpper);
    if(ittr != -1)
	mdlprintf(pst->mdl, "  Low %d %f,  High %d %f\n",
		  nLow,outWtLow.fLow + outWtHigh.fLow, nHigh,
		  outWtLow.fHigh + outWtHigh.fHigh);
    nLow = 0;
    nHigh = 0;
    fLow = 0.0;
    fHigh = 0.0;
	
    /*
    ** Now we see if the TOTAL number of particles in the lower and upper
    ** subsets exceeds the local pStores. If so then we need to find a new
    ** boundary to distribute the INACTIVE particles so that everything 
    ** fits.  
    */
    inWtWrap.iSplitDim = d;
    fl = pst->fSplit + 1e-6*pst->bnd.fMax[dBnd];
    fu = pst->fSplit - 1e-6*pst->bnd.fMax[dBnd];
	
    if (!bDoSplitDimFind) fm = pst->fSplitInactive;
    else {
	fm = 0.5*(fl+fu);
	if (fm < pst->bnd.fCenter[dBnd]) fm = pst->bnd.fCenter[dBnd] + 1.000001*pst->bnd.fMax[dBnd];
	else fm = pst->bnd.fCenter[dBnd] - 1.000001*pst->bnd.fMax[dBnd];
	}
    mdlprintf(pst->mdl, "id: %d (%d) Zeroeth guess reverse split: %f (%f,%f)\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd], 
	      pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd]);
    inWtWrap.fSplit = fm;	
    inWtWrap.fSplit2 = pst->fSplit;
    inWtWrap.ittr = 0;
    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
    inWtWrap.iSplitSide = 1;
    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
    inWtWrap.iSplitSide = 0;
    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
    mdlGetReply(pst->mdl,pst->idUpper,&outWtWrHigh,NULL);
    /*
    ** Add lower and Upper subsets numbers of particles
    */ 
    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    
    nSafeTot = nLowerStore + nUpperStore - (nLowTot + nHighTot);
    if(nSafeTot/pst->nLeaves < NUM_SAFETY) {
	NUM_SAFETY = nSafeTot/pst->nLeaves;
	sprintf(ach,"id: %d tripped inactive NUM_SAFETY %d  Low %d/%d  High %d/%d\n",
		pst->idSelf, NUM_SAFETY, nLowTot, nLowerStore, nHighTot, nUpperStore);
	mdlDiag(pst->mdl,ach);
	mdlprintf(pst->mdl,"id: %d tripped inactive NUM_SAFETY %d  Low %d/%d  High %d/%d\n",
		  pst->idSelf, NUM_SAFETY, nLowTot, nLowerStore, nHighTot, nUpperStore);
	}
	    
    margin = nSafeTot/pst->nLeaves/20;
    if (margin < NUM_SAFETY/2) margin = NUM_SAFETY/2;

    mdlprintf(pst->mdl,"id: %d  %d Low %d/%d   %d High %d/%d  NUM_SAFETY %d margin %d\n",
	      pst->idSelf, pst->nLower,nLowTot, nLowerStore, pst->nUpper,nHighTot, nUpperStore,NUM_SAFETY,margin);
	    
	    
    if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) {
	sprintf(ach,"id: %d: nLowTot > nLowerStore-NUM_SAFETY*pst->nLower %d %d %d %d\n",
		pst->idSelf, nLowTot, nLowerStore, NUM_SAFETY, pst->nLower);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fl = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1;
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;	
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*
	      mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
	      pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
	      pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if(nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    /*
	      if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) fl = fm;
	      else if (nLowTot < nLowerStore-2*NUM_SAFETY*pst->nLower) fu = fm;
	    */
	    /*
	      if (nLowTot/pst->nLower > nHighTot/pst->nUpper) fl = fm;
	      else if (nLowTot/pst->nLower < nHighTot/pst->nUpper-NUM_SAFETY) fu = fm;
	    */
	    if (nLowTot > nLowerStore-margin*pst->nLower) fl = fm;
	    else if (nLowTot < nLowerStore-2*margin*pst->nLower) fu = fm;
	    else {
		fl = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if(nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
	    if(abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nLowTot <= nLowerStore);
	mdlPrintTimer(pst->mdl,"TIME fix lower II _pstRootSplit ",&t);
	}
    else if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) {
	sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %d %d %d %d\n",
		pst->idSelf, nHighTot, nUpperStore, NUM_SAFETY, pst->nUpper);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fu = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1; 
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*				
				mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
				pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
				pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if(nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    /*
	      if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) fu = fm;
	      else if (nHighTot < nUpperStore-2*NUM_SAFETY*pst->nUpper) fl = fm;
	    */
	    /*
	      if (nHighTot/pst->nUpper > nLowTot/pst->nLower) fu = fm;
	      else if (nHighTot/pst->nUpper < nLowTot/pst->nLower-NUM_SAFETY) fl = fm;
	    */
	    if (nHighTot > nUpperStore-margin*pst->nUpper) fu = fm;
	    else if (nHighTot < nUpperStore-2*margin*pst->nUpper) fl = fm;
	    else {
		fu = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if(nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
	    if(abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFETY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nHighTot <= nUpperStore);
	mdlPrintTimer(pst->mdl,"TIME fix upper II _pstRootSplit ",&t);
	}

    if (nLowTot < NUM_SAFETY*pst->nLower) {
	sprintf(ach,"id: %d: nLowTot < NUM_SAFETY*pst->nLower %d %d %d %d\n",
		pst->idSelf, nLowTot, nLowerStore, NUM_SAFETY, pst->nLower);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fu = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1;
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;	
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*
	      mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
	      pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
	      pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if(nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    if (nLowTot > margin*pst->nLower) fl = fm;
	    else if (nLowTot < NUM_SAFETY*pst->nLower) fu = fm;
	    else {
		fl = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix too few Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if(nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
	    if(abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nLowTot <= nLowerStore);
	mdlPrintTimer(pst->mdl,"TIME fix lower II _pstRootSplit ",&t);
	}
    if (nHighTot < NUM_SAFETY*pst->nUpper) {
	sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %d %d %d %d\n",
		pst->idSelf, nHighTot, nUpperStore, NUM_SAFETY, pst->nUpper);
	mdlDiag(pst->mdl,ach);
	if (fm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd];
	if (fm < pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd]) fm=pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd];
	fl = fm;
	if (fu > fl) fmm = 0.5*(fl+fu);
	else {
	    fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
	    if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
	    mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
		      fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
	    }
	ittr = 1; 
	nLast = nLowTot;
	while (ittr < MAX_ITTR) {
	    fm = fmm;
	    inWtWrap.iSplitDim = d;
	    inWtWrap.fSplit = fm;
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = ittr;
	    inWtWrap.iVASplitSide = (bSplitVA)?pst->iVASplitSide:0;
	    inWtWrap.iSplitSide = 1;
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtWrHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    /*				
				mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
				pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
				pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	    */
	    if(nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
	    if (nHighTot > margin*pst->nUpper) fu = fm;
	    else if (nHighTot < NUM_SAFETY*pst->nUpper) fl = fm;
	    else {
		fu = fm;
		break;
		}
	    if (fu == fl)
		break;
	    else if (fu > fl) fmm = 0.5*(fl+fu);
	    else {
		fmm = 0.5*(fl+fu)+pst->bnd.fMax[dBnd];
		if (fmm > pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu)-pst->bnd.fMax[dBnd];
		mdlassert(pst->mdl, fmm >= pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd] && 
			  fmm <= pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
		}
	    ++ittr;
	    }
	mdlprintf(pst->mdl, "id: %d (%d) Fix Too few High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if(nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
	    if(abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFETY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nHighTot <= nUpperStore);
	mdlPrintTimer(pst->mdl,"TIME fix upper II _pstRootSplit ",&t);
	}

    mdlassert(pst->mdl, nLowTot >= pst->nLower);
    mdlassert(pst->mdl, nHighTot >= pst->nUpper);
    mdlassert(pst->mdl, nLowTot <= nLowerStore);
    mdlassert(pst->mdl, nHighTot <= nUpperStore);

    pst->nLowTot = nLowTot;
    pst->nHighTot = nHighTot;

    mdlPrintTimer(pst->mdl,"TIME Total Split _pstRootSplit ",&t);
	
    mdlprintf(pst->mdl, "id: %d (%d) Chose reverse split: %f (%f,%f)\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd], 
	      pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
    pst->fSplitInactive = fm;
	
    /*
    ** First Collect rejects.
    **
    ** Careful, SERVICE PST_COLREJECTS does NOT conform strictly to
    ** the proper use of MDL. This should be fixed in the future.
    ** FIXED -- MDL modified.
    */
    pLowerRej = malloc(pst->nLower*sizeof(OREJ));
    mdlassert(pst->mdl,pLowerRej != NULL);
    pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
    mdlassert(pst->mdl,pUpperRej != NULL);
    pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
    mdlassert(pst->mdl,pidSwap != NULL);
	
    mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,NULL,0);
    pstColRejects(pst->pstLower,NULL,0,pLowerRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
    mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
	
    mdlPrintTimer(pst->mdl,"TIME Collected Rejects _pstRootSplit ",&t);
	
	
    ittr = 0;
    while (1) {
	iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,
			    pst->nUpper,pUpperRej,pidSwap);
	if (!iRet) break;
	mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
		      mdlThreads(pst->mdl)*sizeof(int));
	pstSwapRejects(pst->pstLower,pidSwap,
		       mdlThreads(pst->mdl)*sizeof(int),pLowerRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);

	++ittr;
	}

	
    free(pLowerRej);
    free(pUpperRej);
    free(pidSwap);
	
    mdlPrintTimer(pst->mdl,"TIME (FINISH) Swapped Rejects _pstRootSplit ",&t);
	
    }


#define NEWSPLITDIMCUT 0.707
#define NMINFORROOTFIND 16

void pstDomainDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    int d=0,j,nBndWrapd;
    double dimsize;
    struct inDomainDecomp *in = vin;
    double l,u;
	
    mdlTimer t;

    mdlassert(pst->mdl,nIn == sizeof(struct inDomainDecomp));

    mdlZeroTimer(pst->mdl,&t);
    mdlprintf(pst->mdl,"Starting pstDomainDecomp\n");

    pst->bnd = in->bnd;
    if (pst->nLeaves > 1) {
	if (pst->iSplitDim != -1 && in->nActive < NMINFORROOTFIND) {
	    mdlprintf(pst->mdl,"Aborting RootFind -- too few actives.\n");
	    in->bDoSplitDimFind = 0;
	    in->bDoRootFind = 0;
	    }

	/*
	** Next determine the longest axis based on the bounds.
	** Don't switch dimensions unless the bounds have changed significantly.
	**
	** NB: Standard bnds don't work for Wrapped dimensions
	*/
	if (in->bDoSplitDimFind) {
	    d = pst->iSplitDim;
	    if (d==-1) {
		dimsize = -1;
		nBndWrapd = 0;
		}
	    else {
		dimsize = pst->bnd.fMax[d]*NEWSPLITDIMCUT;
		nBndWrapd = in->nBndWrap[d];
		}

	    for (j=0;j<3;++j) {
		if (in->nBndWrap[j] < nBndWrapd || pst->bnd.fMax[j] > dimsize) {
		    d=j;
		    dimsize = pst->bnd.fMax[d];
		    nBndWrapd = in->nBndWrap[d];
		    }
		}
	    }

	mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
	_pstRootSplit(pst,d,in->bDoRootFind,in->bDoSplitDimFind,in->bSplitVA);
	mdlPrintTimer(pst->mdl,"TIME RootSplit done in pstDomainDecomp",&t);

	mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
	/*
	** Now go on to DD of next levels, but pass correct wrapping bounds.
	*/
	d = pst->iSplitDim;
	nBndWrapd = in->nBndWrap[d];
	if (pst->nUpper > 1) {
	    l = pst->bnd.fCenter[d] - pst->bnd.fMax[d];
	    u = pst->bnd.fCenter[d] + pst->bnd.fMax[d];
	    in->nBndWrap[d] = nBndWrapd;
	    if (pst->fSplitInactive <= l || pst->fSplitInactive >= u) {
		l = pst->fSplit;
		}
	    else if (pst->fSplitInactive > pst->fSplit) {
		l = pst->fSplit;
		u = pst->fSplitInactive;				
		}
	    else 
		in->nBndWrap[d]++;

	    in->bnd.fMax[d] = 0.5*(u - l);
	    in->bnd.fCenter[d] = 0.5*(u + l);
	    mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINDECOMP,
			  vin,sizeof(*in));
	    }
	if (pst->nLower > 1) {
	    l = pst->bnd.fCenter[d] - pst->bnd.fMax[d];
	    u = pst->bnd.fCenter[d] + pst->bnd.fMax[d];
	    in->nBndWrap[d] = nBndWrapd;
	    if (pst->fSplitInactive <= l || pst->fSplitInactive >= u) {
		u = pst->fSplit;
		}
	    else if (pst->fSplitInactive < pst->fSplit) {
		u = pst->fSplit;
		l = pst->fSplitInactive;				
		}
	    else
		in->nBndWrap[d]++;

	    in->bnd.fMax[d] = 0.5*(u - l);
	    in->bnd.fCenter[d] = 0.5*(u + l);
	    pstDomainDecomp(pst->pstLower,vin,sizeof(*in),NULL,NULL);
	    }
	if (pst->nUpper > 1) 
	    mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    if (pnOut) *pnOut = 0;
    }


void pstCalcBound(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct outCalcBound *out = vout;
    struct outCalcBound outBnd;
	
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_CALCBOUND,NULL,0);
	pstCalcBound(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outBnd,NULL);
	BND_COMBINE(out->bnd,out->bnd,outBnd.bnd);
	}
    else {
	pkdCalcBound(plcl->pkd,&out->bnd);
	}
    if (pnOut) *pnOut = sizeof(struct outCalcBound); 
    }

void pstRungDDWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inRungDDWeight *in = vin;
	
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_RUNGDDWEIGHT,vin,nIn);
	pstRungDDWeight(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdRungDDWeight(plcl->pkd,in->iMaxRung,in->dWeight);
	}
    if (pnOut) *pnOut = 0;
    }


/*
** Make sure that the local particles are split into active and inactive
** when passing pFlag != 0.
** pFlag == 0 => weight all particles.
** pFlag > 0 => weight active particles.
** pFlag < 0 => weight inactive particles.
*/
void pstWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inWeight *in = vin;
    struct outWeight *out = vout;
    struct outWeight outWt;
    FLOAT fSplit,fLow,fHigh;
    int iSplitSide;

    mdlassert(pst->mdl,nIn == sizeof(struct inWeight));
    /*
      pkdStartTimer(plcl->pkd,7);
    */
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,in,nIn);
	pstWeight(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outWt,NULL);
	out->nLow += outWt.nLow;
	out->nHigh += outWt.nHigh;
	out->fLow += outWt.fLow;
	out->fHigh += outWt.fHigh;
	}
    else {
	fSplit = in->fSplit;
	iSplitSide = in->iSplitSide;
	if (in->ittr == 0) {
	    /*
	    ** Initialize.
	    */
	    plcl->fSplit = fSplit;
	    if (in->pFlag == 0) {
		plcl->iWtFrom = 0;
		plcl->iWtTo = pkdLocal(plcl->pkd)-1;
		}
	    else if (in->pFlag > 0) {
		/*
		** Particles must be in the active-inactive order here!
		*/
		plcl->iWtFrom = 0;
		plcl->iWtTo = pkdActive(plcl->pkd)-1;
		}
	    else {
		/*
		** Particles must be in the active-inactive order here!
		*/
		plcl->iWtFrom = pkdActive(plcl->pkd);
		plcl->iWtTo = pkdLocal(plcl->pkd)-1;
		}
	    plcl->fWtLow = 0.0;
	    plcl->fWtHigh = 0.0;
	    }
	else {
	    /*
	    ** Update the Weight Sums and use smaller weight region.
	    */
	    if (fSplit < plcl->fSplit) {
		plcl->fWtHigh += plcl->fHigh;
		if (iSplitSide) plcl->iWtFrom = plcl->iPart;
		else plcl->iWtTo = plcl->iPart-1;
		}
	    else {
		plcl->fWtLow += plcl->fLow;
		if (iSplitSide) plcl->iWtTo = plcl->iPart-1;
		else plcl->iWtFrom = plcl->iPart;
		}
	    plcl->fSplit = fSplit;
	    }
	plcl->iPart = pkdWeight(plcl->pkd,in->iSplitDim,fSplit,iSplitSide,
				plcl->iWtFrom,plcl->iWtTo,
				&out->nLow,&out->nHigh,&fLow,&fHigh);
	out->fLow = fLow + plcl->fWtLow;
	out->fHigh = fHigh + plcl->fWtHigh;
	plcl->fLow = fLow;
	plcl->fHigh = fHigh;
	if(in->pFlag > 0) {
	    if(iSplitSide) out->nLow -= pkdInactive(plcl->pkd);
	    else out->nHigh -= pkdInactive(plcl->pkd);
	    }
	if(in->pFlag < 0) {
	    if(iSplitSide) out->nHigh -= pkdActive(plcl->pkd);
	    else out->nLow -= pkdActive(plcl->pkd);
	    }
	}
    if (pnOut) *pnOut = sizeof(struct outWeight); 
    /*
      pkdStopTimer(plcl->pkd,7);
    */
    }



void pstCountVA(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    LCL *plcl = pst->plcl;
    struct inCountVA *in = vin;
    struct outCountVA *out = vout;
    struct outCountVA outCt;

    mdlassert(pst->mdl,nIn == sizeof(struct inCountVA));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_COUNTVA,in,nIn);
	pstCountVA(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outCt,NULL);
	out->nLow += outCt.nLow;
	out->nHigh += outCt.nHigh;
	}
    else {
	pkdCountVA(plcl->pkd,in->iSplitDim,in->fSplit,&out->nLow,&out->nHigh);
	}
    if (pnOut) *pnOut = sizeof(struct outCountVA); 
    }



/*
** Make sure that the local particles are split into active and inactive
** when passing pFlag != 0.
** pFlag == 0 => weight all particles.
** pFlag > 0 => weight active particles.
** pFlag < 0 => weight inactive particles.
*/
void pstWeightWrap(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inWeightWrap *in = vin;
    struct outWeightWrap *out = vout;
    struct outWeightWrap outWt;

    mdlassert(pst->mdl,nIn == sizeof(struct inWeightWrap));
    /*
      pkdStartTimer(plcl->pkd,7);
    */
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,in,nIn);
	pstWeightWrap(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outWt,NULL);
	out->nLow += outWt.nLow;
	out->nHigh += outWt.nHigh;
	}
    else {
	if (in->ittr == 0) {
	    /*
	    ** Initialize.
	    */
	    plcl->fSplit = in->fSplit;
	    plcl->iWtFrom = 0;
	    plcl->iWtTo = pkdLocal(plcl->pkd)-1;
	    }
	else {
	    /*
	    ** Update the Weight Sums and use smaller weight region.
	    */
	    if ((in->fSplit < plcl->fSplit && in->fSplit<in->fSplit2 && plcl->fSplit<in->fSplit2) ||
		(in->fSplit < plcl->fSplit && in->fSplit>in->fSplit2 && plcl->fSplit>in->fSplit2) ||
		(in->fSplit > plcl->fSplit && in->fSplit>in->fSplit2 && plcl->fSplit<in->fSplit2)) {
		if (!in->iSplitSide) plcl->iWtFrom = plcl->iPart;
		else plcl->iWtTo = plcl->iPart-1;
		}
	    else {
		if (!in->iSplitSide) plcl->iWtTo = plcl->iPart-1;
		else plcl->iWtFrom = plcl->iPart;
		}
	    plcl->fSplit = in->fSplit;
	    }
	plcl->iPart = pkdWeightWrap(plcl->pkd,in->iSplitDim,in->fSplit,in->fSplit2,in->iSplitSide,
				    in->iVASplitSide,plcl->iWtFrom,plcl->iWtTo,&out->nLow,&out->nHigh);
	/* For collect rejects */
	plcl->nSplit = plcl->iPart;
	}
    if (pnOut) *pnOut = sizeof(struct outWeightWrap); 
    /*
      pkdStopTimer(plcl->pkd,7);
    */
    }


/*
** Weight request for splitting into iOrder order.
*/
void pstOrdWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inOrdWeight *in = vin;
    struct outOrdWeight *out = vout;
    struct outOrdWeight outWt;

    mdlassert(pst->mdl,nIn == sizeof(struct inOrdWeight));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_ORDWEIGHT,in,nIn);
	pstOrdWeight(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outWt,NULL);
	out->nLow += outWt.nLow;
	out->nHigh += outWt.nHigh;
	}
    else {
	if (in->ittr == 0) {
	    /*
	    ** Initialize.
	    */
	    plcl->iOrdSplit = in->iOrdSplit;
	    plcl->iWtFrom = 0;
	    plcl->iWtTo = pkdLocal(plcl->pkd)-1;
	    }
	else {
	    /*
	    ** Update the Weight Sums and use smaller weight region.
	    */
	    if (in->iOrdSplit < plcl->iOrdSplit) {
		if (in->iSplitSide) plcl->iWtFrom = plcl->iPart;
		else plcl->iWtTo = plcl->iPart-1;
		}
	    else {
		if (in->iSplitSide) plcl->iWtTo = plcl->iPart-1;
		else plcl->iWtFrom = plcl->iPart;
		}
	    plcl->iOrdSplit = in->iOrdSplit;
	    }
	plcl->iPart = pkdOrdWeight(plcl->pkd,in->iOrdSplit,in->iSplitSide,
				   plcl->iWtFrom,plcl->iWtTo,
				   &out->nLow,&out->nHigh);
	}
    if (pnOut) *pnOut = sizeof(struct outOrdWeight); 
    }


void pstFreeStore(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct outFreeStore *out = vout;
    int nLowerStore,nUpperStore;

    mdlassert(pst->mdl,nIn == 0);
    /*
      pkdStartTimer(plcl->pkd,4);
    */
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
	pstFreeStore(pst->pstLower,NULL,0,out,NULL);
	nLowerStore = out->nFreeStore;
	mdlGetReply(pst->mdl,pst->idUpper,out,NULL);
	nUpperStore = out->nFreeStore;
	out->nFreeStore = nLowerStore + nUpperStore;
	}
    else {
	out->nFreeStore = pkdFreeStore(plcl->pkd);
	}
    if (pnOut) *pnOut = sizeof(struct outFreeStore);
    /*
      pkdStopTimer(plcl->pkd,4);
    */
    }


void pstColRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper;
	
    mdlassert(pst->mdl,nIn == 0);
    /*
      pkdStartTimer(plcl->pkd,9);
    */
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,vin,nIn);
	pstColRejects(pst->pstLower,vin,nIn,&pOutRej[0],&nLower);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,pst->idUpper,&pOutRej[iUpper],&nUpper);
	if (pnOut) *pnOut = nLower + nUpper;
	}
    else {
	pOutRej->nRejects = pkdColRejects(plcl->pkd, plcl->nSplit);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	pOutRej->nLocal = pkdLocal(plcl->pkd);
	if (pnOut) *pnOut = sizeof(OREJ);
	}
    /*
      pkdStopTimer(plcl->pkd,9);
    */
    }


void pstColOrdRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inColOrdRejects *in = vin;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper;
	
    mdlassert(pst->mdl,nIn == sizeof(struct inColOrdRejects));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,in,nIn);
	pstColOrdRejects(pst->pstLower,in,nIn,&pOutRej[0],&nLower);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,pst->idUpper,&pOutRej[iUpper],&nUpper);
	if (pnOut) *pnOut = nLower + nUpper;
	}
    else {
	pOutRej->nRejects = pkdColOrdRejects(plcl->pkd,in->iOrdSplit,
					     in->iSplitSide);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	if (pnOut) *pnOut = sizeof(OREJ);
	}
    }


void pstSwapRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    int *pidSwap = vin;
    OREJ *pOutRej = vout;
    int nLower,nUpper,iUpper,idSwap;
	
    /*	pkdStartTimer(plcl->pkd,8);*/
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,vin,nIn);
	pstSwapRejects(pst->pstLower,vin,nIn,&pOutRej[0],&nLower);
	iUpper = nLower/sizeof(OREJ);
	mdlGetReply(pst->mdl,pst->idUpper,&pOutRej[iUpper],&nUpper);
	if (pnOut) *pnOut = nLower + nUpper;
	}
    else {
	idSwap = pidSwap[pst->idSelf];
	pOutRej->nRejects = pkdSwapRejects(plcl->pkd,idSwap);
	pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
	pOutRej->id = pst->idSelf;
	if (pnOut) *pnOut = sizeof(OREJ);
	}
    /*	pkdStopTimer(plcl->pkd,8); */
    }

/*
 * Routine to swap all particles.  Note that this does not walk the pst
 * but simply works with one other processor.
 */
void pstSwapAll(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl;
    int *pidSwap = vin;
    PST lpst;

    mdlassert(pst->mdl,nIn == sizeof(*pidSwap));
    lpst = pst;
    while(lpst->nLeaves > 1)
	lpst = lpst->pstLower;

    plcl = lpst->plcl;
    pkdSwapAll(plcl->pkd, *pidSwap);
    }


void _pstOrdSplit(PST pst,int iMaxOrder)
    {
    struct outFreeStore outFree;
    struct inOrdWeight inWt;
    struct outOrdWeight outWtLow,outWtHigh;
    int im=-1,imm,il,iu;
    int nLowerStore,nUpperStore,nLow=-1,nHigh=-1;
    struct inColOrdRejects inCol;
    OREJ *pLowerRej,*pUpperRej;
    int *pidSwap,iRet,nOut,ittr;

    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
    pstFreeStore(pst->pstLower,NULL,0,&outFree,NULL);
    nLowerStore = outFree.nFreeStore;
    mdlGetReply(pst->mdl,pst->idUpper,&outFree,NULL);
    nUpperStore = outFree.nFreeStore;
    /*
    ** Find the correct iOrdSplit, such that all processors will
    ** have close to the same number of particles in the end.
    ** Start the ROOT finder based on balancing number of particles.
    */
    il = 0;
    iu = iMaxOrder;
    imm = (il + iu)/2;
    ittr = 0;
    while (il < imm && imm < iu && ittr < MAX_ITTR) {
	im = imm;
	inWt.iOrdSplit = im;
	inWt.ittr = ittr;
	inWt.iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_ORDWEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	pstOrdWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	/*
	** Add lower and Upper subsets weights and numbers
	*/
	nLow = outWtLow.nLow + outWtHigh.nLow;
	nHigh = outWtLow.nHigh + outWtHigh.nHigh;
	/*
	  printf("ittr:%d l:%d u:%d\n",ittr,nLow,nHigh);
	*/
	if(nLow == 1 && nHigh == 1) /* break on trivial case */
	    break;
	else {		/* split on number */
	    if (nLow/(double)pst->nLower >
		nHigh/(double)pst->nUpper) iu = im;
	    else if (nLow/(double)pst->nLower <
		     nHigh/(double)pst->nUpper) il = im;
	    else break;
	    }
	imm = (il + iu)/2;
	++ittr;
	}
    mdlassert(pst->mdl,nLow <= nLowerStore);
    mdlassert(pst->mdl,nHigh <= nUpperStore);
    pst->iOrdSplit = im;
    /*
    ** Collect rejects.
    */
    pLowerRej = malloc(pst->nLower*sizeof(OREJ));
    mdlassert(pst->mdl,pLowerRej != NULL);
    pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
    mdlassert(pst->mdl,pUpperRej != NULL);
    pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
    mdlassert(pst->mdl,pidSwap != NULL);
    inCol.iOrdSplit = pst->iOrdSplit;
    inCol.iSplitSide = 1;
    mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,&inCol,
		  sizeof(inCol));
    inCol.iSplitSide = 0;
    pstColOrdRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
    mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
    mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
    while (1) {
	iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,pst->nUpper,
			    pUpperRej,pidSwap);
	if (!iRet) break;
	mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
		      mdlThreads(pst->mdl)*sizeof(int));
	pstSwapRejects(pst->pstLower,pidSwap,mdlThreads(pst->mdl)*sizeof(int),
		       pLowerRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
	}
    free(pLowerRej);
    free(pUpperRej);
    free(pidSwap);
    }


void pstDomainOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inDomainOrder *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inDomainOrder));
    if (pst->nLeaves > 1) {
	_pstOrdSplit(pst,in->iMaxOrder);
	/*
	** Now go on to Domain Order of next levels.
	*/
	if (pst->nUpper > 1) 
	    mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINORDER,in,nIn);
	if (pst->nLower > 1) 
	    pstDomainOrder(pst->pstLower,in,nIn,NULL,NULL);
	if (pst->nUpper > 1) 
	    mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    if (pnOut) *pnOut = 0;
    }


void pstLocalOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_LOCALORDER,NULL,0);
	pstLocalOrder(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdLocalOrder(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }


void pstActiveTypeOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    int *pnActive = vout;
    struct inActiveTypeOrder *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVETYPEORDER,vin,nIn);
	pstActiveTypeOrder(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveTypeOrder(plcl->pkd, in->iTestMask );
	}
    if (pnOut) *pnOut = sizeof(int);
    }


void pstActiveOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    int *pnActive = vout;

    mdlassert(pst->mdl,nIn == 0);
    /*
      pkdStartTimer(pst->plcl->pkd,5);
    */
    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEORDER,NULL,0);
	pstActiveOrder(pst->pstLower,NULL,0,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveOrder(plcl->pkd);
	}
    if (pnOut) *pnOut = sizeof(int);
    /*
      pkdStopTimer(pst->plcl->pkd,5);
    */
    }


void pstOutArray(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inOutArray *in = vin;
    char achOutFile[PST_FILENAME_SIZE];

    mdlassert(pst->mdl,nIn == sizeof(struct inOutArray));
    if (pst->nLeaves > 1) {
	/*
	** Non-Recursive Text output.
	*/
	pstOutArray(pst->pstLower,in,nIn,NULL,NULL);
	mdlReqService(pst->mdl,pst->idUpper,PST_OUTARRAY,in,nIn);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	/*
	** Add the local Data Path to the provided filename.
	*/
	achOutFile[0] = 0;
	if (plcl->pszDataPath) {
	    strcat(achOutFile,plcl->pszDataPath);
	    strcat(achOutFile,"/");
	    }
	strcat(achOutFile,in->achOutFile);
	pkdOutArray(plcl->pkd,achOutFile,in->iType);
	}
    if (pnOut) *pnOut = 0;
    }


void pstOutVector(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inOutVector *in = vin;
    char achOutFile[PST_FILENAME_SIZE];

    mdlassert(pst->mdl,nIn == sizeof(struct inOutVector));
    if (pst->nLeaves > 1) {
	/*
	** Non-Recursive Text output.
	*/
	pstOutVector(pst->pstLower,in,nIn,NULL,NULL);
	mdlReqService(pst->mdl,pst->idUpper,PST_OUTVECTOR,in,nIn);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	/*
	** Add the local Data Path to the provided filename.
	*/
	achOutFile[0] = 0;
	if (plcl->pszDataPath) {
	    strcat(achOutFile,plcl->pszDataPath);
	    strcat(achOutFile,"/");
	    }
	strcat(achOutFile,in->achOutFile);
	pkdOutVector(plcl->pkd,achOutFile,in->iDim,in->iType);
	}
    if (pnOut) *pnOut = 0;
    }


void pstWriteTipsy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inWriteTipsy *in = vin;
    char achOutFile[PST_FILENAME_SIZE];

    mdlassert(pst->mdl,nIn == sizeof(struct inWriteTipsy));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_WRITETIPSY,in,nIn);
	pstWriteTipsy(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	/*
	** Add the local Data Path to the provided filename.
	*/
	achOutFile[0] = 0;
	if (plcl->pszDataPath) {
	    strcat(achOutFile,plcl->pszDataPath);
	    strcat(achOutFile,"/");
	    }
	strcat(achOutFile,in->achOutFile);
	pkdWriteTipsy(plcl->pkd,achOutFile,plcl->nWriteStart,
		      in->bStandard,in->dvFac,in->duTFac,in->bDoublePos);
	}
    if (pnOut) *pnOut = 0;
    }

#ifdef USE_MDL_IO
void pstFindIOS(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inFindIOS *in = vin;
    struct outFindIOS *out = vout;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(struct inFindIOS));
    mdlassert(pst->mdl,mdlIO(pst->mdl)<=MDL_MAX_IO_PROCS);
    if (pst->nLeaves > 1) {
	struct outFindIOS right;
	mdlassert(pst->mdl,pst->nLowTot != 0 && pst->nHighTot != 0);
	in->nLower += pst->nLowTot;
	mdlReqService(pst->mdl,pst->idUpper,PST_FINDIOS,in,nIn);
	in->nLower -= pst->nLowTot;
	pstFindIOS(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&right,pnOut);
	for( i=0; i<MDL_MAX_IO_PROCS; i++ )
	    out->nCount[i] += right.nCount[i];
	}
    else {
	mdlassert(pst->mdl,in->N>mdlIO(pst->mdl));
	pst->ioIndex = in->nLower / (in->N / mdlIO(pst->mdl) );
	mdlassert(pst->mdl,pst->ioIndex >= 0 && pst->ioIndex < mdlIO(pst->mdl));
	for( i=0; i<MDL_MAX_IO_PROCS; i++ )
	    out->nCount[i] = 0;
	out->nCount[pst->ioIndex] = pkdLocal(plcl->pkd);
        }
    if (pnOut) *pnOut = sizeof(struct outFindIOS);
    }

typedef struct ctxIO {
    PST pst;
    int nStart;
} * CTXIO;

static int pstPackIO(void *vctx, int nSize, void *vBuff)
{
    CTXIO ctx = vctx;
    LCL *plcl = ctx->pst->plcl;
    PIO *io = (PIO *)vBuff;
    int nPack;

    nPack = pkdPackIO(plcl->pkd, io, ctx->nStart, nSize/sizeof(PIO));
    ctx->nStart += nPack;
    return nPack * sizeof(PIO);
}



void pstStartIO(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inStartIO *in = vin;
    struct ctxIO ctx;
    char achOutFile[PST_FILENAME_SIZE];

    mdlassert(pst->mdl,nIn == sizeof(struct inStartIO));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_STARTIO,in,nIn);
	pstStartIO(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {

	ctx.pst = pst;
	ctx.nStart = 0;
	mdlSetComm(pst->mdl,1);
	mdlSend(pst->mdl,pst->ioIndex,pstPackIO,&ctx);
	mdlSetComm(pst->mdl,0);

	/*
	** Add the local Data Path to the provided filename.
	*/
	achOutFile[0] = 0;
	if (plcl->pszDataPath) {
	    strcat(achOutFile,plcl->pszDataPath);
	    strcat(achOutFile,"/");
	    }
	strcat(achOutFile,in->achOutFile);
	/*	pkdWriteTipsy(plcl->pkd,achOutFile,plcl->nWriteStart,
		in->bStandard,in->dvFac,in->duTFac,in->bDoublePos);*/
	}
    if (pnOut) *pnOut = 0;
    }
#endif

void pstSetSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inSetSoft *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSetSoft));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SETSOFT,in,nIn);
	pstSetSoft(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdSetSoft(plcl->pkd,in->dSoft);
	}
    if (pnOut) *pnOut = 0;
    }


void pstBuildTree(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inBuildTree *in = vin;
    KDN *pkdn = vout;
    KDN *ptmp;
    int i,iCell,iLower,iNext;

    mdlassert(pst->mdl,nIn == sizeof(struct inBuildTree));
    iCell = in->iCell;
    if (pst->nLeaves > 1) {
	in->iCell = UPPER(iCell);
	mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,in,nIn);
	in->iCell = LOWER(iCell);
	pstBuildTree(pst->pstLower,in,nIn,vout,pnOut);
	in->iCell = iCell;
	ptmp = malloc(in->nCell*sizeof(KDN));
	mdlassert(pst->mdl,ptmp != NULL);
	mdlGetReply(pst->mdl,pst->idUpper,ptmp,pnOut);
	for (i=1;i<in->nCell;++i) {
	    if (ptmp[i].pUpper) {
		pkdn[i] = ptmp[i];
		}
	    }
	free(ptmp);
	/*
	** Combine to find cell CoM, bounds and multipoles.
	** This also computes the opening radius for gravity.
	*/
	iLower = LOWER(iCell);
	iNext = UPPER(iCell);
	pkdCombineCells(&pkdn[iCell],&pkdn[iLower],&pkdn[iNext],1);
	CALCOPEN(&pkdn[iCell],in->diCrit2);
	/*
	** Set all the pointers and flags.
	*/
	pkdn[iCell].iLower = iLower;
	pkdn[iCell].iParent = 0;
	pkdn[iCell].pLower = -1;
	pkdn[iCell].pUpper = 1;
	pkdn[iLower].iParent = iCell;
	pkdn[iNext].iParent = iCell;
	}
    else {
	for (i=1;i<in->nCell;++i) pkdn[i].pUpper = 0; /* used flag = unused */
	/*
	** Passing a null bounds pointer forces pkdTreeBuild to
	** compute the bounds of the root cell before starting.
	*/
	pkdTreeBuild(plcl->pkd,in->nBucket,in->diCrit2,&pkdn[iCell],in->bTreeSqueeze,in->bExcludeVeryActive,in->dTimeStamp);
	pkdn[iCell].iLower = 0;
	pkdn[iCell].pLower = pst->idSelf;
	pkdn[iCell].pUpper = 1;
	}
    /*
    ** Calculated all cell properties, now pass up this cell info.
    */
    if (pnOut) *pnOut = in->nCell*sizeof(KDN);
    }


void pstDistribCells(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    KDN *pkdn = vin;
    int nCell;

    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBCELLS,vin,nIn);
	pstDistribCells(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	nCell = nIn/sizeof(KDN);
	pkdDistribCells(plcl->pkd,nCell,pkdn);
	}
    if (pnOut) *pnOut = 0;
    }


#ifdef GASOLINE
void pstCalcBoundBall(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inCalcBoundBall *in = vin;
    BND *pbnd = vout;
    BND *ptmp;
    int i,iCell,iLower,iNext;

    mdlassert(pst->mdl,nIn == sizeof(struct inCalcBoundBall));
    iCell = in->iCell;
    if (pst->nLeaves > 1) {
	in->iCell = UPPER(iCell);
	mdlReqService(pst->mdl,pst->idUpper,PST_CALCBOUNDBALL,in,nIn);
	in->iCell = LOWER(iCell);
	pstCalcBoundBall(pst->pstLower,in,nIn,vout,pnOut);
	in->iCell = iCell;
	ptmp = malloc(in->nCell*sizeof(BND));
	mdlassert(pst->mdl,ptmp != NULL);
	mdlGetReply(pst->mdl,pst->idUpper,ptmp,pnOut);
	for (i=1;i<in->nCell;++i) {
	    if (ptmp[i].fMax[0] > 0) {
		pbnd[i] = ptmp[i];
		}
	    }
	free(ptmp);
	/*
	** Combine the ball bounds.
	*/
	iLower = LOWER(iCell);
	iNext = UPPER(iCell);
	BND_COMBINE(pbnd[iCell],pbnd[iLower],pbnd[iNext]);
	}
    else {
	for (i=1;i<in->nCell;++i) pbnd[i].fMax[0] = 0; /* used flag = unused */
	/*
	** Passing a null bounds pointer forces pkdTreeBuild to
	** compute the bounds of the root cell before starting.
	*/
	pkdCalcBoundBall(plcl->pkd,in->fBallFactor,&pbnd[iCell]);
	}
    /*
    ** Calculated all cell properties, now pass up this cell info.
    */
    if (pnOut) *pnOut = in->nCell*sizeof(BND);
    }


void pstDistribBoundBall(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    BND *pbnd = vin;
    int nCell;

    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBBOUNDBALL,vin,nIn);
	pstDistribBoundBall(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	nCell = nIn/sizeof(BND);
	pkdDistribBoundBall(plcl->pkd,nCell,pbnd);
	}
    if (pnOut) *pnOut = 0;
    }
#endif /* of GASOLINE */

void pstCalcRoot(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct ioCalcRoot *out = vout;
    struct ioCalcRoot temp;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_CALCROOT,vin,nIn);
	pstCalcRoot(pst->pstLower,vin,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&temp,NULL);
	momAddMomc(&out->momc,&temp.momc);
	}
    else {
	pkdCalcRoot(plcl->pkd,&out->momc);
	}
    if (pnOut) *pnOut = sizeof(struct ioCalcRoot);
    }

void pstDistribRoot(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct ioCalcRoot *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct ioCalcRoot));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBROOT,vin,nIn);
	pstDistribRoot(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdDistribRoot(plcl->pkd,&in->momc);
	}
    if (pnOut) *pnOut = 0;
    }


#ifdef CHANGESOFT
void pstPhysicalSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inPhysicalSoft *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPhysicalSoft));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_PHYSICALSOFT,in,nIn);
	pstPhysicalSoft(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdPhysicalSoft(plcl->pkd,in->dSoftMax,in->dFac,in->bSoftMaxMul);
	}
    if (pnOut) *pnOut = 0;
    }


void pstPreVariableSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_PREVARIABLESOFT,vin,nIn);
	pstPreVariableSoft(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdPreVariableSoft(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }

void pstPostVariableSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inPostVariableSoft *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPostVariableSoft));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_POSTVARIABLESOFT,in,nIn);
	pstPostVariableSoft(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdPostVariableSoft(plcl->pkd,in->dSoftMax,in->bSoftMaxMul);
	}
    if (pnOut) *pnOut = 0;
    }
#endif


void pstSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SMOOTH,in,nIn);
	pstSmooth(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;

	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bGasOnly,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType,in->dfBall2OverSoft2);
	smSmooth(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }


void pstReSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inReSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inReSmooth));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTH,in,nIn);
	pstReSmooth(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;

	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bGasOnly,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType,in->dfBall2OverSoft2);
	smReSmooth(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }

#ifdef GASOLINE
void pstReSmoothWalk(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inReSmooth *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inReSmooth));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTH,in,nIn);
	pstReSmooth(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;

	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bGasOnly,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType,in->dfBall2OverSoft2);
	smReSmoothWalk(smx,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }
#endif

void pstGravity(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    CASTAT cs;
    struct inGravity *in = vin;
    struct outGravity *out = vout;
    struct outGravity outUp;
    int j,bSunLower,bSunUpper;

    mdlassert(pst->mdl,nIn == sizeof(struct inGravity));
    if (pst->nLeaves > 1) {
	/*
	** Here we need to determine which domain is closest, or containing,
	** coordinate (0,0,0), the location of the Sun!
	if (pst->fSplit < 0) { 
	    bSunLower = 0;
	    bSunUpper = in->bDoSun;
	    }
	else {
	    bSunLower = in->bDoSun;
	    bSunUpper = 0;
	    }
	in->bDoSun = bSunUpper; */

	mdlReqService(pst->mdl,pst->idUpper,PST_GRAVITY,in,nIn);
	in->bDoSun = bSunLower;
	pstGravity(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

	/*if (bSunUpper) {
	    for (j=0;j<3;++j) out->aSun[j] = outUp.aSun[j];
	    }*/

	out->dPartSum += outUp.dPartSum;
	out->dCellSum += outUp.dCellSum;
	out->dFlop += outUp.dFlop;
	out->nActive += outUp.nActive;
	out->dWSum += outUp.dWSum;
	out->dISum += outUp.dISum;
	out->dESum += outUp.dESum;
	if (outUp.dWMax > out->dWMax) out->dWMax = outUp.dWMax;
	if (outUp.dIMax > out->dIMax) out->dIMax = outUp.dIMax;
	if (outUp.dEMax > out->dEMax) out->dEMax = outUp.dEMax;
	if (outUp.dWMin < out->dWMin) out->dWMin = outUp.dWMin;
	if (outUp.dIMin < out->dIMin) out->dIMin = outUp.dIMin;
	if (outUp.dEMin < out->dEMin) out->dEMin = outUp.dEMin;
	/*
	** Cache statistics sums.
	*/
	out->dpASum += outUp.dpASum;
	out->dpMSum += outUp.dpMSum;
	out->dpCSum += outUp.dpCSum;
	out->dpTSum += outUp.dpTSum;
	out->dcASum += outUp.dcASum;
	out->dcMSum += outUp.dcMSum;
	out->dcCSum += outUp.dcCSum;
	out->dcTSum += outUp.dcTSum;
    }
    else {
      pkdGravAll(plcl->pkd,in->dTime,in->nReps,in->bPeriodic,4,in->bEwald,
		   in->iEwOrder,in->dEwCut,in->dEwhCut,in->bDoSun,out->aSun,
		   &out->nActive,&out->dPartSum,&out->dCellSum,
		   &cs,&out->dFlop);

	out->dWSum = pkdGetWallClockTimer(plcl->pkd,1);
	out->dISum = pkdGetWallClockTimer(plcl->pkd,2);
	out->dESum = pkdGetWallClockTimer(plcl->pkd,3);
	out->dWSum -= out->dISum + out->dESum;
	out->dWMax = out->dWSum;
	out->dIMax = out->dISum;
	out->dEMax = out->dESum;
	out->dWMin = out->dWSum;
	out->dIMin = out->dISum;
	out->dEMin = out->dESum;
	/*
	** Cache statistics
	*/
	out->dpASum = cs.dpNumAccess;
	out->dpMSum = cs.dpMissRatio;
	out->dpCSum = cs.dpCollRatio;
	out->dpTSum = cs.dpMinRatio;
	out->dcASum = cs.dcNumAccess;
	out->dcMSum = cs.dcMissRatio;
	out->dcCSum = cs.dcCollRatio;
	out->dcTSum = cs.dcMinRatio;
	}
    if (pnOut) *pnOut = sizeof(struct outGravity);
    }


void pstCalcEandL(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct outCalcEandL *out = vout;
    struct outCalcEandL outLcl;
    int k;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_CALCEANDL,NULL,0);
	pstCalcEandL(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outLcl,NULL);
	out->T += outLcl.T;
	out->U += outLcl.U;
	out->Eth += outLcl.Eth;
	for (k=0;k<3;k++) out->L[k] = outLcl.L[k];
	}
    else {
	pkdCalcEandL(plcl->pkd,&out->T,&out->U,&out->Eth,out->L);
	}
    if (pnOut) *pnOut = sizeof(struct outCalcEandL);
    }


void pstDrift(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inDrift *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_DRIFT,in,nIn);
	pstDrift(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdDrift(plcl->pkd,in->dTime,in->dDelta,in->fCenter,in->bPeriodic,in->bFandG,
		 in->fCentMass);
	}
    if (pnOut) *pnOut = 0;
    }

void pstDriftInactive(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inDrift *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_DRIFTINACTIVE,in,nIn);
	pstDriftInactive(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdDriftInactive(plcl->pkd,in->dTime,in->dDelta,in->fCenter,in->bPeriodic,
			 in->bFandG, in->fCentMass);
	}
    if (pnOut) *pnOut = 0;
    }

void pstCacheBarrier(PST pst,void *vin,int nIn,void *vout,int *pnOut) {
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	pstCacheBarrier(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	mdlCacheBarrier(pst->mdl,CID_CELL);
	}
    if (pnOut) *pnOut = 0;
    }

void pstStepVeryActiveKDK(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inStepVeryActive *in = vin;
    struct outStepVeryActive *out = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inStepVeryActive));
    if (pst->nLeaves > 1) {
	if(pst->iVASplitSide > 0) {
	    mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	    pstStepVeryActiveKDK(pst->pstLower,in,nIn,out,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	    }
	else if (pst->iVASplitSide < 0) {
	    mdlReqService(pst->mdl,pst->idUpper,PST_STEPVERYACTIVE,in,nIn);
	    pstCacheBarrier(pst->pstLower,NULL,0,NULL,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,out,NULL);
	    }
	else {
	    mdlassert(pst->mdl,pst->iVASplitSide != 0);
	    }
	}
    else {
	assert(plcl->pkd->nVeryActive > 0);
	
	out->nMaxRung = in->nMaxRung;
	pkdStepVeryActiveKDK(plcl->pkd,in->dStep,in->dTime,in->dDelta,
			     in->iRung, in->iRung, in->iRung, 0, in->diCrit2,
			     &out->nMaxRung, in->aSunInact, in->adSunInact, 
			     in->dSunMass);
	mdlCacheBarrier(pst->mdl,CID_CELL);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

#ifdef HERMITE
/* Hermite */
void pstStepVeryActiveHermite(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inStepVeryActiveH *in = vin;
    struct outStepVeryActiveH *out = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inStepVeryActiveH));
    if (pst->nLeaves > 1) {
	if(pst->iVASplitSide > 0) {
	    mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	    pstStepVeryActiveHermite(pst->pstLower,in,nIn,out,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	    }
	else if (pst->iVASplitSide < 0) {
	    mdlReqService(pst->mdl,pst->idUpper,PST_STEPVERYACTIVEH,in,nIn);
	    pstCacheBarrier(pst->pstLower,NULL,0,NULL,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,out,NULL);
	    }
	else {
	    mdlassert(pst->mdl,pst->iVASplitSide != 0);
	    }
	}
    else {
	assert(plcl->pkd->nVeryActive > 0);
	
	out->nMaxRung = in->nMaxRung;
	pkdStepVeryActiveHermite(plcl->pkd,in->dStep,in->dTime,in->dDelta,
			     in->iRung, in->iRung, in->iRung, 0, in->diCrit2,
			     &out->nMaxRung, in->aSunInact, in->adSunInact, 
			     in->dSunMass);
	mdlCacheBarrier(pst->mdl,CID_CELL);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }


void pstCopy0(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inCopy0 *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inCopy0));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_COPY0,in,nIn);
	pstCopy0(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdCopy0(plcl->pkd,in->dTime);
	}
    if (pnOut) *pnOut = 0;
    }

void pstPredictor(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inPredictor *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPredictor));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_PREDICTOR,in,nIn);
	pstPredictor(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdPredictor(plcl->pkd,in->dTime);
	}
    if (pnOut) *pnOut = 0;
    }

void pstCorrector(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inCorrector *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inCorrector));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_CORRECTOR,in,nIn);
	pstCorrector(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdCorrector(plcl->pkd,in->dTime);
	}
    if (pnOut) *pnOut = 0;
    }

void pstSunCorrector(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inSunCorrector *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inSunCorrector));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SUNCORRECTOR,in,nIn);
	pstSunCorrector(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdSunCorrector(plcl->pkd,in->dTime,in->dSunMass);
	}
    if (pnOut) *pnOut = 0;
    }

void pstPredictorInactive(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inPredictorInactive *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inPredictorInactive));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_PREDICTORINACTIVE,in,nIn);
	pstPredictorInactive(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdPredictorInactive(plcl->pkd,in->dTime);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstAarsethStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inAarsethStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inAarsethStep));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_AARSETHSTEP,vin,nIn);
	pstAarsethStep(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
	}
    else {
	pkdAarsethStep(plcl->pkd,in->dEta);
	}
    if (pnOut) *pnOut = 0;
    }

void pstFirstDt(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_FIRSTDT,vin,nIn);
	pstFirstDt(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
	}
    else {
	pkdFirstDt(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }
/* Hermite end */
#endif

void pstROParticleCache(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_ROPARTICLECACHE,vin,nIn);
	pstROParticleCache(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	/*
	** Start particle caching space.
	*/
	PKD pkd = plcl->pkd;
	mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
		   pkdLocal(pkd));
	
	}
    if (pnOut) *pnOut = 0;
    }

void pstParticleCacheFinish(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_PARTICLECACHEFINISH,vin,nIn);
	pstParticleCacheFinish(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	/*
	** Start particle caching space.
	*/
	PKD pkd = plcl->pkd;
	mdlFinishCache(pkd->mdl,CID_PARTICLE);
	}
    if (pnOut) *pnOut = 0;
    }

void pstKick(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inKick *in = vin;
    struct outKick *out = vout;
    struct outKick outUp;

    mdlassert(pst->mdl,nIn == sizeof(struct inKick));

    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_KICK,in,nIn);
	pstKick(pst->pstLower,in,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

	out->SumTime += outUp.SumTime;
	out->nSum += outUp.nSum;
	if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
	}
    else {
	pkdKick(plcl->pkd,in->dvFacOne,in->dvFacTwo,
		in->dvPredFacOne,in->dvPredFacTwo,in->duDelta,in->duPredDelta,
		in->iGasModel,in->z,in->duDotLimit);
	out->Time = pkdGetTimer(plcl->pkd,1);
	out->MaxTime = out->Time;
	out->SumTime = out->Time;
	out->nSum = 1;
	}
    if (pnOut) *pnOut = sizeof(struct outKick);
    }


void pstSetTotal(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {	
    LCL *plcl = pst->plcl;
    struct outSetTotal *out = vout;
    struct outSetTotal oute;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SETTOTAL,NULL,0);
	pstSetTotal(pst->pstLower,NULL,0,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&oute,NULL);
	out->nTotal += oute.nTotal;
	pst->nTotal = out->nTotal;
	}
    else {
	pst->nTotal = pkdLocal(plcl->pkd);
	out->nTotal = pst->nTotal;
	}
    if (pnOut) *pnOut = sizeof(struct outSetTotal);
    }


void pstSetWriteStart(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {	
    LCL *plcl = pst->plcl;
    struct inSetWriteStart *in = vin;
    int nWriteStart;

    mdlassert(pst->mdl,nIn == sizeof(struct inSetWriteStart));
    nWriteStart = in->nWriteStart;
    if (pst->nLeaves > 1) {
	in->nWriteStart = nWriteStart + pst->pstLower->nTotal;
	mdlReqService(pst->mdl,pst->idUpper,PST_SETWRITESTART,in,nIn);
	in->nWriteStart = nWriteStart;
	pstSetWriteStart(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	plcl->nWriteStart = nWriteStart;
	}
    if (pnOut) *pnOut = 0;
    }


void
pstSetRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inSetRung *in = vin;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SETRUNG,vin,nIn);
	pstSetRung(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdSetRung(plcl->pkd, in->iRung);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstInitStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inInitStep *in = vin;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_INITSTEP,vin,nIn);
	pstInitStep(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
      pkdInitStep(plcl->pkd,&in->param,&in->csm);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstActiveRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveRung *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVERUNG,vin,nIn);
	pstActiveRung(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveRung(plcl->pkd, in->iRung, in->bGreater);
	}
    if (pnOut) *pnOut = sizeof(int);
    }

void
pstCurrRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inCurrRung *in = vin;
    struct outCurrRung *out = vout;
    int iCurrent;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_CURRRUNG,vin,nIn);
	pstCurrRung(pst->pstLower,vin,nIn,vout,pnOut);
	iCurrent = out->iCurrent;
	mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
	if(iCurrent)
	    out->iCurrent = iCurrent;
	}
    else {
	out->iCurrent = pkdCurrRung(plcl->pkd, in->iRung);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

void
pstGravStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inGravStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inGravStep));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_GRAVSTEP,vin,nIn);
	pstGravStep(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
	}
    else {
	pkdGravStep(plcl->pkd,in->dEta,in->dRhoFac);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstAccelStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inAccelStep *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inAccelStep));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_ACCELSTEP,vin,nIn);
	pstAccelStep(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
	}
    else {
	pkdAccelStep(plcl->pkd,in->dEta,in->dVelFac,in->dAccFac,
		     in->bDoGravity,in->bEpsAcc,in->bSqrtPhi,in->dhMinOverSoft);
	}
    if (pnOut) *pnOut = 0;
    }

void
pstDensityStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inDensityStep *in = vin;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_DENSITYSTEP,vin,nIn);
	pstDensityStep(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
	}
    else {
	pkdDensityStep(plcl->pkd,in->dEta,in->dRhoFac);
	}
    if (pnOut) *pnOut = 0;
    }

void pstDtToRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct outDtToRung outTemp;
    struct inDtToRung *in = vin;
    struct outDtToRung *out = vout;
    int i;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_DTTORUNG,vin,nIn);
	pstDtToRung(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&outTemp,pnOut);
	for (i=0;i<in->iMaxRung;++i) {
	    out->nRungCount[i] += outTemp.nRungCount[i];
	    }
	}
    else {
	pkdDtToRung(plcl->pkd,in->iRung,in->dDelta, in->iMaxRung, in->bAll, out->nRungCount);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

void pstInitDt(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inInitDt *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_INITDT,vin,nIn);
	pstInitDt(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
	}
    else {
	pkdInitDt(plcl->pkd,in->dDelta);
	}
    if (pnOut) *pnOut = 0;
    }


void pstActiveExactType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveType *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));

    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEEXACTTYPE,vin,nIn);
	pstActiveExactType(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveExactType(plcl->pkd,in->iFilterMask,in->iTestMask,in->iSetMask);
	}
    if (pnOut) *pnOut = sizeof(int);
    }


void pstSetType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveType *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(struct inActiveType));

    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_SETTYPE,vin,nIn);
	pstSetType(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdSetType(plcl->pkd,in->iTestMask,in->iSetMask);
	}
    if (pnOut) *pnOut = sizeof(int);
    }


void pstResetType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveType *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(struct inActiveType));

    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_RESETTYPE,vin,nIn);
	pstResetType(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdResetType(plcl->pkd,in->iTestMask,in->iSetMask);
	}
    if (pnOut) *pnOut = sizeof(int);
    }

void pstCountType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveType *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));

    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_COUNTTYPE,vin,nIn);
	pstCountType(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdCountType(plcl->pkd,in->iFilterMask,in->iTestMask);
	}
    if (pnOut) *pnOut = sizeof(int);
    }


void pstActiveType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveType *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));

    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVETYPE,vin,nIn);
	pstActiveType(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveType(plcl->pkd,in->iTestMask,in->iSetMask);
	}
    if (pnOut) *pnOut = sizeof(int);
    }


void pstActiveMaskRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveType *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEMASKRUNG,vin,nIn);
	pstActiveMaskRung(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveMaskRung(plcl->pkd,in->iSetMask,in->iRung,in->bGreater);
	}
    if (pnOut) *pnOut = sizeof(int);
    }


void pstActiveTypeRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inActiveType *in = vin;
    int *pnActive = vout;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	int nActiveLeaf;
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVETYPERUNG,vin,nIn);
	pstActiveTypeRung(pst->pstLower,vin,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
	*pnActive += nActiveLeaf;
	}
    else {
	*pnActive = pkdActiveTypeRung(plcl->pkd,in->iTestMask,in->iSetMask,
				      in->iRung,in->bGreater);
	}
    if (pnOut) *pnOut = sizeof(int);
    }


void pstSetParticleTypes(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inSetParticleTypes *in = vin;
	
    mdlassert(pst->mdl,nIn == sizeof(*in));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SETPARTICLETYPES,vin,nIn);
	pstSetParticleTypes(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdSetParticleTypes(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }


void
pstColNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct outColNParts *out = vout;
    struct outColNParts *ptmp;
    int nThreads;
    int i;

    nThreads = mdlThreads(pst->mdl);
    for(i=0;i<nThreads;i++)
	out[i].nNew = -1;
    if(pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_COLNPARTS,vin,nIn);
	pstColNParts(pst->pstLower,vin,nIn,vout,pnOut);
	ptmp = malloc(nThreads*sizeof(*ptmp));
	mdlassert(pst->mdl,ptmp != NULL);
	mdlGetReply(pst->mdl,pst->idUpper,ptmp,pnOut);
	for(i = 0; i < nThreads; i++) {
	    if(ptmp[i].nNew != -1)
		out[i] = ptmp[i];
	    }
	free(ptmp);
	}
    else {
	pkdColNParts(plcl->pkd, &out[pst->idSelf].nNew,
		     &out[pst->idSelf].nDeltaGas,
		     &out[pst->idSelf].nDeltaDark,
		     &out[pst->idSelf].nDeltaStar);
	}
    if(pnOut) *pnOut = nThreads*sizeof(*out);
    }

void
pstNewOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    int *in = vin;
    
    if(pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_NEWORDER,vin,nIn);
	pstNewOrder(pst->pstLower, vin, nIn, vout, pnOut);
	mdlGetReply(pst->mdl, pst->idUpper, vout, pnOut);
	}
    else {
	pkdNewOrder(plcl->pkd, in[pst->idSelf]);
	}
    if(pnOut) *pnOut = 0;
    }

void
pstSetNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inSetNParts *in = vin;
    
    if(pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SETNPARTS,vin,nIn);
	pstSetNParts(pst->pstLower, vin, nIn, vout, pnOut);
	mdlGetReply(pst->mdl, pst->idUpper, vout, pnOut);
	}
    else {
	pkdSetNParts(pst->plcl->pkd, in->nGas, in->nDark, in->nStar,
		     in->nMaxOrderGas, in->nMaxOrderDark);
	}
    if(pnOut) *pnOut = 0;
    }


void 
pstClearTimer(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inClearTimer *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inClearTimer));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_CLEARTIMER,in,nIn);
	pstClearTimer(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdClearTimer(pst->plcl->pkd,in->iTimer);
	}
    if (pnOut) *pnOut = 0;
    }


void pstFof(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inFof *in = vin;
	
    mdlassert(pst->mdl,nIn == sizeof(struct inFof));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_FOF,in,nIn);
	pstFof(pst->pstLower,in,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;
	(&in->smf)->pkd = pst->plcl->pkd;
	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,0,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType,0.0);
	smFof(smx,in->nFOFsDone,&in->smf);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = 0;
    }

void pstGroupMerge(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inGroupMerge *in = vin;
    int *nGroups = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inGroupMerge));
    if (pst->nLeaves > 1) {
	int nGroupsLeaf;	
	mdlReqService(pst->mdl,pst->idUpper,PST_GROUPMERGE,in,nIn);
	pstGroupMerge(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nGroupsLeaf,pnOut);
	*nGroups += nGroupsLeaf;
	}
    else {
	(&in->smf)->pkd = pst->plcl->pkd;
	*nGroups = smGroupMerge(&in->smf,in->bPeriodic);
	}
    if (pnOut) *pnOut = sizeof(int);
    }

void pstGroupProfiles(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    struct inGroupProfiles *in = vin;
    int *nBins = vout;

    mdlassert(pst->mdl,nIn == sizeof(struct inGroupProfiles));
    if (pst->nLeaves > 1) {
	int nBinsLeaf;	
	mdlReqService(pst->mdl,pst->idUpper,PST_GROUPPROFILES,in,nIn);
	pstGroupProfiles(pst->pstLower,in,nIn,vout,pnOut);
	mdlGetReply(pst->mdl,pst->idUpper,&nBinsLeaf,pnOut);
	*nBins = nBinsLeaf;
	}
    else {
	LCL *plcl = pst->plcl;
	SMX smx;
	(&in->smf)->pkd = pst->plcl->pkd;
	smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,0,
		     in->bPeriodic,in->bSymmetric,in->iSmoothType,0.0);
	*nBins = smGroupProfiles(smx, &in->smf,in->bPeriodic,in->nTotalGroups,in->bLogBins,in->nFOFsDone);
	smFinish(smx,&in->smf);
	}
    if (pnOut) *pnOut = sizeof(int);
    }

#ifdef RELAXATION
void pstInitRelaxation(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
	
    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_INITRELAXATION,vin,nIn);
	pstInitRelaxation(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	}
    else {
	pkdInitRelaxation(plcl->pkd);
	}
    if (pnOut) *pnOut = 0;
    }
#endif /* RELAXATION */

/* PLANETS begin */
#ifdef PLANETS 
void
pstReadSS(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inReadSS *in = vin;
	int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
	char achInFile[PST_FILENAME_SIZE];

	mdlassert(pst->mdl,nIn == sizeof(struct inReadSS));
	nFileStart = in->nFileStart;
	nFileEnd = in->nFileEnd;
	nFileTotal = nFileEnd - nFileStart + 1;
	if (pst->nLeaves > 1) {
		nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
		in->nFileStart = nFileSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READSS,vin,nIn);
		in->nFileStart = nFileStart;
		in->nFileEnd = nFileSplit - 1;
		pstReadSS(pst->pstLower,vin,nIn,NULL,NULL);
		in->nFileEnd = nFileEnd;
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achInFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achInFile,plcl->pszDataPath);
			strcat(achInFile,"/");
			}
		strcat(achInFile,in->achInFile);
		/*
		 ** Determine the size of the local particle store.
		 */
		nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,nStore,in->fPeriod,
		      in->nDark,in->nGas,in->nStar);		
		pkdReadSS(plcl->pkd,achInFile,nFileStart,nFileTotal);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstWriteSS(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWriteSS *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	mdlassert(pst->mdl,nIn == sizeof(struct inWriteSS));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WRITESS,vin,nIn);
		pstWriteSS(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achOutFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achOutFile,plcl->pszDataPath);
			strcat(achOutFile,"/");
			}
		strcat(achOutFile,in->achOutFile);
		pkdWriteSS(plcl->pkd,achOutFile,plcl->nWriteStart);
		}
	if (pnOut) *pnOut = 0;
	}

void pstSunIndirect(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    /* calculate acceralation on Sun by direct summation */
    LCL *plcl = pst->plcl;
    struct inSunIndirect *in = vin;
    struct outSunIndirect *out = vout;
    struct outSunIndirect outLcl;
    int k;

    mdlassert(pst->mdl,nIn == sizeof(struct inSunIndirect));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_SUNINDIRECT,vin,nIn);
	pstSunIndirect(pst->pstLower,vin,nIn,out,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outLcl,NULL);
	for (k=0;k<3;k++){ 
         out->aSun[k] += outLcl.aSun[k];
	 out->adSun[k] += outLcl.adSun[k];
	 }
	}
    else {
	pkdSunIndirect(plcl->pkd,out->aSun,out->adSun,in->iFlag);
	}
    if (pnOut) *pnOut = sizeof(struct outSunIndirect);
     }

void pstGravSun(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
    struct inGravSun *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inGravSun));
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_GRAVSUN,vin,nIn);
	pstGravSun(pst->pstLower,vin,nIn,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
    }
    else {
      pkdGravSun(plcl->pkd,in->aSun,in->adSun,in->dSunMass);
    }
    if (pnOut) *pnOut = 0;
    }

void
pstNextCollision(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outNextCollision local,*out = vout;

	mdlassert(pst->mdl,nIn == 0);
	out->dt = DBL_MAX;
	out->iOrder1 = out->iOrder2 = -1;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_NEXTCOLLISION,NULL,0);
		pstNextCollision(pst->pstLower,NULL,0,vout,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		if (local.dt < out->dt) *out = local;
		}
	else {
		pkdNextCollision(plcl->pkd,&out->dt,&out->iOrder1,&out->iOrder2);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstGetColliderInfo(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGetColliderInfo *in = vin;
	struct outGetColliderInfo local,*out = vout;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	local.Collider.id.iOrder = out->Collider.id.iOrder = -1;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GETCOLLIDERINFO,vin,nIn);
		pstGetColliderInfo(pst->pstLower,vin,nIn,vout,NULL);
		if (out->Collider.id.iOrder == in->iOrder) local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		if (local.Collider.id.iOrder == in->iOrder) *out = local;
		}
	else {
		pkdGetColliderInfo(plcl->pkd,in->iOrder,&out->Collider);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstDoCollision(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDoCollision *in = vin;
	struct outDoCollision local,*out = vout;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	local.nOut = out->nOut = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DOCOLLISION,vin,nIn);
		pstDoCollision(pst->pstLower,vin,nIn,vout,NULL);
		if (out->nOut) local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		if (local.nOut) *out = local;
		}
	else {
		PKD pkd = plcl->pkd;
		if (in->Collider1.id.iPid == pkd->idSelf ||
		    in->Collider2.id.iPid == pkd->idSelf) {
		  pkdDoCollision(pkd,in->dt,&in->Collider1,&in->Collider2,
				 in->bPeriodic,&in->CP,&out->iOutcome,&out->dT,
				 out->Out,&out->nOut);
			}
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstResetColliders(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inResetColliders *in = vin;
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RESETCOLLIDERS,vin,nIn);
		pstResetColliders(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdResetColliders(plcl->pkd,in->iOrder1,in->iOrder2);
		}
	if (pnOut) *pnOut = 0;
	}

void pstGetVariableVeryActive(PST pst,void *vin,int nIn,void *vout,int *pnOut)
    {
    LCL *plcl = pst->plcl;
   
    struct outGetVariableVeryActive *out = vout;

    mdlassert(pst->mdl,nIn == 0);
    if (pst->nLeaves > 1) {
	if(pst->iVASplitSide > 0) {
	    mdlReqService(pst->mdl,pst->idUpper,PST_CACHEBARRIER,NULL,0);
	    pstGetVariableVeryActive(pst->pstLower,NULL,0,out,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	    }
	else if (pst->iVASplitSide < 0) {
	    mdlReqService(pst->mdl,pst->idUpper,PST_GETVARIABLEVERYACTIVE,NULL,0);
	    pstCacheBarrier(pst->pstLower,NULL,0,NULL,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,out,NULL);
	    }
	else {
	    mdlassert(pst->mdl,pst->iVASplitSide != 0);
	    }
	}
    else {
	assert(plcl->pkd->nVeryActive > 0);	
	pkdGetVariableVeryActive(plcl->pkd, &out->dDeltaEcoll);
	mdlCacheBarrier(pst->mdl,CID_CELL);
	}
    if (pnOut) *pnOut = sizeof(*out);
    }

#endif
/* PLANETS end */
