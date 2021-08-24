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
#include "olddd.h"
#include "freestore.h"

namespace OldDD {

constexpr double NEWSPLITDIMCUT = 0.707;
constexpr int NMINFORROOTFIND = 16;
constexpr int MAX_ITTR = 64;

/*****************************************************************************\
* ServiceDomainDecomp
\*****************************************************************************/

static_assert(std::is_void<ServiceDomainDecomp::input>()  || std::is_trivial<ServiceDomainDecomp::input>());
static_assert(std::is_void<ServiceDomainDecomp::output>() || std::is_trivial<ServiceDomainDecomp::output>());

int ServiceDomainDecomp::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto in = static_cast<input *>(vin);
    auto out = static_cast<output *>(vout);
    auto pkd = pst->plcl->pkd;
    assert(nIn == sizeof(input));

    pst->bnd = in->bnd;

    int nBndWrapd, d=0;
    mdlTimer t;

#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Domain Decomposition");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
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
	double dimsize;
	if (d==-1) {
	    dimsize = -1;
	    nBndWrapd = 0;
	    }
	else {
	    dimsize = pst->bnd.fMax[d]*NEWSPLITDIMCUT;
	    nBndWrapd = in->nBndWrap[d];
	    }

	for (auto j=0;j<3;++j) {
	    if (in->nBndWrap[j] < nBndWrapd || pst->bnd.fMax[j] > dimsize) {
		d=j;
		dimsize = pst->bnd.fMax[d];
		nBndWrapd = in->nBndWrap[d];
		}
	    }
	}

    mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
    RootSplit(pst,d,in->bDoRootFind,in->bDoSplitDimFind);
    mdlPrintTimer(pst->mdl,"TIME RootSplit done in pstDomainDecomp",&t);

    mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
    /*
    ** Now go on to DD of next levels, but pass correct wrapping bounds.
    */
    d = pst->iSplitDim;
    nBndWrapd = in->nBndWrap[d];

    auto l = pst->bnd.fCenter[d] - pst->bnd.fMax[d];
    auto u = pst->bnd.fCenter[d] + pst->bnd.fMax[d];
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
    auto rID = mdl->ReqService(pst->idUpper,PST_DOMAINDECOMP,vin,sizeof(*in));

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
    Traverse(pst->pstLower,vin,sizeof(*in),NULL,0);

    mdl->GetReply(rID);
#ifdef USE_ITT
    __itt_task_end(domain);
#endif

    return 0;
    }

int ServiceDomainDecomp::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in   = static_cast<input*>(vin);
    auto out  = static_cast<output*>(vout);
    auto pkd = pst->plcl->pkd;

    pst->bnd = in->bnd;

    float offs;
    /*
    ** We always set plcl->pkd->bnd from pst->bnd.
    */
    pkd->bnd = pst->bnd;   /* This resets the local bounding box, but doesn't squeeze! */
    offs= 0.5f / (pkd->nLocal*1.0f - 1.0f);
    for (auto j=0; j < 3; j++) {
	pst->bnd.fMax[j] += offs;
	}


    return 0;
    }

void ServiceDomainDecomp::RootSplit(PST pst,int iSplitDim,int bDoRootFind,int bDoSplitDimFind) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    int NUM_SAFETY = 4;			/* slop space when filling up memory */
    uint64_t nSafeTot;			/* total slop space we have to play with */
    uint64_t margin;			/* more slop */
    int d,ittr,nOut;
    /*
    ** Why are these initialized to -1 here???
    int nLow=-1,nHigh=-1;
    */
    uint64_t nLow,nHigh;
    uint64_t nLowerStore,nUpperStore;
    uint64_t nLowTot,nHighTot;
    uint64_t nLast;					/* number of particles at the last split iteration */
    uint64_t nTotalActive;
    int nDiff=0;				/* Difference between one iteration and the next, seems to only be used to warn. */
    double fLow,fHigh;
    double fl,fu,fm=-1,fmm;
    ServiceFreeStore::output outFree;
    ServiceWeight::output outWtLow, outWtHigh;
    ServiceWeightWrap::output outWtWrLow,outWtWrHigh;
    int iRet;
    char ach[256];				/* Debug */
    mdlTimer t;
    int pFlag;					/* 0 => we are splitting all particles at once. 1 => we first split active, and then inactive. */
    int dBnd;
    int rID;

    mdlZeroTimer(pst->mdl,&t);
    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    rID = mdl->ReqService(pst->idUpper,PST_FREESTORE);
    Traverse(PST_FREESTORE,pst->pstLower,NULL,0,&outFree,sizeof(outFree));
    nLowerStore = outFree;
    mdl->GetReply(rID,&outFree);
    nUpperStore = outFree;

    mdlprintf(pst->mdl,"_pstRootSplit: id %d Level %d\n",pst->idSelf,pst->iLvl);
    mdlPrintTimer(pst->mdl,"TIME START _pstRootSplit ",&t);
    mdlprintf(pst->mdl,"_pstRootSplit: fA0 %f fIA0 %f RS? %d DC? %d\n",
	      pst->fSplit,pst->fSplitInactive,bDoRootFind, bDoSplitDimFind);
    /* Debug */
    /*
      sprintf(ach,"id: %d _pstRootSplit\n", pst->idSelf );
      mdlDiag(pst->mdl,ach);
    */

    if (bDoSplitDimFind || pst->iSplitDim == -1) {
	pst->iSplitDim = iSplitDim;
	}

    d = dBnd = pst->iSplitDim;
    fl = pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd];
    fu = pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd];
    fm = pst->fSplit;
    ittr = -1;

    if (bDoRootFind || fm<fl || fm>fu) {
	/*
	** First order the particles into active/inactive order...
	*/
	uint64_t nActiveOrder;
	pstActiveOrder(pst, NULL, 0, &nActiveOrder, sizeof(nActiveOrder)); /* SOON NO MORE ACTIVE ORDER */

	fmm = (fl + fu)/2;
	/*
	 * First find total number of active particles.
	 */
	ServiceWeight::input inWt;
	inWt.iSplitDim = d;
	inWt.fSplit = fmm;
	inWt.ittr = 0;
	inWt.iSplitSide = 1;
	inWt.pFlag = 1;
	rID = mdl->ReqService(pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	Traverse(PST_WEIGHT,pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
	mdl->GetReply(rID,&outWtHigh);
	nTotalActive = outWtLow.nLow + outWtHigh.nLow
		       + outWtLow.nHigh + outWtHigh.nHigh;
	mdlassert(pst->mdl,nActiveOrder == nTotalActive);
	pFlag = 1;
	if (nTotalActive <=1) {
	    pFlag = 0;			/* Divide them all */
	    rID = mdl->ReqService(pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    Traverse(PST_WEIGHT,pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
	    mdl->GetReply(rID,&outWtHigh);
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
	    rID = mdl->ReqService(pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    Traverse(PST_WEIGHT,pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
	    mdl->GetReply(rID,&outWtHigh);
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
	    if (nLow == 1 && nHigh == 1) /* break on trivial case */
		break;
	    if (pFlag) {			/* split on work */
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
	mdlPrintTimer(pst->mdl,"TIME active split _pstRootSplit ",&t);
	}

    pst->fSplit = fm;

    mdlprintf(pst->mdl, "id: %d (%d) Chose split: %f (%f,%f) %d %d\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd] - pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd] + pst->bnd.fMax[dBnd], pst->nLower, pst->nUpper);
    if (ittr != -1)
	mdlprintf(pst->mdl, "  Low %" PRIu64 " %f,  High %" PRIu64 " %f, ittr=%d\n",
		  nLow,outWtLow.fLow + outWtHigh.fLow, nHigh,
		  outWtLow.fHigh + outWtHigh.fHigh,ittr);
    nLow = 0;
    nHigh = 0;
    fLow = 0.0;
    fHigh = 0.0;

    /*
    ** Now we see if the TOTAL number of particles in the lower and upper
    ** subsets exceeds the local particle stores. If so then we need to
    ** find a new boundary to distribute the INACTIVE particles so that
    ** everything fits.
    */
    ServiceWeightWrap::input inWtWrap;
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
    inWtWrap.iSplitSide = 1;
    rID = mdl->ReqService(pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
    inWtWrap.iSplitSide = 0;
    Traverse(PST_WEIGHTWRAP,pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
    mdl->GetReply(rID,&outWtWrHigh);
    /*
    ** Add lower and Upper subsets numbers of particles
    */
    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;

    nSafeTot = nLowerStore + nUpperStore - (nLowTot + nHighTot);
    if (nSafeTot/pst->nLeaves < NUM_SAFETY) {
	NUM_SAFETY = nSafeTot/pst->nLeaves;
	sprintf(ach,"id: %d tripped inactive NUM_SAFETY %d  Low %" PRIu64 "/%" PRIu64 "  High %" PRIu64 "/%" PRIu64 "\n",
		pst->idSelf, NUM_SAFETY, nLowTot, nLowerStore, nHighTot, nUpperStore);
	mdlDiag(pst->mdl,ach);
	mdlprintf(pst->mdl,"id: %d tripped inactive NUM_SAFETY %d  Low %%" PRIu64 "/%" PRIu64 "  High %" PRIu64 "/%" PRIu64 "\n",
		  pst->idSelf, NUM_SAFETY, nLowTot, nLowerStore, nHighTot, nUpperStore);
	}

    margin = nSafeTot/pst->nLeaves/20;
    if (margin < NUM_SAFETY/2) margin = NUM_SAFETY/2;

    mdlprintf(pst->mdl,"id: %d  %d Low %" PRIu64 "/%" PRIu64 "   %d High %" PRIu64 "/%" PRIu64 "  NUM_SAFETY %d margin %d\n",
	      pst->idSelf, pst->nLower,nLowTot, nLowerStore, pst->nUpper,nHighTot, nUpperStore,NUM_SAFETY,margin);


    if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) {
	sprintf(ach,"id: %d: nLowTot > nLowerStore-NUM_SAFETY*pst->nLower %" PRIu64 " %" PRIu64 " %d %d\n",
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
	    inWtWrap.iSplitSide = 1;
	    rID = mdl->ReqService(pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    Traverse(PST_WEIGHTWRAP,pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
	    mdl->GetReply(rID,&outWtWrHigh);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    if (nLowTot != nLast)
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
	mdlprintf(pst->mdl, "id: %d (%d) Fix Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
	    if (abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nLowTot <= nLowerStore);
	mdlPrintTimer(pst->mdl,"TIME fix lower II _pstRootSplit ",&t);
	}
    else if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) {
	sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %" PRIu64 " %" PRIu64 " %d %d\n",
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
	    inWtWrap.iSplitSide = 1;
	    rID = mdl->ReqService(pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    Traverse(PST_WEIGHTWRAP,pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
	    mdl->GetReply(rID,&outWtWrHigh);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    if (nLowTot != nLast)
		nDiff = nLowTot - nLast;
	    nLast = nLowTot;
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
	mdlprintf(pst->mdl, "id: %d (%d) Fix High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
	    if (abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFETY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nHighTot <= nUpperStore);
	mdlPrintTimer(pst->mdl,"TIME fix upper II _pstRootSplit ",&t);
	}

    if (nLowTot < NUM_SAFETY*pst->nLower) {
	sprintf(ach,"id: %d: nLowTot < NUM_SAFETY*pst->nLower %" PRIu64 " %" PRIu64 " %d %d\n",
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
	    inWtWrap.iSplitSide = 1;
	    rID = mdl->ReqService(pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    Traverse(PST_WEIGHTWRAP,pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
	    mdl->GetReply(rID,&outWtWrHigh);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    if (nLowTot != nLast)
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
	mdlprintf(pst->mdl, "id: %d (%d) Fix too few Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
	    if (abs(nDiff) > 1)
		mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
			  pst->idSelf, nDiff);
	    }
	mdlassert(pst->mdl,nLowTot <= nLowerStore);
	mdlPrintTimer(pst->mdl,"TIME fix lower II _pstRootSplit ",&t);
	}
    if (nHighTot < NUM_SAFETY*pst->nUpper) {
	sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %" PRIu64 " %" PRIu64 " %d %d\n",
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
	    inWtWrap.iSplitSide = 1;
	    rID = mdl->ReqService(pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    Traverse(PST_WEIGHTWRAP,pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtWrLow,sizeof(outWtWrLow));
	    mdl->GetReply(rID,&outWtWrHigh);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */
	    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
	    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
	    if (nLowTot != nLast)
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
	mdlprintf(pst->mdl, "id: %d (%d) Fix Too few High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n",
		  pst->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
		  pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd],nLowTot,nHighTot);
	if (nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
	    if (abs(nDiff) > 1)
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

    mdlPrintTimer(pst->mdl,"TIME Total Split _pstRootSplit ",&t);

    mdlprintf(pst->mdl, "id: %d (%d) Chose reverse split: %f (%f,%f)\n",
	      pst->idSelf, pst->iLvl, fm, pst->bnd.fCenter[dBnd]-pst->bnd.fMax[dBnd],
	      pst->bnd.fCenter[dBnd]+pst->bnd.fMax[dBnd]);
    pst->fSplitInactive = fm;

    auto pLowerRej = std::make_unique<ServiceColRejects::output[]>(pst->nLower);
    auto pUpperRej = std::make_unique<ServiceColRejects::output[]>(pst->nUpper);
    auto pidSwap = std::make_unique<int[]>(mdl->Threads());

    rID = mdl->ReqService(pst->idUpper,PST_COLREJECTS,NULL,0);
    nOut = Traverse(PST_COLREJECTS,pst->pstLower,NULL,0,pLowerRej.get(),mdl->Threads()*sizeof(ServiceColRejects::output));
    assert(nOut/sizeof(ServiceColRejects::output) == pst->nLower);
    nOut = mdl->GetReply(rID,pUpperRej.get());
    assert(nOut/sizeof(ServiceColRejects::output) == pst->nUpper);

    mdlPrintTimer(pst->mdl,"TIME Collected Rejects _pstRootSplit ",&t);


    ittr = 0;
    while (1) {
	iRet = RejMatch(pst,pst->nLower,pLowerRej.get(),pst->nUpper,pUpperRej.get(),pidSwap.get());
	if (!iRet) break;
	rID = mdl->ReqService(pst->idUpper,PST_SWAPREJECTS,pidSwap.get(),
		      mdl->Threads()*sizeof(int));
	nOut = Traverse(PST_SWAPREJECTS,pst->pstLower,
			pidSwap.get(),mdl->Threads()*sizeof(int),
			pLowerRej.get(),mdl->Threads()*sizeof(ServiceSwapRejects::output));
	assert(nOut/sizeof(ServiceSwapRejects::output) == pst->nLower);
	nOut = mdl->GetReply(rID,pUpperRej.get());
	assert(nOut/sizeof(ServiceSwapRejects::output) == pst->nUpper);

	++ittr;
	}

    mdlPrintTimer(pst->mdl,"TIME (FINISH) Swapped Rejects _pstRootSplit ",&t);
    }

/*****************************************************************************\
* ServiceDomain
\*****************************************************************************/

int ServiceDomain::RejMatch(PST pst,int n1,ServiceColRejects::output *p1,int n2,ServiceColRejects::output *p2,int *pidSwap) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    int id,i,i1=-1,i2=-1,id1,id2;
    local_t nLarge;
    total_t s1,s2,r1,r2;

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
    assert(r1 <= s2);
    assert(r2 <= s1);
    /*
    ** First invalidate the pidSwap array.
    */
    for (id=0;id<mdl->Threads();++id) pidSwap[id] = -1;
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
    for (i=0;i<mdl->Threads();++i)
	if (pidSwap[i] != -1) return(1);
    return(0);
    }

/*****************************************************************************\
* ServiceColRejects
\*****************************************************************************/

static_assert(std::is_void<ServiceColRejects::input>()  || std::is_trivial<ServiceColRejects::input>());
static_assert(std::is_void<ServiceColRejects::output>() || std::is_trivial<ServiceColRejects::output>());

int ServiceColRejects::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto pOutRej = static_cast<output*>(vout);
    auto rID = mdl->ReqService(pst->idUpper,PST_COLREJECTS,vin,nIn);
    auto nLower = Traverse(pst->pstLower,vin,nIn,pOutRej,nOut);
    auto iUpper = nLower/sizeof(output);
    int nUpper;
    nUpper = mdl->GetReply(rID,pOutRej+iUpper);
    return nLower + nUpper;
    }

int ServiceColRejects::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto pkd = pst->plcl->pkd;
    auto pOutRej = static_cast<output*>(vout);
    pOutRej->nRejects = pkdColRejects(pkd, pst->plcl->nSplit);
    pOutRej->nSpace = pkdSwapSpace(pkd);
    pOutRej->id = pst->idSelf;
    pOutRej->nLocal = pkdLocal(pkd);
    return sizeof(output);
    }

/*****************************************************************************\
* ServiceSwapRejects
\*****************************************************************************/

static_assert(std::is_void<ServiceSwapRejects::input>()  || std::is_trivial<ServiceSwapRejects::input>());
static_assert(std::is_void<ServiceSwapRejects::output>() || std::is_trivial<ServiceSwapRejects::output>());

int ServiceSwapRejects::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto pOutRej = static_cast<output*>(vout);
    auto rID = mdl->ReqService(pst->idUpper,PST_SWAPREJECTS,vin,nIn);
    auto nLower = Traverse(pst->pstLower,vin,nIn,pOutRej,nOut);
    auto iUpper = nLower/sizeof(ServiceSwapRejects::output);
    int nUpper;
    nUpper = mdl->GetReply(rID,pOutRej+iUpper);
    return nLower + nUpper;
    }

int ServiceSwapRejects::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto pkd = pst->plcl->pkd;
    auto pidSwap = static_cast<input*>(vin);
    auto pOutRej = static_cast<output*>(vout);
    auto idSwap = pidSwap[pst->idSelf];
    pOutRej->nRejects = pkdSwapRejects(pkd,idSwap);
    pOutRej->nSpace = pkdSwapSpace(pkd);
    pOutRej->id = pst->idSelf;
    return sizeof(output);
    }

/*****************************************************************************\
* ServiceColOrdRejects
\*****************************************************************************/

static_assert(std::is_void<ServiceColOrdRejects::input>()  || std::is_trivial<ServiceColOrdRejects::input>());
static_assert(std::is_void<ServiceColOrdRejects::output>() || std::is_trivial<ServiceColOrdRejects::output>());

int ServiceColOrdRejects::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto pOutRej = static_cast<output*>(vout);
    auto rID = mdl->ReqService(pst->idUpper,PST_COLORDREJECTS,vin,nIn);
    auto nLower = Traverse(pst->pstLower,vin,nIn,&pOutRej[0],nOut);
    auto iUpper = nLower/sizeof(output);
    int nUpper;
    mdlGetReply(pst->mdl,rID,&pOutRej[iUpper],&nUpper);
    return nLower + nUpper;
    }

int ServiceColOrdRejects::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto pkd = pst->plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto pOutRej = static_cast<output*>(vout);
    pOutRej->nRejects = pkdColOrdRejects(pkd,in->iOrdSplit,in->iSplitSide);
    pOutRej->nSpace = pkdSwapSpace(pkd);
    pOutRej->id = pst->idSelf;
    return sizeof(output);
    }

/*****************************************************************************\
* ServiceDomainOrder
\*****************************************************************************/

static_assert(std::is_void<ServiceDomainOrder::input>()  || std::is_trivial<ServiceDomainOrder::input>());
static_assert(std::is_void<ServiceDomainOrder::output>() || std::is_trivial<ServiceDomainOrder::output>());

int ServiceDomainOrder::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto in = static_cast<input*>(vin);
    int rID;
    uint64_t iMinOrder = in->iMinOrder;
    uint64_t iMaxOrder = in->iMaxOrder;
    uint64_t iMidOrder;
    iMidOrder = OrdSplit(pst,iMinOrder,iMaxOrder);
    /*
    ** Now go on to Domain Order of next levels.
    */
    in->iMinOrder = iMidOrder;
    if (pst->nUpper > 1) rID = mdl->ReqService(pst->idUpper,PST_DOMAINORDER,in,nIn);
    in->iMinOrder = iMinOrder;
    in->iMaxOrder = iMidOrder-1;
    if (pst->nLower > 1) Traverse(pst->pstLower,in,nIn,NULL,0);
    if (pst->nUpper > 1) mdl->GetReply(rID);

    return 0;
    }

int ServiceDomainOrder::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    return 0;
    }

uint64_t ServiceDomainOrder::OrdSplit(PST pst,uint64_t iMinOrder,uint64_t iMaxOrder) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;

    ServiceOrdWeight::input inWt;
    ServiceOrdWeight::output outWtLow,outWtHigh;
    uint64_t im,imm,il,iu;
    uint64_t nLowerStore,nUpperStore,nLow,nHigh;
    ServiceColOrdRejects::input inCol;
    int iRet,nOut,ittr;
    int rID;

    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    ServiceFreeStore::output outFree;
    rID = mdl->ReqService(pst->idUpper,PST_FREESTORE,NULL,0);
    Traverse(PST_FREESTORE,pst->pstLower,NULL,0,&outFree,sizeof(outFree));
    nLowerStore = outFree;
    mdl->GetReply(rID,&outFree);
    nUpperStore = outFree;
    /*
    ** Find the correct iOrdSplit, such that all processors will
    ** have close to the same number of particles in the end.
    ** Start the ROOT finder based on balancing number of particles.
    */ 
    il = iMinOrder;
    iu = iMaxOrder;
    im = 0;        /* just initialization */
    nLow = 0;      /* just initialization */
    nHigh = iu+1;  /* just initialization */
    imm = (il + iu + 1)/2;
    ittr = 0;
    while (il < imm && imm < iu && ittr < MAX_ITTR) {
	im = imm;
	inWt.iOrdSplit = im;
	inWt.ittr = ittr;
	inWt.iSplitSide = 1;
	rID = mdl->ReqService(pst->idUpper,PST_ORDWEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	Traverse(PST_ORDWEIGHT,pst->pstLower,&inWt,sizeof(inWt),&outWtLow,sizeof(outWtLow));
	mdl->GetReply(rID,&outWtHigh);
	/*
	** Add lower and Upper subsets weights and numbers
	*/
	nLow = outWtLow.nLow + outWtHigh.nLow;
	nHigh = outWtLow.nHigh + outWtHigh.nHigh;
	/*
	  printf("ittr:%d l:%d u:%d : %llu < %llu < %llu\n",
	      ittr,nLow,nHigh, il, imm, iu);
	*/
	if (nLow == 1 && nHigh == 1) /* break on trivial case */
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
    auto pLowerRej = std::make_unique<ServiceColRejects::output[]>(pst->nLower);
    auto pUpperRej = std::make_unique<ServiceColRejects::output[]>(pst->nUpper);
    auto pidSwap = std::make_unique<int[]>(mdl->Threads());
    assert(pLowerRej);
    assert(pUpperRej);
    assert(pidSwap);
    inCol.iOrdSplit = pst->iOrdSplit;
    inCol.iSplitSide = 1;
    rID = mdl->ReqService(pst->idUpper,PST_COLORDREJECTS,&inCol,sizeof(inCol));
    inCol.iSplitSide = 0;
    nOut = Traverse(PST_COLORDREJECTS,pst->pstLower,&inCol,sizeof(inCol),pLowerRej.get(),
    	     	     mdl->Threads()*sizeof(ServiceColRejects::output));
    mdlassert(pst->mdl,nOut/sizeof(ServiceColRejects::output) == pst->nLower);
    nOut = mdl->GetReply(rID,pUpperRej.get());
    mdlassert(pst->mdl,nOut/sizeof(ServiceColRejects::output) == pst->nUpper);
    while (1) {
	iRet = RejMatch(pst,pst->nLower,pLowerRej.get(),pst->nUpper,pUpperRej.get(),pidSwap.get());
	if (!iRet) break;
	rID = mdl->ReqService(pst->idUpper,PST_SWAPREJECTS,pidSwap.get(),
		      mdl->Threads()*sizeof(int));
	nOut = Traverse(PST_SWAPREJECTS,pst->pstLower,pidSwap.get(),mdl->Threads()*sizeof(int),
		       pLowerRej.get(),mdl->Threads()*sizeof(ServiceSwapRejects::output));
	mdlassert(pst->mdl,nOut/sizeof(ServiceSwapRejects::output) == pst->nLower);
	nOut = mdl->GetReply(rID,pUpperRej.get());
	mdlassert(pst->mdl,nOut/sizeof(ServiceSwapRejects::output) == pst->nUpper);
	}
    return pst->iOrdSplit;
    }

/*****************************************************************************\
* ServiceLocalOrder
\*****************************************************************************/

static_assert(std::is_void<ServiceLocalOrder::input>()  || std::is_trivial<ServiceLocalOrder::input>());
static_assert(std::is_void<ServiceLocalOrder::output>() || std::is_trivial<ServiceLocalOrder::output>());

int ServiceLocalOrder::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto iMinOrder = in->iMinOrder;
    auto iMaxOrder = in->iMaxOrder;
    auto iMidOrder = pst->iOrdSplit;
    assert(iMidOrder >= iMinOrder && iMidOrder <= iMaxOrder);
    in->iMinOrder = iMidOrder;
    auto rID = mdl->ReqService(pst->idUpper,PST_LOCALORDER,in,nIn);
    in->iMinOrder = iMinOrder;
    in->iMaxOrder = iMidOrder-1;
    Traverse(pst->pstLower,in,nIn,NULL,0);
    mdl->GetReply(rID);
    return 0;
    }

int ServiceLocalOrder::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto pkd = pst->plcl->pkd;
    auto in = static_cast<input*>(vin);
    pkdLocalOrder(pkd,in->iMinOrder,in->iMaxOrder);
    return 0;
    }

/*****************************************************************************\
* ServiceWeight
* 
* Make sure that the local particles are split into active and inactive
* when passing pFlag != 0.
* pFlag == 0 => weight all particles.
* pFlag > 0 => weight active particles.
* pFlag < 0 => weight inactive particles.
\*****************************************************************************/

static_assert(std::is_void<ServiceWeight::input>()  || std::is_trivial<ServiceWeight::input>());
static_assert(std::is_void<ServiceWeight::output>() || std::is_trivial<ServiceWeight::output>());

int ServiceWeight::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto rID = mdl->ReqService(pst->idUpper,PST_WEIGHT,in,nIn);
    Traverse(pst->pstLower,in,nIn,out,nOut);
    output outWt;
    mdl->GetReply(rID,&outWt);
    out->nLow += outWt.nLow;
    out->nHigh += outWt.nHigh;
    out->fLow += outWt.fLow;
    out->fHigh += outWt.fHigh;
    return sizeof(output);
    }

int ServiceWeight::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto plcl = pst->plcl;
    auto pkd = plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);

    auto fSplit = in->fSplit;
    auto iSplitSide = in->iSplitSide;
    if (in->ittr == 0) {
	/*
	** Initialize.
	*/
	plcl->fSplit = fSplit;
	if (in->pFlag == 0) {
	    plcl->iWtFrom = 0;
	    plcl->iWtTo = pkdLocal(pkd)-1;
	    }
	else if (in->pFlag > 0) {
	    /*
	    ** Particles must be in the active-inactive order here!
	    */
	    plcl->iWtFrom = 0;
	    plcl->iWtTo = pkdActive(pkd)-1;
	    }
	else {
	    /*
	    ** Particles must be in the active-inactive order here!
	    */
	    plcl->iWtFrom = pkdActive(pkd);
	    plcl->iWtTo = pkdLocal(pkd)-1;
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

    int nLow,nHigh;
    double fLow,fHigh;
    plcl->iPart = pkdWeight(pkd,in->iSplitDim,fSplit,iSplitSide,
			    plcl->iWtFrom,plcl->iWtTo,
			    &nLow,&nHigh,&fLow,&fHigh);
    out->nLow = nLow;
    out->nHigh = nHigh;
    out->fLow = fLow + plcl->fWtLow;
    out->fHigh = fHigh + plcl->fWtHigh;
    plcl->fLow = fLow;
    plcl->fHigh = fHigh;
    if (in->pFlag > 0) {
	if (iSplitSide) out->nLow -= pkdInactive(pkd);
	else out->nHigh -= pkdInactive(pkd);
	}
    if (in->pFlag < 0) {
	if (iSplitSide) out->nHigh -= pkdActive(pkd);
	else out->nLow -= pkdActive(pkd);
	}

    return sizeof(output);
    }

/*****************************************************************************\
* ServiceWeightWrap
* 
* Make sure that the local particles are split into active and inactive
* when passing pFlag != 0.
* pFlag == 0 => weight all particles.
* pFlag > 0 => weight active particles.
* pFlag < 0 => weight inactive particles.
\*****************************************************************************/

static_assert(std::is_void<ServiceWeightWrap::input>()  || std::is_trivial<ServiceWeightWrap::input>());
static_assert(std::is_void<ServiceWeightWrap::output>() || std::is_trivial<ServiceWeightWrap::output>());

int ServiceWeightWrap::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto rID = mdl->ReqService(pst->idUpper,PST_WEIGHTWRAP,in,nIn);
    Traverse(pst->pstLower,in,nIn,out,nOut);
    output outWt;
    mdl->GetReply(rID,&outWt);
    out->nLow += outWt.nLow;
    out->nHigh += outWt.nHigh;
    return sizeof(output);
    }

int ServiceWeightWrap::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto plcl = pst->plcl;
    auto pkd = plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);

    if (in->ittr == 0) {
	/*
	** Initialize.
	*/
	plcl->fSplit = in->fSplit;
	plcl->iWtFrom = 0;
	plcl->iWtTo = pkdLocal(pkd)-1;
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
    int nLow,nHigh;
    plcl->iPart = pkdWeightWrap(pkd,in->iSplitDim,in->fSplit,in->fSplit2,in->iSplitSide,
				plcl->iWtFrom,plcl->iWtTo,&nLow,&nHigh);
    out->nLow = nLow;
    out->nHigh = nHigh;
    /* For collect rejects */
    plcl->nSplit = plcl->iPart;

    return sizeof(output);
    }

/*****************************************************************************\
* ServiceOrdWeight
\*****************************************************************************/

static_assert(std::is_void<ServiceOrdWeight::input>()  || std::is_trivial<ServiceOrdWeight::input>());
static_assert(std::is_void<ServiceOrdWeight::output>() || std::is_trivial<ServiceOrdWeight::output>());

int ServiceOrdWeight::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto pkd = pst->plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);
    auto rID = mdl->ReqService(pst->idUpper,PST_ORDWEIGHT,in,nIn);
    Traverse(pst->pstLower,in,nIn,out,nOut);
    output outWt;
    mdl->GetReply(rID,&outWt);
    out->nLow += outWt.nLow;
    out->nHigh += outWt.nHigh;
    return sizeof(output);
    }

int ServiceOrdWeight::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto plcl = pst->plcl;
    auto pkd = plcl->pkd;
    auto in = static_cast<input*>(vin);
    auto out = static_cast<output*>(vout);

    if (in->ittr == 0) {
	/*
	** Initialize.
	*/
	plcl->iOrdSplit = in->iOrdSplit;
	plcl->iWtFrom = 0;
	plcl->iWtTo = pkdLocal(pkd)-1;
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
    int nLow,nHigh;
    plcl->iPart = pkdOrdWeight(pkd,in->iOrdSplit,in->iSplitSide,
				plcl->iWtFrom,plcl->iWtTo,
				&nLow,&nHigh);
    out->nLow = nLow;
    out->nHigh = nHigh;
    return sizeof(output);
    }
} // namespace OldDD