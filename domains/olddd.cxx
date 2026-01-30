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
#include "orderinfo.h"
#include "freestore.h"
#include "pkd.h"
#include <cstdint>
#include <vector>

namespace OldDD {

constexpr double NEWSPLITDIMCUT = 0.707;
constexpr int NMINFORROOTFIND = 16;
constexpr int MAX_ITTR = 64;

// Stored the LCL registry for weight finding
struct WeightInfo {
    WeightInfo() = default;
    int iWtFrom{};
    int iWtTo{};
    int iPart{};
    int nSplit{};
    uint64_t iOrdSplit{};
    double fSplit{};
    double fWtLow{};
    double fWtHigh{};
    double fLow{};
    double fHigh{};
};

struct DomainInfo {
    double fSplit = 0.0;
    double fSplitInactive{};
    int iSplitDim = -1; // -1 means no split choosen yet
};

/*****************************************************************************\
* ServiceDomainDecomp
\*****************************************************************************/

static_assert(std::is_void_v<ServiceDomainDecomp::input> || std::is_standard_layout_v<ServiceDomainDecomp::input>);
static_assert(std::is_void_v<ServiceDomainDecomp::output> || std::is_standard_layout_v<ServiceDomainDecomp::output>);

int ServiceDomainDecomp::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto in = reinterpret_cast<input*>(vin.data());
    static_assert(std::is_void_v<output>);
    assert(vin.size() == sizeof(input));
    auto& dinfo = pst->data_ensure<DomainInfo>();
    auto bnd = pst->data_emplace<Bound>(in->bnd);

    int /*nBndWrapd,*/ d = 0;
    mdlTimer t;

#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Domain Decomposition");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
    if (dinfo.iSplitDim != -1 && in->nActive < NMINFORROOTFIND) {
        mdlprintf(pst->mdl(), "Aborting RootFind -- too few actives.\n");
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
        d = bnd.maxdim(); // Find the longest axis
        // If we already have a split, and the new size isn't significantly different
        // then just keep the existing split dimension (so less data movement)
        if (dinfo.iSplitDim >= 0 && bnd.width(d) < bnd.width(dinfo.iSplitDim) * NEWSPLITDIMCUT) d = dinfo.iSplitDim;
#if 0
        d = dinfo.iSplitDim;
        double dimsize;
        if (d==-1) {
            dimsize = -1;
            nBndWrapd = 0;
        }
        else {
            dimsize = bnd.width(d)*NEWSPLITDIMCUT;
            nBndWrapd = in->nBndWrap[d];
        }

        for (auto j=0; j<3; ++j) {
            if (in->nBndWrap[j] < nBndWrapd || bnd.width(j) > dimsize) {
                d=j;
                dimsize = bnd.width(d);
                nBndWrapd = in->nBndWrap[d];
            }
        }
#endif
    }
    else d = dinfo.iSplitDim;
    assert(d >= 0);

    mdlPrintTimer(pst->mdl(), "TIME Mass Check done in pstDomainDecomp", &t);
    RootSplit(pst, d, in->bDoRootFind, in->bDoSplitDimFind);
    mdlPrintTimer(pst->mdl(), "TIME RootSplit done in pstDomainDecomp", &t);

    mdlPrintTimer(pst->mdl(), "TIME Mass Check done in pstDomainDecomp", &t);
    /*
    ** Now go on to DD of next levels, but pass correct wrapping bounds.
    */
    auto l = bnd.lower();
    auto u = bnd.upper();
    // in->nBndWrap[d] = nBndWrapd;
    if (dinfo.fSplitInactive <= l[d] || dinfo.fSplitInactive >= u[d]) {
        l[d] = dinfo.fSplit;
    }
    else if (dinfo.fSplitInactive > dinfo.fSplit) {
        l[d] = dinfo.fSplit;
        u[d] = dinfo.fSplitInactive;
    }
    else in->nBndWrap[d]++;

    in->bnd = Bound(l, u);
    auto rID = mdl->ReqService<OldDD::ServiceDomainDecomp>(pst->upper(), vin);

    l = bnd.lower();
    u = bnd.upper();
    // in->nBndWrap[d] = nBndWrapd;
    if (dinfo.fSplitInactive <= l[d] || dinfo.fSplitInactive >= u[d]) {
        u[d] = dinfo.fSplit;
    }
    else if (dinfo.fSplitInactive < dinfo.fSplit) {
        u[d] = dinfo.fSplit;
        l[d] = dinfo.fSplitInactive;
    }
    else in->nBndWrap[d]++;

    in->bnd = Bound(l, u);
    Traverse(pst->lower(), vin, vout);

    GetReply(rID);
#ifdef USE_ITT
    __itt_task_end(domain);
#endif

    return 0;
}

int ServiceDomainDecomp::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    static_assert(std::is_void_v<output>);
    auto pkd = pst->lcl()->pkd;

    pkd->bnd = in->bnd;

    float offs = 0.5f / (pkd->Local() * 1.0f - 1.0f);
    pst->data_emplace<Bound>(Bound(in->bnd.lower() - offs, in->bnd.upper() + offs));

    return 0;
}

void ServiceDomainDecomp::RootSplit(PST pst, int iSplitDim, int bDoRootFind, int bDoSplitDimFind) {
    auto mdl = pst->mdl();
    auto& dinfo = pst->data_ensure<DomainInfo>();
    int NUM_SAFETY = 4; /* slop space when filling up memory */
    uint64_t nSafeTot;  /* total slop space we have to play with */
    uint64_t margin;    /* more slop */
    int d, ittr, nOut;
    /*
    ** Why are these initialized to -1 here???
    int nLow=-1,nHigh=-1;
    */
    uint64_t nLow, nHigh;
    uint64_t nLowerStore, nUpperStore;
    uint64_t nLowTot, nHighTot;
    uint64_t nLast; /* number of particles at the last split iteration */
    uint64_t nTotalActive;
    int nDiff = 0; /* Difference between one iteration and the next, seems to only be used to warn. */
    double fLow, fHigh;
    double fl, fu, fm = -1, fmm;
    ServiceFreeStore::output outFree;
    ServiceWeight::output outWtLow, outWtHigh;
    ServiceWeightWrap::output outWtWrLow, outWtWrHigh;
    int iRet;
    char ach[256]; /* Debug */
    mdlTimer t;
    int pFlag; /* 0 => we are splitting all particles at once. 1 => we first split active, and then inactive. */
    int dBnd;
    int rID;

    mdlZeroTimer(pst->mdl(), &t);
    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    rID = mdl->ReqService<ServiceFreeStore>(pst->upper());
    Traverse<ServiceFreeStore>(pst->lower(), {}, outFree);
    nLowerStore = outFree;
    GetReply(rID, outFree);
    nUpperStore = outFree;

    mdlprintf(pst->mdl(), "_pstRootSplit: id %d Level %d\n", pst->Self(), pst->depth());
    mdlPrintTimer(pst->mdl(), "TIME START _pstRootSplit ", &t);
    mdlprintf(pst->mdl(), "_pstRootSplit: fA0 %f fIA0 %f RS? %d DC? %d\n", dinfo.fSplit, dinfo.fSplitInactive, bDoRootFind, bDoSplitDimFind);
    /* Debug */
    /*
      sprintf(ach,"id: %d _pstRootSplit\n", pst->Self() );
      mdlDiag(pst->mdl(),ach);
    */

    if (bDoSplitDimFind || dinfo.iSplitDim == -1) {
        dinfo.iSplitDim = iSplitDim;
    }

    auto bnd = pst->data<Bound>();
    d = dBnd = dinfo.iSplitDim;
    fl = bnd.lower(dBnd);
    fu = bnd.upper(dBnd);
    fm = dinfo.fSplit;
    ittr = -1;

    if (bDoRootFind || fm < fl || fm > fu) {
        /*
        ** First order the particles into active/inactive order...
        */
        uint64_t nActiveOrder;
        Traverse<ServiceActiveOrder>(pst, {}, nActiveOrder);

        fmm = (fl + fu) / 2;
        /*
         * First find total number of active particles.
         */
        ServiceWeight::input inWt;
        inWt.iSplitDim = d;
        inWt.fSplit = fmm;
        inWt.ittr = 0;
        inWt.iSplitSide = 1;
        inWt.pFlag = 1;
        rID = mdl->ReqService<ServiceWeight>(pst->upper(), inWt);
        inWt.iSplitSide = 0;
        Traverse<ServiceWeight>(pst->lower(), inWt, outWtLow);
        GetReply(rID, outWtHigh);
        nTotalActive = outWtLow.nLow + outWtHigh.nLow + outWtLow.nHigh + outWtHigh.nHigh;
        mdlassert(pst->mdl(), nActiveOrder == nTotalActive);
        pFlag = 1;
        if (nTotalActive <= 1) {
            pFlag = 0; /* Divide them all */
            rID = mdl->ReqService<ServiceWeight>(pst->upper(), inWt);
            inWt.iSplitSide = 0;
            Traverse<ServiceWeight>(pst->lower(), inWt, outWtLow);
            GetReply(rID, outWtHigh);
            nTotalActive = outWtLow.nLow + outWtHigh.nLow + outWtLow.nHigh + outWtHigh.nHigh;
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
            rID = mdl->ReqService<ServiceWeight>(pst->upper(), inWt);
            inWt.iSplitSide = 0;
            Traverse<ServiceWeight>(pst->lower(), inWt, outWtLow);
            GetReply(rID, outWtHigh);
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
            if (pFlag) { /* split on work */
                if (fLow / pst->n_lower() > fHigh / pst->n_upper()) fu = fm;
                else if (fLow / pst->n_lower() < fHigh / pst->n_upper()) fl = fm;
                else break;
            }
            else { /* split on number */
                if (nLow / (double)pst->n_lower() > nHigh / (double)pst->n_upper()) fu = fm;
                else if (nLow / (double)pst->n_lower() < nHigh / (double)pst->n_upper()) fl = fm;
                else break;
            }
            fmm = (fl + fu) / 2;
            ++ittr;
        }
        mdlPrintTimer(pst->mdl(), "TIME active split _pstRootSplit ", &t);
    }

    dinfo.fSplit = fm;

    mdlprintf(pst->mdl(), "id: %d (%d) Chose split: %f (%f,%f) %d %d\n", pst->Self(), pst->depth(), fm, bnd.lower(dBnd), bnd.upper(dBnd), pst->n_lower(),
              pst->n_upper());
    if (ittr != -1)
        mdlprintf(pst->mdl(), "  Low %" PRIu64 " %f,  High %" PRIu64 " %f, ittr=%d\n", nLow, outWtLow.fLow + outWtHigh.fLow, nHigh,
                  outWtLow.fHigh + outWtHigh.fHigh, ittr);
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
    fl = dinfo.fSplit + 1e-6 * bnd.width(dBnd);
    fu = dinfo.fSplit - 1e-6 * bnd.width(dBnd);

    if (!bDoSplitDimFind) fm = dinfo.fSplitInactive;
    else {
        fm = 0.5 * (fl + fu);
        if (fm < bnd.center(dBnd)) fm = bnd.upper(dBnd) + 1.000001 * bnd.width(dBnd);
        else fm = bnd.lower(dBnd) - 1.000001 * bnd.width(dBnd);
    }
    mdlprintf(pst->mdl(), "id: %d (%d) Zeroeth guess reverse split: %f (%f,%f)\n", pst->Self(), pst->depth(), fm, bnd.lower(dBnd), bnd.upper(dBnd));
    inWtWrap.fSplit = fm;
    inWtWrap.fSplit2 = dinfo.fSplit;
    inWtWrap.ittr = 0;
    inWtWrap.iSplitSide = 1;
    rID = mdl->ReqService<ServiceWeightWrap>(pst->upper(), inWtWrap);
    inWtWrap.iSplitSide = 0;
    Traverse<ServiceWeightWrap>(pst->lower(), inWtWrap, outWtWrLow);
    GetReply(rID, outWtWrHigh);
    /*
    ** Add lower and Upper subsets numbers of particles
    */
    nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
    nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;

    nSafeTot = nLowerStore + nUpperStore - (nLowTot + nHighTot);
    if (nSafeTot / pst->n_total() < NUM_SAFETY) {
        NUM_SAFETY = nSafeTot / pst->n_total();
        snprintf(ach, sizeof(ach), "id: %d tripped inactive NUM_SAFETY %d  Low %" PRIu64 "/%" PRIu64 "  High %" PRIu64 "/%" PRIu64 "\n", pst->Self(),
                 NUM_SAFETY, nLowTot, nLowerStore, nHighTot, nUpperStore);
        mdlDiag(pst->mdl(), ach);
        mdlprintf(pst->mdl(), "id: %d tripped inactive NUM_SAFETY %d  Low %%" PRIu64 "/%" PRIu64 "  High %" PRIu64 "/%" PRIu64 "\n", pst->Self(), NUM_SAFETY,
                  nLowTot, nLowerStore, nHighTot, nUpperStore);
    }

    margin = nSafeTot / pst->n_total() / 20;
    if (margin < NUM_SAFETY / 2) margin = NUM_SAFETY / 2;

    mdlprintf(pst->mdl(), "id: %d  %d Low %" PRIu64 "/%" PRIu64 "   %d High %" PRIu64 "/%" PRIu64 "  NUM_SAFETY %d margin %d\n", pst->Self(), pst->n_lower(),
              nLowTot, nLowerStore, pst->n_upper(), nHighTot, nUpperStore, NUM_SAFETY, margin);

    if (nLowTot > nLowerStore - NUM_SAFETY * pst->n_lower()) {
        snprintf(ach, sizeof(ach), "id: %d: nLowTot > nLowerStore-NUM_SAFETY*pst->n_lower() %" PRIu64 " %" PRIu64 " %d %d\n", pst->Self(), nLowTot, nLowerStore,
                 NUM_SAFETY, pst->n_lower());
        mdlDiag(pst->mdl(), ach);
        fm = std::max(std::min(fm, bnd.upper(dBnd)), bnd.lower(dBnd));
        fl = fm;

        if (fu > fl) fmm = 0.5 * (fl + fu);
        else {
            fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
            if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
            mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
        }
        ittr = 1;
        nLast = nLowTot;
        while (ittr < MAX_ITTR) {
            fm = fmm;
            inWtWrap.iSplitDim = d;
            inWtWrap.fSplit = fm;
            inWtWrap.fSplit2 = dinfo.fSplit;
            inWtWrap.ittr = ittr;
            inWtWrap.iSplitSide = 1;
            rID = mdl->ReqService<ServiceWeightWrap>(pst->upper(), inWtWrap);
            inWtWrap.iSplitSide = 0;
            Traverse<ServiceWeightWrap>(pst->lower(), inWtWrap, outWtWrLow);
            GetReply(rID, outWtWrHigh);
            /*
            ** Add lower and Upper subsets numbers of particles
            */
            nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
            nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
            if (nLowTot != nLast) nDiff = nLowTot - nLast;
            nLast = nLowTot;
            /*
              if (nLowTot > nLowerStore-NUM_SAFETY*pst->n_lower()) fl = fm;
              else if (nLowTot < nLowerStore-2*NUM_SAFETY*pst->n_lower()) fu = fm;
            */
            /*
              if (nLowTot/pst->n_lower() > nHighTot/pst->n_upper()) fl = fm;
              else if (nLowTot/pst->n_lower() < nHighTot/pst->n_upper()-NUM_SAFETY) fu = fm;
            */
            if (nLowTot > nLowerStore - margin * pst->n_lower()) fl = fm;
            else if (nLowTot < nLowerStore - 2 * margin * pst->n_lower()) fu = fm;
            else {
                fl = fm;
                break;
            }
            if (fu == fl) break;
            else if (fu > fl) fmm = 0.5 * (fl + fu);
            else {
                fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
                if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
                mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
            }
            ++ittr;
        }
        mdlprintf(pst->mdl(), "id: %d (%d) Fix Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n", pst->Self(), pst->depth(),
                  ittr, fm, fl, fu, bnd.lower(dBnd), bnd.upper(dBnd), nLowTot, nHighTot);
        if (nLowTot != nLowerStore - NUM_SAFETY * pst->n_lower()) {
            if (abs(nDiff) > 1) mdlprintf(pst->mdl(), "id: %d delta of %d, check NUM_SAFTEY\n", pst->Self(), nDiff);
        }
        mdlassert(pst->mdl(), nLowTot <= nLowerStore);
        mdlPrintTimer(pst->mdl(), "TIME fix lower II _pstRootSplit ", &t);
    }
    else if (nHighTot > nUpperStore - NUM_SAFETY * pst->n_upper()) {
        snprintf(ach, sizeof(ach), "id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->n_upper() %" PRIu64 " %" PRIu64 " %d %d\n", pst->Self(), nHighTot,
                 nUpperStore, NUM_SAFETY, pst->n_upper());
        mdlDiag(pst->mdl(), ach);
        fm = std::max(std::min(fm, bnd.upper(dBnd)), bnd.lower(dBnd));
        fu = fm;
        if (fu > fl) fmm = 0.5 * (fl + fu);
        else {
            fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
            if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
            mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
        }
        ittr = 1;
        nLast = nLowTot;
        while (ittr < MAX_ITTR) {
            fm = fmm;
            inWtWrap.iSplitDim = d;
            inWtWrap.fSplit = fm;
            inWtWrap.fSplit2 = dinfo.fSplit;
            inWtWrap.ittr = ittr;
            inWtWrap.iSplitSide = 1;
            rID = mdl->ReqService<ServiceWeightWrap>(pst->upper(), inWtWrap);
            inWtWrap.iSplitSide = 0;
            Traverse<ServiceWeightWrap>(pst->lower(), inWtWrap, outWtWrLow);
            GetReply(rID, outWtWrHigh);
            /*
            ** Add lower and Upper subsets numbers of particles
            */
            nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
            nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
            if (nLowTot != nLast) nDiff = nLowTot - nLast;
            nLast = nLowTot;
            if (nHighTot > nUpperStore - margin * pst->n_upper()) fu = fm;
            else if (nHighTot < nUpperStore - 2 * margin * pst->n_upper()) fl = fm;
            else {
                fu = fm;
                break;
            }
            if (fu == fl) break;
            else if (fu > fl) fmm = 0.5 * (fl + fu);
            else {
                fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
                if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
                mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
            }
            ++ittr;
        }
        mdlprintf(pst->mdl(), "id: %d (%d) Fix High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n", pst->Self(),
                  pst->depth(), ittr, fm, fl, fu, bnd.lower(dBnd), bnd.upper(dBnd), nLowTot, nHighTot);
        if (nHighTot != nUpperStore - NUM_SAFETY * pst->n_upper()) {
            if (abs(nDiff) > 1) mdlprintf(pst->mdl(), "id: %d delta of %d, check NUM_SAFETY\n", pst->Self(), nDiff);
        }
        mdlassert(pst->mdl(), nHighTot <= nUpperStore);
        mdlPrintTimer(pst->mdl(), "TIME fix upper II _pstRootSplit ", &t);
    }

    if (nLowTot < NUM_SAFETY * pst->n_lower()) {
        snprintf(ach, sizeof(ach), "id: %d: nLowTot < NUM_SAFETY*pst->n_lower() %" PRIu64 " %" PRIu64 " %d %d\n", pst->Self(), nLowTot, nLowerStore, NUM_SAFETY,
                 pst->n_lower());
        mdlDiag(pst->mdl(), ach);
        fm = std::max(std::min(fm, bnd.upper(dBnd)), bnd.lower(dBnd));
        fu = fm;
        if (fu > fl) fmm = 0.5 * (fl + fu);
        else {
            fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
            if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
            mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
        }
        ittr = 1;
        nLast = nLowTot;
        while (ittr < MAX_ITTR) {
            fm = fmm;
            inWtWrap.iSplitDim = d;
            inWtWrap.fSplit = fm;
            inWtWrap.fSplit2 = dinfo.fSplit;
            inWtWrap.ittr = ittr;
            inWtWrap.iSplitSide = 1;
            rID = mdl->ReqService<ServiceWeightWrap>(pst->upper(), inWtWrap);
            inWtWrap.iSplitSide = 0;
            Traverse<ServiceWeightWrap>(pst->lower(), inWtWrap, outWtWrLow);
            GetReply(rID, outWtWrHigh);
            /*
            ** Add lower and Upper subsets numbers of particles
            */
            nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
            nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
            if (nLowTot != nLast) nDiff = nLowTot - nLast;
            nLast = nLowTot;
            if (nLowTot > margin * pst->n_lower()) fl = fm;
            else if (nLowTot < NUM_SAFETY * pst->n_lower()) fu = fm;
            else {
                fl = fm;
                break;
            }
            if (fu == fl) break;
            else if (fu > fl) fmm = 0.5 * (fl + fu);
            else {
                fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
                if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
                mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
            }
            ++ittr;
        }
        mdlprintf(pst->mdl(), "id: %d (%d) Fix too few Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n", pst->Self(),
                  pst->depth(), ittr, fm, fl, fu, bnd.lower(dBnd), bnd.upper(dBnd), nLowTot, nHighTot);
        if (nLowTot != nLowerStore - NUM_SAFETY * pst->n_lower()) {
            if (abs(nDiff) > 1) mdlprintf(pst->mdl(), "id: %d delta of %d, check NUM_SAFTEY\n", pst->Self(), nDiff);
        }
        mdlassert(pst->mdl(), nLowTot <= nLowerStore);
        mdlPrintTimer(pst->mdl(), "TIME fix lower II _pstRootSplit ", &t);
    }
    if (nHighTot < NUM_SAFETY * pst->n_upper()) {
        snprintf(ach, sizeof(ach), "id: %d: nHighTot < NUM_SAFETY*pst->n_upper() %" PRIu64 " %" PRIu64 " %d %d\n", pst->Self(), nHighTot, nUpperStore,
                 NUM_SAFETY, pst->n_upper());
        mdlDiag(pst->mdl(), ach);
        fm = std::max(std::min(fm, bnd.upper(dBnd)), bnd.lower(dBnd));
        fl = fm;
        if (fu > fl) fmm = 0.5 * (fl + fu);
        else {
            fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
            if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
            mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
        }
        ittr = 1;
        nLast = nLowTot;
        while (ittr < MAX_ITTR) {
            fm = fmm;
            inWtWrap.iSplitDim = d;
            inWtWrap.fSplit = fm;
            inWtWrap.fSplit2 = dinfo.fSplit;
            inWtWrap.ittr = ittr;
            inWtWrap.iSplitSide = 1;
            rID = mdl->ReqService<ServiceWeightWrap>(pst->upper(), inWtWrap);
            inWtWrap.iSplitSide = 0;
            Traverse<ServiceWeightWrap>(pst->lower(), inWtWrap, outWtWrLow);
            GetReply(rID, outWtWrHigh);
            /*
            ** Add lower and Upper subsets numbers of particles
            */
            nLowTot = outWtWrLow.nLow + outWtWrHigh.nLow;
            nHighTot = outWtWrLow.nHigh + outWtWrHigh.nHigh;
            if (nLowTot != nLast) nDiff = nLowTot - nLast;
            nLast = nLowTot;
            if (nHighTot > margin * pst->n_upper()) fu = fm;
            else if (nHighTot < NUM_SAFETY * pst->n_upper()) fl = fm;
            else {
                fu = fm;
                break;
            }
            if (fu == fl) break;
            else if (fu > fl) fmm = 0.5 * (fl + fu);
            else {
                fmm = 0.5 * (fl + fu) + bnd.apothem(dBnd);
                if (fmm > bnd.upper(dBnd)) fmm = 0.5 * (fl + fu) - bnd.apothem(dBnd);
                mdlassert(pst->mdl(), fmm >= bnd.lower(dBnd) && fmm <= bnd.upper(dBnd));
            }
            ++ittr;
        }
        mdlprintf(pst->mdl(), "id: %d (%d) Fix Too few High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %" PRIu64 " High %" PRIu64 "\n", pst->Self(),
                  pst->depth(), ittr, fm, fl, fu, bnd.lower(dBnd), bnd.upper(dBnd), nLowTot, nHighTot);
        if (nHighTot != nUpperStore - NUM_SAFETY * pst->n_upper()) {
            if (abs(nDiff) > 1) mdlprintf(pst->mdl(), "id: %d delta of %d, check NUM_SAFETY\n", pst->Self(), nDiff);
        }
        mdlassert(pst->mdl(), nHighTot <= nUpperStore);
        mdlPrintTimer(pst->mdl(), "TIME fix upper II _pstRootSplit ", &t);
    }

    mdlassert(pst->mdl(), nLowTot >= pst->n_lower());
    mdlassert(pst->mdl(), nHighTot >= pst->n_upper());
    mdlassert(pst->mdl(), nLowTot <= nLowerStore);
    mdlassert(pst->mdl(), nHighTot <= nUpperStore);

    mdlPrintTimer(pst->mdl(), "TIME Total Split _pstRootSplit ", &t);
    mdlprintf(pst->mdl(), "id: %d (%d) Chose reverse split: %f (%f,%f)\n", pst->Self(), pst->depth(), fm, bnd.lower(dBnd), bnd.upper(dBnd));
    dinfo.fSplitInactive = fm;

    std::vector<ServiceColRejects::output> lowerRej(pst->n_lower());
    std::vector<ServiceColRejects::output> upperRej(pst->n_upper());
    std::vector<int> pidSwap(mdl->Threads());

    rID = mdl->ReqService<ServiceColRejects>(pst->upper());
    nOut = Traverse<ServiceColRejects>(pst->lower(), {}, lowerRej);
    assert(nOut / sizeof(ServiceColRejects::output) == pst->n_lower());
    nOut = GetReply(rID, upperRej);
    assert(nOut / sizeof(ServiceColRejects::output) == pst->n_upper());

    mdlPrintTimer(pst->mdl(), "TIME Collected Rejects _pstRootSplit ", &t);

    ittr = 0;
    while (1) {
        iRet = RejMatch(pst, pst->n_lower(), lowerRej.data(), pst->n_upper(), upperRej.data(), pidSwap.data());
        if (!iRet) break;
        rID = mdl->ReqService<ServiceSwapRejects>(pst->upper(), pidSwap);
        nOut = Traverse<ServiceSwapRejects>(pst->lower(), pidSwap, lowerRej);
        assert(nOut / sizeof(ServiceSwapRejects::output) == pst->n_lower());
        nOut = GetReply(rID, upperRej);
        assert(nOut / sizeof(ServiceSwapRejects::output) == pst->n_upper());

        ++ittr;
    }

    mdlPrintTimer(pst->mdl(), "TIME (FINISH) Swapped Rejects _pstRootSplit ", &t);
}

/*****************************************************************************\
* ServiceDomain
\*****************************************************************************/

int ServiceDomain::RejMatch(PST pst, int n1, ServiceColRejects::output* p1, int n2, ServiceColRejects::output* p2, int* pidSwap) {
    auto mdl = pst->mdl();
    int id, i, i1 = -1, i2 = -1, id1, id2;
    local_t nLarge;
    total_t s1, s2, r1, r2;

    /*
    ** Check to see if there is enough space...
    */
    s1 = 0;
    r1 = 0;
    for (i = 0; i < n1; ++i) {
        s1 += p1[i].nSpace;
        r1 += p1[i].nRejects;
    }
    s2 = 0;
    r2 = 0;
    for (i = 0; i < n2; ++i) {
        s2 += p2[i].nSpace;
        r2 += p2[i].nRejects;
    }
    auto pkd = pst->lcl()->pkd;
    // printf("%d\t %d \t %d\n", pkd->Self(), r1, s2);
    // printf("%d\t %d \t %d\n", pkd->Self(), r2, s1);
    assert(r1 <= s2);
    assert(r2 <= s1);
    /*
    ** And also check for decorators
    */
    s1 = 0;
    r1 = 0;
    for (i = 0; i < n1; ++i) {
        s1 += p1[i].nDecoratorSpace;
        r1 += p1[i].nDecoratorRejects;
    }
    s2 = 0;
    r2 = 0;
    for (i = 0; i < n2; ++i) {
        s2 += p2[i].nDecoratorSpace;
        r2 += p2[i].nDecoratorRejects;
    }
    assert(r1 <= s2);
    assert(r2 <= s1);
    /*
    ** First invalidate the pidSwap array.
    */
    for (id = 0; id < mdl->Threads(); ++id) pidSwap[id] = -1;
    /*
    ** Now map largest nReject of p1 to largest nSpace of p2.
    */
    while (1) {
        nLarge = 0;
        for (i = 0; i < n1; ++i) {
            auto maxRejects = std::max(p1[i].nRejects, p1[i].nDecoratorRejects);
            if (maxRejects > nLarge) {
                nLarge = maxRejects;
                i1 = i;
            }
        }
        if (nLarge == 0) break;
        nLarge = 0;
        for (i = 0; i < n2; ++i) {
            if (p2[i].nSpace > nLarge) {
                nLarge = p2[i].nSpace;
                i2 = i;
            }
        }
        if (nLarge == 0) break;
        p1[i1].nRejects = 0;
        p1[i1].nDecoratorRejects = 0;
        p1[i1].nSpace = 0;
        p2[i2].nRejects = 0;
        p2[i2].nDecoratorRejects = 0;
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
        for (i = 0; i < n2; ++i) {
            auto maxRejects = std::max(p2[i].nRejects, p2[i].nDecoratorRejects);
            if (maxRejects > nLarge) {
                nLarge = maxRejects;
                i2 = i;
            }
        }
        if (nLarge == 0) break;
        nLarge = 0;
        for (i = 0; i < n1; ++i) {
            if (p1[i].nSpace > nLarge) {
                nLarge = p1[i].nSpace;
                i1 = i;
            }
        }
        if (nLarge == 0) break;
        p1[i1].nRejects = 0;
        p1[i1].nDecoratorRejects = 0;
        p1[i1].nSpace = 0;
        p2[i2].nRejects = 0;
        p2[i2].nDecoratorRejects = 0;
        p2[i2].nSpace = 0;
        id1 = p1[i1].id;
        id2 = p2[i2].id;
        pidSwap[id1] = id2;
        pidSwap[id2] = id1;
    }
    for (i = 0; i < mdl->Threads(); ++i)
        if (pidSwap[i] != -1) return (1);
    return (0);
}

/*****************************************************************************\
* ServiceColRejects
\*****************************************************************************/

static_assert(std::is_void_v<ServiceColRejects::input> || std::is_trivial_v<ServiceColRejects::input>);
static_assert(std::is_void_v<ServiceColRejects::output> || std::is_trivial_v<ServiceColRejects::output>);

int ServiceColRejects::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto rejects = span_cast<output>(vout);
    auto lowerRej = rejects.subspan(0, pst->n_lower());
    auto upperRej = rejects.subspan(pst->n_lower());
    auto rID = mdl->ReqService<ServiceColRejects>(pst->upper(), vin);
    auto nLower = Traverse<ServiceColRejects>(pst->lower(), vin, lowerRej);
    auto nUpper = GetReply(rID, upperRej);
    return nLower + nUpper;
}

int ServiceColRejects::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto pkd = pst->lcl()->pkd;
    auto pOutRej = reinterpret_cast<output*>(vout.data());
    auto& winfo = pst->lcl()->data<WeightInfo>();
    std::tie(pOutRej->nRejects, pOutRej->nDecoratorRejects) = pkdColRejects(pkd, winfo.nSplit);
    pOutRej->nSpace = pkdSwapSpace(pkd);
    pOutRej->nDecoratorSpace = pkdSwapDecoratorSpace(pkd);
    pOutRej->id = pst->Self();
    pOutRej->nLocal = pkd->Local();
    return sizeof(output);
}

/*****************************************************************************\
* ServiceSwapRejects
\*****************************************************************************/

static_assert(std::is_void_v<ServiceSwapRejects::input> || std::is_trivial_v<ServiceSwapRejects::input>);
static_assert(std::is_void_v<ServiceSwapRejects::output> || std::is_trivial_v<ServiceSwapRejects::output>);

int ServiceSwapRejects::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto rejects = span_cast<output>(vout);
    auto lowerRej = rejects.subspan(0, pst->n_lower());
    auto upperRej = rejects.subspan(pst->n_lower());
    auto rID = mdl->ReqService<ServiceSwapRejects>(pst->upper(), vin);
    auto nLower = Traverse<ServiceSwapRejects>(pst->lower(), vin, lowerRej);
    auto nUpper = GetReply(rID, upperRej);
    return nLower + nUpper;
}

int ServiceSwapRejects::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto pkd = pst->lcl()->pkd;
    auto pidSwap = reinterpret_cast<const input*>(vin.data());
    auto pOutRej = reinterpret_cast<output*>(vout.data());
    auto idSwap = pidSwap[pst->Self()];
    std::tie(pOutRej->nRejects, pOutRej->nDecoratorRejects) = pkdSwapRejects(pkd, idSwap);
    pOutRej->nSpace = pkdSwapSpace(pkd);
    pOutRej->nDecoratorSpace = pkdSwapDecoratorSpace(pkd);
    pOutRej->id = pst->Self();
    return sizeof(output);
}

/*****************************************************************************\
* ServiceColOrdRejects
\*****************************************************************************/

static_assert(std::is_void_v<ServiceColOrdRejects::input> || std::is_trivial_v<ServiceColOrdRejects::input>);
static_assert(std::is_void_v<ServiceColOrdRejects::output> || std::is_trivial_v<ServiceColOrdRejects::output>);

int ServiceColOrdRejects::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto rejects = span_cast<output>(vout);
    auto lowerRej = rejects.subspan(0, pst->n_lower());
    auto upperRej = rejects.subspan(pst->n_lower());
    auto rID = mdl->ReqService<ServiceColOrdRejects>(pst->upper(), vin);
    auto nLower = Traverse<ServiceColOrdRejects>(pst->lower(), vin, lowerRej);
    auto nUpper = GetReply(rID, upperRej);
    return nLower + nUpper;
}

int ServiceColOrdRejects::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto pkd = pst->lcl()->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pOutRej = reinterpret_cast<output*>(vout.data());
    std::tie(pOutRej->nRejects, pOutRej->nDecoratorRejects) = pkdColOrdRejects(pkd, in->iOrdSplit, in->iSplitSide);
    pOutRej->nSpace = pkdSwapSpace(pkd);
    pOutRej->nDecoratorSpace = pkdSwapDecoratorSpace(pkd);
    pOutRej->id = pst->Self();
    return sizeof(output);
}

#ifndef NEW_REORDER
/*****************************************************************************\
* ServiceDomainOrder
\*****************************************************************************/

static_assert(std::is_void_v<ServiceDomainOrder::input> || std::is_trivial_v<ServiceDomainOrder::input>);
static_assert(std::is_void_v<ServiceDomainOrder::output> || std::is_trivial_v<ServiceDomainOrder::output>);

int ServiceDomainOrder::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto in = reinterpret_cast<input*>(vin.data());
    int rID;
    uint64_t iMinOrder = in->iMinOrder;
    uint64_t iMaxOrder = in->iMaxOrder;
    uint64_t iMidOrder;
    iMidOrder = OrdSplit(pst, iMinOrder, iMaxOrder);
    /*
    ** Now go on to Domain Order of next levels.
    */
    in->iMinOrder = iMidOrder;
    rID = mdl->ReqService<ServiceDomainOrder>(pst->upper(), vin);
    in->iMinOrder = iMinOrder;
    in->iMaxOrder = iMidOrder - 1;
    Traverse<ServiceDomainOrder>(pst->lower(), vin, vout);
    GetReply(rID);

    return 0;
}

int ServiceDomainOrder::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto in = reinterpret_cast<const input*>(vin.data());
    pst->data_emplace<OrderInfoPST>(OrderInfoPST{.iOrdSplit = in->iMinOrder});
    return 0;
}

uint64_t ServiceDomainOrder::OrdSplit(PST pst, uint64_t iMinOrder, uint64_t iMaxOrder) {
    auto mdl = pst->mdl();
    ServiceOrdWeight::input inWt;
    ServiceOrdWeight::output outWtLow, outWtHigh;
    uint64_t im, imm, il, iu;
    uint64_t nLowerStore, nUpperStore, nLow, nHigh;
    ServiceColOrdRejects::input inCol;
    int iRet, nOut, ittr;
    int rID;
    auto& winfo = pst->data_ensure<OrderInfoPST>();

    /*
    ** First find out how much free storage there is available for particles
    ** on the lower and upper subset of processors.
    */
    ServiceFreeStore::output outFree;
    rID = mdl->ReqService<ServiceFreeStore>(pst->upper());
    Traverse<ServiceFreeStore>(pst->lower(), {}, outFree);
    nLowerStore = outFree;
    GetReply(rID, outFree);
    nUpperStore = outFree;
    /*
    ** Find the correct iOrdSplit, such that all processors will
    ** have close to the same number of particles in the end.
    ** Start the ROOT finder based on balancing number of particles.
    */
    il = iMinOrder;
    iu = iMaxOrder;
    im = 0;         /* just initialization */
    nLow = 0;       /* just initialization */
    nHigh = iu + 1; /* just initialization */
    imm = (il + iu + 1) / 2;
    ittr = 0;
    while (il < imm && imm < iu && ittr < MAX_ITTR) {
        im = imm;
        inWt.iOrdSplit = im;
        inWt.ittr = ittr;
        inWt.iSplitSide = 1;
        rID = mdl->ReqService<ServiceOrdWeight>(pst->upper(), inWt);
        inWt.iSplitSide = 0;
        Traverse<ServiceOrdWeight>(pst->lower(), inWt, outWtLow);
        GetReply(rID, outWtHigh);
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
        else { /* split on number */
            if (nLow / (double)pst->n_lower() > nHigh / (double)pst->n_upper()) iu = im;
            else if (nLow / (double)pst->n_lower() < nHigh / (double)pst->n_upper()) il = im;
            else break;
        }
        imm = (il + iu) / 2;
        ++ittr;
    }
    mdlassert(pst->mdl(), nLow <= nLowerStore);
    mdlassert(pst->mdl(), nHigh <= nUpperStore);
    winfo.iOrdSplit = im;
    /*
    ** Collect rejects.
    */
    std::vector<ServiceColRejects::output> lowerRej(pst->n_lower());
    std::vector<ServiceColRejects::output> upperRej(pst->n_upper());
    std::vector<int> pidSwap(mdl->Threads());
    inCol.iOrdSplit = winfo.iOrdSplit;
    inCol.iSplitSide = 1;
    rID = mdl->ReqService<ServiceColOrdRejects>(pst->upper(), inCol);
    inCol.iSplitSide = 0;
    nOut = Traverse<ServiceColOrdRejects>(pst->lower(), inCol, lowerRej);
    mdlassert(pst->mdl(), nOut / sizeof(ServiceColRejects::output) == pst->n_lower());
    nOut = GetReply(rID, upperRej);
    mdlassert(pst->mdl(), nOut / sizeof(ServiceColRejects::output) == pst->n_upper());
    while (1) {
        iRet = RejMatch(pst, pst->n_lower(), lowerRej.data(), pst->n_upper(), upperRej.data(), pidSwap.data());
        if (!iRet) break;
        rID = mdl->ReqService<ServiceSwapRejects>(pst->upper(), pidSwap);
        nOut = Traverse<ServiceSwapRejects>(pst->lower(), pidSwap, lowerRej);
        mdlassert(pst->mdl(), nOut / sizeof(ServiceSwapRejects::output) == pst->n_lower());
        nOut = GetReply(rID, upperRej);
        mdlassert(pst->mdl(), nOut / sizeof(ServiceSwapRejects::output) == pst->n_upper());
    }
    return winfo.iOrdSplit;
}

/*****************************************************************************\
* ServiceLocalOrder
\*****************************************************************************/

static_assert(std::is_void_v<ServiceLocalOrder::input> || std::is_trivial_v<ServiceLocalOrder::input>);
static_assert(std::is_void_v<ServiceLocalOrder::output> || std::is_trivial_v<ServiceLocalOrder::output>);

int ServiceLocalOrder::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto in = reinterpret_cast<input*>(vin.data());
    auto iMinOrder = in->iMinOrder;
    auto iMaxOrder = in->iMaxOrder;
    auto iMidOrder = pst->data<OrderInfoPST>().iOrdSplit;
    assert(iMidOrder >= iMinOrder && iMidOrder <= iMaxOrder);
    in->iMinOrder = iMidOrder;
    auto rID = mdl->ReqService<ServiceLocalOrder>(pst->upper(), vin);
    in->iMinOrder = iMinOrder;
    in->iMaxOrder = iMidOrder - 1;
    Traverse<ServiceLocalOrder>(pst->lower(), vin, vout);
    GetReply(rID);
    return 0;
}

int ServiceLocalOrder::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto pkd = pst->lcl()->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    pkdLocalOrder(pkd, in->iMinOrder, in->iMaxOrder);
    return 0;
}
#endif

/*****************************************************************************\
* ServiceWeight
*
* Make sure that the local particles are split into active and inactive
* when passing pFlag != 0.
* pFlag == 0 => weight all particles.
* pFlag > 0 => weight active particles.
* pFlag < 0 => weight inactive particles.
\*****************************************************************************/

static_assert(std::is_void_v<ServiceWeight::input> || std::is_trivial_v<ServiceWeight::input>);
static_assert(std::is_void_v<ServiceWeight::output> || std::is_trivial_v<ServiceWeight::output>);

int ServiceWeight::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto in = reinterpret_cast<input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto rID = mdl->ReqService<ServiceWeight>(pst->upper(), vin);
    Traverse<ServiceWeight>(pst->lower(), vin, vout);
    output outWt;
    GetReply(rID, outWt);
    out->nLow += outWt.nLow;
    out->nHigh += outWt.nHigh;
    out->fLow += outWt.fLow;
    out->fHigh += outWt.fHigh;
    return sizeof(output);
}

int ServiceWeight::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto plcl = pst->lcl();
    auto pkd = plcl->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto fSplit = in->fSplit;
    auto iSplitSide = in->iSplitSide;
    auto& winfo = plcl->data_ensure<WeightInfo>();
    if (in->ittr == 0) {
        /*
        ** Initialize.
        */
        winfo.fSplit = fSplit;
        if (in->pFlag == 0) {
            winfo.iWtFrom = 0;
            winfo.iWtTo = pkd->Local() - 1;
        }
        else if (in->pFlag > 0) {
            /*
            ** Particles must be in the active-inactive order here!
            */
            winfo.iWtFrom = 0;
            winfo.iWtTo = pkdActive(pkd) - 1;
        }
        else {
            /*
            ** Particles must be in the active-inactive order here!
            */
            winfo.iWtFrom = pkdActive(pkd);
            winfo.iWtTo = pkd->Local() - 1;
        }
        winfo.fWtLow = 0.0;
        winfo.fWtHigh = 0.0;
    }
    else {
        /*
        ** Update the Weight Sums and use smaller weight region.
        */
        if (fSplit < winfo.fSplit) {
            winfo.fWtHigh += winfo.fHigh;
            if (iSplitSide) winfo.iWtFrom = winfo.iPart;
            else winfo.iWtTo = winfo.iPart - 1;
        }
        else {
            winfo.fWtLow += winfo.fLow;
            if (iSplitSide) winfo.iWtTo = winfo.iPart - 1;
            else winfo.iWtFrom = winfo.iPart;
        }
        winfo.fSplit = fSplit;
    }

    int nLow, nHigh;
    double fLow, fHigh;
    winfo.iPart = pkdWeight(pkd, in->iSplitDim, fSplit, iSplitSide, winfo.iWtFrom, winfo.iWtTo, &nLow, &nHigh, &fLow, &fHigh);

    out->nLow = nLow;
    out->nHigh = nHigh;
    out->fLow = fLow + winfo.fWtLow;
    out->fHigh = fHigh + winfo.fWtHigh;
    winfo.fLow = fLow;
    winfo.fHigh = fHigh;
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

static_assert(std::is_void_v<ServiceWeightWrap::input> || std::is_trivial_v<ServiceWeightWrap::input>);
static_assert(std::is_void_v<ServiceWeightWrap::output> || std::is_trivial_v<ServiceWeightWrap::output>);

int ServiceWeightWrap::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto in = reinterpret_cast<input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto rID = mdl->ReqService<ServiceWeightWrap>(pst->upper(), vin);
    Traverse<ServiceWeightWrap>(pst->lower(), vin, vout);
    output outWt;
    GetReply(rID, outWt);
    out->nLow += outWt.nLow;
    out->nHigh += outWt.nHigh;
    return sizeof(output);
}

int ServiceWeightWrap::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto plcl = pst->lcl();
    auto pkd = plcl->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto& info = plcl->data_ensure<WeightInfo>();
    if (in->ittr == 0) {
        /*
        ** Initialize.
        */
        info.fSplit = in->fSplit;
        info.iWtFrom = 0;
        info.iWtTo = pkd->Local() - 1;
    }
    else {
        /*
        ** Update the Weight Sums and use smaller weight region.
        */
        if ((in->fSplit < info.fSplit && in->fSplit < in->fSplit2 && info.fSplit < in->fSplit2) ||
            (in->fSplit < info.fSplit && in->fSplit > in->fSplit2 && info.fSplit > in->fSplit2) ||
            (in->fSplit > info.fSplit && in->fSplit > in->fSplit2 && info.fSplit < in->fSplit2)) {
            if (!in->iSplitSide) info.iWtFrom = info.iPart;
            else info.iWtTo = info.iPart - 1;
        }
        else {
            if (!in->iSplitSide) info.iWtTo = info.iPart - 1;
            else info.iWtFrom = info.iPart;
        }
        info.fSplit = in->fSplit;
    }
    int nLow, nHigh;
    info.iPart = pkdWeightWrap(pkd, in->iSplitDim, in->fSplit, in->fSplit2, in->iSplitSide, info.iWtFrom, info.iWtTo, &nLow, &nHigh);

    out->nLow = nLow;
    out->nHigh = nHigh;
    /* For collect rejects */
    info.nSplit = info.iPart;

    return sizeof(output);
}

/*****************************************************************************\
* ServiceOrdWeight
\*****************************************************************************/

static_assert(std::is_void_v<ServiceOrdWeight::input> || std::is_trivial_v<ServiceOrdWeight::input>);
static_assert(std::is_void_v<ServiceOrdWeight::output> || std::is_trivial_v<ServiceOrdWeight::output>);

int ServiceOrdWeight::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto mdl = pst->mdl();
    auto in = reinterpret_cast<input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto rID = mdl->ReqService<ServiceOrdWeight>(pst->upper(), vin);
    Traverse<ServiceOrdWeight>(pst->lower(), vin, vout);
    output outWt;
    GetReply(rID, outWt);
    out->nLow += outWt.nLow;
    out->nHigh += outWt.nHigh;
    return sizeof(output);
}

int ServiceOrdWeight::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    auto plcl = pst->lcl();
    auto pkd = plcl->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto& info = plcl->data_ensure<WeightInfo>();
    if (in->ittr == 0) {
        /*
        ** Initialize.
        */
        info.iOrdSplit = in->iOrdSplit;
        info.iWtFrom = 0;
        info.iWtTo = pkd->Local() - 1;
    }
    else {
        /*
        ** Update the Weight Sums and use smaller weight region.
        */
        if (in->iOrdSplit < info.iOrdSplit) {
            if (in->iSplitSide) info.iWtFrom = info.iPart;
            else info.iWtTo = info.iPart - 1;
        }
        else {
            if (in->iSplitSide) info.iWtTo = info.iPart - 1;
            else info.iWtFrom = info.iPart;
        }
        info.iOrdSplit = in->iOrdSplit;
    }
    int nLow, nHigh;
    info.iPart = pkdOrdWeight(pkd, in->iOrdSplit, in->iSplitSide, info.iWtFrom, info.iWtTo, &nLow, &nHigh);
    out->nLow = nLow;
    out->nHigh = nHigh;
    return sizeof(output);
}

auto ServiceActiveOrder::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->lcl()->pkd;
    *out = pkd->ActiveOrder();
    return sizeof(output);
}
} // namespace OldDD

auto pkdContext::ActiveOrder() -> uint64_t {
    auto i = std::partition(particles.begin(), particles.end(), [](auto& p) { return p.is_active(); });
    return (nActive = i - particles.begin());
}
