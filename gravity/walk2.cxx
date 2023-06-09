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

#ifdef __SSE__
    #include <xmmintrin.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
    #include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include <string.h>
#include "pkd.h"
#include "walk.h"
#include "gravity/grav.h"
#include "group/hop.h"
#include "smooth/smooth.h"
#include "moments.h"
#include "vmoments.h"
#include "gravity/opening.h"
#include "cl.h"
#include "cuda/cudautil.h"
#include "cuda/cudapppc.h"
#include "cuda/cudaewald.h"
#include "../SPH/SPHOptions.h"
#include "../SPH/SPHEOS.h"
#include "../SPH/SPHpredict.h"

static void addChild(PKD pkd, int iCache, clList *cl, int iChild, int id, blitz::TinyVector<float,3> fOffset) {
    auto c = (id == pkd->Self()) ? pkd->tree[iChild] :
             pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,iCache,iChild,id))];
    auto nc = (c->is_remote()|c->is_top_tree()) ? 1000000000 : c->count();
    auto cOpen = c->bMax() * pkd->fiCritTheta;
    auto c_r = c->position();
    auto cbnd = c->bound();
    auto [iLower, idLower, iUpper, idUpper] = c->get_child_cells(id);
    const auto &SPHbob = c->have_BOB() ? c->BOB() : SPHBOB();
    cl->append(iCache,id,iChild,idLower,iLower,idUpper,iUpper,nc,cOpen,
               c->mass(),4.0f*c->fSoft2(),c_r,fOffset,cbnd,
               SPHbob);

}
/*
** Returns total number of active particles for which gravity was calculated.
*/
static int processCheckList(PKD pkd, SMX smx, SMF smf, int iRoot, int iRoot2,
                            struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
                            double dTime,int bEwald,bool bGPU,
                            double dThetaMin, double *pdFlop, double *pdPartSum,double *pdCellSum,SPHOptions *SPHoptions) {
    treeStore::NodePointer c(pkd->tree);
    int id,iCell,iSib,iCheckCell,iCheckLower;
    FMOMR monoPole;
    LOCR L;
    FLOCR Lf;
    float vpred[3];
    blitz::TinyVector<float,3> a;
    blitz::TinyVector<double,3> r, k_r, c_r, dOffset, dx;
    double xParent,yParent,zParent;
    blitz::TinyVector<float,3> fOffset;
    float dirLsum,normLsum;
    float fMass,fSoft;
    int iStack;
    int j,jTile,pj,nActive,nTotActive;
    static const blitz::TinyVector<float,3> fZero3(0,0,0);
    int iCidPart;
    int bReferenceFound;
    double dFlop;

    pkd->dFlop = 0.0; /* Flops are accumulated here! */
    iStack = -1;

    /*
    ** Clear monopole sentinel and local expansion and the timestepping sums.
    */
    momClearFmomr(&monoPole);
    momClearLocr(&L);
    dirLsum = 0;
    normLsum = 0;
    nTotActive = 0;

    a = 0;

    /*
    ** We are now going to work on the local tree.
    ** Make iCell point to the root of the tree again.
    */
    auto k = pkd->tree[iCell = iRoot];
    k_r = k->position();
    while (1) {
#ifdef ILP_ILC_CAN_BE_NON_EMPTY
        ILPTILE tile;
        ILCTILE ctile;
        blitz::TinyVector<double,3> reference;
        double d2c;
        int iLower,iUpper;
        /*
        ** Find the next active particle that will be encountered in the walk algorithm below
        ** in order to set a good initial center for the P-P interaction list.
        */
        auto kFind = k;
        while (kFind->iLower) {
            int iCellDescend;
            kFind = pkd->tree[iCellDescend = kFind->iLower];
            if (kFind->uMinRung>uRungHi || kFind->uMaxRung < uRungLo) {
                /*
                ** Move onto processing the sibling.
                */
                kFind = pkd->tree[++iCellDescend];
            }
        }
        for (auto &p : *kfind) {
            if (!p.is_active()) continue;
            reference = p.position();
            goto found_it;
        }
        printf("%d: TREE ERROR\n", pkd->Self());
        assert(0); /* We didn't find an active particle */
found_it:
        reference -= pkd->ilp->getReference();
        d2c = blitz::dot(reference,reference);
        if ( d2c > 1e-5 ) {
            // Correct all remaining PP/PC interactions to this new center.
            pkd->ilp.recenter(reference);
            pkd->ilc.recenter(reference);
        }
        bReferenceFound = 1;
#else
        bReferenceFound = 0;
#endif
        while (1) {
            /*
            ** Process the Checklist.
            */
#ifdef USE_SIMD_FMM
            a = (ts->bGravStep) ? k->acceleration() : 0.0f;
#else
            float imaga;
            if (ts->bGravStep) {
                a = k->acceleration();
                imaga = 1.0 / sqrtf(blitz::dot(a,a));
            }
#endif
            /*
            ** For cells which will remain on the checklist for a further deeper level of
            ** level of the tree we will be using the stack cl (even if no push happens later).
            */
            pkd->S[iStack+1].cl->clear();
#ifdef USE_SIMD_FMM
            pkd->ill.clear();
#endif
            do {
                for (auto &tile : *pkd->cl) {
                    iOpenOutcomeSIMD(pkd,k,tile,dThetaMin,SPHoptions);
                }
                pkd->clNew->clear();
                for (auto &tile : *pkd->cl) {
                    auto nBlocks = tile.count() / tile.width;
                    for (auto iBlock=0; iBlock<=nBlocks; ++iBlock) {
                        int n = iBlock<nBlocks ? tile.width : tile.count() - nBlocks*tile.width;
                        auto &blk = tile[iBlock];
                        for (jTile=0; jTile<n; ++jTile) {
                            switch (blk.iOpen[jTile]) {
                            case 0:
                                /*
                                ** This checkcell stays on the checklist.
                                */
                                pkd->S[iStack+1].cl->append(blk,jTile);
                                break;
                            case 1:
                                /*
                                ** This checkcell's particles are added to the P-P list.
                                */
                                iCheckCell = blk.iCell[jTile];
                                id = blk.idCell[jTile];
                                if (iCheckCell < 0) {
                                    pj = -1 - iCheckCell;
                                    assert(id >= 0);
                                    iCidPart = blk.iCache[jTile]==CID_CELL ? CID_PARTICLE : CID_PARTICLE2;
                                    if (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) {
                                        auto p = (id == pkd->Self()) ? pkd->particles[pj] : pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,iCidPart,pj,id))];
                                        if (SPHoptions->doSetDensityFlags) p.set_marked(true);
                                        if (SPHoptions->doSetNNflags) p.set_NN_flag(true);
                                        if (id != pkd->Self()) mdlRelease(pkd->mdl,iCidPart,&p);
                                    }
                                    else {
                                        auto p = (id == pkd->Self()) ? pkd->particles[pj]
                                                 : pkd->particles[static_cast<PARTICLE *>(mdlFetch(pkd->mdl,iCidPart,pj,id))];
                                        r = p.position();
                                        if (!bReferenceFound) {
                                            bReferenceFound=1;
                                            pkd->ilp.setReference(r);
                                            pkd->ilc.setReference(r);
                                        }
                                        if (pkd->particles.present(PKD_FIELD::oNewSph)) {
                                            const auto &NewSph = p.newsph();
                                            float Omega = NewSph.Omega;
                                            float P = 0.0f;
                                            float cs = 0.0f;
                                            float T = 0.0f;
                                            float expImb2 = NewSph.expImb2;
                                            SPHpredictOnTheFly(pkd, p, kick, SPHoptions->nPredictRung, vpred, &P, &cs, &T, SPHoptions);
                                            pkd->ilp.append(
                                                r[0] + blk.xOffset[jTile],
                                                r[1] + blk.yOffset[jTile],
                                                r[2] + blk.zOffset[jTile],
                                                blk.m[jTile], blk.fourh2[jTile],
                                                vpred[0], vpred[1], vpred[2],
                                                p.ball(), Omega, p.density(), P, cs, p.species(), p.rung(), p.imaterial(), T, expImb2);
                                        }
                                        else {
                                            auto v = p.velocity();
                                            pkd->ilp.append(
                                                r[0] + blk.xOffset[jTile],
                                                r[1] + blk.yOffset[jTile],
                                                r[2] + blk.zOffset[jTile],
                                                blk.m[jTile], blk.fourh2[jTile],
                                                v[0], v[1], v[2],
                                                0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, p.rung(), 0, 0.0f, 0.0f);
                                        }
                                    }
                                }
                                else {
                                    assert(id >= 0);
                                    c = (id == pkd->Self()) ? pkd->tree[iCheckCell]
                                        : pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,blk.iCache[jTile],iCheckCell,id))];
                                    iCidPart = blk.iCache[jTile]==CID_CELL ? CID_PARTICLE : CID_PARTICLE2;
                                    if (!bReferenceFound) {
                                        bReferenceFound=1;
                                        auto p = (id == pkd->Self()) ? pkd->particles[c->lower()] : ((SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) ? pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,iCidPart,c->lower(),id))] : pkd->particles[static_cast<PARTICLE *>(mdlFetch(pkd->mdl,iCidPart,c->lower(),id))]);
                                        r = p.position();
                                        pkd->ilp.setReference(r);
                                        pkd->ilc.setReference(r);
                                        if ((id != pkd->Self()) && (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags)) mdlRelease(pkd->mdl,iCidPart,&p);
                                    }
                                    for (pj=c->lower(); pj<=c->upper(); ++pj) {
                                        if (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) {
                                            auto p = (id == pkd->Self()) ? pkd->particles[pj] : pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,iCidPart,pj,id))];
                                            if (SPHoptions->doSetDensityFlags) p.set_marked(true);
                                            if (SPHoptions->doSetNNflags) p.set_NN_flag(true);
                                            if (id != pkd->Self()) mdlRelease(pkd->mdl,iCidPart,&p);
                                        }
                                        else {
                                            auto p = (id == pkd->Self()) ? pkd->particles[pj]
                                                     : pkd->particles[static_cast<PARTICLE *>(mdlFetch(pkd->mdl,iCidPart,pj,id))];
                                            fMass = p.mass();
                                            fSoft = p.soft();
                                            r = p.position();
                                            if (p.have_newsph()) {
                                                const auto &NewSph = p.newsph();
                                                float Omega = NewSph.Omega;
                                                float P = 0.0f;
                                                float cs = 0.0f;
                                                float T = 0.0f;
                                                float expImb2 = NewSph.expImb2;
                                                SPHpredictOnTheFly(pkd, p, kick, SPHoptions->nPredictRung, vpred, &P, &cs, &T, SPHoptions);
                                                pkd->ilp.append(
                                                    r[0] + blk.xOffset[jTile],
                                                    r[1] + blk.yOffset[jTile],
                                                    r[2] + blk.zOffset[jTile],
                                                    fMass, 4*fSoft*fSoft,
                                                    vpred[0], vpred[1], vpred[2],
                                                    p.ball(), Omega, p.density(), P, cs, p.species(), p.rung(), p.imaterial(), T, expImb2);
                                            }
                                            else {
                                                auto v = p.velocity();
                                                pkd->ilp.append(
                                                    r[0] + blk.xOffset[jTile],
                                                    r[1] + blk.yOffset[jTile],
                                                    r[2] + blk.zOffset[jTile],
                                                    fMass, 4*fSoft*fSoft,
                                                    v[0], v[1], v[2],
                                                    0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0, p.rung(), 0, 0.0f, 0.0f);
                                            }
                                        }
                                    }
                                }
                                break;
                            case 2:
                                /*
                                ** Now I am trying to open a bucket, which means I add each particle back on the
                                ** checklist with a cell size of zero.
                                */
                                iCheckCell = blk.iCell[jTile];
                                assert(iCheckCell>=0);
                                id = blk.idCell[jTile];
                                fOffset[0] = blk.xOffset[jTile];
                                fOffset[1] = blk.yOffset[jTile];
                                fOffset[2] = blk.zOffset[jTile];
                                assert(id >= 0);
                                c = (id == pkd->Self()) ? pkd->tree[iCheckCell]
                                    : pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,blk.iCache[jTile],iCheckCell,id))];
                                iCidPart = blk.iCache[jTile]==CID_CELL ? CID_PARTICLE : CID_PARTICLE2;
                                for (pj=c->lower(); pj<=c->upper(); ++pj) {
                                    auto p = (id == pkd->Self()) ? pkd->particles[pj] : ((SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags) ? pkd->particles[static_cast<PARTICLE *>(mdlAcquire(pkd->mdl,iCidPart,pj,id))] : pkd->particles[static_cast<PARTICLE *>(mdlFetch(pkd->mdl,iCidPart,pj,id))]);
                                    r = p.position();
                                    fMass = p.mass();
                                    fSoft = p.soft();
                                    if (p.have_newsph()) {
                                        float fBallFactor = (SPHoptions->dofBallFactor) ? SPHoptions->fBallFactor : 1.0f;
                                        float limitedBallSize = std::min(SPHoptions->ballSizeLimit,fBallFactor * p.ball());
                                        pkd->clNew->append(blk.iCache[jTile],id,-1 - pj,0,0,0,0,1,0.0,fMass,4.0f*fSoft*fSoft,
                                                           r,                   // center of mass
                                                           fOffset,             // fOffset
                                                           Bound(r,r),          // zero size box at r
                                                           SPHBOB(r,limitedBallSize));
                                    }
                                    else {
                                        pkd->clNew->append(blk.iCache[jTile],id,-1 - pj,0,0,0,0,1,0.0,fMass,4.0f*fSoft*fSoft,
                                                           r,                   // center of mass
                                                           fOffset,             // fOffset
                                                           Bound(r,r),          // zero size box at r
                                                           SPHBOB(r,0.0));
                                    }
                                    if ((id != pkd->Self()) && (SPHoptions->doSetDensityFlags || SPHoptions->doSetNNflags)) mdlRelease(pkd->mdl,iCidPart,&p);
                                }
                                break;
                            case 3:
                                /*
                                ** Open the cell.
                                ** We could do a prefetch here for non-local
                                ** cells.
                                */
                                iCheckCell = blk.iCell[jTile];                 assert(iCheckCell >= 0);
                                iCheckLower = blk.iLower[jTile];               assert(iCheckLower > 0);

                                fOffset[0] = blk.xOffset[jTile];
                                fOffset[1] = blk.yOffset[jTile];
                                fOffset[2] = blk.zOffset[jTile];

                                addChild(pkd,blk.iCache[jTile],pkd->clNew,blk.iLower[jTile],blk.idLower[jTile],fOffset);
                                addChild(pkd,blk.iCache[jTile],pkd->clNew,blk.iUpper[jTile],blk.idUpper[jTile],fOffset);

                                break;
                            case 4:
                                /*
                                ** Accept multipole!
                                ** Interact += Moment(c);
                                */
                                if (SPHoptions->doGravity) {
                                    iCheckCell = blk.iCell[jTile];
                                    assert(iCheckCell>=0);
                                    id = blk.idCell[jTile];
                                    c = (id == pkd->Self()) ? pkd->tree[iCheckCell]
                                        : pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,blk.iCache[jTile],iCheckCell,id))];
                                    r[0] = blk.x[jTile] + blk.xOffset[jTile];
                                    r[1] = blk.y[jTile] + blk.yOffset[jTile];
                                    r[2] = blk.z[jTile] + blk.zOffset[jTile];
                                    if (!bReferenceFound) {
                                        bReferenceFound=1;
                                        pkd->ilp.setReference(r);
                                        pkd->ilc.setReference(r);
                                    }
                                    pkd->ilc.append(r[0],r[1],r[2],&c->moment(),c->bMax());
                                }
                                break;
                            case 8:
                                /*
                                ** Local expansion accepted!
                                */
                                if (SPHoptions->doGravity) {
                                    iCheckCell = blk.iCell[jTile];
                                    if (iCheckCell<0) {
                                        fOffset[0] = blk.xOffset[jTile];
                                        fOffset[1] = blk.yOffset[jTile];
                                        fOffset[2] = blk.zOffset[jTile];
                                        dx[0] = k_r[0] - (blk.x[jTile] + blk.xOffset[jTile]);
                                        dx[1] = k_r[1] - (blk.y[jTile] + blk.yOffset[jTile]);
                                        dx[2] = k_r[2] - (blk.z[jTile] + blk.zOffset[jTile]);
#ifdef USE_SIMD_FMM
                                        monoPole.m = blk.m[jTile];
                                        pkd->ill.append((float)dx[0],(float)dx[1],(float)dx[2],&monoPole,1.0);
#else
                                        double d2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
                                        double dir = 1.0/sqrt(d2);
                                        /* monoPole.m = blk.m[jTile];*/
                                        /* *pdFlop += momLocrAddFmomr5cm(&L,&monoPole,0.0,dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);*/
                                        double tax,tay,taz;
                                        dFlop = momLocrAddMono5(&L,blk.m[jTile],dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
                                        *pdFlop += dFlop;
                                        pkd->dFlopDoubleCPU += dFlop;
                                        if (ts->bGravStep) {
                                            float adotai = a[0]*tax + a[1]*tay + a[2]*taz;
                                            if (adotai > 0) {
                                                adotai *= imaga;
                                                dirLsum += dir*adotai*adotai;
                                                normLsum += adotai*adotai;
                                            }
                                        }
#endif
                                    }
                                    else {
                                        id = blk.idCell[jTile];
                                        dOffset[0] = blk.xOffset[jTile];
                                        dOffset[1] = blk.yOffset[jTile];
                                        dOffset[2] = blk.zOffset[jTile];
                                        c = (id == pkd->Self()) ? pkd->tree[iCheckCell]
                                            : pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,blk.iCache[jTile],iCheckCell,id))];
                                        c_r = c->position();
#ifdef USE_SIMD_FMM
                                        dx = k_r - (c_r + dOffset);
                                        pkd->ill.append((float)dx[0],(float)dx[1],(float)dx[2],&c->moment(),c->bMax());
#else
                                        dx = k_r - (c_r + dOffset);
                                        double d2 = blitz::dot(dx,dx);
                                        double dir = 1.0/sqrt(d2);
                                        double tax,tay,taz;
                                        dFlop = momLocrAddFmomr5cm(&L,&c->moment(),c->bMax(),dir,dx[0],dx[1],dx[2],&tax,&tay,&taz);
                                        *pdFlop += dFlop;
                                        pkd->dFlopDoubleCPU += dFlop;
                                        if (ts->bGravStep) {
                                            float adotai = a[0]*tax + a[1]*tay + a[2]*taz;
                                            if (adotai > 0) {
                                                adotai *= imaga;
                                                dirLsum += dir*adotai*adotai;
                                                normLsum += adotai*adotai;
                                            }
                                        }
#endif
                                    }
                                }
                                break;
                            case 10:
                                /*
                                ** This checkcell is removed from the checklist since it has zero/negative mass.
                                */
                                break;
                            default:
                                assert(0);
                            }
                        }
                    }
                } /* end of CL_LOOP */
                std::swap(pkd->cl,pkd->clNew);
            } while (pkd->cl->count());
            /*
            ** Now calculate the local expansion.
            */
#ifdef USE_SIMD_FMM
            // Need to do something here with dirLsum and normLsum for GravStep
            // Need to get the scaling factor correct here
            if (pkd->ill.count()) {
                float v = k->bMax();
                dFlop = momFlocrSetVFmomr5cm(&Lf,v,pkd->ill,a.data(),&dirLsum,&normLsum);
                *pdFlop += dFlop;
                pkd->dFlopSingleCPU += dFlop;
                momLocrAddFlocr(&L,&Lf,v);
            }
#endif
            /*
            ** Done processing of the Checklist.
            ** Now prepare to proceed to the next deeper
            ** level of the tree.
            */
            if (k->is_group()) break; /* A bucket is ALWAYS a group */
            xParent = k_r[0];
            yParent = k_r[1];
            zParent = k_r[2];
            for (j=0; j<3; ++j) fOffset[j] = 0.0f;
            iCell = k->lchild();
            k = pkd->tree[iCell];
            k_r = k->position();
            /*
            ** Check iCell is active. We eventually want to just to a
            ** rung check here when we start using tree repair, but
            ** for now this is just as good.
            */
            if ((k->min_rung()<=ts->uRungHi && k->max_rung() >= ts->uRungLo) || (SPHoptions->useDensityFlags && k->is_marked()) || (SPHoptions->useNNflags && k->is_NN())) {
                /*
                ** iCell is active, continue processing it.
                ** Put the sibling onto the checklist.
                */
                iSib = iCell+1;
                c = pkd->tree[iSib];
                if ((c->min_rung()<=ts->uRungHi && c->max_rung() >= ts->uRungLo) || (SPHoptions->useDensityFlags && c->is_marked()) || (SPHoptions->useNNflags && c->is_NN())) {
                    /*
                    ** Sibling is active so we need to clone the checklist!
                    */
                    pkd->cl->clone(*pkd->S[iStack+1].cl);
                }
                else {
                    /*
                    ** Otherwise we can simple grab it.
                    */
                    std::swap(pkd->cl,pkd->S[iStack+1].cl);
                }
                /*
                ** Test whether the sibling is active as well.
                ** If not we don't push it onto the stack, but we
                ** have to be careful to not pop the stack when we
                ** hit the sibling.
                */
                if ((c->min_rung()<=ts->uRungHi && c->max_rung() >= ts->uRungLo) || (SPHoptions->useDensityFlags && c->is_marked()) || (SPHoptions->useNNflags && c->is_NN())) {
                    /*
                    ** Sibling is active as well.
                    ** Push Checklist for the sibling onto the stack
                    ** before proceeding deeper in the tree.
                    */
                    ++iStack;
                    assert(iStack < pkd->nMaxStack);
                    pkd->S[iStack].iNodeIndex = iSib;
                    /*
                    ** At this point, the ILP is normally empty if you never do P-P except when reaching a bucket.
                    ** Softened multi-poles are also an exception.
                    */
                    //ilpCheckPt(pkd->ilp,&pkd->S[iStack].PartChkPt);
                    //ilcCheckPt(pkd->ilc,&pkd->S[iStack].CellChkPt);
                    assert(pkd->ilp.count()==0);
                    assert(pkd->ilc.count()==0);
                    /*
                    ** Note here we already have the correct elements in S[iStack] (iStack+1 was used previously), just need to add one.
                    */
                    pkd->S[iStack].L = L;
                    pkd->S[iStack].dirLsum = dirLsum;
                    pkd->S[iStack].normLsum = normLsum;
                    c_r = c->position();
                    momShiftLocr(&pkd->S[iStack].L,
                                 c_r[0] - xParent,
                                 c_r[1] - yParent,
                                 c_r[2] - zParent);
                }
            }
            else {
                /*
                ** Skip iCell, but add it to the Checklist.
                ** No need to push anything onto the stack here.
                ** We can simply grab it the checklist from the Stack.
                */
                std::swap(pkd->cl,pkd->S[iStack+1].cl);
                /*
                ** Move onto processing the sibling.
                */
                k = pkd->tree[++iCell];
                k_r = k->position();
            }
            dFlop = momShiftLocr(&L,k_r[0] - xParent,
                                 k_r[1] - yParent,
                                 k_r[2] - zParent);
            *pdFlop += dFlop;
            pkd->dFlopDoubleCPU += dFlop;
        }
        /*
        ** Now the interaction list should be complete and the
        ** Checklist should be empty! Calculate gravity on this
        ** Bucket!
        */
        nActive = pkdGravInteract(pkd,kick,lc,ts,
                                  k,&L,pkd->ilp,pkd->ilc,dirLsum,normLsum,bEwald,pdFlop,
                                  smx, &smf, iRoot, iRoot2, SPHoptions, bGPU);
        /*
        ** Update the limit for a shift of the center here based on the opening radius of this
        ** cell (the one we just evaluated).
        */
        if (nActive) {
            /*
            ** Here we used to set the weights of particles based on the work done, but now we just assume that
            ** all particles cost the same in domain decomposition, so we really don't need to set anything here.
            */
            *pdPartSum += nActive * pkd->ilp.count();
            *pdCellSum += nActive * pkd->ilc.count();
            nTotActive += nActive;
        }
        /* Get the next cell to process from the stack */
        if (iStack == -1) goto doneCheckList;
        k = pkd->tree[iCell = pkd->S[iStack].iNodeIndex];
        k_r = k->position();
        /*
        ** Pop the Checklist from the top of the stack,
        ** also getting the state of the interaction list.
        */
        //ilpRestore(pkd->ilp,&pkd->S[iStack].PartChkPt);
        //ilcRestore(pkd->ilc,&pkd->S[iStack].CellChkPt);
        pkd->ilp.clear();
        pkd->ilc.clear();
        /*
        ** Grab the checklist from the stack.
        */
        assert(pkd->cl->count() == 0);
        std::swap(pkd->cl,pkd->S[iStack].cl);
        L = pkd->S[iStack].L;
        dirLsum = pkd->S[iStack].dirLsum;
        normLsum = pkd->S[iStack].normLsum;
        --iStack;
    }
doneCheckList:
#ifdef USE_CUDA
    pkd->cudaClient->flushCUDA();
    if (SPHoptions->doDensity) {
        /*
        ** In the density pass, a wp can requeue itself to perform the Newton iterations.
        ** To make sure that after the thread has reached the end of the tree walk
        ** the work that is still being resent by the last few wps will still be treated,
        ** we not only have to check (as done in mdlCompleteAllWork),
        ** but also flush the queues repeatedly.
        ** Once this while loop exits, ALL density work has been done.
        */
        while (pkd->mdl->gpu.flushCompleted()) pkd->cudaClient->flushCUDA();
    }
#endif
#ifdef USE_METAL
    pkd->metalClient->flushMETAL();
#endif
    mdlCompleteAllWork(pkd->mdl);
    *pdFlop += pkd->dFlop; /* Accumulate work flops (notably Ewald) */
    return (nTotActive);
}

static void doneGravWalk(PKD pkd,SMX smx,SMF *smf) {
    if (smx) {
        smSmoothFinish(smx);
        smFinish(smx,smf);
    }
}

static void initGravWalk(PKD pkd,double dTime,double dThetaMin,int bPeriodic,int bGravStep,int nPartRhoLoc,int iTimeStepCrit,
                         SMX *smx, SMF *smf) {

    pkd->dEnergyU = 0.0;
    pkd->dEnergyT = 0.0;
    pkd->dEnergyW = 0.0;
    pkd->dEnergyF[0] = pkd->dEnergyF[1] = pkd->dEnergyF[2] = 0.0;
    pkd->dEnergyL[0] = pkd->dEnergyL[1] = pkd->dEnergyL[2] = 0.0;

    pkd->fiCritTheta = 1.0f / dThetaMin;

    assert(pkd->tree.present(KDN_FIELD::oNodeMass));
    if (bGravStep) {
        assert(pkd->tree.present(KDN_FIELD::oNodeAcceleration));
        if (iTimeStepCrit == 1) {
            assert(pkd->tree.present(KDN_FIELD::oNodeVelocity));
            assert(pkd->particles.present(PKD_FIELD::oVelocity));
        }
    }

    /*
    ** Setup smooth for calculating local densities when a particle has too few P-P interactions.
    */
    if (bGravStep) {
        smInitializeRO(smx,pkd,smf,nPartRhoLoc,bPeriodic,SMX_DENSITY_F1);
        smSmoothInitialize(*smx);
    }
    else (*smx) = NULL;
}

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalkHop(PKD pkd,double dTime, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    int id,iRoot,iRootSelf;
    float fOffset[3];
    int nActive;
    int i,j,gid;
    SMX smx;
    SMF smf;

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkd->particles,pkd->particles.ParticleSize(), pkd->Local());
    initGravWalk(pkd,dTime,dThetaMin,0,0,0,0,&smx,&smf);
    nActive = 0;
    for (gid=1; gid<pkd->nGroups; ++gid) {
        if (!pkd->hopGroups[gid].bNeedGrav) continue;
        pkd->hopGroups[gid].bNeedGrav = 0;
        pkd->ilp.clear();
        pkd->ilc.clear();
        pkd->cl->clear();
        iRootSelf = pkd->hopGroups[gid].iAllRoots;
        for (i=iRootSelf; i<iRootSelf + pkd->hopGroups[gid].nRemote+1; ++i) {
            for (j=0; j<3; ++j) fOffset[j] = 0.0f;
            id = pkd->hopRoots[i].iPid;
            iRoot = pkd->hopRoots[i].iIndex;
            assert(iRoot>0);
            addChild(pkd,CID_CELL,pkd->cl,iRoot,id,fOffset);
        }
        assert(pkd->hopRoots[iRootSelf].iPid==pkd->Self());
        // nActive += processCheckList(pkd, smx, smf, pkd->hopRoots[iRootSelf].iIndex, 0, 0, MAX_RUNG,
        //     NULL,NULL,1.0,dTime,
        //     0, dThetaMin, 0, 0, pdFlop, pdPartSum, pdCellSum);
    }
    doneGravWalk(pkd,smx,&smf);
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    return nActive;
}


/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalk(PKD pkd,struct pkdKickParameters *kick,struct pkdLightconeParameters *lc,struct pkdTimestepParameters *ts,
                double dTime,int nReps,int bEwald,bool bGPU,
                int iLocalRoot1, int iLocalRoot2,int iVARoot,
                double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum,SPHOptions *SPHoptions) {
    int id;
    float fOffset[3];
    int ix,iy,iz;
    int nActive = 0;
    SMX smx;
    SMF smf;
    int iTop1, iTop2;

    initGravWalk(pkd,dTime,dThetaMin,nReps?1:0,ts->bGravStep,ts->nPartRhoLoc,ts->iTimeStepCrit,&smx,&smf);

    iTop1 = pkd->iTopTree[iLocalRoot1];
    iTop2 = pkd->iTopTree[iLocalRoot2];
    id = pkd->Self();

    /*
    ** Walk tree 1 against trees 1 (and optionally 2) if there are active particles
    */
    auto k = pkd->tree[iLocalRoot1];
    if (k->count() && (k->is_active(ts->uRungLo,ts->uRungHi) || (SPHoptions->useDensityFlags && k->is_marked())  || (SPHoptions->useNNflags && k->is_NN()))) {
        /*
        ** Initially we set our cell pointer to
        ** point to the top tree.
        */
        pkd->ilp.clear();
        pkd->ilc.clear();
        pkd->cl->clear();

        /*
        ** Add all replicas of the entire box to the Checklist.
        ** We add at least one box (0,0,0). The root cell is alway on processor 0.
        */
        for (ix=-nReps; ix<=nReps; ++ix) {
            fOffset[0] = ix*pkd->fPeriod[0];
            for (iy=-nReps; iy<=nReps; ++iy) {
                fOffset[1] = iy*pkd->fPeriod[1];
                for (iz=-nReps; iz<=nReps; ++iz) {
                    fOffset[2] = iz*pkd->fPeriod[2];
                    addChild(pkd,CID_CELL,pkd->cl,iTop1,id,fOffset);
#ifndef SINGLE_CACHES
                    if (iLocalRoot2>0) addChild(pkd,CID_CELL2,pkd->cl,iTop2,id,fOffset);
#else
                    if (iLocalRoot2>0) addChild(pkd,CID_CELL,pkd->cl,iTop2,id,fOffset);
#endif
                }
            }
        }
        nActive += processCheckList(pkd, smx, smf, iLocalRoot1, iLocalRoot2, kick,lc,ts,
                                    dTime,bEwald, bGPU, dThetaMin, pdFlop, pdPartSum, pdCellSum, SPHoptions);
    }
#if 0
    /*
    ** Walk tree 2 against tree 1.
    */
    if (iLocalRoot2>0) {
        /* Check that the iRoot has active particles! */
        if (!pkd->tree[iLocalRoot2]->is_active(uRungLo,uRungHi)) return 0;
        pkd->ilp.clear();
        pkd->ilc.clear();
        pkd->cl->clear();
        /*
        ** Add all replicas of the entire box to the Checklist.
        ** We add at least one box (0,0,0). The root cell is alway on processor 0.
        */
        for (ix=-nReps; ix<=nReps; ++ix) {
            fOffset[0] = ix*pkd->fPeriod[0];
            for (iy=-nReps; iy<=nReps; ++iy) {
                fOffset[1] = iy*pkd->fPeriod[1];
                for (iz=-nReps; iz<=nReps; ++iz) {
                    fOffset[2] = iz*pkd->fPeriod[2];
                    addChild(pkd,CID_CELL,pkd->cl,iTop1,id,fOffset);
                }
            }
        }
        nActive += processCheckList(pkd, smx, smf, iLocalRoot2, iLocalRoot1,
                                    kick,lc,ts,dTime,0, dThetaMin, pdFlop, pdPartSum, pdCellSum, SPHoptions);
    }
#endif
    doneGravWalk(pkd,smx,&smf);
    return nActive;
}

/*
** Returns total number of active particles for which gravity was calculated.
*/
int pkdGravWalkGroups(PKD pkd,double dTime, double dThetaMin,double *pdFlop,double *pdPartSum,double *pdCellSum) {
    int id,iRoot;
    blitz::TinyVector<float,3> fOffset;
    int i,k;
    SMX smx;
    SMF smf;

    initGravWalk(pkd,dTime,dThetaMin,0,0,0,0,&smx,&smf);
    /*
    ** Initially we set our cell pointer to
    ** point to the top tree.
    */

    struct psGroup *gd = pkd->psGroupTable.pGroup;
    int nActive=0;
    for (i=1; i < pkd->psGroupTable.nGroups; i++) {
        if (gd[i].nLocal == 0) continue;
        pkd->ilp.clear();
        pkd->ilc.clear();
        pkd->cl->clear();
#if 1
        for (k=1; k < gd[i].nTreeRoots; k++) {
            fOffset = 0.0f;
            id = gd[i].treeRoots[k].iPid;
            iRoot = gd[i].treeRoots[k].iLocalRootId;
            addChild(pkd,CID_CELL,pkd->cl,iRoot,id,fOffset);
        }
#endif
        // nActive += processCheckList(pkd, smx, smf, gd[i].treeRoots[0].iLocalRootId, 0, 0, MAX_RUNG,
        //     NULL,NULL,1.0,dTime,
        //     0, dThetaMin, 0, 0, pdFlop, pdPartSum, pdCellSum);
    }
    doneGravWalk(pkd,smx,&smf);
    return nActive;
}
