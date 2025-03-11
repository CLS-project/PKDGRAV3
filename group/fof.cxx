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
#include <math.h>
#include "basetype.h"
#include "fof.h"
#include "pkd.h"
#include "group.h"
#include "blitz/array.h"
#include "../SPH/SPHOptions.h"
#include "core/remote.h"

#include <algorithm>

using blitz::TinyVector;
using blitz::sum;
using blitz::dot;
using blitz::any;

uint32_t pkdFofGatherLocal(PKD pkd,int *S,double fBall2,TinyVector<double,3> r,uint32_t iGroup,
                           uint32_t iTail,uint32_t *Fifo,int *pbCurrFofContained,
                           const Bound::coord_type fMinFofContained,const Bound::coord_type fMaxFofContained) {
    double min2,dx,dy,dz,fDist2;
    int sp = 0;
    int iCell,pj,j;
    uint32_t iPartGroup;

    auto kdn = pkd->tree[iCell = ROOT];
    while (1) {
        auto bnd = kdn->bound();
        min2 = bnd.mindist(r);
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
        for (pj=kdn->lower(); pj<=kdn->upper(); ++pj) {
            auto p = pkd->particles[pj];
            iPartGroup = p.group();
            if (iPartGroup) continue;
            auto p_r = p.position();
            dx = r[0] - p_r[0];
            dy = r[1] - p_r[1];
            dz = r[2] - p_r[2];
            fDist2 = dx*dx + dy*dy + dz*dz;
            if (fDist2 <= fBall2) {
                /*
                **  Mark particle and add it to the do-fifo
                */
                p.set_group(iGroup);
                Fifo[iTail++] = pj;
                if (*pbCurrFofContained) {
                    for (j=0; j<3; ++j) {
                        if (p_r[j] < fMinFofContained[j]) {
                            *pbCurrFofContained = 0;
                            break;
                        }
                        else if (p_r[j] > fMaxFofContained[j]) {
                            *pbCurrFofContained = 0;
                            break;
                        }
                    }
                }
            }
        }
NoIntersect:
        if (sp) kdn = pkd->tree[iCell = S[--sp]];
        else break;
    }
    return (iTail);
}

static void iOpenRemoteFof(PKD pkd,treeStore::NodePointer k,clTile &tile,float dTau2) {
    float dx,minbnd2,kOpen;
    int i,iOpen;

    auto kbnd = k->bound();
    kOpen = sum(kbnd.apothem()); /* Manhatten metric */

    auto nBlocks = tile.count() / tile.width;
    for (auto iBlock=0; iBlock<=nBlocks; ++iBlock) {
        int n = iBlock<nBlocks ? tile.width : tile.count() - nBlocks*tile.width;
        clBlock &blk = tile[iBlock];
        for (i=0; i<n; ++i) {
            if (blk.idCell[i] > pkd->Self()) iOpen = 10;  /* ignore this cell, but this never ignores the top tree */
            else {
                minbnd2 = 0;
                dx = kbnd.lower(0) -  blk.xCenter[i] - blk.xOffset[i] - blk.xMax[i];
                if (dx > 0) minbnd2 += dx*dx;
                dx = blk.xCenter[i] + blk.xOffset[i] - blk.xMax[i] - kbnd.upper(0);
                if (dx > 0) minbnd2 += dx*dx;
                dx = kbnd.lower(1) - blk.yCenter[i] - blk.yOffset[i] - blk.yMax[i];
                if (dx > 0) minbnd2 += dx*dx;
                dx = blk.yCenter[i] + blk.yOffset[i] - blk.yMax[i] - kbnd.upper(1);
                if (dx > 0) minbnd2 += dx*dx;
                dx = kbnd.lower(2) - blk.zCenter[i] - blk.zOffset[i] - blk.zMax[i];
                if (dx > 0) minbnd2 += dx*dx;
                dx = blk.zCenter[i] + blk.zOffset[i] - blk.zMax[i] - kbnd.upper(2);
                if (dx > 0) minbnd2 += dx*dx;
                if (minbnd2 > dTau2) iOpen = 10;  /* ignore this cell */
                else if (k->is_bucket()) {
                    if (blk.iLower[i] == 0) iOpen = 1;
                    else iOpen = 3;
                }
                else if (kOpen > blk.cOpen[i] || blk.iLower[i] == 0) iOpen = 0;
                else iOpen = 3;
            }
            blk.iOpen[i] = iOpen;
        }
    }
}

static void addChildFof(PKD pkd, remoteTree &tree, clList *cl, int iChild, int id, TinyVector<float,3> fOffset) {
    auto c = tree(iChild,id);
    auto c_r = c->position();
    auto cbnd = c->bound();
    auto iCache = 0;
    auto cOpen = blitz::sum(cbnd.apothem()); /* Manhatten metric */
    auto [iLower, idLower, iUpper, idUpper] = c->get_child_cells(id);
    const auto &SPHbob = c->have_BOB() ? c->BOB() : SPHBOB();
    cl->append(iCache,id,iChild,idLower,iLower,idUpper,iUpper,c->count(),cOpen,
               c->mass(),4.0f*c->fSoft2(),c_r,fOffset,cbnd,SPHbob);
}

void pkdFofRemoteSearch(PKD pkd,double dTau2,int bPeriodic,int nReplicas,int nBucket) {
    double d2;
    int npi;
    uint32_t pjGroup;
    uint32_t pi,pj;
    int iRemote;
    int i,j,ix,iy,iz,bRep;
    int idSelf,iTop,iCell,id,iSib,iCheckCell,iCheckLower;
    int jTile,M,iStack;
    TinyVector<float,3> fOffset;

    remoteTree tree(pkd->mdl,pkd->tree,CID_CELL);
    remoteParticles particles(pkd->mdl,pkd->particles,CID_PARTICLE);

    /*
    ** Allocate the vectors to be large enough to handle all particles in a bucket.
    ** M should be aligned to 4*sizeof(double)
    */
    M = nBucket;
    auto ri = std::make_unique<TinyVector<double,3>[]>(M);
    auto piGroup = std::make_unique<uint32_t[]>(M);

    iStack = 0;
    pkd->cl->clear();

    auto bndSelf = pkd->tree[ROOT]->bound();
    idSelf = mdlSelf(pkd->mdl);
    iTop = pkd->iTopTree[ROOT];
    id = idSelf;
    for (j=0; j<3; ++j) fOffset[j] = 0.0;
    /*
    ** Add all siblings of the top tree down to local root (but not including it) to
    ** the checklist.
    */
    iCell = iTop;
    auto c = pkd->tree[iCell];
    while (c->is_top_tree()) {
        auto [iCellLo,idLo,iCellUp,idUp] = c->get_child_cells(id);
        if (idLo == idSelf) {
            c = pkd->tree[iCellLo];
            auto bnd = c->bound();
            if (any(fabs(bndSelf.center() - bnd.center()) > bnd.apothem())) {
                addChildFof(pkd,tree,pkd->cl,iCellLo,idLo,fOffset);
                id = idUp;
                assert(id == idSelf);
                iCell = iCellUp;
                c = pkd->tree[iCell];
                goto NextCell;
            }
        }
        addChildFof(pkd,tree,pkd->cl,iCellUp,idUp,fOffset);
        iCell = iCellLo;
        id = idLo;
        assert(id == idSelf);
NextCell:
        ;
    }
    /*
    ** Add all replica global roots to the checklist for periodic BCs.
    */
    if (bPeriodic) {
        int nReps = nReplicas;
        for (ix=-nReps; ix<=nReps; ++ix) {
            fOffset[0] = ix*pkd->fPeriod[0];
            for (iy=-nReps; iy<=nReps; ++iy) {
                fOffset[1] = iy*pkd->fPeriod[1];
                for (iz=-nReps; iz<=nReps; ++iz) {
                    fOffset[2] = iz*pkd->fPeriod[2];
                    bRep = ix || iy || iz;
                    if (bRep) addChildFof(pkd,tree,pkd->cl,iTop,idSelf,fOffset);
                }
            }
        }
    }

    iCell = ROOT;
    id = idSelf;
    auto k = pkd->tree[iCell];
    /*
    ** The checklist is ready for a bound-bound walk of the remote particles.
    */
    while (1) {
        while (1) {
            /*
            ** Process the Checklist for the cell pointed to by k.
            */
            pkd->S[iStack+1].cl->clear();
            do {
                for (auto &tile : *pkd->cl) {
                    iOpenRemoteFof(pkd,k,tile,dTau2);
                }
                pkd->clNew->clear();
                for (auto &tile : *pkd->cl) {

                    auto nBlocks = tile.count() / tile.width;
                    for (auto iBlock=0; iBlock<=nBlocks; ++iBlock) {
                        int n = iBlock<nBlocks ? tile.width : tile.count() - nBlocks*tile.width;
                        clBlock &blk = tile[iBlock];
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
                                ** We check individual particles against each other here.
                                */
                                iCheckCell = blk.iCell[jTile];
                                id = blk.idCell[jTile];
                                c = tree(iCheckCell,id);
                                /*
                                ** Convert all the coordinates in the k-cell and store them in vectors.
                                */
                                for (pi=k->lower(),npi=0; pi<=k->upper(); ++pi,++npi) {
                                    auto p = pkd->particles[pi];
                                    ri[npi] = p.position();
                                    piGroup[npi] = p.group();
                                }
                                for (pj=c->lower(); pj<=c->upper(); ++pj) {
                                    auto p = particles(pj,id);
                                    pjGroup = p.group();
                                    auto rj = p.position();
                                    rj[0] += blk.xOffset[jTile];
                                    rj[1] += blk.yOffset[jTile];
                                    rj[2] += blk.zOffset[jTile];
                                    /*
                                    ** The following could be vectorized over the vectors xi,yi and zi!
                                    */
                                    for (i=0; i<npi; ++i) {
                                        auto r = rj - ri[i];
                                        d2 = dot(r,r);
                                        if (d2 < dTau2) {
                                            /*
                                            ** We have found a remote group link!
                                            ** Check if it is already linked to this group and if not add the link.
                                            */
                                            if (pjGroup == 0) printf("UNGROUPED PARTICLE FOUND at %.15g < %.15g\n",d2,dTau2);
                                            assert(pjGroup > 0);
                                            iRemote = pkd->ga[piGroup[i]].iLink;
                                            while (iRemote) {
                                                if (pkd->tmpFofRemote[iRemote].key.iIndex == pjGroup &&
                                                        pkd->tmpFofRemote[iRemote].key.iPid == id) break;
                                                iRemote = pkd->tmpFofRemote[iRemote].iLink;
                                            }
                                            if (!iRemote) {
                                                /*
                                                ** Add this remote group to the list.
                                                */
                                                iRemote = pkd->iRemoteGroup++;
                                                /*
                                                ** Make sure we don't run out of ephemeral storage!
                                                ** (we really should be ok as long as our nMinMembers is greater than 1 or 2)
                                                */
                                                assert(pkd->iRemoteGroup < pkd->nMaxRemoteGroups);
                                                pkd->tmpFofRemote[iRemote].key.iIndex = pjGroup;
                                                assert(pjGroup > 0);
                                                pkd->tmpFofRemote[iRemote].key.iPid = id;
                                                pkd->tmpFofRemote[iRemote].iLink = pkd->ga[piGroup[i]].iLink;
                                                pkd->ga[piGroup[i]].iLink = iRemote;
                                            }
                                        }
                                    }
                                }
                                break;
                            case 3:
                                /*
                                ** Open the cell.
                                ** We could do a prefetch here for non-local
                                §              ** cells.
                                */
                                iCheckCell = blk.iCell[jTile];    assert(iCheckCell >= 0);
                                iCheckLower = blk.iLower[jTile];  assert(iCheckLower > 0);

                                fOffset[0] = blk.xOffset[jTile];
                                fOffset[1] = blk.yOffset[jTile];
                                fOffset[2] = blk.zOffset[jTile];

                                addChildFof(pkd,tree,pkd->clNew,blk.iLower[jTile],blk.idLower[jTile],fOffset);
                                addChildFof(pkd,tree,pkd->clNew,blk.iUpper[jTile],blk.idUpper[jTile],fOffset);
                                break;
                            case 10:
                                /*
                                ** This checkcell is removed from the checklist since it has no overlap with the current cell.
                                */
                                break;
                            default:
                                assert(0);
                            } /* end of switch */
                        } /* end of for (jTile) */
                    } /* end of for (nLeft) */
                } /* end of CL_LOOP */
                std::swap(pkd->cl,pkd->clNew);
            } while (pkd->cl->count());
            /*
            ** Done processing of the Checklist.
            ** Now prepare to proceed to the next deeper
            ** level of the tree.
            */
            if (k->is_bucket()) break;
            iCell = k->lchild();
            iSib = k->rchild();
            k = pkd->tree[iCell];
            /*
            ** Push the sibling onto the stack.
            */
            c = pkd->tree[iSib];
            pkd->cl->clone(*pkd->S[iStack+1].cl);
            ++iStack;
            assert(iStack < pkd->nMaxStack);
            pkd->S[iStack].iNodeIndex = iSib;
        }
        /*
        ** Now the checklist should be empty and we should have dealt with all
        ** links going across the processor domains for this bucket.
        */
        assert(pkd->S[iStack+1].cl->count()==0);
        /*
        ** Get the next cell to process from the stack.
        */
        if (iStack == -1) break;  /* we are done! */
        k = pkd->tree[iCell = pkd->S[iStack].iNodeIndex];
        /*
        ** Grab the checklist from the stack.
        */
        std::swap(pkd->cl,pkd->S[iStack].cl);

        --iStack;
    }
}

#if 0
void updateInteriorBound(Bound &ibnd,const Bound &bnd) {
    int j;
    for (j=0; j<3; ++j) {
        double bmin = bnd.lower(j);
        double bmax = bnd.upper(j);
        double fmin = ibnd.lower(j);
        double fmax = ibnd.upper(j);
        assert(bmin <= bmax);
        if (bmin > fmax || bmax < fmin) continue;
        else if (bmax < fmax && bmin > fmin) {
            if (bmin - fmin > fmax - bmax)
                fmax = bmin;
            else
                fmin = bmax;
        }
        else {
            if (bmin < fmax) fmax = bmin;
            if (bmax > fmin) fmin = bmax;
        }
        ibnd.set(j,fmin,fmax);
    }
}
#endif

void pkdNewFof(PKD pkd,double dTau2,int nMinMembers,int bPeriodic,int nReplicas,int nBucket) {
    uint32_t iGroup;
    int pn,i,j;
    uint32_t iHead;
    uint32_t iTail;
    const uint32_t uGroupMax = (pkd->bNoParticleOrder)?IGROUPMAX:0xffffffff;
    uint32_t *Fifo;
    int bCurrFofContained;

    assert(pkd->particles.present(PKD_FIELD::oGroup) || pkd->bNoParticleOrder); /* Validate memory model */
    auto S = new int[1024]; assert(S);
    /*
    ** Set up the bounds for the FOF groups that are certainly contained in the domain.
    ** This is a little trickier for domains which could potentially overlap a bit.
    ** For now I assume that the domains do NOT overlap, but the calculation for overlapping
    ** domains just involves a tree walk.
    */
    auto bndSelf = pkd->tree[ROOT]->bound();
    pkd->bndInterior = bndSelf;
#if 0
    /*
    ** The following code would allow doing fof on a substep, which isn't
    ** forseen in the near future. We can test it at a later stage.
    ** Check bounds against all siblings of the top tree down to local root.
    */
    iCell = iTop;
    c = pkd->tree[iCell];
    auto bndTop = c->bound();
    while (c->bTopTree) {
        auto [iCellLo,idLo,iCellUp,idUp] = c->get_child_cells(id);
        if (idLo == idSelf) {
            c = pkd->tree[iCellLo];
            auto bnd = c->bound();
            for (j=0; j<3; ++j) {
                if (fabs(bndSelf.fCenter[j]-bnd.fCenter[j]) > bnd.fMax[j]) {
                    /*
                    ** Check bounds against this sibling.
                    */
                    updateInteriorBound(pkd->bndInterior,bnd);
                    id = idUp;
                    assert(id == idSelf);
                    iCell = iCellUp;
                    c = pkd->tree[iCell];
                    goto NextCell;
                }
            }
        }
        assert(idUp == idSelf);
        c = pkd->tree[iCellUp];
        bnd = c->bound();
        /*
        ** Check bounds against this sibling.
        */
        updateInteriorBound(pkd->bndInterior,bnd);
        iCell = iCellLo;
        id = idLo;
NextCell:
        ;
    }
    /*
    ** Check bounds against first replica global roots for periodic BCs.
    */
    if (bPeriodic) {
        Bound rbnd;
        for (j=0; j<3; ++j) rbnd.fMax[j] = bndTop->fMax[j];
        for (ix=-1; ix<=1; ++ix) {
            fOffset[0] = ix*pkd->fPeriod[0];
            for (iy=-1; iy<=1; ++iy) {
                fOffset[1] = iy*pkd->fPeriod[1];
                for (iz=-1; iz<=1; ++iz) {
                    fOffset[2] = iz*pkd->fPeriod[2];
                    bRep = ix || iy || iz;
                    if (bRep) {
                        /*
                        ** Check bounds against this replica.
                        */
                        for (j=0; j<3; ++j)
                            rbnd.fCenter[j] = bndTop.fCenter[j] + fOffset[j];
                        updateInteriorBound(pkd->bndInterior,rbnd);
                    }
                }
            }
        }
    }
#endif
    bndSelf = pkd->bndInterior;
    /*
    ** Finally make the contained region be dTau smaller on each side.
    */
    TinyVector<double,3> fMinFofContained = bndSelf.lower() + sqrt(dTau2);
    TinyVector<double,3> fMaxFofContained = bndSelf.upper() - sqrt(dTau2);
    /*
    ** Clear the group numbers!
    */
    for (auto &p : pkd->particles) {
        p.set_group(0);
    }
    /*
    ** The following *just* fits into ephemeral storage of 4 bytes/particle.
    */
    assert(pkd->EphemeralBytes() >= 4);
    Fifo = (uint32_t *)(pkd->pLite);
    iGroup = 1;
    for (pn=0; pn<pkd->Local(); ++pn) {
        auto p = pkd->particles[pn];
        if (p.group()) continue;
        /*
        ** Mark particle and add it to the do-fifo
        */
        iHead = iTail = 0;
        Fifo[iTail++] = pn;
        assert(iGroup < uGroupMax);
        p.set_group(iGroup);
        bCurrFofContained = 1;
        auto p_r = p.position();
        for (j=0; j<3; ++j) {
            if (p_r[j] < fMinFofContained[j]) {
                bCurrFofContained = 0;
                break;
            }
            else if (p_r[j] > fMaxFofContained[j]) {
                bCurrFofContained = 0;
                break;
            }
        }
        while (iHead < iTail) {
            p_r = pkd->particles[Fifo[iHead++]].position();
            iTail = pkdFofGatherLocal(pkd,S,dTau2,p_r,iGroup,iTail,Fifo,
                                      &bCurrFofContained,fMinFofContained,fMaxFofContained);
        }
        assert(iTail <= pkd->Local());
        /*
        ** Now check if this fof group is contained and has fewer than nMinFof particles.
        */
        if (bCurrFofContained && iTail < nMinMembers) {
            /*
            ** In this case mark the group particles as belonging to a removed group.
            */
            for (iHead=0; iHead<iTail; ++iHead) {
                pkd->particles[Fifo[iHead]].set_group(uGroupMax);
            }
        }
        else {
            ++iGroup;
        }
    }
    /*
    ** Clear group ids for removed small groups.
    */
    for (auto &p : pkd->particles) {
        if (p.group() == uGroupMax) p.set_group(0);
    }
    pkd->nLocalGroups = iGroup-1;
    pkd->nGroups = pkd->nLocalGroups + 1;
    delete [] S;  /* this stack is no longer needed */
    Fifo = NULL;  /* done with the Fifo, can use the storage for other stuff now */
    /*
    ** Create initial group table. The assert below is a very minimal requirement as it doesn't account for remote
    ** links (tmpFofRemote). However, we check this again everytime we add a new remote link.
    */
    assert(sizeof(*pkd->ga)*pkd->nGroups+sizeof(*pkd->tmpFofRemote) <= 1ul*pkd->EphemeralBytes()*pkd->FreeStore());
    pkd->nMaxRemoteGroups = (1ul*pkd->EphemeralBytes()*pkd->FreeStore() - sizeof(*pkd->ga)*pkd->nGroups) / sizeof(*pkd->tmpFofRemote);
    pkd->ga = (struct smGroupArray *)(pkd->pLite);
    pkd->tmpFofRemote = (FOFRemote *)&pkd->ga[pkd->nGroups];
    for (i=0; i<pkd->nGroups; ++i) {
        pkd->ga[i].id.iIndex = i;
        pkd->ga[i].id.iPid = pkd->Self();
        pkd->ga[i].iGid = i;
        pkd->ga[i].iLink = 0;   /* this is a linked list of remote groups linked to this local group */
        pkd->ga[i].minPot = FLOAT_MAXVAL;
        pkd->ga[i].iMinPart = 0xffffffff;
    }
    /*
    ** Set a local reference point for each group.
    */
    for (pn=0; pn<pkd->Local(); ++pn) {
        auto p = pkd->particles[pn];
        if ( (i = p.group()) != 0 ) {
            if (pkd->ga[i].iMinPart == 0xffffffff) {
                pkd->ga[i].iMinPart = pn;
                pkd->ga[i].minPot = (float)pkd->Self(); /* this makes the reference particle be in the lowest processor number */
            }
        }
        /*
        ** Note that IF we calculate gravity this reference gets overwritten with
        ** the true minimum potential particle (since the potentials are always
        ** negative). Otherwise it remains a useful reference point for the group.
        ** This is useful in case where a group straddles a periodic boundary.
        */
    }
    pkd->iRemoteGroup = 1;  /* The first entry is a dummy one for a null index */
    /*
    ** Now lets go looking for local particles which have a remote neighbor that is part of
    ** a group.
    */
    pkd->mdl->CacheInitialize(CID_PARTICLE,NULL,pkd->particles,pkd->Local(),pkd->particles.ParticleSize());
    pkdFofRemoteSearch(pkd,dTau2,bPeriodic,nReplicas,nBucket);
    pkd->mdl->FinishCache(CID_PARTICLE);
}

class FetchNames : public mdl::CACHEhelper {
protected:
    virtual void pack(void *dst, const void *src) override {
        auto g1 = static_cast<remoteID *>(dst);           // Packed value
        auto g2 = static_cast<smGroupArray const *>(src); // Regular element
        g1->iPid = g2->id.iPid;
        g1->iIndex = g2->id.iIndex;
    }
    virtual void  unpack(void *dst, const void *src, const void *key) override {
        auto g1 = static_cast<smGroupArray *>(dst);   // Regular element
        auto g2 = static_cast<remoteID const *>(src); // Packed value
        g1->id.iPid = g2->iPid;
        g1->id.iIndex = g2->iIndex;
    }
    virtual uint32_t pack_size()  override {return sizeof(remoteID);}
public:
    explicit FetchNames() : CACHEhelper(sizeof(struct smGroupArray)) {}
};

/*
** When we virtual fetch a name of one of the groups we may already have fetched the
** same key on the same processor. For this reason we need to initialize the
** virtual fetch to something that will for sure be updated on the first fetch.
** the second fetch will only update it if the new name is actually smaller than the
** one set by the first fetch.
*/
class PropagateNames : public mdl::CACHEhelper {
protected:
    virtual void pack(void *dst, const void *src) override {
        assert(0); abort(); // We use virtual fetch only, so pack() is not used
    }
    virtual void  unpack(void *dst, const void *src, const void *key) override {
        assert(0); abort(); // We use virtual fetch only, so unpack() is not used
    }
    // virtual uint32_t pack_size()  override {return sizeof(remoteID);}
    virtual void init(void *dst) override {
        auto g = static_cast<smGroupArray *>(dst);
        g->id.iPid = INT32_MAX;
        g->id.iIndex = INT32_MAX;
    }
    virtual void combine(void *dst, const void *src,const void *key) override {
        auto g1 = static_cast<smGroupArray *>(dst);   // Regular element
        auto g2 = static_cast<remoteID const *>(src); // Packed value
        if ( g1->id.iPid>g2->iPid || (g1->id.iPid==g2->iPid && g1->id.iIndex>g2->iIndex) ) {
            g1->id.iPid = g2->iPid;
            g1->id.iIndex = g2->iIndex;
        }
    }
    virtual void flush(void *dst, const void *src) override {
        auto g1 = static_cast<remoteID *>(dst);           // Packed value
        auto g2 = static_cast<smGroupArray const *>(src); // Regular element
        g1->iPid = g2->id.iPid;
        g1->iIndex = g2->id.iIndex;
    }
    virtual uint32_t flush_size() override {return sizeof(remoteID);}

public:
    explicit PropagateNames() : CACHEhelper(sizeof(struct smGroupArray),true) {}
};

int pkdFofPhases(PKD pkd) {
    MDL mdl = pkd->mdl;
    int bMadeProgress=0;
    int iIndex,iPid,iRemote,iLink,i;
    remoteID name;
    struct smGroupArray *pRemote;

    /*
    ** Phase 1: fetch remote names.
    */
    pkd->mdl->CacheInitialize(CID_GROUP,NULL,pkd->ga,pkd->nGroups,std::make_shared<FetchNames>());
    for (iRemote=1; iRemote<pkd->iRemoteGroup; ++iRemote) {
        iIndex = pkd->tmpFofRemote[iRemote].key.iIndex;
        assert(iIndex > 0);
        iPid = pkd->tmpFofRemote[iRemote].key.iPid;
        pRemote = static_cast<struct smGroupArray *>(mdlFetch(mdl,CID_GROUP,iIndex,iPid));
        /*
        ** Now update the name in the local table of remote group links.
        */
        assert(pRemote->id.iIndex > 0);
        pkd->tmpFofRemote[iRemote].name.iIndex = pRemote->id.iIndex;
        pkd->tmpFofRemote[iRemote].name.iPid = pRemote->id.iPid;
    }
    pkd->mdl->FinishCache(CID_GROUP);
    /*
    ** Phase 2: update to unique names.
    */
    for (i=1; i<pkd->nGroups; ++i) {
        iLink = pkd->ga[i].iLink;
        if (iLink) {
            name.iIndex = pkd->ga[i].id.iIndex;
            name.iPid = pkd->ga[i].id.iPid;
            /*
            ** Find current master name (this is the lowest iPid,iIndex pair found).
            */
            while (iLink) {
                iPid = pkd->tmpFofRemote[iLink].name.iPid;
                iIndex = pkd->tmpFofRemote[iLink].name.iIndex;
                iLink = pkd->tmpFofRemote[iLink].iLink;
                if (iPid < name.iPid) {
                    name.iPid = iPid;
                    name.iIndex = iIndex;
                    bMadeProgress = 1;
                }
                else if (iPid == name.iPid && iIndex < name.iIndex) {
                    name.iIndex = iIndex;
                    bMadeProgress = 1;
                }
            }
            assert(name.iIndex > 0);
            pkd->ga[i].id.iIndex = name.iIndex;
            pkd->ga[i].id.iPid = name.iPid;
            iLink = pkd->ga[i].iLink;
            while (iLink) {
                if (pkd->tmpFofRemote[iLink].name.iPid == name.iPid &&
                        pkd->tmpFofRemote[iLink].name.iIndex == name.iIndex) {
                    /*
                    ** There is no update to be made. We mark this with a
                    ** a dummy iPid! It is ok to destroy it here since it
                    ** will be refetched in phase 1 (all remote links names
                    ** are refetched).
                    */
                    pkd->tmpFofRemote[iLink].name.iPid = -1;
                }
                else {
                    pkd->tmpFofRemote[iLink].name.iPid = name.iPid;
                    pkd->tmpFofRemote[iLink].name.iIndex = name.iIndex;
                }
                iLink = pkd->tmpFofRemote[iLink].iLink;
            }
        }
    }
    /*
    ** Phase 3: propagate.
    */
    pkd->mdl->CacheInitialize(CID_GROUP,NULL,pkd->ga,pkd->nGroups,std::make_shared<PropagateNames>());
    for (i=1; i<pkd->nGroups; ++i) {
        iLink = pkd->ga[i].iLink;
        if (iLink) {
            while (iLink) {
                name.iPid = pkd->tmpFofRemote[iLink].name.iPid;
                if (name.iPid >= 0) {
                    name.iIndex = pkd->tmpFofRemote[iLink].name.iIndex;
                    iPid = pkd->tmpFofRemote[iLink].key.iPid;
                    iIndex = pkd->tmpFofRemote[iLink].key.iIndex;
                    pRemote = static_cast<struct smGroupArray *>(mdlVirtualFetch(mdl,CID_GROUP,iIndex,iPid));
                    assert(pRemote);
                    if (name.iPid < pRemote->id.iPid) {
                        pRemote->id.iPid = name.iPid;
                        pRemote->id.iIndex = name.iIndex;
                    }
                    else if (name.iPid == pRemote->id.iPid && name.iIndex < pRemote->id.iIndex) {
                        pRemote->id.iIndex = name.iIndex;
                    }
                }
                iLink = pkd->tmpFofRemote[iLink].iLink;
            }
        }
    }
    pkd->mdl->FinishCache(CID_GROUP);
    return (bMadeProgress);
}

uint64_t pkdFofFinishUp(PKD pkd,int nMinGroupSize) {
    int i;

    /*
    ** Merge local groups in order to have only one group per process informing the master of its existence.
    ** For a given name in the table, we need to know if some other entry has the same name and remap that
    ** index.
    */
    for (i=1; i<pkd->nGroups; ++i) {
        assert(pkd->ga[i].id.iIndex > 0);
    }
    pkd->nGroups = pkdGroupCombineDuplicateIds(pkd,pkd->nGroups,pkd->ga,1);
    pkd->nGroups = pkdPurgeSmallGroups(pkd,pkd->nGroups,pkd->ga,nMinGroupSize);
    /*
    ** This last call to pkdGroupCounts just gets the final group counts in pkd->ga restored
    ** which is optional if group stats are calculated later (which does the same thing again).
    ** pkd->nLocalGroups is kept consistent within pkdPurgeSmallGroups as well.
    */
    pkdGroupCounts(pkd,pkd->nGroups,pkd->ga);
    return pkd->nLocalGroups;
}
