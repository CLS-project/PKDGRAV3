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
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
    #include <malloc.h>
#endif
#include <assert.h>
#include "pkd.h"
#include "gravity/moments.h"
#include "../SPH/SPHOptions.h"

#ifdef HAVE_SYS_TIME_H
    #include <sys/time.h>
#endif
#ifdef USE_ITT
    #include "ittnotify.h"
#endif
#include <vector>
#include <numeric>
#include <algorithm>

uint32_t pkdDistribTopTree(PKD pkd, uint32_t uRoot, uint32_t nTop, KDN *pTop, int allocateMemory) {
    int i, iTop;
    auto LocalRoot = pkd->tree[uRoot];

    if (allocateMemory) {
        iTop = pkd->iTopTree[uRoot] = pkd->tree.AllocNode(nTop);
    }
    else {
        iTop = pkd->iTopTree[uRoot];
    }
    for (i=0; i<nTop; ++i) {
        treeStore::NodeReference node(pkd->tree,pTop,i);
        auto Local = pkd->tree[iTop+i];
        *Local = node;
        // Local->copy(node);
        assert(Local->is_top_tree());
        Local->set_group(false);
        if (Local->is_local()) { /* Fixup the links if necessary */
            (*Local)->iLower += iTop + i;
            (*Local)->pUpper += iTop + i;
        }
        else if (Local->remote_id() == pkd->Self()) {
            Local->set_local(LocalRoot->lower(),LocalRoot->upper(),LocalRoot->lchild());
        }
    }
    return iTop;
}

static auto InitializeRootCommon(PKD pkd,uint32_t uRoot) {
    auto pRoot = pkd->tree[uRoot];
    pRoot->set_local();
    pRoot->set_bound(pkd->bnd);
    return pRoot;
}

/*
** Creates a single root at ROOT with only marked particles
*/
void pkdTreeInitMarked(PKD pkd) {
    auto pRoot    = InitializeRootCommon(pkd,ROOT);
    auto pRootFix = InitializeRootCommon(pkd,FIXROOT);

    if (mdlCacheStatus(pkd->mdl,CID_CELL)) mdlFinishCache(pkd->mdl,CID_CELL);

    auto i = std::partition(pkd->particles.begin(),pkd->particles.end(),[](auto &p) {return p.marked();}) - pkd->particles.begin();
    pRoot->set_local(0,i-1);
    pRootFix->set_local(i,pkd->Local() - 1);
    pkd->tree.SetNodeCount(NRESERVED_NODES);
    pRoot->set_group(true);
    pRootFix->set_group(true);
}

void pkdDumpTrees(PKD pkd,int bOnlyVA, uint8_t uRungDD) {
    auto pRoot    = InitializeRootCommon(pkd,ROOT);

    if (pkd->Nodes() == 0) {
        assert(bOnlyVA==0);
        pkd->tree.SetNodeCount(NRESERVED_NODES);
        pkd->tree[ROOT]->set_local();
        pkd->tree[FIXROOT]->set_local();
    }

    /* We can always close the "active" caches */
    if (mdlCacheStatus(pkd->mdl,CID_CELL)) mdlFinishCache(pkd->mdl,CID_CELL);
// done elsewhere   if (mdlCacheStatus(pkd->mdl,CID_PARTICLE)) mdlFinishCache(pkd->mdl,CID_PARTICLE);

    /* Full Normal tree build. Only ROOT will be used */
    if (uRungDD == IRUNGMAX) {
        pRoot->set_local(0,pkd->Local() - 1);
        pRoot->set_group(true);
        pkd->tree.SetNodeCount(NRESERVED_NODES);
    }
    /* Just rebuilding (active) ROOT. Truncate it. pLower and pUpper are still valid. */
    else if (bOnlyVA) {
        if (pRoot->is_cell()) {
            assert(pRoot->lchild() >= NRESERVED_NODES);
            pkd->tree.SetNodeCount(pRoot->lchild()); /* This effectively truncates the nodes used by this tree */
            pRoot->set_local(pRoot->lower(),pRoot->upper());
            pRoot->set_group(true);
        }
    }

    /* Need to build two trees which is more complicated. */
    else {
        auto pRootFix = InitializeRootCommon(pkd,FIXROOT);

#ifndef SINGLE_CACHES
        /* Here we also have to close the "fixed" caches */
        if (mdlCacheStatus(pkd->mdl,CID_CELL2))     mdlFinishCache(pkd->mdl,CID_CELL2);
        if (mdlCacheStatus(pkd->mdl,CID_PARTICLE2)) mdlFinishCache(pkd->mdl,CID_PARTICLE2);
#endif
        auto iLast = std::partition(pkd->particles.begin(),pkd->particles.end(),
        [uRungDD](auto &p) {return p.rung() < uRungDD;})
        - pkd->particles.begin();

        pkd->tree.SetNodeCount(NRESERVED_NODES);
        pRoot->set_local(iLast+1,pkd->Local() - 1);
        pRoot->set_group(true);
        pRootFix->set_local(0,iLast);
        pRootFix->set_group(true);
    }
}

#define MIN_SRATIO    0.05


/*
** Partition the particles between pLower and pUpper (inclusive)
** around "Split" in dimension "d". Return the index of the split.
*/
template<class split_t>
static int PartPart(PKD pkd,int pLower,int pUpper,int d,split_t Split) {
    auto pi = pkd->particles.begin() + pLower;
    auto pj = pkd->particles.begin() + pUpper;
    return std::partition(pi,pj+1,[Split,d](auto &p) {return p.template raw_position<split_t>(d) < Split;}) - pkd->particles.begin();
}

/*
** M is the bucket size.
** This function assumes that the root node is correctly set up (particularly the bounds).
*/
void BuildTemp(PKD pkd,int iNode,int M,int nGroup,double dMaxMax) {
    auto pNode = pkd->tree[iNode];
    double lrMax;
    std::vector<int> S;     /* this is the stack */
    int d;
    int i;
    int nr,nl;
    int lc,rc;
    int nBucket = 0;

    pNode->set_depth(0);
    pNode->set_split_dim(3);

    // Single bucket? We are done.
    if (pNode->count() <= M) return;

    S.reserve(100); // Start with a sensible stack size

    auto bnd = pNode->bound();
    while (1) {
        // Begin new stage! Calculate the appropriate fSplit.
        // Pick longest dimension and split it in half.
        // Calculate the new left and right cells. We will use
        // either the left, the right or (usually) both.
        pNode->set_split_dim(d = bnd.maxdim());
        auto [lbnd,rbnd] = bnd.split(d);
        lrMax = 0.5*lbnd.maxside();
        // Now start the partitioning of the particles about
        // fSplit on dimension given by d.
        if (pkd->bIntegerPosition) {
            int32_t Split = pkdDblToIntPos(pkd,bnd.center(d));
            i = PartPart(pkd,pNode->lower(),pNode->upper(),d,Split);
        }
        else i = PartPart(pkd,pNode->lower(),pNode->upper(),d,bnd.center(d));
        nl = i - pNode->lower();
        nr = pNode->upper() + 1 - i;
        if (nl > 0 && nr > 0) {
            // Split this node into two children
            auto [pLeft,pRight] = pNode->split(i);

            pNode->set_group(pNode->count() <= nGroup);
            pLeft->set_bound(lbnd);
            pRight->set_bound(rbnd);

            /*
            ** Now figure out which subfile to process next.
            */
            if (lrMax > dMaxMax) {
                lc = (nl > 0); /* this condition means the left child is not a bucket */
                rc = (nr > 0);
            }
            else {
                lc = (nl > M); /* this condition means the left child is not a bucket */
                rc = (nr > M);
            }
            if (rc && lc) {
                if (nr > nl) {
                    S.push_back(pNode->rchild());   // push tr
                    iNode = pNode->lchild();        // process lower subfile
                }
                else {
                    S.push_back(pNode->lchild());   // push tl
                    iNode = pNode->rchild();        // process upper subfile
                }
            }
            else if (lc) {
                /*
                ** Right must be a bucket in this case!
                */
                iNode = pNode->lchild();   /* process lower subfile */
                pRight->set_group(true);
                ++nBucket;
            }
            else if (rc) {
                /*
                ** Left must be a bucket in this case!
                */
                iNode = pNode->rchild();   /* process upper subfile */
                pLeft->set_group(true);
                ++nBucket;
            }
            else {
                /*
                ** Both are buckets (we need to pop from the stack to get the next subfile.
                */
                pLeft->set_group(true);
                ++nBucket;
                pRight->set_group(true);
                ++nBucket;
                if (S.empty()) break;
                iNode = S.back();
                S.pop_back();
            }
        }
        else {
            pNode->set_depth(pNode->depth()+1);
            // No nodes allocated: just change the bounds
            if (nl > 0) {
                pNode->set_bound(lbnd);
                lc = (lrMax>dMaxMax || nl > M); /* this condition means the node is not a bucket */
                if (!lc) {
                    pNode->set_group(true);
                    ++nBucket;
                    if (S.empty()) break;
                    iNode = S.back();
                    S.pop_back();
                }
            }
            else {
                pNode->set_bound(rbnd);
                rc = (lrMax>dMaxMax || nr > M);
                if (!rc) {
                    pNode->set_group(true);
                    ++nBucket;
                    if (S.empty()) break;
                    iNode = S.back();
                    S.pop_back();
                }
            }
        }
        pNode = pkd->tree[iNode];
        bnd = pNode->bound();
    }
}


/*
** With more than a single tree, we must be careful to make sure
** that they match or the P-P, P-C and hence checklists will explode.
** This routine takes a template tree (usually the large "fixed" tree)
** and contructs a tree with buckets no larger than the corresponding
** buckets in the template tree. Once we reach a bucket in the template
** tree, then we can either create a bucket, or continue to create the
** new tree in the usual way.
*/
void BuildFromTemplate(PKD pkd,int iNode,int M,int nGroup,int iTemplate) {
    int d;
    int i, nr, nl;
    double dMaxMax;
    struct buildStack {
        int iTemp; // Template
        int iNode; // New tree
        buildStack() = default;
        buildStack(int iTemp, int iNode) : iTemp(iTemp), iNode(iNode) {}
    };
    std::vector<buildStack> S;
    // If we have an empty tree root then we are "done"
    auto pNode = pkd->tree[iNode];
    pNode->set_depth(0);
    pNode->set_leaf();
    pNode->set_group(true);
    if (pNode->count() == 0) return;

    // Setup stack and push the two root cells of each tree.
    S.reserve(100); // Start with a sensible stack size

    S.emplace_back(iTemplate,iNode);
    auto pTemp = pkd->tree[iTemplate];
    dMaxMax = pTemp->bMax();

    while (!S.empty()) {
        // Pop the next cells to process
        pTemp = pkd->tree[S.back().iTemp];
        iNode = S.back().iNode;
        S.pop_back();

        pNode = pkd->tree[iNode];
        auto bnd = pNode->bound();

        // Follow the template tree to wherever it leads.
        while (pTemp->is_cell()) {
            d = pTemp->split_dim();
            if (pTemp->depth() != pNode->depth() || d>2) break;

            // Split is between left and right child nodes in the given dimension
            auto ptLeft = pkd->tree[pTemp->lchild()];
            auto ptRight = pkd->tree[pTemp->rchild()];
            auto ltbnd = ptLeft->bound();
            auto rtbnd = ptRight->bound();
            auto dSplit = 0.5 * (ltbnd.upper(d) + rtbnd.lower(d));

            // Partition the particles on either side of the split
            if (pkd->bIntegerPosition) {
                int32_t Split = pkdDblToIntPos(pkd,dSplit);
                i = PartPart(pkd,pNode->lower(),pNode->upper(),d,Split);
            }
            else i = PartPart(pkd,pNode->lower(),pNode->upper(),d,dSplit);
            nl = i - pNode->lower();
            nr = pNode->upper() + 1 - i;
            assert(nl>0 || nr>0);

            // Calculate bounding regions
            auto [lbnd,rbnd] = bnd.split(d,dSplit);
            assert(rbnd.width(d) > 0.0);
            assert(lbnd.width(d) > 0.0);

            if (nl==0) { // Particles on the right only
                pNode->set_depth(pNode->depth()+1);
                pNode->set_bound(rbnd);
                pTemp = ptRight;
                continue;
            }
            else if (nr==0) { // Particles on the left only
                pNode->set_depth(pNode->depth()+1);
                pNode->set_bound(lbnd);
                pTemp = ptLeft;
                continue;
            }

            // Good. We have particles on both sides, so we need to split this cell
            // and setup the appropriate bounds.
            auto [pLeft,pRight] = pNode->split(i);
            pLeft->set_bound(lbnd);
            pRight->set_bound(rbnd);

            pNode->set_group(pNode->count() <= nGroup && pTemp->is_group());
            pLeft->set_group(true);
            pRight->set_group(true);

            S.emplace_back(pTemp->rchild(),pNode->rchild());
            pTemp = ptLeft;
            iNode = pNode->lchild();
            pNode = pkd->tree[iNode];
            bnd = pNode->bound();
        }

        // Bucket in the template tree: Now just build, but set a sensible maximum cell size
        BuildTemp(pkd,iNode,M,nGroup,dMaxMax * pow(2.0,-pNode->depth()/3.0));
    }
}

void Create(PKD pkd,int iRoot,double ddHonHLimit) {
    int iNode = iRoot;
    blitz::TinyVector<double,3> kdn_r;
    double fSoft,d2,dih2,b;
    SPHBOB fBoB(0.0);
    float m, fMass, fBall;
    int pj,nDepth,ism;
    const int nMaxStackIncrease = 1;

    /* If the tree is empty, we just create a sensible moment and we are done. */
    auto pkdn = pkd->tree[iNode];
    if (pkdn->count()==0) {
        pkdn->bMax() = 1.0;
        pkdn->set_min_rung(MAX_RUNG);
        pkdn->set_max_rung(0);
        pkdn->set_marked(false);
        pkdn->set_NN(false);
        if (pkdn->have_moment()) momClearFmomr(&pkdn->moment());
        else pkdn->mass() = 0;
        return;
    }

    nDepth = 1;
    while (1) {
        while ((pkdn=pkd->tree[iNode])->is_cell()) {
            pkd->S[nDepth-1].iNodeIndex = iNode;
            iNode = pkdn->lchild();
            ++nDepth;
            /*
            ** Is this the deepest in the tree so far? We might need to have more stack
            ** elements for the tree walk!
            ** nMaxStack == nDepth guarantees that there is at least one deeper
            ** stack entry available than what is needed to walk the tree.
            */
            if (nDepth > pkd->nMaxStack) {
                pkd->S = static_cast<CSTACK *>(realloc(pkd->S,(pkd->nMaxStack+nMaxStackIncrease)*sizeof(CSTACK)));
                assert(pkd->S != NULL);
                for (ism=pkd->nMaxStack; ism<(pkd->nMaxStack+nMaxStackIncrease); ++ism) {
                    pkd->S[ism].cl = new clList(pkd->clFreeList);
                }
                pkd->nMaxStack += nMaxStackIncrease;
            }
        }

        // Now calculate all bucket quantities!
        // This includes M,CoM,Moments and special
        // bounds and iMaxRung.
        pkdn = pkd->tree[iNode];
        // Before squeezing the bounds, calculate a minimum b value based on the splitting bounds alone.
        // This gives us a better feel for the "size" of a bucket with only a single particle.
        auto bmin = pkdn->bound().minside();
        pkdn->update_bound();       // Update the bound to shrink wrap the particles
        auto bnd = pkdn->bound();

        pj = pkdn->lower();
        auto p = pkd->particles[pj];
        const auto a = p.have_acceleration() ? p.acceleration() : blitz::TinyVector<acc_t,3>(0.0);
        const auto v = p.have_velocity()     ? p.velocity()     : blitz::TinyVector<vel_t,3>(0.0);
        m = p.mass();
        fSoft = p.soft();
        fMass = m;
        dih2 = fSoft;
        auto r = p.position();

        if (pkd->particles.present(PKD_FIELD::oBall)) {
            /* initialize ball or box of balls */
            float limitedBallSize = std::min(pkd->SPHoptions.ballSizeLimit,pkd->SPHoptions.fBallFactor * p.ball());
            fBoB = SPHBOB(r,limitedBallSize);
        }
        /* initialize marked flag */
        pkdn->set_marked(p.marked());
        pkdn->set_NN(p.NN_flag());

        r *= m;

        blitz::TinyVector<vel_t,3>  vr = m * v;
        blitz::TinyVector<double,3> ar = m * a;
        pkdn->set_rung(p.rung());
        for (++pj; pj<=pkdn->upper(); ++pj) {
            auto p = pkd->particles[pj];
            const auto a = p.have_acceleration() ? p.acceleration() : blitz::TinyVector<acc_t,3>(0.0);
            const auto v = p.have_velocity()     ? p.velocity()     : blitz::TinyVector<vel_t,3>(0.0);
            m = p.mass();
            fSoft = p.soft();
            fMass += m;
            if (fSoft>dih2) dih2=fSoft;
            auto ft = p.position();

            if (pkd->particles.present(PKD_FIELD::oBall)) {
                float limitedBallSize = std::min(pkd->SPHoptions.ballSizeLimit,pkd->SPHoptions.fBallFactor * p.ball());
                fBoB = fBoB.combine(ft,limitedBallSize);
            }
            if (p.marked()) pkdn->set_marked(true);
            if (p.NN_flag()) pkdn->set_NN(true);

            r += m*ft;
            vr += m*v;
            ar += m*a;
            pkdn->set_max_rung(std::max(pkdn->max_rung(),p.rung()));
            pkdn->set_min_rung(std::min(pkdn->min_rung(),p.rung()));
        }
        m = 1.0f / fMass;
        kdn_r = m*r;
        pkdn->set_position(kdn_r);
        if (pkdn->have_velocity())     pkdn->velocity() = m*vr;
        if (pkdn->have_acceleration()) pkdn->acceleration() = m*ar;
        pkdn->fSoft2() = dih2*dih2;
        if (pkdn->have_BOB()) pkdn->BOB() = fBoB;

#ifdef USE_MAXSIDE
        b = bnd.maxside();
#else
        double d2Max = 0.0;
        for (auto &p : *pkdn) {
            r = p.position() - kdn_r;
            d2 = blitz::dot(r,r);
            d2Max = std::max(d2,d2Max);
        }
        b = sqrt(d2Max);
#endif
        if (b==0.0) b = 1.0f; /* FIXME: Single particle. Perhaps momMakeFmomr should be more robust. */
        else if (b < bmin) b = bmin;
        pkdn->bMax() = b;
        assert(b>=0);
        /*
        ** Now calculate the reduced multipole moment.
        ** Note that we use the cell's openening radius as the scaling factor!
        */
        if (pkdn->have_moment()) {
            auto &moment = pkdn->moment();
            momClearFmomr(&moment);
            for (auto &p : *pkdn) {
                FMOMR mom;
                auto ft = p.position();
                r = ft - kdn_r;
                m = p.mass();
                momMakeFmomr(&mom,m,pkdn->bMax(),r[0],r[1],r[2]);
                momAddFmomr(&moment,&mom);
            }
        }
        else pkdn->mass() = std::accumulate(pkdn->begin(),pkdn->end(),mass_t(0),[](mass_t a,auto &p) {return a + p.mass(); });
        /*
        ** Calculate bucket fast gas bounds.
        */
        if (pkdn->have_sph_bound()) {
            auto &bn = pkdn->sphbounds();
            /*
            ** Default bounds always makes the cell look infinitely far away, regardless from where.
            */
            bn.A.min = HUGE_VAL;
            bn.A.max = -HUGE_VAL;
            bn.B.min = HUGE_VAL;
            bn.B.max = -HUGE_VAL;
            bn.BI.min = HUGE_VAL;
            bn.BI.max = -HUGE_VAL;
            for (auto &p : *pkdn) {
                fBall = p.ball();
                if (p.is_gas()) {
                    auto r = p.position();
                    /*
                    ** This first ball bound over all gas particles is only used for remote searching.
                    */
                    bn.B.min = blitz::min(bn.B.min,r - (1+ddHonHLimit)*fBall);
                    bn.B.max = blitz::max(bn.B.max,r + (1+ddHonHLimit)*fBall);
                    if (p.is_active()) {
                        bn.A.min = blitz::min(bn.A.min,r);
                        bn.A.max = blitz::max(bn.A.max,r);
                    }
                    else {
                        bn.BI.min = blitz::min(bn.BI.min,r - (1+ddHonHLimit)*fBall);
                        bn.BI.max = blitz::max(bn.BI.max,r + (1+ddHonHLimit)*fBall);
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
        }
        /*
        ** Finished with the bucket, move onto the next one,
        ** or to the parent.
        */
        while ((iNode & 1) || iNode==iRoot ) {
            if ( --nDepth == 0) return; /* exit point!!! */
            iNode = pkd->S[nDepth-1].iNodeIndex;
            /*
            ** Now combine quantities from each of the children (2) of
            ** this cell to form the quantities for this cell.
            ** First find the CoM, just like for the bucket.
            */
            pkdn = pkd->tree[iNode];
            auto bnd = pkdn->bound();
            /*
            ** Before squeezing the bounds, calculate a minimum b value based on the splitting bounds alone.
            ** This gives us a better feel for the "size" of a bucket with only a single particle.
            */
            bmin = bnd.minside();
            pj = pkdn->lower();
            auto pkdl = pkd->tree[pkdn->lchild()];
            auto pkdu = pkd->tree[pkdn->rchild()];
            pkdCombineCells1(pkd,pkdn,pkdl,pkdu);
            kdn_r = pkdn->position();
            if (pkdn->count() <= NMAX_OPENCALC) {
                assert(pj<=pkdn->upper());
                double d2Max = 0;
                for (; pj<=pkdn->upper(); ++pj) {
                    auto p = pkd->particles[pj];
// #if defined(__AVX__) && defined(USE_SIMD)
//                     if (pkd->bIntegerPosition) {
//                         __m256d v = _mm256_sub_pd(pkdGetPos(pkd,p),_mm256_setr_pd(kdn_r[0],kdn_r[1],kdn_r[2],0.0));
//                         v = _mm256_mul_pd(v,v);
//                         __m128d t0 = _mm256_extractf128_pd(v,0);
//                         __m128d t2 = _mm256_extractf128_pd(v,1);
//                         __m128d t1 = _mm_unpackhi_pd(t0,t0);
//                         t0 = _mm_add_sd(t0,t2);
//                         t0 = _mm_add_sd(t0,t1);
//                         d2Max = _mm_cvtsd_f64(_mm_max_sd(t0,_mm_set_sd(d2Max)));
//                     }
//                     else
// #endif
                    {
                        r = p.position() - kdn_r;
                        d2 = blitz::dot(r,r);
                        d2Max = (d2 > d2Max)?d2:d2Max;
                    }
                }
                assert(d2Max>0);
                /*
                ** Now determine the opening radius for gravity.
                */
#ifdef USE_MAXSIDE
                b = bnd.maxside();
                if (b < bmin) b = bmin;
                if (d2Max>b) b = d2Max;
                pkdn->bMax = b;
#else
                pkdn->bMax() = sqrt(d2Max);
                if (pkdn->bMax() < bmin) pkdn->bMax() = bmin;
#endif
                assert(pkdn->bMax() >= 0);
            }
            else {
                pkdn->calc_open(bmin); // set bMax
            }
            pkdCombineCells2(pkd,pkdn,pkdl,pkdu);
        }
        ++iNode;
    }
}


void pkdCombineCells1(PKD pkd,treeStore::NodePointer pkdn,treeStore::NodePointer p1,treeStore::NodePointer p2) {
    auto p1bnd = p1->bound();
    auto p2bnd = p2->bound();
    auto bnd = p1bnd.combine(p2bnd);
    pkdn->set_bound(bnd);
    double m1,m2,ifMass;

    m1 = p1->mass();
    m2 = p2->mass();
    // In the case where a cell has all its particles source inactive mass == 0, which is ok, but we
    // still need a reasonable center in order to define opening balls in the tree code.
    if ( m1==0.0 || m2 == 0.0 ) {
        ifMass = 1.0;
        m1 = m2 = 0.5;
    }
    else ifMass = 1/(m1 + m2);
    blitz::TinyVector<double,3> p1_r, p2_r, kdn_r;
    p1_r = p1->position();
    p2_r = p2->position();
    kdn_r = ifMass*(m1*p1_r + m2*p2_r);
    pkdn->set_position(kdn_r);
    if (pkdn->have_velocity()) {
        pkdn->velocity() = ifMass*(m1*p1->velocity() + m2*p2->velocity());
    }
    if (pkdn->have_acceleration()) {
        pkdn->acceleration() = ifMass*(m1*p1->acceleration() + m2*p2->acceleration());
    }
    pkdn->fSoft2() = std::max(p1->fSoft2(),p2->fSoft2());
    pkdn->set_min_rung(std::min(p1->min_rung(),p2->min_rung()));
    pkdn->set_max_rung(std::max(p1->max_rung(),p2->max_rung()));

    if (pkdn->have_BOB()) {
        auto &SPHbob = pkdn->BOB();
        const auto &SPHbob1 = p1->BOB();
        const auto &SPHbob2 = p2->BOB();

        if (pkd->particles.present(PKD_FIELD::oNewSph)) {
            /* Combine ball or box of balls */
            SPHbob = SPHbob1.combine(SPHbob2);
        }
        else {
            SPHbob = SPHBOB(0.0);
        }
    }
    /* Combine marked flag */
    pkdn->set_marked(p1->is_marked() || p2->is_marked());
    pkdn->set_NN(p1->is_NN() || p2->is_NN());
}


void pkdCombineCells2(PKD pkd,treeStore::NodePointer pkdn,treeStore::NodePointer p1,treeStore::NodePointer p2) {
    /*
    ** Now calculate the reduced multipole moment.
    ** Shift the multipoles of each of the children
    ** to the CoM of this cell and add them up.
    */
    if (pkdn->have_moment()) {
        blitz::TinyVector<double,3> dr;
        auto &pkdn_moment = pkdn->moment();
        auto kdn_r = pkdn->position();
        auto p1_r  = p1->position();
        auto p2_r  = p2->position();
        pkdn_moment = p1->moment();
        dr = p1_r - kdn_r;
        momShiftFmomr(&pkdn_moment,p1->bMax(),dr[0],dr[1],dr[2]);

        momRescaleFmomr(&pkdn_moment,pkdn->bMax(),p1->bMax());

        auto mom = p2->moment();
        dr = p2_r - kdn_r;
        momShiftFmomr(&mom,p2->bMax(),dr[0],dr[1],dr[2]);
        momScaledAddFmomr(&pkdn_moment,pkdn->bMax(),&mom,p2->bMax());

    }
    else pkdn->mass() = p1->mass() + p2->mass();
    /*
    ** Combine the special fast gas ball bounds for SPH.
    */
    if (pkd->tree.present(KDN_FIELD::oNodeSphBounds)) {
        const auto &b1 = p1->sphbounds();
        const auto &b2 = p2->sphbounds();
        auto &bn = pkdn->sphbounds();
        bn.A.min = blitz::min(b1.A.min,b2.A.min);
        bn.A.max = blitz::max(b1.A.max,b2.A.max);
        bn.B.min = blitz::min(b1.B.min,b2.B.min);
        bn.B.max = blitz::max(b1.B.max,b2.B.max);
        bn.BI.min = blitz::min(b1.BI.min,b2.BI.min);
        bn.BI.max = blitz::max(b1.BI.max,b2.BI.max);
    }
}

void pkdTreeBuild(PKD pkd,int nBucket, int nGroup, uint32_t uRoot,uint32_t uTemp, double ddHonHLimit) {
#ifdef USE_ITT
    __itt_domain *domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle *shMyTask = __itt_string_handle_create("Tree Build");
    __itt_string_handle *shMySubtask = __itt_string_handle_create("My SubTask");
#endif
#ifdef USE_ITT
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif

    /*
    ** The KDN at "uRoot" (e.g., ROOT) is already setup (pLower and pUpper are correct)
    ** For more information look a pkdDumpTrees and the Initialize*() routines above.
    */

    if (uTemp==0) BuildTemp(pkd,uRoot,nBucket,nGroup,HUGE_VAL);
    else  BuildFromTemplate(pkd,uRoot,nBucket,nGroup,uTemp);
    Create(pkd,uRoot,ddHonHLimit);

    if (uRoot == FIXROOT) {
#ifndef SINGLE_CACHES
        mdlROcache(pkd->mdl,CID_CELL2,pkdTreeNodeGetElement,pkd,pkd->NodeSize(),pkd->Nodes());
        mdlROcache(pkd->mdl,CID_PARTICLE2,NULL,pkd->particles,pkd->particles.ParticleSize(),pkd->Local());
#endif
    }
    else {
        mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,pkd->NodeSize(),pkd->Nodes());
    }

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
}

/*
** The array iGrpOffset[i] passed in here must have 2*pkd->nGroups entries!
*/
void pkdGroupOrder(PKD pkd,uint32_t *iGrpOffset) {
    uint32_t i,iTree;
    uint32_t *iGrpEnd = &iGrpOffset[pkd->nGroups]; /* tricky because the 0th element is not used! */

    /* Count the number of particles in each group */
    for (i=0; i<=pkd->nGroups; ++i) iGrpOffset[i] = 0;
    for (i=0; i<pkd->Local(); ++i) {
        auto P = pkd->particles[i];
        auto gid = P.group();
        ++iGrpOffset[gid+1];
    }
    iGrpOffset[1] = 0;
    /* Calculate starting offsets for particles in a group */
    for (i=2; i<=pkd->nGroups; ++i) {
        iGrpOffset[i] += iGrpOffset[i-1];
        iGrpEnd[i-1] = iGrpOffset[i];
    }
    /* Reorder the particles into group order */
    for (iTree=1; iTree<pkd->nGroups; ++iTree) {
        for (i=iGrpOffset[iTree]; i<iGrpEnd[iTree];) {
            auto P = pkd->particles[i];
            auto gid = P.group();
            if (!gid) gid = pkd->nGroups;
            if (gid == iTree) ++i;
            else {
                auto p2 = pkd->particles[iGrpOffset[gid]++];
                swap(P,p2);
            }
        }
    }
}


void pkdTreeBuildByGroup(PKD pkd, int nBucket, int nGroup) {
    int i,k,n,iRoot;
    int iTree;

    /*
    ** Should use the pkdGroupOrder function to just reorder the particles
    ** without building the trees first. This is useful for some of the
    ** groupstats functions.
    */
    assert(0); /* pLite is gone -- this code path needs to be tested */
    if (pkd->Nodes() > 0) {
        /*
        ** Close cell caching space and free up nodes.
        */
        mdlFinishCache(pkd->mdl,CID_CELL);
    }

    /*
    ** It is only forseen that there are 4 reserved nodes at present 0-NULL, 1-ROOT, 2-UNUSED, 3-VAROOT.
    */
    pkd->tree.SetNodeCount(NRESERVED_NODES);

    if (pkd->hopSavedRoots == 0) {
        /* Sort particle by group, but with group 0 at the end */
        auto iGrpOffset = new int[pkd->nGroups+1];
        auto iGrpEnd = new int[pkd->nGroups+1];

        /* Count the number of particles in each group */
        for (i=0; i<=pkd->nGroups; ++i) iGrpOffset[i] = 0;
        for (i=0; i<pkd->Local(); ++i) {
            auto P = pkd->particles[i];
            auto gid = P.group();
            ++iGrpOffset[gid+1];
        }
        iGrpOffset[0] = iGrpOffset[1];
        iGrpOffset[1] = 0;
        /* Calculate starting offsets for particles in a group */
        for (i=2; i<=pkd->nGroups; ++i) {
            iGrpOffset[i] += iGrpOffset[i-1];
            iGrpEnd[i-1] = iGrpOffset[i];
        }

        /* Now construct the top tree node for each group */
        i = 0;
        for (auto gid=1; gid<pkd->nGroups; ++gid) {
            iRoot = pkd->TreeAllocRootNode();
            pkd->hopGroups[gid].iTreeRoot = iRoot;
            auto pNode = pkd->tree[iRoot];
            i = iGrpOffset[gid];
            pNode->set_local(i,i-1);
        }

        /* Reorder the particles into group order */
        for (iTree=1; iTree<pkd->nGroups; ++iTree) {
            for (i=iGrpOffset[iTree]; i<iGrpEnd[iTree]; ) {
                auto P = pkd->particles[i];
                auto gid = P.group();
                if (gid==0) gid = pkd->nGroups;
                if (gid == iTree) ++i;
                else {
                    auto p2 = pkd->particles[iGrpOffset[gid]++];
                    swap(P,p2);
                }
            }
        }

        delete [] iGrpOffset;
        delete [] iGrpEnd;

        /* Calculate the bounds for each group */
        for (i=0; i<pkd->Local();) {
            auto P = pkd->particles[i];
            auto r = P.position();
            auto gid = P.group();
            if (gid==0) break;
            iRoot = pkd->hopGroups[gid].iTreeRoot;
            auto pNode = pkd->tree[iRoot];

            // pNode->iLower = 0;
            pNode->set_group(true);
            // assert(pNode->pLower == i);

            auto dMin=r, dMax=r;
            for (auto p = pkd->particles[++i]; i<pkd->Local() && p.group()==gid; ++i) {
                r = p.position();
                dMin = blitz::min(dMin,r);
                dMax = blitz::max(dMax,r);
            }
            pNode->set_bound(Bound(dMin,dMax));
            // assert(pNode->pUpper == i-1);
        }
        /* We can use this to quickly rebuild the trees */
        pkd->hopSavedRoots = pkd->Nodes();
    }
    else {
        pkd->tree.SetNodeCount(pkd->hopSavedRoots);
        for (auto gid2=1; gid2<pkd->nGroups; ++gid2) {
            if (!pkd->hopGroups[gid2].bNeedGrav) continue;
            iRoot = pkd->hopGroups[gid2].iTreeRoot;
            auto pNode = pkd->tree[iRoot];
            // pNode->iLower = 0;
            pNode->set_group(true);
            n = pNode->upper();
            blitz::TinyVector<double,3> dMin, dMax;
            for (i=pNode->lower(); i<=n; ) {
                auto p = pkd->particles[i];
                auto r = p.position();
                auto gid = p.group();
                if (gid) {
                    assert(gid==gid2);
                    if (i==pNode->lower()) dMin = dMax = r;
                    else {
                        dMin = blitz::min(dMin,r);
                        dMax = blitz::max(dMax,r);
                    }
                    ++i;
                }
                else {
                    auto p2 = pkd->particles[n--];
                    swap(p,p2);
                }
            }
            assert(k==n+1);
            pNode->set_local(pNode->lower(),n);
            pNode->set_bound(Bound(dMin,dMax));
        }
    }


    for (auto gid=1; gid<pkd->nGroups; ++gid)
        if (pkd->hopGroups[gid].bNeedGrav)
            BuildTemp(pkd,pkd->hopGroups[gid].iTreeRoot,nBucket,nGroup,HUGE_VAL);
    for (auto gid=1; gid<pkd->nGroups; ++gid)
        if (pkd->hopGroups[gid].bNeedGrav)
            Create(pkd,pkd->hopGroups[gid].iTreeRoot,0.0);
    /*
    ** Finally activate a read only cache for remote access.
    */
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,
               pkd->NodeSize(),pkd->Nodes());


}

void pkdDistribRoot(PKD pkd,double *r,MOMC *pmom) {
    pkd->ew.r[0] = r[0];
    pkd->ew.r[1] = r[1];
    pkd->ew.r[2] = r[2];
    pkd->ew.mom = *pmom;
}

void pkdTreeUpdateFlagBoundsRecurse(PKD pkd,uint32_t uRoot,SPHOptions *SPHoptions) {
    int id = 0;
    int pj;

    if (!pkd->tree.present(KDN_FIELD::oNodeBOB)) return; // We aren't keeping track of the bounds
    float fBallFactor = (SPHoptions->dofBallFactor) ? pkd->SPHoptions.fBallFactor : 1.0f;

    auto c = pkd->tree[uRoot];
    auto &SPHbob = pkd->NodeBOB(c);

    if (c->is_bucket()) {
        pj = c->lower();
        auto p = pkd->particles[pj];
        auto ft = p.position();
        float limitedBallSize = std::min(pkd->SPHoptions.ballSizeLimit,fBallFactor * p.ball());
        SPHbob = SPHBOB(ft,limitedBallSize);
        if (p.marked()) c->set_marked(true);
        if (p.NN_flag()) c->set_NN(true);
        for (++pj; pj<=c->upper(); ++pj) {
            auto p = pkd->particles[pj];
            auto ft = p.position();
            float limitedBallSize = std::min(pkd->SPHoptions.ballSizeLimit,fBallFactor * p.ball());
            SPHbob = SPHbob.combine(ft,limitedBallSize);
            if (p.marked()) c->set_marked(true);
            if (p.NN_flag()) c->set_NN(true);
        }
    }
    else {
        auto [iCellLo,idLo,iCellUp,idUp] = c->get_child_cells(id);
        pkdTreeUpdateFlagBoundsRecurse(pkd,iCellLo,SPHoptions);
        pkdTreeUpdateFlagBoundsRecurse(pkd,iCellUp,SPHoptions);
        auto cLow = pkd->tree[iCellLo];
        auto cUp = pkd->tree[iCellUp];
        auto &SPHbobLow = cLow->BOB();
        auto &SPHbobUp = cUp->BOB();
        SPHbob = SPHbobLow.combine(SPHbobUp);
        c->set_marked(cLow->is_marked() || cUp->is_marked());
        c->set_NN(cLow->is_NN() || cUp->is_NN());
    }
}

void pkdTreeUpdateFlagBounds(PKD pkd,uint32_t uRoot,SPHOptions *SPHoptions) {
    if (mdlCacheStatus(pkd->mdl,CID_CELL)) mdlFinishCache(pkd->mdl,CID_CELL);
    pkdTreeUpdateFlagBoundsRecurse(pkd, uRoot,SPHoptions);
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd,pkd->NodeSize(),pkd->Nodes());
}
