#ifdef  BLACKHOLES
#include <algorithm>
#include "blackhole/merger.h"
#include "smooth/smooth.h"
#include "master.h"

using blitz::TinyVector;
using blitz::any;
using blitz::all;
using blitz::abs;
using blitz::max;
using blitz::floor;
using blitz::dot;

void MSR::BHMerger(double dTime) {
    struct inSmooth in;
    struct outSmooth out;
    struct outGetNParts Nout;

    TimerStart(TIMER_BHS);

#ifdef DEBUG_BH_ONLY
    // This was not performed in the main loop because there is no gas
    ReorderWithinNodes();
#endif

    in.nSmooth = param.nSmooth;
    in.bPeriodic = param.bPeriodic;
    in.bSymmetric = 0;
    in.iSmoothType = SMX_BH_MERGER;
    SmoothSetSMF(&(in.smf), dTime, 0.0, in.nSmooth);

    Nout.n = 0;
    Nout.nDark = 0;
    Nout.nGas = 0;
    Nout.nStar = 0;
    Nout.nBH = 0;

    if (parameters.get_bVStep()) {
        double sec,dsec;
        sec = Time();
        pstReSmoothNode(pst,&in,sizeof(in), &out,sizeof(struct outSmooth));
        pstMoveDeletedParticles(pst, NULL, 0, &Nout, sizeof(struct outGetNParts));
        dsec = Time() - sec;
    }
    else {
        pstReSmoothNode(pst,&in,sizeof(in),&out,sizeof(struct outSmooth));
        pstMoveDeletedParticles(pst, NULL, 0, &Nout, sizeof(struct outGetNParts));
    }

    TimerStop(TIMER_BHS);

    N = Nout.n;
    nDark = Nout.nDark;
    nGas = Nout.nGas;
    nStar = Nout.nStar;
    nBH = Nout.nBH;
}

void pkdRepositionBH(PKD pkd) {
#ifndef DEBUG_BH_ONLY
    for (auto &p : pkd->particles) {
        if (p.is_bh()) {
            auto &bh = p.BH();
            if (bh.doReposition) {
                p.set_position(bh.newPos);
                bh.doReposition = 0;
            }
        }
    }
#endif
}


/* IA: This should be called when we have a merger across different domains...
 * think about it!
 * TODO
 **/
void combBHmerger(void *vpkd, void *p1,const void *p2) {

    assert(0);


}



void smBHmerger(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];

    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
        /* IA: Merging across domain boundaries is a tricky thing,
         * This would require careful atomic operations in order to maintain
         * conservation *and* a deterministic behaviour independent of
         * the number of cores.
         *
         * The easiest solution is to relax the deterministic and invariant
         * behaviour and just merge BH that are in the same domain.
         * This should be conservative, and in standard cosmological simulations
         * this situation will be very rare, and eventually both BH will lie
         * in the same computational domain, so I think this won't be a
         * critical problem.
         *
         * Anyway, the most important thing is *conservation*, so that is what
         * I should check for in the test cases
         */
        if (pkd->Self() != nnList[i].iPid) continue;


        if (&p!=&q && p.is_bh() && q.is_bh() && p.marked() && q.marked()) {
            // At this point, we are sure that both particles are black holes and
            //  the distance between them is, at most, p.soft()
            //
            // So, the following should not give any problem
            // (remove for performance)
            assert( nnList[i].fDist2 < 4.*p.soft()*p.soft() );
            assert( p.have_mass() );


            if ( p.position(0) >= q.position(0) ) {
                const auto &pmass = p.mass();
                auto &pv = p.velocity();
                const auto &qv = q.velocity();
                const auto dv2 = dot(pv - qv, pv - qv);

                const auto &dist2 = nnList[i].fDist2;
                const float inv_dist = dist2 > 0.0 ? sqrt(1./dist2) : FLT_MAX;
                if (dv2 < pmass*inv_dist) {
                    // We have a merger!!
                    const auto &qmass = q.mass();
                    printf("Merger!!! %e %e \n", pmass, qmass);

                    const float newmass = pmass + qmass;
                    const float inv_newmass = 1./newmass;

                    /* We can not update the position while in this loop,
                     * because we are still using information of the tree,
                     * which may be compromised.
                     *
                     * Instead, we save the new position and then update the
                     * merged BHs in a separated loop.
                     * In the case that they are following the minimum potential of
                     * the nearby particles, this update may not be needed.
                     */
                    auto &pbh = p.BH();
                    pbh.newPos = p.position() - qmass*inv_newmass*nnList[i].dr;
                    pbh.doReposition = 2;
                    pbh.dInternalMass += q.BH().dInternalMass;
                    pbh.dAccEnergy += q.BH().dAccEnergy;
                    pv = (pmass*pv + qmass*qv) * inv_newmass;
                    p.set_mass(newmass);
                    q.set_mass(0.);

                    // We do this to check that the deleted particles are not
                    // entering the gravity computation, as this will activate
                    // an assert in the gravity loop

                    pkdDeleteParticle(pkd, q);

                    // We only allow for one merger per particle per timestep
                    p.set_marked(false);
                    return;
                }
            }
        }
    }
}




int smReSmoothBHNode(SMX smx,SMF *smf, int iSmoothType) {
    PKD pkd = smx->pkd;
    int nSmoothed=0;

    smx->nnListSize = 0;
    int nnListMax_p = NNLIST_INCREMENT;

    NN *nnList_p;
    nnList_p = static_cast<NN *>(malloc(sizeof(NN)*nnListMax_p));

    // Here we store the pointers to the particles whose interaction need
    // to be computed
    std::vector<PARTICLE *> sinks;
    sinks.reserve(64); // At most, the size of the bucket



    for (int i = NRESERVED_NODES; i < pkd->Nodes() - 1; ++i) {
        auto node = pkd->tree[i];
        bool bBHinBucket;
#ifdef OPTIM_REORDER_IN_NODES
        bBHinBucket = node->is_bucket() && node->Nbh()>0;
#else
        bBHinBucket = node->is_bucket() && std::any_of(node->begin(), node->end(), [](auto &p) {return p.is_bh();});
#endif
        if (bBHinBucket) {
            // We are in a bucket containing a BH

            // Prepare the interaction list
            auto bnd_node = node->bound();
            TinyVector<double,3> r{bnd_node.center()};

            // First, we add all the particles whose interactions
            // need to be computed
            sinks.clear();
            auto pStart = node->begin();
#ifdef OPTIM_REORDER_IN_NODES
            pStart += node->Ngas();
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
            pStart += node->Nstar();
#endif
            auto pEnd = pStart + node->Nbh();
#else
            auto pEnd = node->end();
#endif

            TinyVector<double,3> fMax_shrink{0.0};
            for (auto pj = pStart; pj < pEnd; ++pj) {
#ifndef OPTIM_REORDER_IN_NODES
                if (!pj->is_bh()) continue;
#endif

#ifdef OPTIM_AVOID_IS_ACTIVE
                if (pj->marked()) {
#else
                if (pj->is_active()) {
#endif
                    TinyVector<double,3> disp {abs(pj->position()-bnd_node.center()) + pj->soft()*2.};
                    fMax_shrink = max(fMax_shrink,disp);

                    sinks.push_back(&*pj);
                }
            }
            // There is no elligible particle in this bucket, go to the next
            if (sinks.empty()) continue;

            bnd_node.shrink(fMax_shrink);

            int nCnt = 0;
            if (smx->bPeriodic) {
                TinyVector<int,3> iStart, iEnd;
                iStart = floor((r - bnd_node.apothem()) / pkd->fPeriod + 0.5);
                iEnd   = floor((r + bnd_node.apothem()) / pkd->fPeriod + 0.5);
                for (int ix=iStart[0]; ix<=iEnd[0]; ++ix) {
                    r[0] = bnd_node.center(0) - ix*pkd->fPeriod[0];
                    for (int iy=iStart[1]; iy<=iEnd[1]; ++iy) {
                        r[1] = bnd_node.center(1) - iy*pkd->fPeriod[1];
                        for (int iz=iStart[2]; iz<=iEnd[2]; ++iz) {
                            r[2] = bnd_node.center(2) - iz*pkd->fPeriod[2];
                            buildCandidateMergerList(smx, smf, node, bnd_node,
                                                     &nCnt, r, ix, iy, iz);
                        }
                    }
                }
            }
            else {
                buildCandidateMergerList(smx, smf, node, bnd_node,
                                         &nCnt, r, 0, 0, 0);
            }



            // IA: Now we should have inside nnList all the particles in the bucket (sinks) and those of which can
            //  interact with them from other buckets (smx->nnList)
            //


            for (auto &P : sinks) {
                auto partj = pkd->particles[P];
                // We need to double check this, as the particle may have
                // been marked for deletion just before this
                if (partj.is_marked()) {
                    const float fBall2_p = 4.*partj.soft()*partj.soft();
                    TinyVector<double,3> dr_node{bnd_node.center() - partj.position()};
                    int nCnt_p = 0;
                    for (auto pk = 0; pk < nCnt; ++pk) {
                        if (P == smx->nnList[pk].pPart) continue;
                        TinyVector<double,3> dr{smx->nnList[pk].dr - dr_node};
                        const auto fDist2 = dot(dr,dr);
                        if (fDist2 <= fBall2_p) {
                            if (nCnt_p >= nnListMax_p) {
                                nnListMax_p += NNLIST_INCREMENT;
                                nnList_p = static_cast<NN *>(realloc(nnList_p,nnListMax_p*sizeof(NN)));
                                assert(nnList_p != NULL);
                            }

                            nnList_p[nCnt_p].fDist2 = fDist2;
                            nnList_p[nCnt_p].dr = dr;
                            nnList_p[nCnt_p].pPart = smx->nnList[pk].pPart;
                            nnList_p[nCnt_p].iIndex = smx->nnList[pk].iIndex;
                            nnList_p[nCnt_p].iPid = smx->nnList[pk].iPid;

                            ++nCnt_p;
                        }

                    }

                    smx->fcnSmooth(&partj,partj.soft(),nCnt_p,nnList_p,smf);

                }
            }

            nSmoothed += sinks.size();

            for (auto pk = 0; pk < nCnt; ++pk) {
                if (smx->nnList[pk].iPid != pkd->Self()) {
                    mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[pk].pPart);
                }
            }
        }
    }
    free(nnList_p);
    return nSmoothed;
}

inline static auto getCell(PKD pkd, int iCell, int id) {
    return (id==pkd->Self()) ? pkd->tree[iCell]
           : pkd->tree[static_cast<KDN *>(mdlFetch(pkd->mdl,CID_CELL,iCell,id))];
}

/* This is mostly a version of buildInteractionList, but that only looks
 * for BH particles.
 *
 * Probably, using C++ templates this could be done avoiding code duplications
 * or extra overheads...
 */
void buildCandidateMergerList(SMX smx, SMF *smf, KDN *node, Bound bnd_node, int *nCnt_tot,
                              TinyVector<double,3> r, int ix, int iy, int iz) {
    PKD pkd = smx->pkd;
    int id, sp, iCell;
    smContext::stStack *S = smx->ST;
    int nCnt = *nCnt_tot;

    // We look for the biggest node that encloses the needed domain
    // We can only take advantage of this if we are are in the original cell
    auto kdn = getCell(pkd, iCell=pkd->iTopTree[ROOT], id=pkd->Self());

    //  Now we start the walk as usual
    sp = 0;
    while (1) {
        auto bnd = kdn->bound();
        //min2 = bnd.mindist(r);
        //if (min2 > fBall2) {
        //    goto NoIntersect;
        //}
        if (any(abs(bnd.center()-r) - bnd.apothem() - bnd_node.apothem() > 0)) goto NoIntersect;
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
                auto pStart = kdn->lower();
#ifdef OPTIM_REORDER_IN_NODES
                pStart += kdn->Ngas();
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
                pStart += kdn->Nstar();
#endif
                auto pEnd = pStart + kdn->Nbh();
#else
                auto pEnd = kdn->upper() + 1;
#endif
                for (auto pj = pStart; pj < pEnd; ++pj) {
                    auto p = pkd->particles[pj];
#ifndef OPTIM_REORDER_IN_NODES
                    if (!p.is_bh()) continue;
#endif
                    TinyVector<double,3> dr {r - p.position()};
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
                        smx->nnList[nCnt].iIndex = pj;
                        smx->nnList[nCnt].iPid = pkd->Self();
                        ++nCnt;
                    }
                }
            }
            else {
                auto pStart = kdn->lower();
#ifdef OPTIM_REORDER_IN_NODES
                pStart += kdn->Ngas();
#if (defined(STAR_FORMATION) && defined(FEEDBACK)) || defined(STELLAR_EVOLUTION)
                pStart += kdn->Nstar();
#endif
                auto pEnd = pStart + kdn->Nbh();
#else
                auto pEnd = kdn->upper() + 1;
#endif
                for (auto pj = pStart; pj < pEnd; ++pj) {
                    auto p = pkd->particles[static_cast<PARTICLE *>(mdlFetch(pkd->mdl,CID_PARTICLE,pj,id))];
#ifndef OPTIM_REORDER_IN_NODES
                    if (!p.is_bh()) continue;
#endif
                    TinyVector<double,3> dr {r - p.position()};
                    if (all(abs(dr) <= bnd_node.apothem())) {
                        if (nCnt >= smx->nnListMax) {
                            smx->nnListMax += NNLIST_INCREMENT;
                            smx->nnList = static_cast<NN *>(realloc(smx->nnList,smx->nnListMax*sizeof(NN)));
                            //printf("realloc \n");
                            assert(smx->nnList != NULL);
                        }

                        smx->nnList[nCnt].fDist2 = dot(dr,dr);
                        smx->nnList[nCnt].dr = dr;
                        smx->nnList[nCnt].pPart = static_cast<PARTICLE *>(mdlAcquire(pkd->mdl, CID_PARTICLE, pj, id));
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

#endif
