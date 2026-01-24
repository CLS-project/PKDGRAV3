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

#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include <memory>
#include <sys/types.h>
#include <vector>
#include <boost/unordered/unordered_flat_set.hpp>
#include "pkd.h"
#include "smooth/smoothfcn.h"
#include "group/group.h"
#include "blitz/array.h"

#define NNLIST_INCREMENT 200 /* number of extra neighbor elements added to nnList */

struct smContext {
    PKD pkd;
    PARTICLE* pSentinel = nullptr;
    std::shared_ptr<SmoothFunction> function;
    int nSmooth = 0;
    bool bPeriodic = false;
    int bOwnCache = 0;
    blitz::TinyVector<double, 3> rLast{0.0, 0.0, 0.0}; /* For the snake */
    std::vector<PQ> vpq;
    /*
    ** Hash table to indicate whether a remote particle is already present in the
    ** priority queue.
    */
    boost::unordered_flat_set<void*> hash;
    std::vector<NN> nnList;
    /*
     ** Also need the stacks for the tree search
     */
    std::vector<uint32_t> S;
    struct stStack {
        int id;
        int iCell;
        double min;
    };
    std::vector<stStack> ST;

    int iSmoothType = 0;
    static auto makeHelper(PKD pkd, SMF* smf, int bSymmetric, int iSmoothType) -> std::shared_ptr<SmoothHelper>;

    void initialize();
    void finish();

    auto pqSearch(PQ* pq, blitz::TinyVector<double, 3> r, int iRoot) -> PQ*;
    void Gather(double fBall2, blitz::TinyVector<double, 3> r);

  protected:
    // These calls (Re)SmoothSingle for all particles which ultimately calls the Smooth() method of the SmoothFunction
    // and can be overridden to provide efficient particle filtering.
    virtual void SmoothParticles();
    virtual auto ReSmoothParticles() -> uint64_t;

  public:
    smContext(PKD pkdIn, int iSmoothType, SMF* smf, bool bSymmetric, int nSmooth, bool bPeriodic, bool bMakeCache = true);
    smContext(PKD pkdIn, std::shared_ptr<SmoothHelper> helper, int nSmooth, bool bPeriodic, bool bMakeCache = true);
    virtual ~smContext();
    auto SmoothSingle(particleStore::ParticleReference& p, int iRoot1 = ROOT, int iRoot2 = 0) -> float;
    auto ReSmoothSingle(particleStore::ParticleReference& p, double fBall) -> bool;
    void Smooth();
    auto ReSmooth() -> uint64_t { return ReSmoothParticles(); }
};
using SMX = smContext*;

void smSmooth(SMX, SMF*);
auto smReSmooth(SMX, int) -> uint64_t;
#ifdef OPTIM_SMOOTH_NODE
auto smReSmoothNode(SMX, SMF*, int) -> uint64_t;
void buildInteractionList(SMX smx, KDN* node, Bound bnd_node, blitz::TinyVector<double, 3> r, int ix, int iy, int iz);
#endif

#endif
