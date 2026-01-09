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
#include "pkd.h"
#include "smooth/smoothfcn.h"
#include "hydro/hydro.h"
#include "group/group.h"
#include "blitz/array.h"

#define NNLIST_INCREMENT 200 /* number of extra neighbor elements added to nnList */

struct hashElement {
    void* p;
    struct hashElement* coll;
};

struct smContext {
    smContext(PKD pkdIn, std::shared_ptr<SmoothHelper> helper, int nSmooth, bool bPeriodic, bool bMakeCache = true);
    virtual ~smContext();
    PKD pkd;
    PARTICLE* pSentinel = nullptr;
    std::shared_ptr<SmoothFunction> function;
    int nSmooth = 0;
    int nQueue = 0;
    bool bPeriodic = false;
    int bOwnCache = 0;
    int iSmoothType = 0;
    blitz::TinyVector<double, 3> rLast{0.0, 0.0, 0.0}; /* For the snake */
    std::vector<PQ> pq;
    /*
    ** Flags to mark local particles which are finished in the processing
    ** of the smFastGas routine. They have updated their densities.
    */
    char* bDone;
    /*
    ** Hash table to indicate whether a remote particle is already present in the
    ** priority queue.
    */
    int nHash = 0; /* should be a prime number > nSmooth */
    struct hashElement* pHash = nullptr;
    struct hashElement* pFreeHash = nullptr;
    std::vector<NN> nnList;
    int* S = nullptr;
    double* Smin = nullptr;
    /*
     ** Also need the stacks for the tree search
     */
    struct stStack {
        int id;
        int iCell;
        double min;
    }* ST;
    /*
    ** Context for nearest neighbor lists.
    */
    // LCODE lcmp;
    /*
    ** Some variables needed for smNewFof().
    */
    uint32_t iHead = 0;
    uint32_t iTail = 0;
    int* Fifo = nullptr;
};
using SMX = smContext*;

auto smInitialize(SMX* psmx, PKD pkd, SMF* smf, int nSmooth, int bPeriodic, int bSymmetric, int iSmoothType, bool bMakeCache = true) -> int;
void smFinish(SMX, SMF*);
void smSmoothInitialize(SMX smx);
void smSmoothFinish(SMX smx);
float smSmoothSingle(SMX smx, SMF* smf, particleStore::ParticleReference& p, int iRoot1, int iRoot2);
void smSmooth(SMX, SMF*);
auto smReSmoothSingle(SMX smx, SMF* smf, particleStore::ParticleReference& p, double fBall) -> bool;
uint64_t smReSmooth(SMX, SMF*, int);
#ifdef OPTIM_SMOOTH_NODE
uint64_t smReSmoothNode(SMX, SMF*, int);
void buildInteractionList(SMX smx, KDN* node, Bound bnd_node, blitz::TinyVector<double, 3> r, int ix, int iy, int iz);
#endif

void smGather(SMX smx, double fBall2, blitz::TinyVector<double, 3> r);

#endif
