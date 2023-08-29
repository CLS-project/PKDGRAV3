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

#ifndef METAL_H
#define METAL_H
#ifndef xxUSE_METALxx

#include "mdl.h"
#include "mdlmetal.h"
#include "pppc.h"

class MetalClient {
protected:
    mdl::mdlClass &mdl;
    mdl::messageQueue<MessagePP> freePP;
    MessagePP *pp;
    mdl::messageQueue<MessagePC> freePC;
    MessagePC *pc;
protected:
    template<class MESSAGE>
    void flush(MESSAGE *&M) {
        if (M) {
            M->prepare();
            mdl.enqueue(*M);
            M = nullptr;
        }
    }
    template<class MESSAGE,class QUEUE,class TILE>
    int queue(MESSAGE *&M,QUEUE &Q, workParticle *work, TILE &tile, bool bGravStep) {
        if (M) { // If we are in the middle of building data for a kernel
            if (M->queue(work,tile,bGravStep)) return work->nP; // Successfully queued
            flush(M); // The buffer is full, so send it
        }
        mdl.gpu.flushCompleted();
        if (Q.empty()) return 0; // No buffers so the CPU has to do this part
        M = & Q.dequeue();
        if (M->queue(work,tile,bGravStep)) return work->nP; // Successfully queued
        return 0; // Not sure how this would happen, but okay.
    }
public:
    explicit MetalClient(mdl::mdlClass &mdl);
    void flushMETAL();
    int queuePP(workParticle *work, ilpTile &tile, bool bGravStep) {
        return queue(pp,freePP,work,tile,bGravStep);
    }
    int queuePC(workParticle *work, ilcTile &tile, bool bGravStep) {
        return queue(pc,freePC,work,tile,bGravStep);
    }
};
#endif
#endif
