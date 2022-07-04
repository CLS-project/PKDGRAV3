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
#ifndef METAL_PPPC_H
#define METAL_PPPC_H
#include "mdl.h"
#include "gpu/pppcdata.h"

template<class TILE,int WIDTH=32>
class MessagePPPC : public mdl::metal::metalMessage, public gpu::pppcData<TILE,WIDTH> {
protected:
    mdl::messageQueue<MessagePPPC> &freeQueue;
    bool bGravStep = false;
    virtual void launch(mdl::metal::Stream &stream,MTL::CommandBuffer *cbuf) override;
    virtual void finish() override;
public:
    explicit MessagePPPC(mdl::messageQueue<MessagePPPC> &freeQueue)
        : freeQueue(freeQueue) {
    }
    virtual ~MessagePPPC() {
    }

};

typedef MessagePPPC<ilpTile> MessagePP;
typedef MessagePPPC<ilcTile> MessagePC;
#endif
