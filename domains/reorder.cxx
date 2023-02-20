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
#include "reorder.h"

namespace NewDD {

int ServiceReorder::OffNode(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto in = static_cast<input *>(vin);
    auto thread = mdl->Self() + pst->nLower;
    auto proc = mdl->ThreadToProc(thread);
    pst->iOrdSplit = in->nPerProc * proc;
    auto rID = ReqService(pst,in,nIn);
    Traverse(pst->pstLower,in,nIn,NULL,0);
    mdl->GetReply(rID);
    return 0;
    // return Recurse(pst,vin,nIn,vout,nOut);
}

int ServiceReorder::AtNode(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto in = static_cast<input *>(vin);
    in->nPerCore = (in->nPerProc + mdl->Cores() - 1) / mdl->Cores();
    return Recurse(pst,vin,nIn,vout,nOut);
}

int ServiceReorder::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(pst->mdl);
    auto in = static_cast<input *>(vin);
    pst->iOrdSplit = in->nPerProc * mdl->Proc() + in->nPerCore * (pst->nLower + mdl->Core());
    auto rID = ReqService(pst,in,nIn);
    Traverse(pst->pstLower,in,nIn,NULL,0);
    mdl->GetReply(rID);
    return 0;
}

int ServiceReorder::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto in = static_cast<input *>(vin);
    auto pkd = pst->plcl->pkd;
    auto mdl = pst->mdl;

    pst->iOrdSplit = in->nPerProc * mdl->Proc() + in->nPerCore * mdl->Core();

    auto reorder = [](auto &particles, auto base) {
        auto get_bin = [base](auto &p) { return p.order() - base;};
        // Move particles to the correct position
        for (auto p=particles.begin(); p!=particles.end(); ++p) {
            auto offset = p - particles.begin();
            dd_offset_type bin;
            while ((bin = get_bin(*p)) != offset) {
                auto q = particles[bin];
                swap(*p,q);
            }
        }
        // assert(std::is_sorted(particles.begin(),particles.end(),[](auto &p1,auto &p2) {return p1.order() < p2.order();}));
    };

    auto shuffle = [](auto &counts,auto &particles, auto nPer, auto base) {
        using container = typename std::remove_reference<decltype(counts)>::type;
        auto get_bin = [base,nPer](auto &p) { return (p.order()-base) / nPer;};
        // Count the number of particles in each bin (process or core)
        for (auto &p : particles) {
            auto bin = get_bin(p);
            assert(bin<counts.size());
            ++counts[bin];
        }
        // Turn this into offsets for the swap
        container offsets;
        typename container::value_type offset = 0;
        for (auto c : counts) {
            offsets.emplace_back(offset);
            offset += c;
        }
        offsets.emplace_back(offset);
        // Now swap particles to the correct position (so they are ordered by bin)
        container targets(offsets);
        int iBin = 0;
        for (auto p=particles.begin(); p!=particles.end(); ++p) {
            auto offset = p - particles.begin();
            while (offset>=offsets[iBin+1]) ++iBin;
            int bin;
            while ((bin = get_bin(*p)) != iBin) {
                auto q = particles[targets[bin]++];
                swap(*p,q);
            }
        }
        // assert(std::is_sorted(particles.begin(),particles.end(),[get_bin](auto &p1,auto &p2) {return get_bin(p1) < get_bin(p2);}));
    };

    // Phase 1: exchange particles with the correct process
    std::vector<dd_offset_type> counts(mdl->Procs(),0);
    shuffle(counts,pkd->particles,in->nPerProc,0);
    pkd->SetLocal(mdl->swapglobal(pkd->particles,pkd->FreeStore(),pkd->particles.ElementSize(),counts.data()));

    // Phase 2: exchange particles with the correct core on each process
    counts.resize(mdl->Cores());
    std::fill(counts.begin(),counts.end(),0);
    shuffle(counts,pkd->particles,in->nPerCore,mdl->Proc() * in->nPerProc);
    pkd->SetLocal(mdl->swaplocal(pkd->particles,pkd->FreeStore(),pkd->particles.ElementSize(),counts.data()));

    // Phase 3: reorder particles locally
    reorder(pkd->particles,mdl->Proc() * in->nPerProc + mdl->Core() * in->nPerCore);

    return 0;
}

} // namespace NewDD