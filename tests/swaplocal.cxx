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

#include "gtest/gtest.h"
#include "mdl.h"

#include <assert.h>
#include <cstdint>
#include <vector>

namespace {

namespace worker {
class Context {
    MDL mdl;
    Context *pLower;
    int     idUpper;
    int nLeaves;
public:
    Context(MDL mdl, int idLast=-1);

    MDL getMDL() const {return mdl;}
    int getLeaves() const {return nLeaves; }
    int nLower() const {return idUpper - mdlSelf(mdl); }
    int nUpper() const {return nLeaves - nLower(); }
    Context *getLower() const {return pLower;}
    int      getUpper() const {return idUpper;}

    Context *split(int idLast);
};

Context::Context(MDL mdl,int iLast)
    : mdl(mdl), pLower(0), idUpper(iLast) {
}

// Split this node as optimally as possible
Context *Context::split(int idLast) {
    assert(pLower==0);
    auto idLower = mdlSelf(mdl);
    nLeaves = idLast - idLower;
    if (nLeaves>1) {
        auto idMiddle = (idLast + idLower) / 2;
        /* Make sure that the ctx lands on core zero */
        auto iProcLower = mdlThreadToProc(mdl,idLower);
        auto iProcUpper = mdlThreadToProc(mdl,idLast-1);
        if (iProcLower!=iProcUpper) {
            idMiddle = mdlProcToThread(mdl,mdlThreadToProc(mdl,idMiddle));
        }
        idUpper = idMiddle;
        pLower = new Context(mdl,idMiddle);
    }
    else idUpper = 0;
    return pLower;
}

enum services {
    STOP,
    SET_ADD,
    TEST_SWAPLOCAL,
};
} // namespace worker

namespace SetAdd {

struct inSetAdd {
    int idUpper;
};

int serviceSetAdd(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto in = reinterpret_cast<struct inSetAdd *>(vin);
    assert(nIn == sizeof(struct inSetAdd));
    auto ctxNew = ctx->split(in->idUpper);
    if (ctxNew) {
        auto rID = mdlReqService(ctx->getMDL(),ctx->getUpper(),worker::SET_ADD,in,nIn);
        in->idUpper = ctx->getUpper();
        serviceSetAdd(ctxNew,in,nIn,NULL,0);
        mdlGetReply(ctx->getMDL(),rID,NULL,NULL);
    }
    return 0;
}
} // namespace

namespace test {

namespace swaplocal {

constexpr int SERVICE = worker::TEST_SWAPLOCAL;

struct DATA {
    int iDomain;
    int iCore;
    DATA() = default;
    DATA(int iDomain, int iCore) : iDomain(iDomain), iCore(iCore) {}
    bool operator<(const DATA &rhs) const {
        return iDomain < rhs.iDomain || (iDomain==rhs.iDomain && iCore<rhs.iCore);
    }
};

int test(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(ctx->getMDL());
    auto pnBAD = reinterpret_cast<std::uint64_t *>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdlReqService(ctx->getMDL(),ctx->getUpper(),SERVICE,NULL,0);
        test(ctx->getLower(),vin,nIn,vout,nOut);
        mdlGetReply(ctx->getMDL(),rID,&nBAD,&nOut);
        *pnBAD += nBAD;
    }
    else {
        nBAD = 0;

        int iCore = mdl->Core();
        int nCore = mdl->Cores();

        int m = 10;
        int n = nCore * (nCore-1) / 2 * m;
        int N = n + n/100 + 128;
        std::vector<DATA> data;

        data.reserve(N);

        // printf("Running test on %d of %d with %d elements\n",iCore,nCore,n);

        auto counts = new mdl::mdlClass::dd_offset_type[nCore];

        int j=0;
        for (auto i=0; i<nCore; ++i) {
            auto nn = (nCore-1-iCore+i) % nCore * m;
            counts[i] = nn;
            for (; nn--; ++j) data.emplace_back(i,iCore);
        }
        assert(j==n);
        mdl->swaplocal(data.data(),N,sizeof(DATA),counts);
        std::sort(data.begin(),data.end());

        j=0;
        for (auto i=0; i<nCore; ++i) {
            auto nn = (2*nCore - i + iCore - 1) % nCore * m;
            for (; nn--; ++j) {
                if (data[j].iDomain != iCore) ++nBAD;
                if (data[j].iCore != i) ++nBAD;
            }
        }
        assert(j==n);
        *pnBAD = nBAD;
    }

    return sizeof(*pnBAD);
}

class SwapLocalTest : public ::testing::Test {};

TEST_F(SwapLocalTest, SwapLocalBasic) {
    auto ctx = reinterpret_cast<worker::Context *>(mdlWORKER());
    std::uint64_t nBAD;
    test::swaplocal::test(ctx,NULL,0,&nBAD,sizeof(nBAD));
    EXPECT_EQ(nBAD,0);
}


} // namespace swaplocal
} // namespace test

/*
** This function is called at the very start by every thread.
** It returns the "worker context"; in this case the ctx.
*/
void *worker_init(MDL mdl) {
    auto ctx = new worker::Context(mdl);
    mdlAddService(mdl,worker::SET_ADD,ctx,(fcnService_t *)SetAdd::serviceSetAdd, sizeof(SetAdd::inSetAdd),0);
    mdlAddService(mdl,test::swaplocal::SERVICE,ctx,(fcnService_t *)test::swaplocal::test, 0,sizeof(std::uint64_t));
    return ctx;
}

/*
** This function is called at the very end for every thread.
** It needs to destroy the worker context (ctx).
*/
void worker_done(MDL mdl, void *vctx) {
    auto ctx = reinterpret_cast<worker::Context *>(vctx);

    delete ctx;
}

}  // namespace


/*
** This is invoked for the "master" process after the worker has been setup.
*/
int master(MDL mdl,void *vctx) {
    auto ctx = reinterpret_cast<worker::Context *>(vctx);
    int argc = mdlGetArgc(mdl);
    char **argv = mdlGetArgv(mdl);

    SetAdd::inSetAdd in;
    in.idUpper = mdlThreads(ctx->getMDL());
    SetAdd::serviceSetAdd(ctx,&in,sizeof(in),NULL,0);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

int main(int argc,char **argv) {
    return mdlLaunch(argc,argv,master,worker_init,worker_done);
}
