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

static const int cacheSize = 64 * 1024;

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

    Context *createLower();
    Context *split(int idLast);
    };

Context::Context(MDL mdl,int iLast)
    : mdl(mdl), pLower(0), idUpper(iLast) {
    }

Context *Context::createLower() {
    assert(pLower==0);
    pLower = new Context(mdl);
    return pLower;
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
    TEST_RO,
    TEST_FLUSH,
    };
} // namespace worker

namespace SetAdd {

struct inSetAdd {
    int idUpper;
    };

int serviceSetAdd(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto in = reinterpret_cast<struct inSetAdd*>(vin);
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

int test_ro(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto pnBAD = reinterpret_cast<std::uint64_t*>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdlReqService(ctx->getMDL(),ctx->getUpper(),worker::TEST_RO,NULL,0);
        test_ro(ctx->getLower(),vin,nIn,vout,nOut);
        mdlGetReply(ctx->getMDL(),rID,&nBAD,&nOut);
	*pnBAD += nBAD;
        }
    else {
    	int idSelf = mdlSelf(ctx->getMDL());
	int nData = cacheSize;
	auto pData = new std::uint64_t[nData];
	for(auto i=0; i<nData; ++i) pData[i] = ((1UL*idSelf)<<33) + 10 + i;
	mdlROcache(ctx->getMDL(),0,NULL,pData,sizeof(pData[0]),nData);

	nBAD = 0;
	for(auto iProc=0; iProc<mdlThreads(ctx->getMDL()); ++iProc) {
	    for(auto i=0; i<nData; ++i) {
		auto pRemote = reinterpret_cast<std::uint64_t*>(mdlAcquire(ctx->getMDL(),0,i,iProc));
		auto expect = ((1UL*iProc)<<33) + 10 + i;
		if (*pRemote != expect) ++nBAD;
		mdlRelease(ctx->getMDL(),0,pRemote);
		}
	    }
	mdlFinishCache(ctx->getMDL(),0);
	delete[] pData;
	*pnBAD = nBAD;
        }

    return sizeof(*pnBAD);
    }

static void initFlush(void *vctx, void *g) {
    //auto ctx = reinterpret_cast<worker::Context*>(vctx);
    auto p1 = reinterpret_cast<std::uint64_t*>(g);
    *p1 = 0;
    }

static void combFlush(void *vctx, void *g1, void *g2) {
    //auto ctx = reinterpret_cast<worker::Context*>(vctx);
    auto p1 = reinterpret_cast<std::uint64_t*>(g1);
    auto p2 = reinterpret_cast<std::uint64_t*>(g2);
    *p1 += *p2;
    }

int test_flush(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto pnBAD = reinterpret_cast<std::uint64_t*>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdlReqService(ctx->getMDL(),ctx->getUpper(),worker::TEST_FLUSH,NULL,0);
        test_flush(ctx->getLower(),vin,nIn,vout,nOut);
        mdlGetReply(ctx->getMDL(),rID,&nBAD,&nOut);
	*pnBAD += nBAD;
        }
    else {
    	int idSelf = mdlSelf(ctx->getMDL());
	int nData = cacheSize;
	auto pData = new std::uint64_t[nData];
	for(auto i=0; i<nData; ++i) pData[i] = ((1UL*idSelf)<<33) + 10 + i;
	mdlCOcache(ctx->getMDL(),0,NULL,pData,sizeof(pData[0]),nData,ctx,initFlush,combFlush);

	nBAD = 0;
	for(auto iProc=0; iProc<mdlThreads(ctx->getMDL()); ++iProc) {
	    for(auto i=0; i<nData; ++i) {
		auto pRemote = reinterpret_cast<std::uint64_t*>(mdlVirtualFetch(ctx->getMDL(),0,i,iProc));
		++*pRemote;
		}
	    }
	mdlFinishCache(ctx->getMDL(),0);
	auto nThreads = mdlThreads(ctx->getMDL());
	for(auto i=0; i<nData; ++i) {
	    if (pData[i] != nThreads + ((1UL*idSelf)<<33) + 10 + i) ++nBAD;
	    }
	delete[] pData;
	*pnBAD = nBAD;
        }

    return sizeof(*pnBAD);
    }

/*
** This function is called at the very start by every thread.
** It returns the "worker context"; in this case the ctx.
*/
void *worker_init(MDL mdl) {
    auto ctx = new worker::Context(mdl);
    mdlAddService(mdl,worker::SET_ADD,ctx,(fcnService_t*)SetAdd::serviceSetAdd,
                  sizeof(SetAdd::inSetAdd),0);
    mdlAddService(mdl,worker::TEST_RO,ctx,(fcnService_t*)test_ro,
                  0,sizeof(std::uint64_t));
    mdlAddService(mdl,worker::TEST_FLUSH,ctx,(fcnService_t*)test_flush,
                  0,sizeof(std::uint64_t));
    return ctx;
    }

/*
** This function is called at the very end for every thread.
** It needs to destroy the worker context (ctx).
*/
void worker_done(MDL mdl, void *vctx) {
    auto ctx = reinterpret_cast<worker::Context*>(vctx);

    delete ctx;
    }

namespace {

class CacheTest : public ::testing::Test {
    };


TEST_F(CacheTest, CacheReadWorks) {
    auto ctx = reinterpret_cast<worker::Context*>(mdlWORKER());
    std::uint64_t nBAD;
    test_ro(ctx,NULL,0,&nBAD,sizeof(nBAD));
    EXPECT_EQ(nBAD,0);
    }

TEST_F(CacheTest, CacheFlushWorks) {
    auto ctx = reinterpret_cast<worker::Context*>(mdlWORKER());
    std::uint64_t nBAD;
    test_flush(ctx,NULL,0,&nBAD,sizeof(nBAD));
    EXPECT_EQ(nBAD,0);
    }

}  // namespace


/*
** This is invoked for the "master" process after the worker has been setup.
*/
int master(MDL mdl,void *vctx) {
    auto ctx = reinterpret_cast<worker::Context*>(vctx);
    int argc = mdlGetArgc(mdl);
    char **argv = mdlGetArgv(mdl);

    mdlSetCacheSize(ctx->getMDL(),cacheSize); // Small cache to make it faster to test


    SetAdd::inSetAdd in;
    in.idUpper = mdlThreads(ctx->getMDL());
    SetAdd::serviceSetAdd(ctx,&in,sizeof(in),NULL,0);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
    }

int main(int argc,char **argv) {
    return mdlLaunch(argc,argv,master,worker_init,worker_done);
    }
