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
#include <array>

mdl::mdlClass *master_mdl;

namespace {

constexpr int cacheSize = 64 * 1024;

namespace worker {
class Context {
    mdl::mdlClass *mdl;
    Context *pLower;
    int     idUpper;
    int nLeaves;
public:
    Context(mdl::mdlClass * mdl, int idLast=-1);

    auto getMDL() const {return mdl;}
    int getLeaves() const {return nLeaves; }
    int nLower() const {return idUpper - mdl->Self(); }
    int nUpper() const {return nLeaves - nLower(); }
    Context *getLower() const {return pLower;}
    int      getUpper() const {return idUpper;}

    Context *split(int idLast);
};

Context::Context(mdl::mdlClass * mdl,int iLast)
    : mdl(mdl), pLower(0), idUpper(iLast) {
}

// Split this node as optimally as possible
Context *Context::split(int idLast) {
    assert(pLower==0);
    auto idLower = mdl->Self();
    nLeaves = idLast - idLower;
    if (nLeaves>1) {
        auto idMiddle = (idLast + idLower) / 2;
        /* Make sure that the ctx lands on core zero */
        auto iProcLower = mdl->ThreadToProc(idLower);
        auto iProcUpper = mdl->ThreadToProc(idLast-1);
        if (iProcLower!=iProcUpper) {
            idMiddle = mdl->ProcToThread(mdl->ThreadToProc(idMiddle));
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
    TEST_HASH,
    TEST_RO,
    TEST_FLUSH,
    TEST_FLUSH_AFTER_READ,
    TEST_ADVANCED_RO,
    TEST_ADVANCED_FLUSH,
};
} // namespace worker

namespace SetAdd {

struct inSetAdd {
    int idUpper;
};

int serviceSetAdd(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = ctx->getMDL();
    auto in = reinterpret_cast<struct inSetAdd *>(vin);
    assert(nIn == sizeof(struct inSetAdd));
    auto ctxNew = ctx->split(in->idUpper);
    if (ctxNew) {
        auto rID = mdl->ReqService(ctx->getUpper(),worker::SET_ADD,in,nIn);
        in->idUpper = ctx->getUpper();
        serviceSetAdd(ctxNew,in,nIn,NULL,0);
        mdl->GetReply(rID);
    }
    return 0;
}
} // namespace

namespace test {

namespace hash {

constexpr int SERVICE = worker::TEST_HASH;
typedef std::array<uint64_t,4> KEY;
typedef mdl::hash::HASH<KEY> HASH;

int test(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = ctx->getMDL();
    auto pnBAD = reinterpret_cast<std::uint64_t *>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdl->ReqService(ctx->getUpper(),SERVICE,NULL,0);
        test(ctx->getLower(),vin,nIn,vout,nOut);
        nOut = mdl->GetReply(rID,nBAD);
        *pnBAD += nBAD;
    }
    else {
        uint32_t idSelf = mdl->Self();
        uint32_t nData = cacheSize;

        nBAD = 0;

        // Construct a hash table filled with data
        HASH hash(nData);
        std::vector<uint64_t> values;
        values.reserve(nData);
        for (auto i=0; i<nData; ++i) {
            values.push_back(0xdeadbeef ^ idSelf ^ i);
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(values.back()),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            hash.insert(uHash,key,&values.back());
        }

        // Check the even keys
        for (auto i=0; i<nData; i+=2) {
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(values[i]),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = hash.lookup(uHash,key);
            if (data == nullptr || data!=&values[i]) ++nBAD;
        }

        // Check the odd keys (and remove)
        for (auto i=1; i<nData; i+=2) {
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(values[i]),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = hash.lookup(uHash,key);
            if (data == nullptr || data!=&values[i]) ++nBAD;
            hash.remove(uHash,key);
        }

        // Check the even keys again
        for (auto i=0; i<nData; i+=2) {
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(values[i]),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = hash.lookup(uHash,key);
            if (data == nullptr || data!=&values[i]) ++nBAD;
        }

        // Check that the odd keys are gone
        for (auto i=1; i<nData; i+=2) {
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(values[i]),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = hash.lookup(uHash,key);
            if (data != nullptr) ++nBAD;
        }

        *pnBAD = nBAD;
    }

    return sizeof(*pnBAD);
}

class HashTest : public ::testing::Test {};

TEST_F(HashTest, HashTestKeys) {
    std::uint64_t nBAD;
    master_mdl->RunService(test::hash::SERVICE,&nBAD);
    EXPECT_EQ(nBAD,0);
}

} // namespace hash

class CacheTest : public ::testing::Test {};

namespace ro {
constexpr int SERVICE = worker::TEST_RO;
int test(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(ctx->getMDL());
    auto pnBAD = reinterpret_cast<std::uint64_t *>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdl->ReqService(ctx->getUpper(),SERVICE,NULL,0);
        test(ctx->getLower(),vin,nIn,vout,nOut);
        nOut = mdl->GetReply(rID,nBAD);
        *pnBAD += nBAD;
    }
    else {
        int idSelf = mdl->Self();
        int nData = cacheSize;
        auto pData = new std::uint64_t[nData];
        for (auto i=0; i<nData; ++i) pData[i] = ((1UL*idSelf)<<33) + 10 + i;
        mdl->CacheInitialize(0,std::make_shared<mdl::arrayCACHEhelper>(sizeof(pData[0]),pData,nData));

        nBAD = 0;
        for (auto iProc=0; iProc<mdl->Threads(); ++iProc) {
            for (auto i=0; i<nData; ++i) {
                auto pRemote = mdl->Acquire<std::uint64_t>(0,i,iProc);
                auto expect = ((1UL*iProc)<<33) + 10 + i;
                if (*pRemote != expect) ++nBAD;
                mdl->Release(0,pRemote);
            }
        }
        mdl->FinishCache(0);
        delete[] pData;
        *pnBAD = nBAD;
    }

    return sizeof(*pnBAD);
}
TEST_F(CacheTest, CacheReadWorks) {
    std::uint64_t nBAD;
    master_mdl->RunService(test::ro::SERVICE,&nBAD);
    EXPECT_EQ(nBAD,0);
}
} // namespace test_ro

namespace flush {
class FlushHelper : public mdl::arrayCACHEhelper {
    virtual void    init(void *dst) override;
    virtual void combine(void *dst, const void *src, const void *key) override;
public:
    explicit FlushHelper(std::uint64_t *pArray,uint32_t nArray) : arrayCACHEhelper(pArray,nArray,true) {}
};

void FlushHelper::init(void *dst) {
    auto p1 = reinterpret_cast<std::uint64_t *>(dst);
    *p1 = 0;
}

void FlushHelper::combine(void *dst, const void *src, const void *key) {
    auto p1 = reinterpret_cast<std::uint64_t *>(dst);
    auto p2 = reinterpret_cast<const std::uint64_t *>(src);
    *p1 += *p2;
}

constexpr int SERVICE = worker::TEST_FLUSH;

int test(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(ctx->getMDL());
    auto pnBAD = reinterpret_cast<std::uint64_t *>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdl->ReqService(ctx->getUpper(),SERVICE,NULL,0);
        test(ctx->getLower(),vin,nIn,vout,nOut);
        nOut = mdl->GetReply(rID,nBAD);
        *pnBAD += nBAD;
    }
    else {
        int idSelf = mdl->Self();
        int nData = cacheSize;
        auto pData = new std::uint64_t[nData];
        for (auto i=0; i<nData; ++i) pData[i] = ((1UL*idSelf)<<33) + 10 + i;
        mdl->CacheInitialize(0,std::make_shared<FlushHelper>(pData,nData));

        nBAD = 0;
        for (auto iProc=0; iProc<mdl->Threads(); ++iProc) {
            for (auto i=0; i<nData; ++i) {
                auto pRemote = mdl->VirtualFetch<std::uint64_t>(0,i,iProc);
                ++*pRemote;
            }
        }
        mdl->FinishCache(0);
        auto nThreads = mdl->Threads();
        for (auto i=0; i<nData; ++i) {
            if (pData[i] != nThreads + ((1UL*idSelf)<<33) + 10 + i) ++nBAD;
        }
        delete[] pData;
        *pnBAD = nBAD;
    }

    return sizeof(*pnBAD);
}
TEST_F(CacheTest, CacheFlushWorks) {
    std::uint64_t nBAD;
    master_mdl->RunService(test::flush::SERVICE,&nBAD);
    EXPECT_EQ(nBAD,0);
}

namespace after_read {
constexpr int SERVICE = worker::TEST_FLUSH_AFTER_READ;

int test(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(ctx->getMDL());
    auto pnBAD = reinterpret_cast<std::uint64_t *>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdl->ReqService(ctx->getUpper(),SERVICE,NULL,0);
        test(ctx->getLower(),vin,nIn,vout,nOut);
        nOut = mdl->GetReply(rID,nBAD);
        *pnBAD += nBAD;
    }
    else {
        int idSelf = mdl->Self();
        int nData = cacheSize / 10; // Smaller is fine for this
        auto pData = new std::uint64_t[nData];
        for (auto i=0; i<nData; ++i) pData[i] = ((1UL*idSelf)<<33) + 10 + i;
        mdl->CacheInitialize(0,std::make_shared<FlushHelper>(pData,nData));

        for (auto iProc=0; iProc<mdl->Threads(); ++iProc) {
            for (auto i=0; i<nData; ++i) (void)mdl->Fetch(0,i,iProc);
        }

        nBAD = 0;
        for (auto iProc=0; iProc<mdl->Threads(); ++iProc) {
            for (auto i=0; i<nData; ++i) {
                auto pRemote = mdl->Acquire<std::uint64_t>(0,i,iProc);
                ++*pRemote;
                mdl->Release(0,pRemote);
            }
        }
        mdl->FinishCache(0);
        auto nThreads = mdl->Threads();
        for (auto i=0; i<nData; ++i) {
            if (pData[i] != nThreads + ((1UL*idSelf)<<33) + 10 + i) ++nBAD;
        }
        delete[] pData;
        *pnBAD = nBAD;
    }

    return sizeof(*pnBAD);
}
TEST_F(CacheTest, CacheFlushAfterReadWorks) {
    std::uint64_t nBAD;
    master_mdl->RunService(test::flush::after_read::SERVICE,&nBAD);
    EXPECT_EQ(nBAD,0);
}

} // namespace after_read
} // namespace flush

namespace advanced {
class AdvancedCache : public ::testing::Test {};

namespace ro {

class test_helper : public mdl::CACHEhelper {
protected:
    virtual uint32_t getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey) override {
        return * static_cast<const uint64_t *>(pKey); // TESTING: we put the remote id here
    }
public:
    explicit test_helper(uint32_t nData, bool bModify=false) : CACHEhelper(sizeof(uint64_t),nData,bModify) {}
};

constexpr int SERVICE = worker::TEST_ADVANCED_RO;
typedef std::array<uint64_t,4> KEY;
typedef mdl::hash::HASH<KEY> HASH;

int test(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(ctx->getMDL());
    auto pnBAD = reinterpret_cast<std::uint64_t *>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdl->ReqService(ctx->getUpper(),SERVICE,NULL,0);
        test(ctx->getLower(),vin,nIn,vout,nOut);
        nOut = mdl->GetReply(rID,nBAD);
        *pnBAD += nBAD;
    }
    else {
        uint32_t idSelf = mdl->Self();
        uint32_t nThreads = mdl->Threads();
        uint32_t nData = cacheSize;

        nBAD = 0;

        // Construct a hash table filled with data
        HASH hash(nData);
        std::vector<uint64_t> values;
        values.reserve(nData);
        for (auto i=0; i<nData; ++i) {
            values.push_back(0xdeadbeef ^ idSelf ^ i);
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(values.back()),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            hash.insert(uHash,key,&values.back());
        }

        mdl->AdvancedCacheInitialize(0,&hash,sizeof(values.front()),std::make_shared<test_helper>(sizeof(values.front()),false));
        auto idNext = (idSelf+1) % nThreads;
        // These should be cache hits
        for (auto i=idSelf; i<nData; i += nThreads) {
            auto expected = (0xdeadbeef ^ idNext ^ i);
            KEY key = {static_cast<uint64_t>(idNext),static_cast<uint64_t>(i),static_cast<uint64_t>(expected),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = mdl->KeyFetch<uint64_t>(0,uHash,&key,0,0,0);
            if (data == nullptr || *data != expected) ++nBAD;
        }

        // These should be cache misses
        for (auto i=idSelf; i<nData; i += nThreads) {
            auto expected = (0xfacefeed ^ idNext ^ i);
            KEY key = {static_cast<uint64_t>(idNext),static_cast<uint64_t>(i),static_cast<uint64_t>(expected),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = mdl->KeyFetch<uint64_t>(0,uHash,&key,0,0,0);
            if (data != nullptr) ++nBAD;
        }

        mdl->FinishCache(0);

        *pnBAD = nBAD;
    }

    return sizeof(*pnBAD);
}

TEST_F(AdvancedCache, ReadOnly) {
    std::uint64_t nBAD;
    master_mdl->RunService(SERVICE,&nBAD);
    EXPECT_EQ(nBAD,0);
}
} // namespace ro

namespace flush {

class test_helper : public mdl::CACHEhelper {
protected:
    virtual void init(void *dst) override {
        auto p1 = reinterpret_cast<std::uint64_t *>(dst);
        *p1 = 0;
    }
    virtual void combine(void *dst, const void *src, const void *key) override {
        auto p1 = reinterpret_cast<std::uint64_t *>(dst);
        auto p2 = reinterpret_cast<const std::uint64_t *>(src);
        *p1 += *p2;
    }
    virtual uint32_t getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey) override {
        return * static_cast<const uint64_t *>(pKey); // TESTING: we put the remote id here
    }
public:
    explicit test_helper(uint32_t nData, bool bModify=true) : CACHEhelper(sizeof(std::uint64_t),nData,bModify) {}
};

constexpr int SERVICE = worker::TEST_ADVANCED_FLUSH;
typedef std::array<uint64_t,4> KEY;
typedef mdl::hash::HASH<KEY> HASH;

int test(worker::Context *ctx,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = static_cast<mdl::mdlClass *>(ctx->getMDL());
    auto pnBAD = reinterpret_cast<std::uint64_t *>(vout);
    std::uint64_t nBAD;
    if (ctx->getLeaves() > 1) {
        int rID = mdl->ReqService(ctx->getUpper(),SERVICE,NULL,0);
        test(ctx->getLower(),vin,nIn,vout,nOut);
        nOut = mdl->GetReply(rID,nBAD);
        *pnBAD += nBAD;
    }
    else {
        uint32_t idSelf = mdl->Self();
        uint32_t nThreads = mdl->Threads();
        uint32_t nData = cacheSize / 10; // Smaller is fine for this

        nBAD = 0;

        // Construct a hash table filled with data
        HASH hash(nData);
        std::vector<uint64_t> values;
        values.reserve(nData);
        for (auto i=0; i<nData; ++i) {
            values.push_back(0xdeadbeef ^ idSelf ^ i);
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(values.back()),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            hash.insert(uHash,key,&values.back());
        }

        mdl->AdvancedCacheInitialize(0,&hash,sizeof(values.front()),std::make_shared<test_helper>(sizeof(values.front()),true));
        auto idNext = (idSelf+1) % nThreads;

        // These should be cache hits
        for (auto i=idSelf; i<nData; i += nThreads) {
            auto expected = (0xdeadbeef ^ idNext ^ i);
            KEY key = {static_cast<uint64_t>(idNext),static_cast<uint64_t>(i),static_cast<uint64_t>(expected),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = mdl->KeyFetch<uint64_t>(0,uHash,&key,0,1,0);
            if (data == nullptr || *data != 0) ++nBAD;
        }
        mdl->CacheBarrier(0);

        for (auto i=0; i<nData; ++i) {
            for (auto iProc=0; iProc<mdl->Threads(); ++iProc) {
                auto expected = (0xdeadbeef ^ iProc ^ i);
                KEY key = {static_cast<uint64_t>(iProc),static_cast<uint64_t>(i),static_cast<uint64_t>(expected),static_cast<uint64_t>(nData-i)};
                auto uHash = mdl::murmur::murmur<4>(key.data());
                auto data = mdl->KeyAcquire<uint64_t>(0,uHash,&key);
                if (data == nullptr) ++nBAD;
                else ++*data;
                mdl->Release(0,data);
            }
        }
        mdl->FinishCache(0);

        for (auto i=0; i<nData; ++i) {
            auto expected = (0xdeadbeef ^ idSelf ^ i);
            KEY key = {static_cast<uint64_t>(idSelf),static_cast<uint64_t>(i),static_cast<uint64_t>(expected),static_cast<uint64_t>(nData-i)};
            auto uHash = mdl::murmur::murmur<4>(key.data());
            auto data = static_cast<uint64_t *>(hash.lookup(uHash,key));
            if (data == nullptr || *data != expected+nThreads) ++nBAD;
        }

        *pnBAD = nBAD;
    }

    return sizeof(*pnBAD);
}

TEST_F(AdvancedCache, Flush) {
    std::uint64_t nBAD;
    master_mdl->RunService(SERVICE,&nBAD);
    EXPECT_EQ(nBAD,0);
}
} // namespace flush

} // namespace advanced
} // namespace test

/*
** This function is called at the very start by every thread.
** It returns the "worker context"; in this case the ctx.
*/
void *worker_init(MDL vmdl) {
    auto mdl = static_cast<mdl::mdlClass *>(vmdl);
    auto ctx = new worker::Context(mdl);
    mdl->SetCacheSize(cacheSize); // Small cache to make it faster to test
    mdl->AddService(worker::SET_ADD,ctx,(fcnService_t *)SetAdd::serviceSetAdd, sizeof(SetAdd::inSetAdd),0);
    mdl->AddService(test::hash::SERVICE,ctx,(fcnService_t *)test::hash::test, 0,sizeof(std::uint64_t));
    mdl->AddService(test::ro::SERVICE,ctx,(fcnService_t *)test::ro::test,0,sizeof(std::uint64_t));
    mdl->AddService(test::flush::SERVICE,ctx,(fcnService_t *)test::flush::test, 0,sizeof(std::uint64_t));
    mdl->AddService(test::flush::after_read::SERVICE,ctx,(fcnService_t *)test::flush::after_read::test,0,sizeof(std::uint64_t));
    mdl->AddService(test::advanced::ro::SERVICE,ctx,(fcnService_t *)test::advanced::ro::test,0,sizeof(std::uint64_t));
    mdl->AddService(test::advanced::flush::SERVICE,ctx,(fcnService_t *)test::advanced::flush::test,0,sizeof(std::uint64_t));
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
int master(MDL vmdl,void *vctx) {
    master_mdl = static_cast<mdl::mdlClass *>(vmdl);
    auto ctx = reinterpret_cast<worker::Context *>(vctx);
    int argc = master_mdl->GetArgc();
    char **argv = master_mdl->GetArgv();

    SetAdd::inSetAdd in;
    in.idUpper = master_mdl->Threads();
    SetAdd::serviceSetAdd(ctx,&in,sizeof(in),NULL,0);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

int main(int argc,char **argv) {
    return mdl::Launch(master,worker_init,worker_done,argc,argv);
}
