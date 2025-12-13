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
#include "mdl.h"
using namespace mdl;

#include <boost/preprocessor.hpp>
#include <boost/preprocessor/facilities/is_empty.hpp>

// Helper macro to generate a custom parameter name: p1, p2, etc.
#define CUSTOM_PARAM_NAME(z, n, data) BOOST_PP_CAT(p, BOOST_PP_INC(n))

// Generate a parameter declaration from a type.
// Produces: <type> p<index>
#define PARAM_DECL(z, n, seq) BOOST_PP_SEQ_ELEM(n, seq) BOOST_PP_CAT(p, BOOST_PP_INC(n))

// Generate a comma-separated list of parameter names for the call.
#define ENUM_PARAM_NAMES(count) BOOST_PP_ENUM(count, CUSTOM_PARAM_NAME, ~)

// Macro for the no‑parameter case.
#define MDL_STUB_IMPL0(name)                                                                                                                                   \
    extern "C" auto mdl##name(MDL cmdl) -> decltype(static_cast<mdlClass*>(cmdl)->name()) {                                                                    \
        return static_cast<mdlClass*>(cmdl)->name();                                                                                                           \
    }

// Macro for the case with one or more parameters.
// __VA_ARGS__ should be a comma‑separated list of types.
#define MDL_STUB_IMPL(name, ...)                                                                                                                               \
    extern "C" auto mdl##name(MDL cmdl, BOOST_PP_ENUM(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__), PARAM_DECL, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))                 \
        -> decltype(static_cast<mdlClass*>(cmdl)->name(ENUM_PARAM_NAMES(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__)))) {                                               \
        return static_cast<mdlClass*>(cmdl)->name(ENUM_PARAM_NAMES(BOOST_PP_VARIADIC_SIZE(__VA_ARGS__)));                                                      \
    }

// The master macro that selects the proper implementation based on
// whether __VA_ARGS__ is empty. Generates something like:
// extern "C" auto mdlFunction(MDL cmdl, int p1, float p2, double p3)
//   -> decltype(static_cast<mdlClass*>(cmdl)->Function(p1, p2, p3)) {
//     return static_cast<mdlClass*>(cmdl)->Function(p1, p2, p3);
// }
#define MDL_STUB(name, ...) BOOST_PP_IF(BOOST_PP_IS_EMPTY(__VA_ARGS__), MDL_STUB_IMPL0(name), MDL_STUB_IMPL(name, __VA_ARGS__))

extern "C" int mdlLaunch(int (*fcnMaster)(MDL, void*), void* (*fcnWorkerInit)(MDL), void (*fcnWorkerDone)(MDL, void*), int argc, char** argv) {
    return Launch(fcnMaster, fcnWorkerInit, fcnWorkerDone, argc, argv);
}

MDL_STUB(Threads)
MDL_STUB(Self)
MDL_STUB(Core)
MDL_STUB(Cores)
MDL_STUB(Proc)
MDL_STUB(Procs)
MDL_STUB(ProcToThread, int)
MDL_STUB(Malloc, size_t)
MDL_STUB(Free, void*)
MDL_STUB(SetCacheSize, int)
MDL_STUB(CacheCheck)
MDL_STUB(CacheStatus, int)
MDL_STUB(CacheBarrier, int)
MDL_STUB(FlushCache, int)
MDL_STUB(ThreadBarrier)
MDL_STUB(CompleteAllWork)
MDL_STUB(Fetch, int, int, int)
MDL_STUB(KeyFetch, int, uint32_t, void*, int, int, int)
MDL_STUB(Acquire, int, int, int)
MDL_STUB(KeyAcquire, int, uint32_t, void*)
MDL_STUB(VirtualFetch, int, int, int)
MDL_STUB(Release, int, void*)
MDL_STUB(AcquireWrite, int, int)
MDL_STUB(ReleaseWrite, int, void*)
MDL_STUB(FinishCache, int)
MDL_STUB(CudaActive)
MDL_STUB(NumAccess, int)
MDL_STUB(MissRatio, int)

class legacyCACHEhelper : public arrayCACHEhelper {
    static void* getArrayElement(void* vData, int i, int iDataSize) {
        char* pData = static_cast<char*>(vData);
        return pData + (size_t)i * (size_t)iDataSize;
    }

  protected:
    void* ctx;
    void* (*get_function)(void* pData, int i, int iDataSize);
    void (*init_function)(void*, void*);
    void (*combine_function)(void*, void*, const void*);
    virtual void* get_element(uint32_t i) override { return (*get_function)(pArray, i, nDataSize); }
    virtual void init(void* dst) override {
        if (init_function) (*init_function)(ctx, dst);
    }
    virtual void combine(void* dst, const void* src, const void* key) override {
        if (combine_function) (*combine_function)(ctx, dst, src);
    }

  public:
    explicit legacyCACHEhelper(void* pArray, uint32_t nArray, uint32_t nDataSize, void* ctx, void* (*get_function)(void*, int, int),
                               void (*init_function)(void*, void*) = nullptr, void (*combine_function)(void*, void*, const void*) = nullptr)
        : arrayCACHEhelper(nDataSize, pArray, nArray, combine_function != nullptr), ctx(ctx),
          get_function(get_function == NULL ? getArrayElement : get_function), init_function(init_function), combine_function(combine_function) {}
};

// This opens a read-only cache. Called from a worker outside of MDL
extern "C" void mdlROcache(MDL mdl, int cid, void* (*getElt)(void* pData, int i, int iDataSize), void* pData, int iDataSize, int nDataSize) {
    static_cast<mdlClass*>(mdl)->CacheInitialize(cid, std::make_shared<legacyCACHEhelper>(pData, nDataSize, iDataSize, nullptr, getElt, nullptr, nullptr));
}

// This opens a combiner (read/write) cache. Called from a worker outside of MDL
extern "C" void mdlCOcache(MDL mdl, int cid, void* (*getElt)(void* pData, int i, int iDataSize), void* pData, int iDataSize, int nData, void* ctx,
                           void (*init)(void*, void*), void (*combine)(void*, void*, const void*)) {
    assert(init);
    assert(combine);
    static_cast<mdlClass*>(mdl)->CacheInitialize(cid, std::make_shared<legacyCACHEhelper>(pData, nData, iDataSize, ctx, getElt, init, combine));
}
