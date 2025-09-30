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

#ifndef MDLBASE_H
#define MDLBASE_H
#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "mdl_config.h"
#endif
#ifdef INSTRUMENT
#include "cycle.h"

/* Last effort to find a clock */
#if defined(HAVE_CLOCK_GETTIME) && !defined(HAVE_TICK_COUNTER)
typedef struct timespec ticks;

static inline ticks getticks(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t;
}

static inline double elapsed(ticks t1, ticks t0) {
    return ((double)t1.tv_sec - (double)t0.tv_sec) * 1.0E9 +
           ((double)t1.tv_nsec - (double)t0.tv_nsec);
}
#define HAVE_TICK_COUNTER
#endif

#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/*
* Compile time mdl debugging options
*
* mdl asserts: define MDLASSERT
* Probably should always be on unless you want no mdlDiag output at all
*
* NB: defining NDEBUG turns off all asserts so MDLASSERT will not assert
* however it will output uding mdlDiag and the code continues.
*/
#define MDLASSERT

typedef int (fcnService_t)(void *p1, void *vin, int nIn, void *vout, int nOut);

#ifdef __cplusplus
#include "mdlbt.h"
#include <vector>
#include <string>
#include <string_view>
#include <unordered_map>
#include <memory>
#include <initializer_list>
#define MAX_NODE_NAME_LENGTH      256

namespace mdl {

// This struct is use to keep track of a variable length buffer for service requests.
// List all data types and the number of times they are repeated in the constructor.
//     ServiceBuffer msg {
//         ServiceBuffer::Field<ServiceWhatever::input>(),
//         ServiceBuffer::Field<int>(100)
//         };
// This would create a message with a single "input", followed by 100 integers.
// Access the individual fields with "data()".
//    auto input = static_cast<ServiceWhatever::input*>(msg.data(0));
//    auto start = static_cast<int*>(msg.data(1));
// The RunService() method has a overload that takes a ServiceBuffer for input.
class ServiceBuffer {
    std::vector<size_t> offsets {0};
    std::unique_ptr<char[]> buffer;
protected:
    class BasicField {
        int nBytes;
    public:
        BasicField(int size) : nBytes(size) {}
        size_t size() const {return nBytes; }
    };
public:
    template<typename T>
    class Field : public BasicField {
        int n;
    public:
        Field(int n=1) : BasicField(sizeof(T)*n) {}
    };
public:
    ServiceBuffer() = delete;
    ServiceBuffer(const ServiceBuffer &) = delete;
    ServiceBuffer &operator=(const ServiceBuffer &) = delete;

    // Move constructor
    ServiceBuffer(ServiceBuffer &&other) noexcept
        : offsets(std::move(other.offsets)), buffer(std::move(other.buffer))
    {}

    // Move assignment
    ServiceBuffer &operator=(ServiceBuffer &&other) noexcept {
        if (&other != this) {
            offsets = std::move(other.offsets);
            buffer = std::move(other.buffer);
        }
        return *this;
    }

    ServiceBuffer(std::initializer_list<BasicField> t) {
        for (auto &i : t) offsets.push_back(offsets.back() + i.size());
        buffer = std::make_unique<char[]>(offsets.back());
    }
    void *data(int i=0) {return static_cast<void *>(buffer.get() + offsets[i]);}
    auto size() {return offsets.back();}
    auto size(int i) {return offsets[i+1]-offsets[i];}
};

class ServiceBufferOut : public ServiceBuffer {
public:
    using ServiceBuffer::ServiceBuffer;
};

class BasicService {
    friend class mdlBASE;
private:
    int nInBytes;
    int nOutBytes;
    int service_id;
    std::string service_name;
protected:
    void setServiceId(int id) { service_id = id; }
public:
    explicit BasicService(int service_id, std::string_view service_name, int nInBytes, int nOutBytes=0)
        : nInBytes(nInBytes), nOutBytes(nOutBytes), service_id(service_id),service_name(service_name) {}
    virtual ~BasicService() = default;
    int getServiceID()  {return service_id;}
    std::string_view getServiceName() const {return service_name;}
    int getMaxBytesIn() {return nInBytes;}
    int getMaxBytesOut() {return nOutBytes;}
protected:
    virtual int operator()(int nIn, void *pIn, void *pOut) = 0;
};

class mdlBASE : protected mdlbt {
public:
#ifdef USE_BT
    using mdlbt::ignore_SIGBUS;
    using mdlbt::show_backtrace;
    using mdlbt::show_backtrace_with_snippets;
#endif
    struct {
        int32_t nThreads; /* Global number of threads (total) */
        int32_t idSelf;   /* Global index of this thread */
        int32_t nProcs;   /* Number of global processes (e.g., MPI ranks) */
        int32_t iProc;    /* Index of this process (MPI rank) */
        int32_t nCores;   /* Number of threads in this process */
        int32_t iCore;    /* Local core id */
    } layout;
    int bDiag;        /* When true, debug output is enabled */
    int argc;

    FILE *fpDiag;
    char **argv;

    /* Services information */
    int nMaxInBytes;
    int nMaxOutBytes;
    std::vector< std::unique_ptr<BasicService>> services;
    std::unordered_map<std::string_view,int> service_map;

    /* Maps a give process (Proc) to the first global thread ID */
    std::vector<int> iProcToThread; /* [0,nProcs] (note inclusive extra element) */

    char nodeName[MAX_NODE_NAME_LENGTH];

    typedef enum {
        TIME_WAITING,
        TIME_COMPUTING,
        TIME_SYNCHRONIZING,
        TIME_COUNT
    } TICK_TIMER;

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
protected:
    ticks nTicks;
#endif
    double dTimer[TIME_COUNT];
private:
    double TimeFraction() const;
public:
    void TimeReset();
    void TimeAddComputing();
    void TimeAddSynchronizing();
    void TimeAddWaiting();
    double TimeGet(TICK_TIMER which) const { return dTimer[which] * TimeFraction(); }
    double TimeComputing() const;
    double TimeSynchronizing() const;
    double TimeWaiting() const;

    void mdl_vprintf(const char *format, va_list ap);
    void mdl_printf(const char *format, ...);

public:
    explicit mdlBASE(int argc,char **argv);
    virtual ~mdlBASE();
    int32_t Threads() const { return layout.nThreads; }
    int32_t Self()    const { return layout.idSelf; }
    int32_t Core()    const { return layout.iCore; }
    int32_t Cores()   const { return layout.nCores; }
    int32_t Proc()    const { return layout.iProc; }
    int32_t Procs()   const { return layout.nProcs; }
    int32_t ProcToThread(int32_t iProc) const;
    int32_t ThreadToProc(int32_t iThread) const;
    int     GetArgc() const { return argc; }
    char **GetArgv() const { return argv; }

    void yield();
    void AddService(int sid, void *p1, fcnService_t *fcnService, int nInBytes, int nOutBytes, std::string_view service_name);
    void AddService(std::unique_ptr<BasicService> &&service);
    int  RunService(int sid, int nIn, void *pIn, void *pOut=nullptr);
    int  RunService(int sid, ServiceBuffer &b, void *pOut=nullptr) { return RunService(sid,b.size(),b.data(),pOut);}
    int  RunService(int sid, int nIn, void *pIn, ServiceBufferOut &out) { return RunService(sid,nIn,pIn,out.data());}
    int  RunService(int sid, ServiceBuffer &b, ServiceBufferOut &out) { return RunService(sid,b.size(),b.data(),out.data());}
    int  RunService(int sid, ServiceBufferOut &out) { return RunService(sid,0,nullptr,out.data());}
    int  RunService(int sid, void *pOut=nullptr) { return RunService(sid,0,nullptr,pOut); }
public:
    int  RunService(std::string_view service_name, int nIn, void *pIn, void *pOut=nullptr);
    int  RunService(std::string_view service_name, ServiceBuffer &b, void *pOut=nullptr) { return RunService(service_name,b.size(),b.data(),pOut);}
    int  RunService(std::string_view service_name, int nIn, void *pIn, ServiceBufferOut &out) { return RunService(service_name,nIn,pIn,out.data());}
    int  RunService(std::string_view service_name, ServiceBuffer &b, ServiceBufferOut &out) { return RunService(service_name,b.size(),b.data(),out.data());}
    int  RunService(std::string_view service_name, ServiceBufferOut &out) { return RunService(service_name,0,nullptr,out.data());}
    int  RunService(std::string_view service_name, void *pOut=nullptr) { return RunService(service_name,0,nullptr,pOut); }
public:
    template<typename T>
    int  RunService(int nIn, void *pIn, void *pOut=nullptr) {
        return RunService(T::name, nIn, pIn, pOut);
    }
    template<typename T>
    int  RunService(ServiceBuffer &b, void *pOut=nullptr) {
        return RunService<T>(b.size(), b.data(), pOut);
    }
    template<typename T>
    int  RunService(int nIn, void *pIn, ServiceBufferOut &out) {
        return RunService<T>(nIn, pIn, out.data());
    }
    template<typename T>
    int  RunService(ServiceBuffer &b, ServiceBufferOut &out) {
        return RunService<T>(b.size(), b.data(), out.data());
    }
    template<typename T>
    int  RunService(ServiceBufferOut &out) {
        return RunService<T>(0, nullptr, out.data());
    }
    template<typename T>
    int  RunService(void *pOut=nullptr) {
        return RunService<T>(0, nullptr, pOut);
    }
    BasicService *GetService(unsigned sid) {return sid<services.size() ? services[sid].get() : nullptr; }
};
int mdlBaseProcToThread(mdlBASE *base, int iProc);
int mdlBaseThreadToProc(mdlBASE *base, int iThread);
} // namespace mdl

#endif

#ifdef __cplusplus
extern "C" {
#endif

const char *mdlName(void *mdl);

void mdlDiag(void *mdl, const char *psz);
#ifdef MDLASSERT
#ifndef __STRING
#define __STRING( arg )   (("arg"))
#endif
#define mdlassert(mdl,expr) \
    { \
    if (!(expr)) { \
    mdlprintf( mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, __STRING(expr) ); \
    assert( expr ); \
        } \
    }
#else
#define mdlassert(mdl,expr)  assert(expr)
#endif

double mdlCpuTimer(void *mdl);
/*
* Timer functions active: define MDLTIMER
* Makes mdl timer functions active
*/
#ifndef _CRAYMPP
#define MDLTIMER
#endif

typedef struct {
    double wallclock;
    double cpu;
    double system;
} mdlTimer;

#ifdef MDLTIMER
void mdlZeroTimer(void *mdl, mdlTimer *);
void mdlGetTimer(void *mdl, mdlTimer *, mdlTimer *);
void mdlPrintTimer(void *mdl, const char *message, mdlTimer *);
#else
#define mdlZeroTimer
#define mdlGetTimer
#define mdlPrintTimer
#endif

#ifdef __cplusplus
}
#endif
#endif
