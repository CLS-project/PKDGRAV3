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
#include "mdl_config.h"
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

static inline uint64_t elapsed(ticks t1, ticks t0) {
    return (t1.tv_sec - t0.tv_sec) * uint64_t(1e9) + (t1.tv_nsec - t0.tv_nsec);
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

typedef int(fcnService_t)(void* p1, void* vin, int nIn, void* vout, int nOut);

#include "mdlbt.h"
#include <span>
#include <concepts>
#include <cstddef>
#include <vector>
#include <string>
#include <string_view>
#include <unordered_map>
#include <memory>
#include <initializer_list>
#define MAX_NODE_NAME_LENGTH 256

namespace mdl {

template <typename C>
concept ByteBufferLike = requires(C& c) {
    { c.data() } -> std::convertible_to<const void*>;
    { c.size() } -> std::convertible_to<std::size_t>;
} && std::is_pointer_v<decltype(std::declval<C&>().data())>;
// && std::is_trivially_copyable_v<std::remove_cv_t<std::remove_pointer_t<decltype(std::declval<C&>().data())>>>;

template <class C>
concept WritableByteBufferLike = ByteBufferLike<C> && !std::is_const_v<std::remove_pointer_t<decltype(std::declval<C&>().data())>>;

// Helper functions to convert various types to byte spans
static inline std::span<std::byte> as_writable_bytes(std::span<std::byte> s) noexcept {
    return s;
}
static inline std::span<const std::byte> as_bytes(std::span<const std::byte> s) noexcept {
    return s;
}
template <WritableByteBufferLike C>
static inline std::span<std::byte> as_writable_bytes(C& c) {
    using T = std::remove_pointer_t<decltype(c.data())>;
    return std::as_writable_bytes(std::span<T>{c.data(), c.size()});
}
template <ByteBufferLike C>
static inline std::span<const std::byte> as_bytes(const C& c) {
    using T = std::remove_pointer_t<decltype(c.data())>;
    return std::as_bytes(std::span<const T>{c.data(), c.size()});
}
template <class T>
    requires(!ByteBufferLike<T>) /* && std::is_trivially_copyable_v<T> */
static inline std::span<std::byte> as_writable_bytes(T& obj) noexcept {
    return {reinterpret_cast<std::byte*>(std::addressof(obj)), sizeof(T)};
}
template <class T>
    requires(!ByteBufferLike<T>) /* && std::is_trivially_copyable_v<T> */
static inline std::span<const std::byte> as_bytes(const T& obj) noexcept {
    return {reinterpret_cast<const std::byte*>(std::addressof(obj)), sizeof(T)};
}

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
    std::vector<std::size_t> offsets{0};
    std::unique_ptr<std::byte[]> buffer;

  protected:
    class BasicField {
        int nBytes;

      public:
        BasicField(int size) : nBytes(size) {}
        std::size_t size() const { return nBytes; }
    };

  public:
    template <typename T>
    class Field : public BasicField {
        int n;

      public:
        Field(int n = 1) : BasicField(sizeof(T) * n) {}
    };

  public:
    ServiceBuffer() = delete;
    ServiceBuffer(const ServiceBuffer&) = delete;
    ServiceBuffer& operator=(const ServiceBuffer&) = delete;

    // Move constructor
    ServiceBuffer(ServiceBuffer&& other) noexcept : offsets(std::move(other.offsets)), buffer(std::move(other.buffer)) {}

    // Move assignment
    ServiceBuffer& operator=(ServiceBuffer&& other) noexcept {
        if (&other != this) {
            offsets = std::move(other.offsets);
            buffer = std::move(other.buffer);
        }
        return *this;
    }

    ServiceBuffer(std::initializer_list<BasicField> t) {
        for (auto& i : t) offsets.push_back(offsets.back() + i.size());
        buffer = std::make_unique<std::byte[]>(offsets.back());
    }
    std::byte* data() { return buffer.get(); }
    const std::byte* data() const { return buffer.get(); }
    auto size() const { return offsets.back(); }
    operator std::span<std::byte>() { return {data(), size()}; }
    operator std::span<const std::byte>() const { return {data(), size()}; }

    void* data(int i) { return static_cast<void*>(buffer.get() + offsets[i]); }
    auto size(int i) { return offsets[i + 1] - offsets[i]; }
};

class ServiceBufferOut : public ServiceBuffer {
  public:
    using ServiceBuffer::ServiceBuffer;
};

class BasicService {
    friend class mdlBASE;
    friend class LegacyService;
    friend class NullService;

  private:
    int nInBytes;
    int nOutBytes;
    int service_id;
    std::string service_name;

  protected:
    void setServiceId(int id) { service_id = id; }
    explicit BasicService(int service_id, std::string_view service_name, int nInBytes, int nOutBytes = 0)
        : nInBytes(nInBytes), nOutBytes(nOutBytes), service_id(service_id), service_name(service_name) {}

  public:
    explicit BasicService(std::string_view service_name, int nInBytes, int nOutBytes = 0)
        : nInBytes(nInBytes), nOutBytes(nOutBytes), service_id(-1), service_name(service_name) {}
    virtual ~BasicService() = default;
    int getServiceID() { return service_id; }
    std::string_view getServiceName() const { return service_name; }
    int getMaxBytesIn() { return nInBytes; }
    int getMaxBytesOut() { return nOutBytes; }

  protected:
    virtual int operator()(std::span<std::byte> In, std::span<std::byte> Out) = 0;
};

class mdlBASE : protected mdlbt {
  protected:
    void CommitServices(size_t nExtraInBytes = 0, size_t nExtraOutBytes = 0);
    std::vector<std::byte> input_buffer;
    std::vector<std::byte> reply_buffer;

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
    int bDiag; /* When true, debug output is enabled */
    int argc;

    FILE* fpDiag;
    char** argv;

    /* Services information */
    int nMaxInBytes;
    int nMaxOutBytes;
    std::vector<std::unique_ptr<BasicService>> services;
    std::unordered_map<std::string_view, int> service_map;

    /* Maps a give process (Proc) to the first global thread ID */
    std::vector<int> iProcToThread; /* [0,nProcs] (note inclusive extra element) */

    char nodeName[MAX_NODE_NAME_LENGTH];

    typedef enum { TIME_WAITING, TIME_COMPUTING, TIME_SYNCHRONIZING, TIME_COUNT } TICK_TIMER;

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

    void mdl_vprintf(const char* format, va_list ap);
    void mdl_printf(const char* format, ...);

  public:
    explicit mdlBASE(int argc, char** argv);
    virtual ~mdlBASE();
    int32_t Threads() const { return layout.nThreads; }
    int32_t Self() const { return layout.idSelf; }
    int32_t Core() const { return layout.iCore; }
    int32_t Cores() const { return layout.nCores; }
    int32_t Proc() const { return layout.iProc; }
    int32_t Procs() const { return layout.nProcs; }
    int32_t ProcToThread(int32_t iProc) const;
    int32_t ThreadToProc(int32_t iThread) const;
    int GetArgc() const { return argc; }
    char** GetArgv() const { return argv; }

    void yield();
    int AddService(int sid, void* p1, fcnService_t* fcnService, int nInBytes, int nOutBytes, std::string_view service_name);
    int AddService(std::unique_ptr<BasicService> service);
    template <class T, class... Args>
    int AddService(Args&&... args) {
        static_assert(std::is_base_of_v<BasicService, T>, "T must derive from BasicService");
        return AddService(std::make_unique<T>(std::forward<Args>(args)...));
    }

  public:
    // Look up the service by ID and run it (base version) -- soon private
    /// @brief Run a service by ID with input and output as byte spans
    /// @param sid Service ID
    /// @param in Input data as a span of bytes
    /// @param out Output data as a span of bytes
    /// @return Number of bytes written to output
    int RunService(int sid, std::span<std::byte> in = {}, std::span<std::byte> out = {});

    // Convert service name to ID and run it
    /// @brief Run a service by name with byte spans for input and output
    /// @tparam I Input container type
    /// @tparam O Output container type
    /// @param service_name Name of the service to run
    /// @param in Input data as a span (converted to a span of bytes)
    /// @param out Output data as a span of bytes (converted if necessary)
    /// @return Number of bytes written to output
    template <typename I = std::byte, typename O = std::byte>
    int RunService(std::string_view service_name, std::span<I> in = {}, std::span<O> out = {}) {
        return RunService(service_map.at(service_name), as_writable_bytes(in), as_writable_bytes(out));
    }

    /// @brief Run a service by name with input and output data/containers by reference
    /// @tparam I Input type
    /// @tparam O Output type
    /// @param service_name Name of the service to run
    /// @param in Input data/container
    /// @param out Output data/container
    /// @return Number of bytes written to output
    template <typename I, typename O>
    int RunService(std::string_view service_name, I& in, O& out) {
        return RunService(service_name, as_writable_bytes(in), as_writable_bytes(out));
    }

    /// @brief Run a service by name with input data/container by reference and output as a byte span
    /// @tparam I Input type
    /// @tparam O Output container type
    /// @param service_name Name of the service to run
    /// @param in Input data/container
    /// @param out Output data as a span of bytes (converted if necessary)
    /// @return Number of bytes written to output
    template <typename I, typename O = std::byte>
    int RunService(std::string_view service_name, I& in, std::span<O> out = {}) {
        return RunService(service_name, as_writable_bytes(in), out);
    }

    /// @brief Run a service by name with input as a byte span and output data/container by reference
    /// @tparam I Input container type
    /// @tparam O Output type
    /// @param service_name Name of the service to run
    /// @param in Input data as a span (converted to a span of bytes)
    /// @param out Output data/container
    /// @return Number of bytes written to output
    template <typename I = std::byte, typename O>
    int RunService(std::string_view service_name, std::span<I> in, O& out) {
        return RunService(service_name, in, as_writable_bytes(out));
    }

    /// @brief Run a service by type with input and output as byte spans
    /// @tparam T Service type with a static name member
    /// @tparam I Input container type
    /// @tparam O Output container type
    /// @param in Input data as a span of bytes
    /// @param out Output data as a span of bytes
    /// @return Number of bytes written to output
    template <typename T, typename I = std::byte, typename O = std::byte>
    int RunService(std::span<I> in = {}, std::span<O> out = {}) {
        return RunService(T::name, in, out);
    }

    /// @brief Run a service by type with input and output data/containers by reference
    /// @tparam T Service type with a static name member
    /// @tparam I Input type
    /// @tparam O Output type
    /// @param in Input data/container
    /// @param out Output data/container
    /// @return Number of bytes written to output
    template <typename T, typename I, typename O>
    int RunService(I& in, O& out) {
        return RunService(T::name, as_writable_bytes(in), as_writable_bytes(out));
    }

    /// @brief Run a service by type with input data/container by reference and output as a byte span
    /// @tparam T Service type with a static name member
    /// @tparam I Input type
    /// @tparam O Output container type
    /// @param in Input data/container
    /// @param out Output data as a span of bytes (converted if necessary)
    /// @return Number of bytes written to output
    template <typename T, typename I, typename O = std::byte>
    int RunService(I& in, std::span<O> out = {}) {
        return RunService(T::name, as_writable_bytes(in), out);
    }

    /// @brief Run a service by type with input as a byte span and output data/container by reference
    /// @tparam T Service type with a static name member
    /// @tparam I Input container type
    /// @tparam O Output type
    /// @param in Input data as a span (converted to a span of bytes)
    /// @param out Output data/container
    /// @return Number of bytes written to output
    template <typename T, typename I = std::byte, typename O>
    int RunService(std::span<I> in, O& out) {
        return RunService(T::name, in, as_writable_bytes(out));
    }

    // Legacy overloads for compatibility
    int RunService(int sid, int nIn, void* pIn, void* pOut, int nOut);
    int RunService(int sid, ServiceBuffer& b) { return RunService(sid, b.size(), b.data(), nullptr, 0); }
    int RunService(int sid, ServiceBuffer& b, void* pOut, int nOut) { return RunService(sid, b.size(), b.data(), pOut, nOut); }
    int RunService(int sid, int nIn, void* pIn, ServiceBufferOut& out) { return RunService(sid, nIn, pIn, out.data(), out.size()); }
    int RunService(int sid, ServiceBuffer& b, ServiceBufferOut& out) { return RunService(sid, b.size(), b.data(), out.data(), out.size()); }
    int RunService(int sid, ServiceBufferOut& out) { return RunService(sid, 0, nullptr, out.data(), out.size()); }
    int RunService(int sid, void* pOut, int nOut) { return RunService(sid, 0, nullptr, pOut, nOut); }

  public:
    BasicService* GetService(unsigned sid) { return sid < services.size() ? services[sid].get() : nullptr; }
    BasicService* GetService(std::string_view service_name) {
        auto it = service_map.find(service_name);
        if (it != service_map.end()) return services[it->second].get();
        return nullptr;
    }
};
int mdlBaseProcToThread(mdlBASE* base, int iProc);
int mdlBaseThreadToProc(mdlBASE* base, int iThread);
} // namespace mdl

const char* mdlName(void* mdl);

void mdlDiag(void* mdl, const char* psz);
#ifdef MDLASSERT
    #ifndef __STRING
        #define __STRING(arg) (("arg"))
    #endif
    #define mdlassert(mdl, expr)                                                                                                                               \
        {                                                                                                                                                      \
            if (!(expr)) {                                                                                                                                     \
                mdlprintf(mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, __STRING(expr));                                                          \
                assert(expr);                                                                                                                                  \
            }                                                                                                                                                  \
        }
#else
    #define mdlassert(mdl, expr) assert(expr)
#endif

double mdlCpuTimer(void* mdl);
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
void mdlZeroTimer(void* mdl, mdlTimer*);
void mdlGetTimer(void* mdl, mdlTimer*, mdlTimer*);
void mdlPrintTimer(void* mdl, const char* message, mdlTimer*);
#else
    #define mdlZeroTimer
    #define mdlGetTimer
    #define mdlPrintTimer
#endif

#endif
