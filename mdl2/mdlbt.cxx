#include "mdl_config.h"
#include "mdlbt.h"
#include <iostream>
#include <cassert>
#include <unistd.h>
#include <stdexcept>
#include <execinfo.h>
#include <pthread.h>
#include <mutex>
#include <condition_variable>
#include <sys/wait.h>
#if defined(USE_CPPTRACE)
#include <cpptrace/cpptrace.hpp>
#elif defined(USE_ELFUTILS)
#include <elfutils/libdwfl.h>
#include <boost/core/demangle.hpp>

struct DebugInfoSession {
    Dwfl_Callbacks callbacks = {};
    char *debuginfo_path = nullptr;
    Dwfl *dwfl = nullptr;

    DebugInfoSession() {
        callbacks.find_elf = dwfl_linux_proc_find_elf;
        callbacks.find_debuginfo = dwfl_standard_find_debuginfo;
        callbacks.debuginfo_path = &debuginfo_path;

        dwfl = dwfl_begin(&callbacks);
        assert(dwfl);

        int r;
        r = dwfl_linux_proc_report(dwfl, getpid());
        assert(!r);
        r = dwfl_report_end(dwfl, nullptr, nullptr);
        assert(!r);
        static_cast<void>(r);
    }

    ~DebugInfoSession() {
        dwfl_end(dwfl);
    }

    DebugInfoSession(DebugInfoSession const &) = delete;
    DebugInfoSession &operator=(DebugInfoSession const &) = delete;
};

struct DebugInfo {
    void *ip;
    std::string function;
    char const *file;
    int line;

    DebugInfo(DebugInfoSession const &dis, void *ip)
        : ip(ip)
        , file()
        , line(-1) {
        // Get function name.
        uintptr_t ip2 = reinterpret_cast<uintptr_t>(ip);
        Dwfl_Module *module = dwfl_addrmodule(dis.dwfl, ip2);
        char const *name = dwfl_module_addrname(module, ip2);
        function = name ? boost::core::demangle(name) : "<unknown>";

        // Get source filename and line number.
        if (Dwfl_Line *dwfl_line = dwfl_module_getsrc(module, ip2)) {
            Dwarf_Addr addr;
            file = dwfl_lineinfo(dwfl_line, &addr, &line, nullptr, nullptr, nullptr);
        }
        else {
            Dwarf_Addr bias = 0;
            if (Dwarf *dwarf = dwfl_module_getdwarf(module, &bias)) {
                const uintptr_t adjusted = ip2 - bias;
                size_t headerSize = 0;
                Dwarf_Off nextOffset = 0;
                for (Dwarf_Off offset = 0;
                dwarf_nextcu(dwarf, offset, &nextOffset, &headerSize, nullptr, nullptr, nullptr) == 0;
                        offset = nextOffset) {
                    Dwarf_Die cudieMemory;
                    Dwarf_Die *cudie = dwarf_offdie(dwarf, offset + headerSize, &cudieMemory);
                    if (!cudie || !dwarf_haspc(cudie, adjusted)) continue;
                    if (Dwarf_Line *lineinfo = dwarf_getsrc_die(cudie, adjusted)) {
                        file = dwarf_linesrc(lineinfo, nullptr, nullptr);
                        dwarf_lineno(lineinfo, &line);
                        //dwarf_linecol(lineinfo, &column);
                    }
                    break;
                }
            }
        }
    }
};

std::ostream &operator<<(std::ostream &s, DebugInfo const &di) {
    s << di.ip << ' ' << di.function;
    if (di.file)
        s << " at " << di.file << ':' << di.line;
    return s;
}

#endif

// We allow only a single thread to generate a backtrace at a time
static bool backtrace_running = false;
static std::mutex backtrace_mutex;
static std::condition_variable backtrace_ready_cv;
static std::condition_variable backgrace_done_cv;

// When we want to do a backtrace, the thread will acquire the backtrace mutext,
// then fill in the backtrace_frames array with the backtrace, and then notify
// the backtrace_cv condition variable. The main thread will wait on this
// condition variable and when it is notified, it will print the backtrace.
// If we do not have cpptrace, we will use the backtrace_stack array.
constexpr std::size_t backtrace_N_MAX = 100;
#ifdef USE_CPPTRACE
    static cpptrace::frame_ptr backtrace_frames[backtrace_N_MAX];
#endif
static void *backtrace_stack[backtrace_N_MAX]; // Used for the backtrace() call
static std::size_t backtrace_count = 0;

// Here we register the function that will generate a backtrace when an
// abnormal termination occurs.
void mdlbt::register_backtrace() {
    std::unique_lock<std::mutex> lock(backtrace_mutex);
    backtrace_running = true;
    backtrace_thread = std::thread(&mdlbt::BacktraceThread, this);

    // std::set_terminate(terminate_handler);

    struct sigaction sa;
    sa.sa_flags = SA_SIGINFO; // | SA_RESETHAND;
    sigfillset(&sa.sa_mask);
    sa.sa_sigaction = signal_handler;
    if (sigaction(SIGBUS, &sa, NULL) == -1) abort();
    if (sigaction(SIGILL, &sa, NULL) == -1) abort();
    if (sigaction(SIGSEGV, &sa, NULL) == -1) abort();
    if (sigaction(SIGABRT, &sa, NULL) == -1) abort();
    if (sigaction(SIGFPE, &sa, NULL) == -1) abort();
#if defined(__APPLE__) && (defined(__arm64__) || defined(__aarch64__))
    // NULL pointer dereference on arm64 is SIGTRAP
    if (sigaction(SIGTRAP, &sa, NULL) == -1) abort();
#endif
}

void mdlbt::ignore_SIGBUS() {
    struct sigaction sa;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_SIGINFO;
    sa.sa_sigaction = signal_sigbus;
    if (sigaction(SIGBUS, &sa, NULL) == -1) abort();
}

void *mdlbt::BacktraceThread() {
    while (true) {
        std::unique_lock<std::mutex> lock(backtrace_mutex);

        // Wait until signaled that the buffer is filled
        backtrace_ready_cv.wait(lock, [] { return backtrace_count || !backtrace_running; });
        if (!backtrace_running) break;

        // Simulate processing the backtrace-filled buffer
#ifdef USE_CPPTRACE
        cpptrace::raw_trace trace;
        for (auto i=0; i<backtrace_count; i++) {
            trace.frames.push_back(backtrace_frames[i]);
        }
        trace.resolve().print_with_snippets();
#else
        print_backtrace(backtrace_stack, backtrace_count);
#endif

        std::cerr.flush();
        // Clear the buffer and reset flags
        backtrace_count = 0;
        // Notify the thread that filled the buffer that processing is done
        backgrace_done_cv.notify_one();
    }
    return NULL;
}

mdlbt::~mdlbt() {
    if (backtrace_thread.joinable()) {
        backtrace_running = false;
        backtrace_ready_cv.notify_one();
        backtrace_thread.join();
    }
}

void mdlbt::terminate_handler() {
    std::unique_lock<std::mutex> lck(backtrace_mutex);
    show_backtrace_with_snippets();
    std::_Exit(EXIT_FAILURE);
}

void mdlbt::signal_sigbus(int signo, siginfo_t *si, void *unused) {
    std::unique_lock<std::mutex> lck(backtrace_mutex);
    char nodeName[200];

    if (gethostname(nodeName, sizeof(nodeName)))
        nodeName[0] = 0;
    else
        nodeName[sizeof(nodeName) - 1] = 0;

    std::cerr << "SIGNAL " << signo << "." << si->si_code << " at " << si->si_addr << " on " << nodeName << '\n';
}

void mdlbt::signal_handler(int signo, siginfo_t *si, void *context) {
    std::unique_lock<std::mutex> lock(backtrace_mutex); // One thread at a time
    std::cerr << "SIGNAL " << signo << " at " << si->si_addr << '\n';
#ifdef USE_CPPTRACE
    backtrace_count = cpptrace::safe_generate_raw_trace(backtrace_frames, backtrace_N_MAX, 1);
#else
    backtrace_count = 0;
#endif
    if (backtrace_count == 0) {
        backtrace_count = ::backtrace(backtrace_stack, backtrace_N_MAX);
        if (backtrace_count > 0) {
            backtrace_count -= 1; // Don't report ourself
            for (auto i=0; i<backtrace_count; i++) {
                backtrace_stack[i] = backtrace_stack[i+1];
            }
        }
#ifdef USE_CPPTRACE
        for (auto i=0; i<backtrace_count; i++) {
            backtrace_frames[i] = (cpptrace::frame_ptr)backtrace_stack[i] - 1;
        }
#endif
    }
    if (backtrace_count>0) {
        backtrace_ready_cv.notify_one();

        // Wait for the processing thread to finish processing
        backgrace_done_cv.wait(lock, [] { return backtrace_count==0; });
    }

    std::_Exit(EXIT_FAILURE);
}

void mdlbt::show_backtrace_with_snippets(int max_frames,int skip_frames) {
#if defined(USE_CPPTRACE)
    cpptrace::generate_trace(skip_frames+1,max_frames).print_with_snippets();
#else
    show_backtrace(max_frames,skip_frames);
#endif
}
void mdlbt::show_backtrace(int max_frames,int skip_frames) {
#ifdef USE_BT
    ++skip_frames; // Don't report ourself
#if defined(USE_CPPTRACE)
    cpptrace::generate_trace(skip_frames,max_frames).print();
#else
    void *stack[512];
    int stack_size = ::backtrace(stack, sizeof stack / sizeof *stack);
    if (stack_size > max_frames) stack_size = max_frames;
    print_backtrace(stack+skip_frames, stack_size);
}

void mdlbt::print_backtrace(void *stack[], int stack_size) {
    std::cerr << "Stacktrace of " << stack_size << " frames:\n";
#ifdef USE_ELFUTILS
    DebugInfoSession dis;
    for (int i = skip_frames; i < stack_size; ++i) {
        std::cerr << i << ": " << DebugInfo(dis, stack[i]) << '\n';
    }
#else
    auto functions = backtrace_symbols(stack, stack_size);
    for (int i=0; i < stack_size; i++) {
        std::cerr << i << ":" << functions[i] << '\n';
    }
    free(functions);
#endif
#endif
#endif
    std::cerr.flush();
}