#include "mdl_config.h"
#include "mdlbt.h"
#include <iostream>
#include <cassert>
#include <unistd.h>
#include <stdexcept>
#include <execinfo.h>
#include <mutex>
#ifdef USE_ELFUTILS
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

static std::mutex backtrace_mutex;

// Here we register the function that will generate a backtrace when an
// abnormal termination occurs.
void mdlbt::register_backtrace() {
    std::set_terminate(terminate_handler);

    struct sigaction sa;
    sa.sa_flags = SA_SIGINFO | SA_RESETHAND;
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = signal_handler;
    if (sigaction(SIGBUS,  &sa, NULL) == -1) abort();
    if (sigaction(SIGILL,  &sa, NULL) == -1) abort();
    if (sigaction(SIGSEGV, &sa, NULL) == -1) abort();
    if (sigaction(SIGABRT, &sa, NULL) == -1) abort();
    if (sigaction(SIGFPE,  &sa, NULL) == -1) abort();

//    signal(SIGBUS, signal_handler);
//    signal(SIGILL, signal_handler);
//    signal(SIGSEGV, signal_handler);
//    signal(SIGABRT, signal_handler);
//    signal(SIGFPE, signal_handler);
}

void mdlbt::ignore_SIGBUS() {
    struct sigaction sa;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_SIGINFO;
    sa.sa_sigaction = signal_sigbus;
    if (sigaction(SIGBUS,  &sa, NULL) == -1) abort();
}

void mdlbt::terminate_handler() {
    std::unique_lock<std::mutex> lck(backtrace_mutex);
    show_backtrace();
    std::_Exit(EXIT_FAILURE);
}

void mdlbt::signal_sigbus(int signo, siginfo_t *si, void *unused) {
    std::unique_lock<std::mutex> lck(backtrace_mutex);
    fflush(stdout);
    fflush(stderr);
    fsync(fileno(stdout));
    fsync(fileno(stderr));
    char nodeName[200];

    if (gethostname(nodeName, sizeof(nodeName)))
        nodeName[0] = 0;
    else
        nodeName[sizeof(nodeName) - 1] = 0;

    std::cerr << "SIGNAL " << signo << "." << si->si_code << " at " << si->si_addr << " on " << nodeName << '\n';
}

void mdlbt::signal_handler(int signo, siginfo_t *si, void *unused) {
    std::unique_lock<std::mutex> lck(backtrace_mutex);
    fflush(stdout);
    fflush(stderr);
    fsync(fileno(stdout));
    fsync(fileno(stderr));
    std::cerr << "SIGNAL " << signo << " at " << si->si_addr << '\n';
    show_backtrace();
    //std::_Exit(EXIT_FAILURE);
}

void mdlbt::show_backtrace(int max_frames,int skip_frames) {
#ifdef USE_BT
    void *stack[512];
    int stack_size = ::backtrace(stack, sizeof stack / sizeof *stack);
    if (stack_size > max_frames) stack_size = max_frames;

    std::cerr << "Stacktrace of " << stack_size-skip_frames << " frames:\n";
#ifdef USE_ELFUTILS
    DebugInfoSession dis;
    for (int i = skip_frames; i < stack_size; ++i) {
        std::cerr << i << ": " << DebugInfo(dis, stack[i]) << '\n';
    }
#else
    auto functions = backtrace_symbols(stack, stack_size);
    for (int i=skip_frames; i < stack_size; i++) {
        std::cerr << i << ":" << functions[i] << '\n';
    }
    free(functions);
#endif
#endif
    std::cerr.flush();
}
