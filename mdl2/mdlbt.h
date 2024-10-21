#ifndef MDLBT_H
#define MDLBT_H
#include <thread>
#include <csignal>
class mdlbt {
private:
    std::thread backtrace_thread;
    void *BacktraceThread();
    static void terminate_handler();
    static void print_backtrace(void *stack[], int stack_size);

    static void signal_handler(int signo, siginfo_t *si, void *unused);
    static void signal_sigbus(int signo, siginfo_t *si, void *unused);
public:
    virtual ~mdlbt();
    static void show_backtrace(int max_frames=63,int skip_frames=0);
    static void show_backtrace_with_snippets(int max_frames=63,int skip_frames=0);
    void register_backtrace();
    void ignore_SIGBUS();
};
#endif
