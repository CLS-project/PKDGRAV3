#ifndef MDLBT_H
#define MDLBT_H
#include <csignal>
class mdlbt {
private:
    static void terminate_handler();

    static void signal_handler(int signo, siginfo_t *si, void *unused);
    static void signal_sigbus(int signo, siginfo_t *si, void *unused);
public:
    static void show_backtrace(int max_frames=63,int skip_frames=0);
    void register_backtrace();
    void ignore_SIGBUS();
};
#endif
