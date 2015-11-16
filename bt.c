#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <execinfo.h>

#include "bt.h"

#define USE_SIGACTION
#ifdef USE_SIGACTION
static void signal_handler(int signo, siginfo_t *si, void *unused) {
#else
static void signal_handler(int signo) {
#endif
    void *stack[20];
    char **functions;
    int count, i;

    /* Shouldn't use printf . . . oh well*/
#ifdef USE_SIGACTION
    printf("Caught signal %d at address %p\n",signo,si->si_addr);
#else
    fprintf(stderr,"Caught signal %d\n", signo);
#endif
    count = backtrace(stack, 20);
    functions = backtrace_symbols(stack, count);
    for (i=0; i < count; i++) {
	fprintf(stderr,"Frame %2d: %s\n", i, functions[i]);
    }
    fflush(stderr);
    free(functions);
    exit(signo);
    }

#ifdef USE_SIGACTION
void bt_initialize(void) {
    struct sigaction sa;
    sa.sa_flags = SA_SIGINFO;
    sigemptyset(&sa.sa_mask);
    sa.sa_sigaction = signal_handler;
    if (sigaction(SIGBUS,  &sa, NULL) == -1) abort();
    if (sigaction(SIGILL,  &sa, NULL) == -1) abort();
    if (sigaction(SIGSEGV, &sa, NULL) == -1) abort();
    if (sigaction(SIGABRT, &sa, NULL) == -1) abort();
    if (sigaction(SIGFPE,  &sa, NULL) == -1) abort();
    }
#else
    void *stack[20];
    char **functions;
    int count, i;
    /* Shouldn't use printf . . . oh well*/
    fprintf(stderr,"Caught signal %d\n", signo);
    count = backtrace(stack, 20);
    functions = backtrace_symbols(stack, count);
    for (i=0; i < count; i++) {
	fprintf(stderr,"Frame %2d: %s\n", i, functions[i]);
    }
    fflush(stderr);
    free(functions);
    exit(signo);
}

void bt_initialize(void) {
    signal(SIGBUS, signal_handler);
    signal(SIGILL, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGABRT, signal_handler);
    signal(SIGFPE, signal_handler);
}
#endif
