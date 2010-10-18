#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <execinfo.h>

#include "bt.h"

static void signal_handler(int signo) {
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

#ifdef USE_PSD
int pst_idSelf=0;

static void exit_handler() {
    fprintf(stderr, "exit_handler %i\n", pst_idSelf);
    if (0);
    {
        void *stack[20];
        char **functions;
        int count, i;
        fprintf(stderr,"Exiting [id=%d]\n", pst_idSelf);
        count = backtrace(stack, 20);
        functions = backtrace_symbols(stack, count);
        for (i=0; i < count; i++) {
            fprintf(stderr,"Frame %2d: %s\n", i, functions[i]);
        }
        fflush(stderr);
        free(functions);
    }
}
#endif


void bt_initialize(void) {
    signal(SIGBUS, signal_handler);
    signal(SIGILL, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGABRT, signal_handler);
#ifdef USE_PSD
    atexit(exit_handler);
#endif
}

