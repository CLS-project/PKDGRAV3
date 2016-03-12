#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <execinfo.h>
#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h> /* for MAXHOSTNAMELEN, if available */
#endif

#include "bt.h"

static void print_traceback(FILE *fp) {
    void *stack[20];
    char **functions;
    int count, i;

    count = backtrace(stack, 20);
    functions = backtrace_symbols(stack, count);
    for (i=0; i < count; i++) {
        fprintf(fp,"Frame %2d: %s\n", i, functions[i]);
    }
    fflush(fp);
    free(functions);
}

#define USE_SIGACTION

#ifdef USE_SIGACTION
static void alarm_handler(int signo, siginfo_t *si, void *unused) {
#else
static void alarm_handler(int signo) {
#endif
    FILE *fp;
#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
    char hostname[MAXHOSTNAMELEN+12];
    strcpy(hostname,"dump/");
    if (gethostname(hostname+5,MAXHOSTNAMELEN)) fp = stderr;
    else {
	strcpy(hostname+strlen(hostname),".XXXXXX");
	mkstemp(hostname);
	fp = fopen(hostname,"a");
	if (fp==NULL) fp = stderr;
	}
#else
    fp = stderr;
#endif
    time_t t;
    time(&t);
    fprintf(fp,"Traceback on : %s",ctime(&t));

    print_traceback(fp);
    if (fp != stderr) fclose(fp);
    }

#ifdef USE_SIGACTION
static void signal_handler(int signo, siginfo_t *si, void *unused) {
#else
static void signal_handler(int signo) {
#endif
    /* Shouldn't use printf . . . oh well*/
#ifdef USE_SIGACTION
    printf("Caught signal %d at address %p\n",signo,si->si_addr);
#else
    fprintf(stderr,"Caught signal %d\n", signo);
#endif
    print_traceback(stderr);
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
    sa.sa_sigaction = alarm_handler;
    if (sigaction(SIGWINCH, &sa, NULL) == -1) abort();
    }
#else
void bt_initialize(void) {
    signal(SIGBUS, signal_handler);
    signal(SIGILL, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGABRT, signal_handler);
    signal(SIGFPE, signal_handler);
    signal(SIGWINCH, alarm_handler);
}
#endif
