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

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
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

#define STACK_DEPTH 50

static void print_traceback(FILE *fp) {
    void *stack[STACK_DEPTH];
    char **functions;
    int count, i;

    count = backtrace(stack, STACK_DEPTH);
    functions = backtrace_symbols(stack, count);
    for (i=0; i < count; i++) {
        fprintf(fp,"Frame %2d: %s\n", i, functions[i]);
    }
    fflush(fp);
    free(functions);
}

#define USE_SIGACTION

#ifdef USE_ALARM_SIGNAL
#ifdef USE_SIGACTION
static void alarm_handler(int signo, siginfo_t *si, void *unused) {
#else
static void alarm_handler(int signo) {
#endif
    FILE *fp;
    int fd;
#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
    char hostname[MAXHOSTNAMELEN+12];
    strcpy(hostname,"dump/");
    if (gethostname(hostname+5,MAXHOSTNAMELEN)) fp = stderr;
    else {
	strcpy(hostname+strlen(hostname),".XXXXXX");
	fd = mkstemp(hostname);
	if (fd>0) {
	    fp = fdopen(fd,"a");
	    if (fp==NULL) fp = stderr;
	    }
	else fp = stderr;
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
#endif

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
#ifdef USE_ALARM_SIGNAL
    sa.sa_sigaction = alarm_handler;
    if (sigaction(SIGWINCH, &sa, NULL) == -1) abort();
#endif
    }
#else
void bt_initialize(void) {
    signal(SIGBUS, signal_handler);
    signal(SIGILL, signal_handler);
    signal(SIGSEGV, signal_handler);
    signal(SIGABRT, signal_handler);
    signal(SIGFPE, signal_handler);
    //signal(SIGWINCH, alarm_handler);
}
#endif
