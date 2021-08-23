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
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "opa_queue.h"

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
#include <memory>
#define MAX_NODE_NAME_LENGTH      256

namespace mdl {
class BasicService {
    friend class mdlBASE;
private:
    int nInBytes;
    int nOutBytes;
    int service_id;
    std::string service_name;
public:
    explicit BasicService(int service_id, int nInBytes, int nOutBytes, const char *service_name="")
	: nInBytes(nInBytes), nOutBytes(nOutBytes), service_id(service_id),service_name(service_name) {}
    explicit BasicService(int service_id, int nInBytes, const char *service_name="")
	: nInBytes(nInBytes), nOutBytes(0), service_id(service_id),service_name(service_name) {}
    explicit BasicService(int service_id, const char *service_name="")
	: nInBytes(0), nOutBytes(0), service_id(service_id),service_name(service_name) {}
    virtual ~BasicService() = default;
    int getServiceID()  {return service_id;}
    int getMaxBytesIn() {return nInBytes;}
    int getMaxBytesOut(){return nOutBytes;}
protected:
    virtual int operator()(int nIn, void *pIn, void *pOut) = 0;
    };

class mdlBASE : protected mdlbt {
public:
    int32_t nThreads; /* Global number of threads (total) */
    int32_t idSelf;   /* Global index of this thread */
    int32_t nProcs;   /* Number of global processes (e.g., MPI ranks) */
    int32_t iProc;    /* Index of this process (MPI rank) */
    int16_t nCores;   /* Number of threads in this process */
    int16_t iCore;    /* Local core id */
    int bDiag;        /* When true, debug output is enabled */
    int argc;

    FILE *fpDiag;
    char **argv;

    /* Services information */
    int nMaxInBytes;
    int nMaxOutBytes;
    std::vector< std::unique_ptr<BasicService> > services;

    /* Maps a give process (Proc) to the first global thread ID */
    std::vector<int> iProcToThread; /* [0,nProcs] (note inclusive extra element) */

    char nodeName[MAX_NODE_NAME_LENGTH];

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
protected:
    ticks nTicks;
#endif
    double dWaiting;
    double dComputing;
    double dSynchronizing;
private:
    double TimeFraction() const;
public:
    void TimeReset();
    void TimeAddComputing();
    void TimeAddSynchronizing();
    void TimeAddWaiting();
    double TimeComputing() const;
    double TimeSynchronizing() const;
    double TimeWaiting() const;

    void mdl_vprintf(const char *format, va_list ap);
    void mdl_printf(const char *format, ...);

public:
    explicit mdlBASE(int argc,char **argv);
    virtual ~mdlBASE();
    int32_t Threads() const { return nThreads; }
    int32_t Self()    const { return idSelf; }
    int16_t Core()    const { return iCore; }
    int16_t Cores()   const { return nCores; }
    int32_t Proc()    const { return iProc; }
    int32_t Procs()   const { return nProcs; }
    int32_t ProcToThread(int32_t iProc) const;
    int32_t ThreadToProc(int32_t iThread) const;
    void yield();
    void AddService(int sid, void *p1, fcnService_t *fcnService, int nInBytes, int nOutBytes, const char *name="");
    void AddService(std::unique_ptr<BasicService> && service);
    int  RunService(int sid, int nIn, void *pIn, void *pOut=nullptr);
    int  RunService(int sid, void *pOut) { return RunService(sid,0,nullptr,pOut); }
    };
int mdlBaseProcToThread(mdlBASE *base, int iProc);
int mdlBaseThreadToProc(mdlBASE *base, int iThread);
} // namespace mdl

#endif

#ifdef __cplusplus
extern "C" {
#endif


int mdlThreads(void *mdl);
int mdlSelf(void *mdl);
int mdlCore(void *mdl);
int mdlCores(void *mdl);
int mdlProc(void *mdl);
int mdlProcs(void *mdl);
const char *mdlName(void *mdl);
int mdlGetArgc(void *mdl);
char **mdlGetArgv(void *mdl);


void mdlDiag(void *mdl, char *psz);
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

double mdlCpuTimer(void * mdl);
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
void mdlZeroTimer(void * mdl, mdlTimer *);
void mdlGetTimer(void * mdl, mdlTimer *, mdlTimer *);
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
