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
#ifndef MDLSTUBS_H
#define MDLSTUBS_H

typedef void* MDL; // This is actually a pointer to a mdlClass object.

// These stubs are used to provide a C interface to the mdlClass class.
// Do not allow them in C++ code.
#ifndef __cplusplus
int mdlLaunch(int (*)(MDL, void*), void* (*)(MDL), void (*)(MDL, void*), int, char**);

int mdlThreads(void* mdl);
int mdlSelf(void* mdl);
int mdlCore(void* mdl);
int mdlCores(void* mdl);
int mdlProc(void* mdl);
int mdlProcs(void* mdl);

int mdlProcToThread(MDL mdl, int iProc);
int mdlThreadToProc(MDL mdl, int iThread);

void mdlAbort(MDL);
int mdlSwap(MDL, int, size_t, void*, size_t, size_t*, size_t*);
void mdlSend(MDL mdl, int id, mdlPack pack, void* ctx);
void mdlRecv(MDL mdl, int id, mdlPack unpack, void* ctx);
void mdlAddService(MDL, int, void*, fcnService_t* fcnService, int, int);

void* mdlMalloc(MDL, size_t);
void mdlFree(MDL, void*);
void mdlSetCacheSize(MDL, int);
void mdlCacheCheck(MDL);
int mdlCacheStatus(MDL mdl, int cid); /* zero means not open */
void mdlCacheBarrier(MDL, int);
void mdlFlushCache(MDL, int);
void mdlThreadBarrier(MDL);
void mdlCompleteAllWork(MDL);

void* mdlAcquire(MDL mdl, int cid, int iIndex, int id);
void* mdlFetch(MDL mdl, int cid, int iIndex, int id);
void* mdlVirtualFetch(MDL mdl, int cid, int iIndex, int id);
const void* mdlKeyFetch(MDL mdl, int cid, uint32_t uHash, void* pKey, int lock, int modify, int virt);
void* mdlKeyAcquire(MDL mdl, int cid, uint32_t uHash, void* pKey);
void mdlRelease(MDL, int, void*);
void* mdlAcquireWrite(MDL mdl, int cid, int iIndex);
void mdlReleaseWrite(MDL mdl, int cid, void* p);

void mdlFinishCache(MDL, int);

int mdlCudaActive(MDL mdl);

double mdlNumAccess(MDL, int);
double mdlMissRatio(MDL, int);

void mdlROcache(MDL mdl, int cid, void* (*getElt)(void* pData, int i, int iDataSize), void* pData, int iDataSize, int nData);

void mdlCOcache(MDL mdl, int cid, void* (*getElt)(void* pData, int i, int iDataSize), void* pData, int iDataSize, int nData, void* ctx,
                void (*init)(void*, void*), void (*combine)(void*, void*, const void*));

#endif // __cplusplus
#endif // MDLSTUBS_H
