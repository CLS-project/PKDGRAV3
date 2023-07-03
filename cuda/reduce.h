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

#ifndef CUDA_REDUCE_H
#define CUDA_REDUCE_H

#include "sm_30_intrinsics.h"

template <typename T,unsigned int blockSize>
inline __device__ T warpReduceAdd(/*volatile T * data, int tid,*/ T t) {
#if CUDART_VERSION >= 9000
    if (blockSize >= 32) t += __shfl_down_sync(0xffffffff,t,16);
    if (blockSize >= 16) t += __shfl_down_sync(0xffffffff,t,8);
    if (blockSize >= 8)  t += __shfl_down_sync(0xffffffff,t,4);
    if (blockSize >= 4)  t += __shfl_down_sync(0xffffffff,t,2);
    if (blockSize >= 2)  t += __shfl_down_sync(0xffffffff,t,1);
#else
    if (blockSize >= 32) t += __shfl_xor(t,16);
    if (blockSize >= 16) t += __shfl_xor(t,8);
    if (blockSize >= 8)  t += __shfl_xor(t,4);
    if (blockSize >= 4)  t += __shfl_xor(t,2);
    if (blockSize >= 2)  t += __shfl_xor(t,1);
#endif
    return t;
}

template <typename T,unsigned int blockSize>
inline __device__ T warpReduceMin(/*volatile T * data, int tid,*/ T t) {
#if CUDART_VERSION >= 9000
    if (blockSize >= 32) t = min(t, __shfl_down_sync(0xffffffff,t,16));
    if (blockSize >= 16) t = min(t, __shfl_down_sync(0xffffffff,t,8));
    if (blockSize >= 8)  t = min(t, __shfl_down_sync(0xffffffff,t,4));
    if (blockSize >= 4)  t = min(t, __shfl_down_sync(0xffffffff,t,2));
    if (blockSize >= 2)  t = min(t, __shfl_down_sync(0xffffffff,t,1));
#else
    if (blockSize >= 32) t = min(t, __shfl_xor(t,16));
    if (blockSize >= 16) t = min(t, __shfl_xor(t,8));
    if (blockSize >= 8)  t = min(t, __shfl_xor(t,4));
    if (blockSize >= 4)  t = min(t, __shfl_xor(t,2));
    if (blockSize >= 2)  t = min(t, __shfl_xor(t,1));
#endif
    return t;
}

template <typename T,unsigned int blockSize>
inline __device__ T warpReduceMax(/*volatile T * data, int tid,*/ T t) {
#if CUDART_VERSION >= 9000
    if (blockSize >= 32) t = max(t, __shfl_down_sync(0xffffffff,t,16));
    if (blockSize >= 16) t = max(t, __shfl_down_sync(0xffffffff,t,8));
    if (blockSize >= 8)  t = max(t, __shfl_down_sync(0xffffffff,t,4));
    if (blockSize >= 4)  t = max(t, __shfl_down_sync(0xffffffff,t,2));
    if (blockSize >= 2)  t = max(t, __shfl_down_sync(0xffffffff,t,1));
#else
    if (blockSize >= 32) t = max(t, __shfl_xor(t,16));
    if (blockSize >= 16) t = max(t, __shfl_xor(t,8));
    if (blockSize >= 8)  t = max(t, __shfl_xor(t,4));
    if (blockSize >= 4)  t = max(t, __shfl_xor(t,2));
    if (blockSize >= 2)  t = max(t, __shfl_xor(t,1));
#endif
    return t;
}

/*
** Overload atomicMin and atomicMax for floats
** Courtesy of Xiaojing An and timothygiraffe
** https://stackoverflow.com/questions/17399119/how-do-i-use-atomicmax-on-floating-point-values-in-cuda/51549250#51549250
** Xiaojing An provided the initial implementation
** timothygiraffe modified it to ensure it works for any values of the inputs
*/
__device__ __forceinline__ float atomicMin(float *addr, float value) {
    float old;
    old = !signbit(value) ? __int_as_float(atomicMin((int *)addr, __float_as_int(value))) :
          __uint_as_float(atomicMax((unsigned int *)addr, __float_as_uint(value)));
    return old;
}

__device__ __forceinline__ float atomicMax(float *addr, float value) {
    float old;
    old = !signbit(value) ? __int_as_float(atomicMax((int *)addr, __float_as_int(value))) :
          __uint_as_float(atomicMin((unsigned int *)addr, __float_as_uint(value)));
    return old;
}

template <typename T,unsigned int blockSize>
__device__ void warpReduceAndStoreAtomicAdd(int tid,T t,T *result) {
    t = warpReduceAdd<T,blockSize>(t);
    if (tid==0) atomicAdd(result,t);
}

template <typename T,unsigned int blockSize>
__device__ void warpReduceAndStoreAtomicMin(int tid,T t,T *result) {
    t = warpReduceMin<T,blockSize>(t);
    if (tid==0) atomicMin(result,t);
}

template <typename T,unsigned int blockSize>
__device__ void warpReduceAndStoreAtomicMax(int tid,T t,T *result) {
    t = warpReduceMax<T,blockSize>(t);
    if (tid==0) atomicMax(result,t);
}
#endif