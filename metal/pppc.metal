#include <metal_compute>
#include <metal_math>
#include <metal_simdgroup>
#include <metal_atomic>
using namespace metal::fast;
using metal::atomic;
using metal::isnan;

inline float maskz_mov(bool p,float a) { return metal::select(0.0f,a,p); /*p ? a : 0.0f;*/ }
inline bool testz(bool p) { return !p; }

#include "gravity/pp.h"
#include "gravity/pc.h"
#include "gpu/workunit.h"

// Perform per-SIMD partial reduction.
inline float shuffle_down(float val,uint simd_size) {
    for (auto offset=simd_size/2; offset>0; offset/=2)
        val += metal::simd_shuffle_down(val, offset);
    return val;
}
template<typename RESULT>
inline void shuffle_down(thread RESULT &result,uint simd_size) {
    result.ax = shuffle_down(result.ax,simd_size);
    result.ay = shuffle_down(result.ay,simd_size);
    result.az = shuffle_down(result.az,simd_size);
    result.pot = shuffle_down(result.pot,simd_size);
    result.ir = shuffle_down(result.ir,simd_size);
    result.norm = shuffle_down(result.norm,simd_size);
}

// It would be nice if this worked, but the barrier is tricky. This would only apply to low power devices.
// template<typename T>
// inline void reduce(thread T &val,threadgroup T *ldata,uint lid,uint lsize,uint simd_size,uint simd_lane_id,uint simd_group_id) {
//     while (lsize>simd_size) {
//         lsize /= simd_size;
//         shuffle_down(val, simd_size);
//         if (simd_lane_id == 0) ldata[simd_group_id] = val;
//         // threadgroup_barrier(metal::mem_flags::mem_threadgroup);
//         if (lid<lsize) val = ldata[lid];
//     }
//     shuffle_down(val, simd_size);
// }

inline int punned_add(int a,float v) {
    float f = * reinterpret_cast<thread float *>(&a) + v;
    return * reinterpret_cast<thread int *>(&f);
}

inline void atomic_add(volatile device float &out,float val) {
    auto data = reinterpret_cast<volatile device atomic<int> *>(&out);
    // atomic_fetch_add_explicit(data,val,metal::memory_order_relaxed);
    int expected = out;
    while (!atomic_compare_exchange_weak_explicit(data,(thread int *)&expected,punned_add(expected,val),
            metal::memory_order_relaxed,metal::memory_order_relaxed)) {}
}


inline auto evalInteraction(const gpu::ppInput pp,const threadgroup gpu::ppInteract *__restrict__ blk,int i) {
    return EvalPP<float,bool>(pp.dx, pp.dy, pp.dz, pp.fSoft2,
                              blk->dx[i], blk->dy[i], blk->dz[i], blk->fourh2[i], blk->m[i],
                              pp.ax, pp.ay, pp.az, pp.dImaga);
}

inline auto evalInteraction(const gpu::ppInput pp,const threadgroup gpu::pcInteract *__restrict__ blk,int i) {
    return EvalPC<float,bool,true>(pp.dx, pp.dy, pp.dz, pp.fSoft2,
                                   blk->dx[i],blk->dy[i],blk->dz[i],blk->m[i],blk->u[i],
                                   blk->xxxx[i],blk->xxxy[i],blk->xxxz[i],blk->xxyz[i],blk->xxyy[i],
                                   blk->yyyz[i],blk->xyyz[i],blk->xyyy[i],blk->yyyy[i],
                                   blk->xxx[i],blk->xyy[i],blk->xxy[i],blk->yyy[i],blk->xxz[i],blk->yyz[i],blk->xyz[i],
                                   blk->xx[i],blk->xy[i],blk->xz[i],blk->yy[i],blk->yz[i],
#ifdef USE_DIAPOLE
                                   blk->x[i],blk->y[i],blk->z[i],
#endif
                                   pp.ax, pp.ay, pp.az, pp.dImaga);
}

// For each thread block
// gid: index of thread block
//   .x: index of work unit
//   .y: always 0
//   .z: always 0
// lid: index within the thread block
//   .x: 0-31, index within interaction block
//   .y: index of the particle group batch
//   .z: index of the interaction block (+ block)
template<typename BLK,typename RESULT>
void evaluateKernel(
    const device gpu::ppWorkUnit *units,
    const device gpu::ppInput    *input,
    const device BLK *inter,
    volatile device gpu::ppResult  *output,
    threadgroup BLK *blocks,
    threadgroup RESULT *ldata,
    uint3 gid,uint3 lid,uint3 lsize,
    uint simd_size,
    uint simd_lane_id,
    uint simd_group_id) {
    // The "work unit" and matching interaction block
    auto iI = gid.z;
    auto iP = units[iI].iP; // Index of the first particle
    auto nP = units[iI].nP; // Number of particles to calculate (0-)
    auto nI = units[iI].nI; // Total number of interactions (0-32)
    // Advance to the correct set of interactions
    inter  += iI;       // Global interaction block
    blocks += lid.z;    // Thread local cache
    input += iP;        // First particle to process
    output += iP;       // Where to store the results
    // Calculation for reduction
    // auto nSIMD = lsize.x / simd_size;
    // ldata += nSIMD;
    // simd_group_id %= nSIMD;

    // copy our interactions into shared memory
    {
        auto nBlocks = sizeof(gpu::pcInteract) / sizeof(gpu::fvector<32>);
        auto dst = reinterpret_cast<threadgroup  gpu::fvector<32>*>(blocks);
        auto src = reinterpret_cast<device const gpu::fvector<32>*>(inter);
        for (auto i=lid.y; i<nBlocks; i+=lsize.y) {
            dst[i][lid.x] = src[i][lid.x];
        }
        threadgroup_barrier(metal::mem_flags::mem_threadgroup);
    }

    decltype(evalInteraction(input[0],blocks,0)) result {0,0,0,0,0,0};
    for (uint i=lid.y; i<nP; i+=lsize.y) {
        if (lid.x<nI) result = evalInteraction(input[i],blocks,lid.x);
        shuffle_down(result,simd_size);
        if (simd_lane_id==0) {
            atomic_add(output[i].ax,result.ax);
            atomic_add(output[i].ay,result.ay);
            atomic_add(output[i].az,result.az);
            atomic_add(output[i].fPot,result.pot);
        }
    }
}

kernel void pp(
    const device gpu::ppWorkUnit *units [[buffer(0)]],
    const device gpu::ppInput    *input [[buffer(1)]],
    const device gpu::ppInteract *inter [[buffer(2)]],
    volatile device gpu::ppResult  *output [[buffer(3)]],
    threadgroup gpu::ppInteract *blocks [[threadgroup(0)]],
    threadgroup ResultPP<float> *ldata [[threadgroup(1)]],
    uint3 gid [[thread_position_in_grid]],
    uint3 lid [[thread_position_in_threadgroup]],
    uint3 lsize [[threads_per_threadgroup]],
    uint simd_size [[threads_per_simdgroup]],
    uint simd_lane_id [[thread_index_in_simdgroup]],
    uint simd_group_id [[simdgroup_index_in_threadgroup]]) {
    evaluateKernel(units,input,inter,output,blocks,ldata,gid,lid,lsize,simd_size,simd_lane_id,simd_group_id);
}

kernel void pc(
    const device gpu::ppWorkUnit *units [[buffer(0)]],
    const device gpu::ppInput    *input [[buffer(1)]],
    const device gpu::pcInteract *inter [[buffer(2)]],
    volatile device gpu::ppResult  *output [[buffer(3)]],
    threadgroup gpu::pcInteract *blocks [[threadgroup(0)]],
    threadgroup ResultPC<float> *ldata [[threadgroup(1)]],
    uint3 gid [[thread_position_in_grid]],
    uint3 lid [[thread_position_in_threadgroup]],
    uint3 lsize [[threads_per_threadgroup]],
    uint simd_size [[threads_per_simdgroup]],
    uint simd_lane_id [[thread_index_in_simdgroup]],
    uint simd_group_id [[simdgroup_index_in_threadgroup]]) {
    evaluateKernel(units,input,inter,output,blocks,ldata,gid,lid,lsize,simd_size,simd_lane_id,simd_group_id);
}
