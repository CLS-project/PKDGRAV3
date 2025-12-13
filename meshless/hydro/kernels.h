#ifndef KERNELS_H
#define KERNELS_H

#ifdef __CUDACC__
    #define PP_CUDA_BOTH __host__ __device__
#else
    #define PP_CUDA_BOTH
#endif

// Cubic spline kernel as defined in Dehnen & Aly (2012)
//
// (1 − q)^3 − 4(1/2 − q)^3  for   0 < q <= 1/2
//                (1 − q)^3  for 1/2 < q <= 1
//                        0  for       q > 1
//
// where q = r/H.
template <typename dtype=dvec, typename mtype=dmask>
PP_CUDA_BOTH inline dtype cubicSplineKernel(dtype r, dtype H) {
    const dtype dNormFac = 16 * M_1_PI;
    const dtype dHInv = dtype(1.0) / H;
    const dtype q = r * dHInv;
    const dtype dNorm = dNormFac * dHInv * dHInv * dHInv;
    const dtype oneMq = 1.0 - q;
    dtype out = 0.0;
    dtype t;
    t = dNorm * oneMq * oneMq * oneMq;
    mtype m = q<1.0;
    out = mask_mov(out, m, t);
    t = dNorm * (0.5 - 3 * q * q * oneMq);
    m = q<0.5;
    out = mask_mov(out, m, t);
    return out;
}

#if defined(SIMD_H) & !defined(__CUDACC__)
template <>
PP_CUDA_BOTH inline double cubicSplineKernel(double r, double H) {
    const double dNormFac = 16 * M_1_PI;
    const double dHInv = 1.0 / H;
    const double q = r * dHInv;
    const double dNorm = dNormFac * dHInv * dHInv * dHInv;

    if (q < 0.5) {
        return dNorm * (0.5 - 3 * q * q * (1.0 - q));
    }
    else if (q < 1.0) {
        return dNorm * (1.0 - q) * (1.0 - q) * (1.0 - q);
    }
    else {
        return 0.0;
    }
}
#endif

// First derivative (i.e. dW / dq) of the cubic spline kernel
template <typename dtype=double, typename mtype=bool>
PP_CUDA_BOTH inline dtype cubicSplineKernelDeriv(dtype r, dtype H) {
    // const dtype dHInv = 1.0 / H;
    // const dtype q = r * dHInv;
    // const dtype dNorm = 48. * M_1_PI * dHInv * dHInv * dHInv;

    // if (q < 0.5) {
    //     return dNorm * q * (3.*q - 2.);
    // }
    // else if (q < 1.0) {
    //     return - dNorm * (1.-q) * (1.-q);
    // }
    // else {
    //     return 0.0;
    // }
    return 0.0;
}

#if defined(SIMD_H) & !defined(__CUDACC__)
template <>
PP_CUDA_BOTH inline dvec cubicSplineKernelDeriv(dvec r, dvec H) {
    // const dvec dHInv = dvec(1.0) / H;
    // const dvec q = r * dHInv;
    // const dvec dNorm = 48. * M_1_PI * dHInv * dHInv * dHInv;
    dvec out{0.}, t;

    // t = - dNorm * (1.-q) * (1.-q);
    // out = mask_mov(out, q<1.0, t);
    // t = dNorm * q * (3.*q - 2.);
    // out = mask_mov(out, q<0.5, t);
    return out;
}
#endif
#endif
