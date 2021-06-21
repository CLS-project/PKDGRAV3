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

#ifndef SPHOPTIONS_HINCLUDED
#define SPHOPTIONS_HINCLUDED

#include <stdint.h>
#include "parameters.h"

typedef struct {
    uint64_t doGravity : 1;
    uint64_t doDensity : 1;
    uint64_t doSPHforces : 1;
    uint64_t useNumDen : 1;
    uint64_t kernelType : 2;
    } SPHOptions;

#ifdef __cplusplus
extern "C"
#endif
SPHOptions initializeSPHOptions(struct parameters param);

#endif

/* Definition of the kernel:
 * W(r) = C  * w(r)
 * r = sqrt(d2)/fball
 * fball is the radius at which the kernel becomes zero, w(r>=1) = 0
 * SPHKERNEL calculates w(r)
 * DSPHKERNEL_DR calculates dw/dr
 * DSPHKERNEL_DFBALL calculates dW/dfball
 */

/* Kernel types: kernelType =
 * 0: M4
 * 1:
 * 2:
 * 3:
 */

/* Initializes the SPH kernel, gives back all masks needed, can calculate temporary variables
 * needed in SPHKERNEL and DSPHKERNEL_DR and calculates the kernel normalization C
 */
#define SPHKERNEL_INIT(r, ifBall, C, t1, mask1, kernelType) { \
    switch(kernelType) { \
    case 0: { \
        mask1 = r < 0.5f; \
        t1 = r - 1.0f; \
        C = 8.0f * M_1_PI * ifBall * ifBall * ifBall; \
        break; } \
    default: assert(0);\
    }\
    }

/* Evaluates the non-normalized kernel function, using the masks and temporary variables
 * initialized in SPHKERNEL_INIT
 * has to be normalized with C at the end
 */
#define SPHKERNEL(r, w, t1, t2, t3, r_lt_one, mask1, kernelType) { \
    switch(kernelType) { \
    case 0: { \
        t2 = 1.0f + 6.0f * r * r * t1; \
        t3 = - 2.0f * t1 * t1 * t1; \
        w = maskz_mov(r_lt_one,t3); \
        w = mask_mov(w,mask1,t2); \
        break; } \
    default: assert(0);\
    }\
    }

/* Evaluates the derivative of the non-normalized kernel function with respect to r,
 * using the masks and temporary variables initialized in SPHKERNEL_INIT
 * has to be normalized with C at the end
 */
#define DSPHKERNEL_DR(r, dwdr, t1, t2, t3, r_lt_one, mask1, kernelType) { \
    switch(kernelType) { \
    case 0: { \
        t2 = 6.0f * r * (3.0f * r - 2.0f); \
        t3 = - 6.0f * t1 * t1; \
        dwdr = maskz_mov(r_lt_one,t3); \
        dwdr = mask_mov(dwdr,mask1,t2); \
        break; } \
    default: assert(0);\
    }\
    }

/* Calculates the derivative of the normalized kernel function with respect to fball
 * Do not normalize, is already normalized.
 */
#define DSPHKERNEL_DFBALL(r, ifBall, w, dwdr, C, dWdfball, kernelType) { \
    switch(kernelType) { \
    case 0: { \
        dWdfball = - C * ifBall * (3.0f * w + dwdr * r); \
        break; } \
    default: assert(0);\
    }\
    }
