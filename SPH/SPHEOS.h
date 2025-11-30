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

#ifndef SPHEOS_HINCLUDED
#define SPHEOS_HINCLUDED

#include "SPHOptions.h"
#include "../pkd.h"
#include "../core/abort_error.h"

#ifdef __cplusplus
extern "C" {
#endif
float SPHEOSPCTofRhoU(PKD pkd, float rho, float u, float *c, float *T, int iMat, SPHOptions *SPHoptions);
float SPHEOSUofRhoT(PKD pkd, float rho, float T, int iMat, SPHOptions *SPHoptions);
float SPHEOSTofRhoU(PKD pkd, float rho, float u, int iMat, SPHOptions *SPHoptions);
float SPHEOSPofRhoT(PKD pkd, float rho, float T, int iMat, SPHOptions *SPHoptions);
float SPHEOSRhoofPT(PKD pkd, float P, float T, int iMat, SPHOptions *SPHoptions);
float SPHEOSIsentropic(PKD pkd, float rho1, float u1, float rho2, int iMat, SPHOptions *SPHoptions);
float SPHEOSGammaofRhoT(PKD pkd, float rho, float T, int iMat, SPHOptions *SPHoptions);
void SPHEOSApplyStrengthLimiter(PKD pkd, float rho, float u, int iMat, float *Sxx, float *Syy, float *Sxy, float *Sxz, float *Syz, SPHOptions *SPHoptions);
static inline float SPHEOSconvertEntropicFunctiontoInternalEnergy(float rho, float u, int iMat, SPHOptions *SPHoptions) {
    return u / (SPHoptions->gamma - 1.0f) * pow(rho, SPHoptions->gamma - 1.0f);
};
static inline float SPHEOSconvertInternalEnergytoEntropicFunction(float rho, float u, int iMat, SPHOptions *SPHoptions) {
    return u * (SPHoptions->gamma - 1.0f) / pow(rho, SPHoptions->gamma - 1.0f);
};
static inline void checkEOSlibVersion() {
    int result = 1;
#ifdef HAVE_EOSLIB_H
#define EOS_VERSION_MAJOR_REQUIRED 1
#define EOS_VERSION_MINOR_REQUIRED 0
#define EOS_VERSION_PATCH_REQUIRED 2
    if ((EOS_VERSION_MAJOR * 10000 + EOS_VERSION_MINOR * 100 + EOS_VERSION_PATCH) < (EOS_VERSION_MAJOR_REQUIRED * 10000 + EOS_VERSION_MINOR_REQUIRED * 100 + EOS_VERSION_PATCH_REQUIRED)) {
        printf("EOSlib version (%s) is outdated. At least version %d.%d.%d is needed.\n",EOS_VERSION_TEXT,EOS_VERSION_MAJOR_REQUIRED,EOS_VERSION_MINOR_REQUIRED,EOS_VERSION_PATCH_REQUIRED);
        result = 0;
    }
    if (!EOSCheckVersions()) {
        result = 0;
    }
#endif
    abort_assert(1004,result);
}
#ifdef __cplusplus
}
#endif
#endif