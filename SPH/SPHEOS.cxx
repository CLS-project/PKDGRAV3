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

#include "SPHEOS.h"

float SPHEOSPCofRhoU(PKD pkd, float rho, float u, float *c, int iMat, SPHOptions *SPHoptions) {
    if (SPHoptions->useIsentropic) {
        u = u / (SPHoptions->gamma - 1.0f) * pow(rho, SPHoptions->gamma - 1.0f);
    }
    *c = sqrtf(SPHoptions->gamma * (SPHoptions->gamma - 1.0f) * u);
    return (SPHoptions->gamma - 1.0f) * rho * u;
}

float SPHEOSUofRhoT(PKD pkd, float rho, float T, int iMat, SPHOptions *SPHoptions) {
    float u = T * SPHoptions->TuFac;
    if (SPHoptions->useIsentropic) {
        u = u * (SPHoptions->gamma - 1.0f) / pow(rho,SPHoptions->gamma - 1.0f);
    }
    return u;
}

float SPHEOSTofRhoU(PKD pkd, float rho, float u, int iMat, SPHOptions *SPHoptions) {
    if (SPHoptions->useIsentropic) {
        u = u / (SPHoptions->gamma - 1.0f) * pow(rho,SPHoptions->gamma - 1.0f);
    }
    float T = u / SPHoptions->TuFac;
    return T;
}