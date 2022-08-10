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

float SPHEOSPCTofRhoU(PKD pkd, float rho, float u, float *c, float *T, int iMat, SPHOptions *SPHoptions) {
    float P = 0.0f;
    if (iMat == 0 && SPHoptions->useBuiltinIdeal) {
        if (SPHoptions->useIsentropic) {
            u = u / (SPHoptions->gamma - 1.0f) * pow(rho, SPHoptions->gamma - 1.0f);
        }
        if (T) *T = u / SPHoptions->TuFac;
        *c = sqrtf(SPHoptions->gamma * (SPHoptions->gamma - 1.0f) * u);
        P = (SPHoptions->gamma - 1.0f) * rho * u;
    }
    else {
#ifdef HAVE_EOSLIB_H
        double ctmp = 0.0;
        double Ttmp = 0.0;
        P = (float)EOSPCTofRhoU(pkd->materials[iMat],rho,u,&ctmp,&Ttmp);
        *c = (float)ctmp;
        if (T) *T = (float)Ttmp;
        if (P < 0.0f) P = 0.0f;
        if (*c < pkd->materials[iMat]->minSoundSpeed) *c = (float)pkd->materials[iMat]->minSoundSpeed;
#endif
    }
    return P;
}

float SPHEOSUofRhoT(PKD pkd, float rho, float T, int iMat, SPHOptions *SPHoptions) {
    float u = 0.0f;
    if (iMat == 0 && SPHoptions->useBuiltinIdeal) {
        u = T * SPHoptions->TuFac;
        if (SPHoptions->useIsentropic) {
            u = u * (SPHoptions->gamma - 1.0f) / pow(rho,SPHoptions->gamma - 1.0f);
        }
    }
    else {
#ifdef HAVE_EOSLIB_H
        u = (float)EOSUofRhoT(pkd->materials[iMat],rho,T);
#endif
    }
    return u;
}

float SPHEOSTofRhoU(PKD pkd, float rho, float u, int iMat, SPHOptions *SPHoptions) {
    float T = 0.0f;
    if (iMat == 0 && SPHoptions->useBuiltinIdeal) {
        if (SPHoptions->useIsentropic) {
            u = u / (SPHoptions->gamma - 1.0f) * pow(rho,SPHoptions->gamma - 1.0f);
        }
        T = u / SPHoptions->TuFac;
    }
    else {
#ifdef HAVE_EOSLIB_H
        T = (float)EOSTofRhoU(pkd->materials[iMat],rho,u);
#endif
    }
    return T;
}

float SPHEOSPofRhoT(PKD pkd, float rho, float T, int iMat, SPHOptions *SPHoptions) {
    float P = 0.0f;
    if (iMat == 0 && SPHoptions->useBuiltinIdeal) {
        float u = T * SPHoptions->TuFac;
        if (SPHoptions->useIsentropic) {
            u = u * (SPHoptions->gamma - 1.0f) / pow(rho,SPHoptions->gamma - 1.0f);
        }
        P = (SPHoptions->gamma - 1.0f) * rho * u;
    }
    else {
#ifdef HAVE_EOSLIB_H
        P = (float)EOSPofRhoT(pkd->materials[iMat],rho,T);
#endif
    }
    return P;
}

float SPHEOSIsentropic(PKD pkd, float rho1, float u1, float rho2, int iMat, SPHOptions *SPHoptions) {
    float u2 = 0.0f;
    if (iMat == 0 && SPHoptions->useBuiltinIdeal) {
        u2 = u1;
    }
    else {
#ifdef HAVE_EOSLIB_H
        u2 = (float)EOSIsentropic(pkd->materials[iMat], rho1, u1, rho2);
#endif
    }
    return u2;
}
