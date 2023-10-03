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

#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

/*
** Don't even think about putting a pointer in here!!
*/
struct parameters {
    /*
    ** Parameters for PKDGRAV.
    */


#ifdef BLACKHOLES
    int bBHMerger;
    int bBHPlaceSeed;
    int bBHAccretion;
    int bBHFeedback;
    double dBHAccretionAlpha;
    double dBHAccretionCvisc;
    double dBHAccretionEddFac;
    double dBHRadiativeEff;
    double dBHFBEff;
    double dBHFBDT;
    double dBHFBEcrit;
    double dBHSeedMass;
    double dBHMhaloMin;
#endif
#ifdef STELLAR_EVOLUTION
    char achStelEvolPath[256];
    char achSNIaDTDType[32];
    int bChemEnrich;
    int nSmoothEnrich;
    double dWindSpecificEkin;
    double dSNIaNorm;
    double dSNIaScale;
    double dSNIaExpScale;
    double dSNIaPLScale;
    double dSNIaPLInitTime;
    double dSNIaPLFinalTime;
    double dSNIaMaxMass;
    double dStellarWindSpeed;
#endif
#if defined(FEEDBACK) || defined(STELLAR_EVOLUTION)
    double dSNIaNumPerMass;
#endif
#ifdef STELLAR_IMF
    char achIMFType[32];
    double dIMFMinMass;
    double dIMFMaxMass;
    double dCCSNMinMass;
    double dCCSNMaxMass;
#endif
};

#endif
