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

    /*
     * IA: Parameters for the meshless hydrodynamics
     */
    double dCFLacc;
    int bMeshlessHydro;
    int bIterativeSmoothingLength;
    int bWakeUpParticles;
    double dNeighborsStd;
    int bOutFineStatistics;

    /*
     * IA: Parameter for fixing all particles to the same rung
     */
    int bGlobalDt;

#if defined(COOLING) || defined(GRACKLE)
    char achCoolingTables[256];
#endif
#ifdef COOLING
    /*
     * IA: Cooling parameters
     */
    double fH_reion_z;
    double fH_reion_eV_p_H;
    double fHe_reion_z_centre;
    double fHe_reion_z_sigma;
    double fHe_reion_eV_p_H;
    double fCa_over_Si_in_Solar;
    double fS_over_Si_in_Solar;
    double fT_CMB_0;
    double dCoolingMinu;
    double dCoolingMinTemp;
#endif

    /*
     * Internal energy floor parameters
     */
    double dEOSFloorDen;
    double dEOSFloorMinOD;
    double dEOSFloorMinBaryonOD;
    double dEOSFloornH;
    double dEOSFlooru;
    double dEOSFloorTemp;
#ifdef EEOS_POLYTROPE
    double dEOSPolyFloorMinOD;
    double dEOSPolyFloorMinBaryonOD;
    double dEOSPolyFloorIndex;
    double dEOSPolyFloorExponent;
    double dEOSPolyFloorDen;
    double dEOSPolyFloornH;
    double dEOSPolyFlooru;
    double dEOSPolyFloorTemp;
#endif
#ifdef EEOS_JEANS
    double dEOSNJeans;
#endif
    /*
     * IA: Initial abundances
     */
    double dInitialH;
#ifdef HAVE_HELIUM
    double dInitialHe;
#endif
#ifdef HAVE_CARBON
    double dInitialC;
#endif
#ifdef HAVE_NITROGEN
    double dInitialN;
#endif
#ifdef HAVE_OXYGEN
    double dInitialO;
#endif
#ifdef HAVE_NEON
    double dInitialNe;
#endif
#ifdef HAVE_MAGNESIUM
    double dInitialMg;
#endif
#ifdef HAVE_SILICON
    double dInitialSi;
#endif
#ifdef HAVE_IRON
    double dInitialFe;
#endif
#ifdef HAVE_METALLICITY
    double dInitialMetallicity;
#endif
#ifdef STAR_FORMATION
    /* IA: Star formation */
#ifdef HAVE_METALLICITY
    int bSFThresholdDenSchaye2004;
#endif
    double dSFMinOverDensity;
    double dSFGasFraction;
    double dSFThresholdDen;
    double dSFThresholdOD;
    double dSFThresholdu;
    double dSFThresholdT;
    double dSFindexKS;
    double dSFnormalizationKS;
    double dSFEfficiency;
#endif
#ifdef FEEDBACK
    int bCCSNFeedback;
    int bSNIaFeedback;
    double dSNFBEfficiency;
    double dSNFBMaxEff;
    double dSNFBEffIndex;
    double dSNFBEffnH0;
    double dSNFBDT;
    double dSNFBDu;
    double dCCSNFBDelay;
    double dCCSNFBSpecEnergy;
    double dCCSNEnergy;
    double dSNIaFBDelay;
    double dSNIaFBSpecEnergy;
    double dSNIaEnergy;
#endif
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
