#include "eEOS/eEOS_struct.h"
#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)

inline static void eEOSFill(const struct parameters param, struct eEOSparam *eEOS) {
#ifdef EEOS_POLYTROPE
    eEOS->dPolyFloorMinOD = param.dEOSPolyFloorMinBaryonOD;
    eEOS->dPolyFloorExponent = param.dEOSPolyFloorExponent;
    eEOS->dPolyFloorDen = param.dEOSPolyFloorDen;
    eEOS->dPolyFlooru = param.dEOSPolyFlooru;
#endif
#ifdef EEOS_JEANS
    eEOS->dNJeans = param.dEOSNJeans;
#endif
    eEOS->dFlooru = param.dEOSFlooru;
    eEOS->dFloorDen = param.dEOSFloorDen;
    eEOS->dFloorMinOD = param.dEOSFloorMinBaryonOD;
}


/*
 * General interface to the eEOS module
 * The density is in comoving coordinates.
 */
inline static double eEOSEnergyFloor(double a_inv3, float fDens, float fBall,
                                     double dConstGamma, struct eEOSparam eEOS);
inline static double eEOSPressureFloor(double a_inv3, float fDens, float fBall,
                                       double dConstGamma, struct eEOSparam eEOS);



/*
 * Simple constant temperature floor
 */
inline static double constantEnergyFloor(double a_inv3, float fDens,
        double dConstGamma, struct eEOSparam eEOS) {
    if ( fDens*a_inv3 > eEOS.dFloorDen && fDens > eEOS.dFloorMinOD )
        return eEOS.dFlooru;
    else
        return NOT_IN_EEOS;
}

inline static double constantPressureFloor(double a_inv3, float fDens,
        double dConstGamma, struct eEOSparam eEOS) {
    return fDens * (dConstGamma - 1. ) *
           constantEnergyFloor(a_inv3, fDens, dConstGamma, eEOS);
}


#ifdef EEOS_POLYTROPE
/*
 * General polytropic eEOS floor.
 * It is composed of a polytrope of arbitrary index and a constant temperature floor
 */
inline static double polytropicEnergyFloor(double a_inv3, float fDens, double dConstGamma,
        struct eEOSparam eEOS) {
    if (fDens > eEOS.dPolyFloorMinOD)
        return eEOS.dPolyFlooru * pow( fDens*a_inv3/eEOS.dPolyFloorDen, eEOS.dPolyFloorExponent );
    else
        return NOT_IN_EEOS;
}

inline static double polytropicPressureFloor(double a_inv3, float fDens, double dConstGamma,
        struct eEOSparam eEOS) {
    return fDens * (dConstGamma - 1. ) * polytropicEnergyFloor(a_inv3, fDens, dConstGamma, eEOS);
}


inline static double eEOSEnergyFloor(double a_inv3, float fDens, float fBall,
                                     double dConstGamma, struct eEOSparam eEOS) {
    const double polyEnergy = polytropicEnergyFloor(a_inv3, fDens, dConstGamma, eEOS);
    const double constantEnergy = constantEnergyFloor(a_inv3, fDens, dConstGamma, eEOS);

    if (constantEnergy != NOT_IN_EEOS)
        return (polyEnergy > constantEnergy) ? polyEnergy : constantEnergy;
    else
        return constantEnergy;
}
inline static double eEOSPressureFloor(double a_inv3, float fDens, float fBall,
                                       double dConstGamma, struct eEOSparam eEOS) {
    const double polyPressure = polytropicPressureFloor(a_inv3, fDens, dConstGamma, eEOS);
    const double constantPressure = constantPressureFloor(a_inv3, fDens, dConstGamma, eEOS);

    if (constantPressure != NOT_IN_EEOS)
        return (polyPressure > constantPressure) ? polyPressure : constantPressure;
    else
        return constantPressure;
}

#endif




#ifdef EEOS_JEANS
/*
 * Pressure floor based on the Jeans criterion, where the only parameter is
 *  the minimum "number" of resolution elements per Jeans length
 *
 * It is similar to polytropicPressureFloor with gamma=4/3, but it takes into
 *  account explicitly the smoothing length
 */
inline static double eEOSPressureFloor(double a_inv3, float fDens, float fBall,
                                       double dConstGamma, struct eEOSparam eEOS) {
    return 1.2/dConstGamma * pow(eEOS.dNJeans,2./3.) * fDens * fDens * fBall * fBall;

}

inline static double eEOSEnergyFloor(double a_inv3, float fDens, float fBall,
                                     double dConstGamma, struct eEOSparam eEOS) {
    return eEOSPressureFloor(a_inv3, fDens, fBall, eEOS) /
           (fDens * (dConstGamma - 1.));

}
#endif
#endif

