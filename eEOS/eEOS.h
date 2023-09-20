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
template <typename dtype=vec<double,double>, typename mtype=mmask<bool>>
inline static dtype eEOSEnergyFloor(dtype a_inv3, dtype fDens, dtype fBall,
                                    dtype dConstGamma, struct eEOSparam eEOS);
template <typename dtype=vec<double,double>, typename mtype=mmask<bool>>
inline static dtype eEOSPressureFloor(dtype a_inv3, dtype fDens, dtype fBall,
                                      dtype dConstGamma, struct eEOSparam eEOS);

/*
 * Simple constant temperature floor
 */
template <typename dtype, typename mtype>
inline static dtype constantEnergyFloor(dtype a_inv3, dtype fDens,
                                        dtype dConstGamma, struct eEOSparam eEOS) {
    mtype mask1 = fDens*a_inv3 > eEOS.dFloorDen;
    mtype mask2 = fDens > eEOS.dFloorMinOD;
    mtype mask = mask1 && mask2;
    return mask_mov( static_cast<dtype>(NOT_IN_EEOS), mask, eEOS.dFlooru);
}

template <typename dtype, typename mtype>
inline static dtype constantPressureFloor(dtype a_inv3, dtype fDens,
        dtype dConstGamma, struct eEOSparam eEOS) {
    return fDens * (dConstGamma - 1. ) *
           constantEnergyFloor(a_inv3, fDens, dConstGamma, eEOS);
}

#ifdef EEOS_POLYTROPE
/*
 * General polytropic eEOS floor.
 * It is composed of a polytrope of arbitrary index and a constant temperature floor
 */
template <typename dtype, typename mtype>
inline static dtype polytropicEnergyFloor(dtype a_inv3, dtype fDens, dtype dConstGamma,
        struct eEOSparam eEOS) {
    mtype mask = fDens > eEOS.dPolyFloorMinOD;
    if (movemask(mask)) {
        dtype eeos = eEOS.dPolyFlooru *
                     pow( fDens*a_inv3/static_cast<dtype>(eEOS.dPolyFloorDen),
                          static_cast<dtype>(eEOS.dPolyFloorExponent) );
        return mask_mov( NOT_IN_EEOS, mask, eeos);
    }
    else {
        return NOT_IN_EEOS;
    }
}

template <typename dtype, typename mtype>
inline static dtype polytropicPressureFloor(dtype a_inv3, dtype fDens, dtype dConstGamma,
        struct eEOSparam eEOS) {
    return fDens * (dConstGamma - 1. ) * polytropicEnergyFloor<dtype,mtype>(a_inv3, fDens, dConstGamma, eEOS);
}

template <typename dtype, typename mtype>
inline static dtype eEOSEnergyFloor(dtype a_inv3, dtype fDens, dtype fBall,
                                    dtype dConstGamma, struct eEOSparam eEOS) {
    const dtype polyEnergy = polytropicEnergyFloor<dtype,mtype>(a_inv3, fDens, dConstGamma, eEOS);
    const dtype constantEnergy = constantEnergyFloor<dtype,mtype>(a_inv3, fDens, dConstGamma, eEOS);

    mtype mask = constantEnergy != NOT_IN_EEOS;
    return mask_mov( constantEnergy, mask, max(polyEnergy, constantEnergy) );
}
template <typename dtype, typename mtype>
inline static dtype eEOSPressureFloor(dtype a_inv3, dtype fDens, dtype fBall,
                                      dtype dConstGamma, struct eEOSparam eEOS) {
    const dtype polyPressure = polytropicPressureFloor(a_inv3, fDens, dConstGamma, eEOS);
    const dtype constantPressure = constantPressureFloor(a_inv3, fDens, dConstGamma, eEOS);

    mtype mask = constantPressure != NOT_IN_EEOS;
    return mask_mov( constantPressure, mask, max(polyPressure, constantPressure) );
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
template <typename dtype, typename mtype>
inline static dtype eEOSPressureFloor(dtype a_inv3, dtype fDens, dtype fBall,
                                      dtype dConstGamma, struct eEOSparam eEOS) {
    return 1.2/dConstGamma * pow(eEOS.dNJeans,2./3.) * fDens * fDens * fBall * fBall;
}

template <typename dtype, typename mtype>
inline static dtype eEOSEnergyFloor(dtype a_inv3, dtype fDens, dtype fBall,
                                    dtype dConstGamma, struct eEOSparam eEOS) {
    return eEOSPressureFloor(a_inv3, fDens, fBall, eEOS) /
           (fDens * (dConstGamma - 1.));
}
#endif
#endif
