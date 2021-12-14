#if defined(EEOS_POLYTROPE) || defined(EEOS_JEANS)

#ifdef EEOS_POLYTROPE
/*
 * General polytropic eEOS floor
 */
inline static double polytropicEnergyFloor(float a_inv3, float fDens,
                     double dEOSPolyFloorIndex, double dEOSPolyFloorDen, double dEOSPolyFlooru){
   return dEOSPolyFlooru * pow( fDens*a_inv3/dEOSPolyFloorDen, dEOSPolyFloorIndex );
}

inline static double polytropicPressureFloor(float a_inv3, float fDens, double dConstGamma,
                     double dEOSPolyFloorIndex, double dEOSPolyFloorDen, double dEOSPolyFlooru){
   return fDens * (dConstGamma - 1. ) *
      polytropicEnergyFloor(a_inv3, fDens, dEOSPolyFloorIndex, dEOSPolyFloorDen,  dEOSPolyFlooru);
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
inline static double jeansPressureFloor(float fDens, float fBall
                                        double dConstGamma, double dEOSNJeans){
   return 1.2/dConstGamma * dEOSNJeans * fDens * fDens * fBall * fBall;

}

inline static double jeansEnergyFloor(float fDens, float fBall,
                                      double dConstGamma, double dEOSNJeans){
   return jeansPressureFloor(fDens, fBall, dConstGamma, dEOSNJeans) /
                        (fDens * (dConstGamma - 1.));

}
#endif
#endif

