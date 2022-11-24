#ifndef EEOS_STRUCT_H
#define EEOS_STRUCT_H

#define NOT_IN_EEOS 0.0
struct eEOSparam {
#ifdef EEOS_POLYTROPE
    double dPolyFloorMinOD;    // Minimum baryon overdensity to apply the polytrope eos
    double dPolyFloorExponent; // Polytropic index - 1
    double dPolyFloorDen;      // Density normalization (in physical units)
    double dPolyFlooru;        // Internal energy normalization
#endif
#ifdef EEOS_JEANS
    double dNJeans;         // "Number" of resolution elements per Jeans length
#endif
    double dFlooru;      // Constant energy floor
    double dFloorDen;    // Minimum physical density at which the floor is active
    double dFloorMinOD;  // Minimum baryon overdensity to apply the floor
};
#endif
