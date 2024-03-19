#ifndef POTENTIAL_HINCLUDED
#define POTENTIAL_HINCLUDED

#if defined(HERNQUIST_POTENTIAL) and defined(NFW_POTENTIAL)
    static_assert(0, "Can not activate two external potentials at the same time");
#endif

typedef std::tuple<blitz::TinyVector<double,3>, double> out_potential;

#define POT_ACC 0
#define POT_DT  1

#ifdef HERNQUIST_POTENTIAL
    #include "potential/hernquist.h"
#endif

#ifdef NFW_POTENTIAL
    #include "potential/nfw.h"
#endif

#endif
