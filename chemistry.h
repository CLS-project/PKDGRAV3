#ifndef CHEMISTRY_HINCLUDED
#define CHEMISTRY_HINCLUDED

enum chemical_elements {
    ELEMENT_H = 0,
#ifdef HAVE_HELIUM
    ELEMENT_He,
#endif
#ifdef HAVE_CARBON
    ELEMENT_C,
#endif
#ifdef HAVE_NITROGEN
    ELEMENT_N,
#endif
#ifdef HAVE_OXYGEN
    ELEMENT_O,
#endif
#ifdef HAVE_NEON
    ELEMENT_Ne,
#endif
#ifdef HAVE_MAGNESIUM
    ELEMENT_Mg,
#endif
#ifdef HAVE_SILICON
    ELEMENT_Si,
#endif
#ifdef HAVE_IRON
    ELEMENT_Fe,
#endif

    ELEMENT_COUNT   /* Should always be the last */
};

#endif /* CHEMISTRY_HINCLUDED */

