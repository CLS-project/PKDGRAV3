#ifndef FLOATTYPE_INCLUDED
#define FLOATTYPE_INCLUDED

#include <limits.h>

#ifndef SINGLE

#define FLOAT			double
#define FLOAT_MAXVAL	HUGE

#else

#define FLOAT			float
#define FLOAT_MAXVAL	HUGE

#endif

#endif
