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

#include "SPHOptions.h"
#include "pkd_parameters.h"
#include <stdio.h>
#include "basetype.h"
#include "cosmo.h"

auto calculateInterfaceCorrectionPrefactor(float nSmooth, int kernelType) -> float {
    float alpha = 0.0f;
    switch (kernelType) {
    case 0: {
        if (nSmooth < 10.0f) nSmooth = 10.0f; // Below 10, the function decreases
        alpha = -7.44731686e+02f / (nSmooth * nSmooth) + 1.42956727e+02f / nSmooth + 4.46685213e+00f;
        break;
    }
    case 1: {
        if (nSmooth < 3.0f) nSmooth = 3.0f; // Below 3, the function decreases
        alpha = -2.27180786e+02f / (nSmooth * nSmooth) + 1.55942641e+02f / nSmooth + 4.53619157e+00f;
        break;
    }
    case 2: {
        alpha = 1.81941313e+03f / (nSmooth * nSmooth) + 1.73546428e+02f / nSmooth + 4.71450606e+00f;
        break;
    }
    case 3: {
        alpha = 6.00324255e+03f / (nSmooth * nSmooth) + 1.71766754e+02f / nSmooth + 4.92489934e+00f;
        break;
    }
    default:
        assert(0);
    }
    return alpha;
}
