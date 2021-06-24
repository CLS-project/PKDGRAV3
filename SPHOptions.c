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
#include "parameters.h"
#include <stdio.h>

SPHOptions initializeSPHOptions(struct parameters param){
    SPHOptions SPHoptions;
    SPHoptions.fKernelTarget = param.fKernelTarget;
    SPHoptions.epsilon = 0.01f;
    SPHoptions.alpha = param.dConstAlpha;
    SPHoptions.beta = param.dConstBeta;
    SPHoptions.EtaCourant = param.dEtaCourant;
    SPHoptions.doGravity = 0;
    SPHoptions.doDensity = 0;
    SPHoptions.doSPHForces = 0;
    SPHoptions.useNumDen = 0;
    SPHoptions.kernelType = 0;
    return SPHoptions;
}