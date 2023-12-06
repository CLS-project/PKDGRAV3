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

#ifndef SPHPREDICT_HINCLUDED
#define SPHPREDICT_HINCLUDED

#include "SPHOptions.h"
#include "SPHEOS.h"
#include "../pkd.h"

#ifdef __cplusplus
extern "C" {
#endif
float getDtPredDrift(struct pkdKickParameters *kick, int bMarked, int uRungLo, int uRung);
void SPHpredictOnTheFly(PKD pkd, particleStore::ParticleReference &p, struct pkdKickParameters *kick, int uRungLo, float *vpred, float *P, float *cs, float *T, SPHOptions *SPHoptions);
void SPHpredictInDensity(PKD pkd, particleStore::ParticleReference &p, struct pkdKickParameters *kick, int uRungLo, SPHOptions *SPHoptions);
#ifdef __cplusplus
}
#endif
#endif