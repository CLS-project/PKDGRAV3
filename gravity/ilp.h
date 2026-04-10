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

#ifndef ILP_H
#define ILP_H

#include "lst.h"

#ifndef ILP_PART_PER_BLK
    #define ILP_PART_PER_BLK 32 /* Don't mess with this: see CUDA */
#endif

// clang-format off
#define ILP_FIELDS(X)                                                                                                                                          \
    X(dx, float) X(dy, float) X(dz, float) X(m, float) X(fourh2, float) X(vx, float) X(vy, float) X(vz, float) X(fBall, float) X(Omega, float) X(rho, float)   \
    X(P, float) X(c, float) X(species, int32_t) X(uRung, float) X(iMat, float) X(T, float) X(expImb2, float) X(isGas, float) X(Sxx, float) X(Syy, float)       \
    X(Sxy, float) X(Sxz, float) X(Syz, float) X(spinx, float) X(spiny, float) X(spinz, float) X(lastPosx, float) X(lastPosy, float) X(lastPosz, float)
// clang-format on

ILIST_DEFINE(PP, ILP_FIELDS)

#define ILP_ASSIGN_BASIC(name, ...) b->name[i] = name;
#define ILP_INIT_BASIC(name, ...) .name = name,

class ilpList : public ListPP<ILP_PART_PER_BLK, 8>, public ilist::ilCenterReference {
  public:
    void append(float dx, float dy, float dz, float m, float fourh2, float vx, float vy, float vz, float fBall, float Omega, float rho, float P, float c,
                int32_t species, float uRung, float iMat, float T, float expImb2, float isGas, float Sxx, float Syy, float Sxy, float Sxz, float Syz,
                float spinx, float spiny, float spinz, float lastPosx, float lastPosy, float lastPosz) {
        ListPP<ILP_PART_PER_BLK, 8>::append(items_type{ILP_FIELDS(ILP_INIT_BASIC)});
    }
    void append(double x, double y, double z, float m, float fourh2, float vx, float vy, float vz, float fBall, float Omega, float rho, float P, float c,
                int32_t species, int uRung, int iMat, float T, float expImb2, float isGas, float Sxx, float Syy, float Sxy, float Sxz, float Syz, float spinx,
                float spiny, float spinz, double lastPosx, double lastPosy, double lastPosz) {
        append((float)(getReference(0) - x), (float)(getReference(1) - y), (float)(getReference(2) - z), m, fourh2, vx, vy, vz, fBall, Omega, rho, P, c,
               species, (float)uRung, (float)iMat, T, expImb2, isGas, Sxx, Syy, Sxy, Sxz, Syz, spinx, spiny, spinz, (float)(getReference(0) - lastPosx),
               (float)(getReference(1) - lastPosy), (float)(getReference(2) - lastPosz));
    }
};
using ilpTile = ilpList::tile_type;
using ilpBlock = ilpList::block_type;

#endif
