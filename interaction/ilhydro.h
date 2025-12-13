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

#ifndef ILHYDRO_H
#define ILHYDRO_H

#include "gravity/lst.h"

#ifndef ILHYDRO_PART_PER_BLK
    #define ILHYDRO_PART_PER_BLK 32 /* Don't mess with this: see CUDA */
#endif

// (float, qgradRhox))((float, qgradRhoy))((float, qgradRhoz))(             \
// (float, qgradPx))((float, qgradPy))((float, qgradPz))(                   \
// (float, qgradVxx))((float, qgradVxy))((float, qgradVxz))(                \
// (float, qgradVyx))((float, qgradVyy))((float, qgradVyz))(                \
// (float, qgradVzx))((float, qgradVzy))((float, qgradVzz))(                \
// (float, qlastAccx))((float, qlastAccy))((float, qlastAccz))(             \


// These are the fields and their types found in the interaction list
#define ILHYDRO_FIELDS_SEQ\
    ((float,dx))((float,dy))((float,dz))\
    ((float,qh))((float,qlast))((float,qDt))((float,qomega))\
    ((float,qspheom))((float,qBXX))((float,qBXY))((float,qBXZ))\
    ((float,qBYY))((float,qBYZ))((float,qBZZ))((float,qDens))\
    ((float,qvx))((float,qvy))((float,qvz))((float,qP))((float,qS))\
    ((float,isGas))((float,isActive))((float,isMarked))


ILIST_DECLARE(HYDRO,ILHYDRO_FIELDS_SEQ)
#define ILHYDRO_FIELD_TYPES ILIST_FIELD_VALUES(ILHYDRO_FIELDS_SEQ,0)
#define ILHYDRO_FIELD_NAMES ILIST_FIELD_VALUES(ILHYDRO_FIELDS_SEQ,1)
#define ILHYDRO_FIELD_PROTOS ILIST_FIELD_PROTOS(ILHYDRO_FIELDS_SEQ)

using ilhydroBlock = BlockHYDRO<ILHYDRO_PART_PER_BLK>;
using ilhydroTile = TileHYDRO<ILHYDRO_PART_PER_BLK,8>;
class ilhydroList : public ListHYDRO<ILHYDRO_PART_PER_BLK,8>,
    public ilist::ilCenterReference {
public:
    void append(float dx, float dy, float dz, float qh,
                float qlast, float qDt, float qomega, float qspheom, float qBXX,
                float qBXY, float qBXZ, float qBYY, float qBYZ, float qBZZ,
                float qDens, float qvx, float qvy, float qvz, float qP, float qS,
                // float qgradRhox, float qgradRhoy, float qgradRhoz, float qgradPx,
                // float qgradPy, float qgradPz, float qgradVxx, float qgradVxy,
                // float qgradVxz, float qgradVyx, float qgradVyy, float qgradVyz,
                // float qgradVzx, float qgradVzy, float qgradVzz, float qlastAccx,
                // float qlastAccy, float qlastAccz,
                float isGas, float isActive, float isMarked) {
        BlockHYDRO<ILHYDRO_PART_PER_BLK> *b;
        int i;
        std::tie(b, i) = ListHYDRO<ILHYDRO_PART_PER_BLK, 8>::create();
        ILIST_ASSIGN_FIELDS(ILHYDRO_FIELDS_SEQ, HYDRO, b, i,)
    }
    void append(double x, double y, double z, float qh,
                float qlast, float qDt, float qomega, float qspheom, float qBXX,
                float qBXY, float qBXZ, float qBYY, float qBYZ, float qBZZ,
                float qDens, float qvx, float qvy, float qvz, float qP, float qS,
                // float qgradRhox, float qgradRhoy, float qgradRhoz, float qgradPx,
                // float qgradPy, float qgradPz, float qgradVxx, float qgradVxy,
                // float qgradVxz, float qgradVyx, float qgradVyy, float qgradVyz,
                // float qgradVzx, float qgradVzy, float qgradVzz, float qlastAccx,
                // float qlastAccy, float qlastAccz,
                float isGas, float isActive, float isMarked) {
        append((float)(getReference(0) - x), (float)(getReference(1) - y),
               (float)(getReference(2) - z), qh, qlast, qDt, qomega, qspheom,
               qBXX, qBXY, qBXZ, qBYY, qBYZ, qBZZ, qDens, qvx, qvy, qvz, qP, qS,
               // qgradRhox, qgradRhoy, qgradRhoz, qgradPx, qgradPy, qgradPz, qgradVxx,
               // qgradVxy, qgradVxz, qgradVyx, qgradVyy, qgradVyz, qgradVzx, qgradVzy,
               // qgradVzz, qlastAccx, qlastAccy, qlastAccz,
               isGas, isActive, isMarked);
    }
};
#endif
