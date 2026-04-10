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

#ifndef CL_H
#define CL_H
#include "lst.h"
#include "core/bound.h"
#include "../SPH/SPHOptions.h"

#ifndef CL_PART_PER_BLK
    #define CL_PART_PER_BLK 16
#endif

// clang-format off
#define CL_FIELDS_BASIC(X)                                                                                                                                     \
    X(iCache, int32_t) X(idCell, int32_t) X(iCell, int32_t) X(idLower, int32_t) X(iLower, int32_t) X(idUpper, int32_t) X(iUpper, int32_t) X(nc, int32_t)       \
    X(cOpen, float) X(m, float) X(fourh2, float) X(x, float) X(y, float) X(z, float) X(xOffset, float) X(yOffset, float) X(zOffset, float)                     \
    X(xCenter, float) X(yCenter, float) X(zCenter, float) X(xMax, float) X(yMax, float) X(zMax, float) X(iOpen, int32_t)
#define CL_FIELDS_BALLS(X)                                                                                                                                     \
    X(fBoBr, float) X(fBoBxCenter, float) X(fBoByCenter, float) X(fBoBzCenter, float) X(fBoBxMin, float) X(fBoByMin, float) X(fBoBzMin, float)                 \
    X(fBoBxMax, float) X(fBoByMax, float) X(fBoBzMax, float)
// clang-format on
#define CL_FIELDS(X) CL_FIELDS_BASIC(X) CL_FIELDS_BALLS(X)
ILIST_DEFINE(CL, CL_FIELDS)

class clList : public ListCL<CL_PART_PER_BLK, 64> {
  public:
    using ListCL<CL_PART_PER_BLK, 64>::append;
    using free_list = ListCL<CL_PART_PER_BLK, 64>;
    clList() = default;
    clList(List::list_type& freeList) : ListCL<CL_PART_PER_BLK, 64>(freeList) {}

    void append(int32_t iCache, int32_t idCell, int32_t iCell, int32_t idLower, int32_t iLower, int32_t idUpper, int32_t iUpper, int32_t nc, float cOpen,
                float m, float fourh2, blitz::TinyVector<float, 3> r, blitz::TinyVector<float, 3> fOffset, Bound bnd, SPHBOB bob) {
        auto fCenter = bnd.center();
        auto fMax = bnd.apothem();
        append(items_type{.iCache = iCache,
                          .idCell = idCell,
                          .iCell = iCell,
                          .idLower = idLower,
                          .iLower = iLower,
                          .idUpper = idUpper,
                          .iUpper = iUpper,
                          .nc = nc,
                          .cOpen = cOpen,
                          .m = m,
                          .fourh2 = fourh2,
                          .x = r[0],
                          .y = r[1],
                          .z = r[2],
                          .xOffset = fOffset[0],
                          .yOffset = fOffset[1],
                          .zOffset = fOffset[2],
                          .xCenter = static_cast<float>(fCenter[0]),
                          .yCenter = static_cast<float>(fCenter[1]),
                          .zCenter = static_cast<float>(fCenter[2]),
                          .xMax = static_cast<float>(fMax[0]),
                          .yMax = static_cast<float>(fMax[1]),
                          .zMax = static_cast<float>(fMax[2]),
                          .iOpen = 0,
                          .fBoBr = bob.fBoBr,
                          .fBoBxCenter = bob.fBoBCenter[0],
                          .fBoByCenter = bob.fBoBCenter[1],
                          .fBoBzCenter = bob.fBoBCenter[2],
                          .fBoBxMin = bob.fBoBMin[0],
                          .fBoByMin = bob.fBoBMin[1],
                          .fBoBzMin = bob.fBoBMin[2],
                          .fBoBxMax = bob.fBoBMax[0],
                          .fBoByMax = bob.fBoBMax[1],
                          .fBoBzMax = bob.fBoBMax[2]});
    }
    void append(block_type& B, int Bi) {

        append(items_type{.iCache = B.iCache[Bi],
                          .idCell = B.idCell[Bi],
                          .iCell = B.iCell[Bi],
                          .idLower = B.idLower[Bi],
                          .iLower = B.iLower[Bi],
                          .idUpper = B.idUpper[Bi],
                          .iUpper = B.iUpper[Bi],
                          .nc = B.nc[Bi],
                          .cOpen = B.cOpen[Bi],
                          .m = B.m[Bi],
                          .fourh2 = B.fourh2[Bi],
                          .x = B.x[Bi],
                          .y = B.y[Bi],
                          .z = B.z[Bi],
                          .xOffset = B.xOffset[Bi],
                          .yOffset = B.yOffset[Bi],
                          .zOffset = B.zOffset[Bi],
                          .xCenter = B.xCenter[Bi],
                          .yCenter = B.yCenter[Bi],
                          .zCenter = B.zCenter[Bi],
                          .xMax = B.xMax[Bi],
                          .yMax = B.yMax[Bi],
                          .zMax = B.zMax[Bi],
                          .iOpen = B.iOpen[Bi],
                          .fBoBr = B.fBoBr[Bi],
                          .fBoBxCenter = B.fBoBxCenter[Bi],
                          .fBoByCenter = B.fBoByCenter[Bi],
                          .fBoBzCenter = B.fBoBzCenter[Bi],
                          .fBoBxMin = B.fBoBxMin[Bi],
                          .fBoByMin = B.fBoByMin[Bi],
                          .fBoBzMin = B.fBoBzMin[Bi],
                          .fBoBxMax = B.fBoBxMax[Bi],
                          .fBoByMax = B.fBoByMax[Bi],
                          .fBoBzMax = B.fBoBzMax[Bi]});
    }
};

using clBlock = clList::block_type;
using clTile = clList::tile_type;

#endif
