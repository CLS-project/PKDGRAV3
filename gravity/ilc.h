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

#ifndef ILC_H
#define ILC_H

#include "lst.h"
#include "gravity/moments.h"

#ifndef ILC_PART_PER_BLK
    #define ILC_PART_PER_BLK 32 /* Don't mess with this: see CUDA */
#endif

// Moments are expanded around the center of mass and thus do not have a dipole term.
// Exanding around the cell center was considered but was found to be less accurate.
#ifdef USE_DIAPOLE
    #define ILC_FIELDS_DIAPOLE(X) X(x, float) X(y, float) X(z, float)
#else
    #define ILC_FIELDS_DIAPOLE(X)
#endif
#define ILC_FIELDS_BASIC(X) X(dx, float) X(dy, float) X(dz, float) ILC_FIELDS_DIAPOLE(X)
// clang-format off
#define ILC_FIELDS_MOMENT(X)                                                                          \
    X(xxxx, float) X(xxxy, float) X(xxxz, float) X(xxyz, float) X(xxyy, float)                        \
    X(yyyz, float) X(xyyz, float) X(xyyy, float) X(yyyy, float)                                       \
    X(xxx, float) X(xyy, float) X(xxy, float) X(yyy, float) X(xxz, float) X(yyz, float) X(xyz, float) \
    X(xx, float) X(xy, float) X(xz, float) X(yy, float) X(yz, float) X(m, float)
// clang-format on

#define ILC_FIELDS(X) ILC_FIELDS_BASIC(X) ILC_FIELDS_MOMENT(X) ILC_FIELDS_DIAPOLE(X) X(u, float)

ILIST_DEFINE(PC, ILC_FIELDS)

#define ILC_PROTO_TYPE(name, ...) __VA_ARGS__ name,
#define ILC_ASSIGN_BASIC(name, ...) b->name[i] = name;
#define ILC_ASSIGN_MOMENT(name, ...) b->name[i] = M->name;
class ilcList : public ListPC<ILC_PART_PER_BLK, 8>, public ilist::ilCenterReference {
  public:
    void append(ILC_FIELDS_BASIC(ILC_PROTO_TYPE) const FMOMR* M, float u) {
        auto [b, i] = this->create();
        ILC_FIELDS_BASIC(ILC_ASSIGN_BASIC)
        ILC_FIELDS_MOMENT(ILC_ASSIGN_MOMENT)
        b->u[i] = u;
    }
    void append(double x, double y, double z, const FMOMR* M, float u) {
        append((float)(getReference(0) - x), (float)(getReference(1) - y), (float)(getReference(2) - z), M, u);
    }
};
using ilcTile = ilcList::tile_type;
using ilcBlock = ilcList::block_type;

#endif
