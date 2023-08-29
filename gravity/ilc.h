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

// These are the fields and their types found in the interaction list
#define ILC_FIELDS_PARAMS_SEQ ((float,dx))((float,dy))((float,dz))((float,u))
#ifdef USE_DIAPOLE
    #define ILC_FIELDS_DIAPOLE_SEQ ((float,x))((float,y))((float,z))
#else
    #define ILC_FIELDS_DIAPOLE_SEQ
#endif
#define ILC_FIELDS_MOMENT_SEQ\
 ((float,xxxx))((float,xxxy))((float,xxxz))((float,xxyz))((float,xxyy))((float,yyyz))((float,xyyz))((float,xyyy))((float,yyyy))\
 ((float,xxx))((float,xyy))((float,xxy))((float,yyy))((float,xxz))((float,yyz))((float,xyz))\
 ((float,xx))((float,xy))((float,xz))((float,yy))((float,yz))\
 ILC_FIELDS_DIAPOLE_SEQ ((float,m))
#define ILC_FIELDS_SEQ ILC_FIELDS_PARAMS_SEQ ILC_FIELDS_MOMENT_SEQ

ILIST_DECLARE(PC,ILC_FIELDS_SEQ)
#define ILC_FIELD_TYPES ILIST_FIELD_VALUES(ILC_FIELDS_SEQ,0)
#define ILC_FIELD_NAMES ILIST_FIELD_VALUES(ILC_FIELDS_SEQ,1)
#define ILC_FIELD_PROTOS ILIST_FIELD_PROTOS(ILC_FIELDS_SEQ)

using ilcBlock = BlockPC<ILC_PART_PER_BLK>;
using ilcTile = TilePC<ILC_PART_PER_BLK,8>;
class ilcList : public ListPC<ILC_PART_PER_BLK,8>, public ilist::ilCenterReference {
public:
    void append(float dx,float dy,float dz,const FMOMR *M,float u) {
        BlockPC<ILC_PART_PER_BLK> *b;
        int i;
        std::tie(b,i) = ListPC<ILC_PART_PER_BLK,8>::create();
        ILIST_ASSIGN_FIELDS(ILC_FIELDS_PARAMS_SEQ,PC,b,i,)
        ILIST_ASSIGN_FIELDS(ILC_FIELDS_MOMENT_SEQ,PC,b,i,M->)
    }
    void append(double x,double y,double z,const FMOMR *M,float u) {
        append((float)(getReference(0)-x),(float)(getReference(1)-y),(float)(getReference(2)-z),M,u);
    }
};
#endif
