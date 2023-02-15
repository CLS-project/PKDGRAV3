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

// These are the fields and their types found in the interaction list
#define CL_FIELDS_PARAMS_SEQ\
    ((int32_t,iCache))((int32_t,idCell))((int32_t,iCell))((int32_t,idLower))((int32_t,iLower))((int32_t,idUpper))((int32_t,iUpper))\
    ((int32_t,nc))((float,cOpen))((float,m))((float,fourh2))((float,x))((float,y))((float,z))\
    ((float,xOffset))((float,yOffset))((float,zOffset))((float,xCenter))((float,yCenter))((float,zCenter))\
    ((float,xMax))((float,yMax))((float,zMax))((int32_t,iOpen))

#if SPHBALLOFBALLS
    #define CL_FIELDS_BALLS_SEQ ((float,fBoBr))((float,fBoBxCenter))((float,fBoByCenter))((float,fBoBzCenter))
#elif SPHBOXOFBALLS
    #define CL_FIELDS_BALLS_SEQ ((float,fBoBxMin))((float,fBoByMin))((float,fBoBzMin))((float,fBoBxMax))((float,fBoByMax))((float,fBoBzMax))
#else
    #define CL_FIELDS_BALLS_SEQ
#endif

#define CL_FIELDS_SEQ CL_FIELDS_PARAMS_SEQ CL_FIELDS_BALLS_SEQ

ILIST_DECLARE(CL,CL_FIELDS_SEQ)
#define CL_FIELD_TYPES ILIST_FIELD_VALUES(CL_FIELDS_SEQ,0)
#define CL_FIELD_NAMES ILIST_FIELD_VALUES(CL_FIELDS_SEQ,1)
#define CL_FIELD_PROTOS ILIST_FIELD_PROTOS(CL_FIELDS_SEQ)

using clBlock = BlockCL<CL_PART_PER_BLK>;
using clTile = TileCL<CL_PART_PER_BLK,64>;
class clList : public ListCL<CL_PART_PER_BLK,64> {
public:
    using ListCL<CL_PART_PER_BLK,64>::append;
    typedef ListCL<CL_PART_PER_BLK,64> free_list;
    clList(free_list &freeList) {setFreeList(freeList);}

    void append(uint32_t iCache,uint32_t idCell,uint32_t iCell,
                uint32_t idLower,uint32_t iLower,uint32_t idUpper,uint32_t iUpper,uint32_t nc,float cOpen,float m,float fourh2,
                blitz::TinyVector<float,3> r,blitz::TinyVector<float,3> fOffset,Bound bnd,
                SPHBOB bob) {
        auto fCenter = bnd.center();
        auto fMax = bnd.apothem();
        append(iCache,idCell,iCell,idLower,iLower,idUpper,iUpper,nc,cOpen,m,fourh2,r[0],r[1],r[2],
               fOffset[0],fOffset[1],fOffset[2],fCenter[0],fCenter[1],fCenter[2],fMax[0],fMax[1],fMax[2],0
#if SPHBALLOFBALLS
               ,bob.fBoBr,bob.fBoBCenter[0],bob.fBoBCenter[1],bob.fBoBCenter[2]
#elif SPHBOXOFBALLS
               ,bob.fBoBMin[0],bob.fBoBMin[1],bob.fBoBMin[2],bob.fBoBMax[0],bob.fBoBMax[1],bob.fBoBMax[2]
#endif
              );
    }
    void append(clBlock &B, int Bi) {
        append(B.iCache[Bi],B.idCell[Bi],B.iCell[Bi],
               B.idLower[Bi],B.iLower[Bi],B.idUpper[Bi],B.iUpper[Bi],
               B.nc[Bi], B.cOpen[Bi], B.m[Bi], B.fourh2[Bi],
               B.x[Bi], B.y[Bi], B.z[Bi],B.xOffset[Bi],B.yOffset[Bi],B.zOffset[Bi],
               B.xCenter[Bi],B.yCenter[Bi],B.zCenter[Bi],B.xMax[Bi],B.yMax[Bi],B.zMax[Bi],B.iOpen[Bi]
#if SPHBALLOFBALLS
               ,B.fBoBr[Bi],B.fBoBxCenter[Bi],B.fBoByCenter[Bi],B.fBoBzCenter[Bi]
#elif SPHBOXOFBALLS
               ,B.fBoBxMin[Bi],B.fBoByMin[Bi],B.fBoBzMin[Bi],B.fBoBxMax[Bi],B.fBoByMax[Bi],B.fBoBzMax[Bi]
#endif
              );
    }
};

#endif
