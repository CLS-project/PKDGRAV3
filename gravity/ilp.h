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

// These are the fields and their types found in the interaction list
#define ILP_FIELDS_SEQ\
    ((float,dx))((float,dy))((float,dz))((float,m))((float,fourh2))\
    ((float,vx))((float,vy))((float,vz))((float,fBall))((float,Omega))\
    ((float,rho))((float,P))((float,c))((int32_t,species))((float,uRung))((float,iMat))\
    ((float,T))((float,expImb2))\
    ((float,Sxx))((float,Syy))((float,Sxy))((float,Sxz))((float,Syz))

ILIST_DECLARE(PP,ILP_FIELDS_SEQ)
#define ILP_FIELD_TYPES ILIST_FIELD_VALUES(ILP_FIELDS_SEQ,0)
#define ILP_FIELD_NAMES ILIST_FIELD_VALUES(ILP_FIELDS_SEQ,1)
#define ILP_FIELD_PROTOS ILIST_FIELD_PROTOS(ILP_FIELDS_SEQ)

using ilpBlock = BlockPP<ILP_PART_PER_BLK>;
using ilpTile = TilePP<ILP_PART_PER_BLK,8>;
class ilpList : public ListPP<ILP_PART_PER_BLK,8>, public ilist::ilCenterReference {
public:
    void append(float dx,float dy,float dz,float m,float fourh2,
                float vx,float vy,float vz,float fBall,float Omega,
                float rho,float P,float c,int32_t species,float uRung,float iMat,
                float T, float expImb2,
                float Sxx, float Syy, float Sxy, float Sxz, float Syz) {
        BlockPP<ILP_PART_PER_BLK> *b;
        int i;
        std::tie(b,i) = ListPP<ILP_PART_PER_BLK,8>::create();
        ILIST_ASSIGN_FIELDS(ILP_FIELDS_SEQ,PP,b,i,)
    }
    void append(double x,double y,double z,float m,float fourh2,
                float vx,float vy,float vz,float fBall,float Omega,
                float rho,float P,float c,int32_t species, int uRung, int iMat,
                float T, float expImb2,
                float Sxx, float Syy, float Sxy, float Sxz, float Syz) {
        append((float)(getReference(0)-x),(float)(getReference(1)-y),(float)(getReference(2)-z),
               m,fourh2,vx,vy,vz,fBall,Omega,rho,P,c,species,(float)uRung,(float)iMat,T,expImb2,Sxx,Syy,Sxy,Sxz,Syz);
    }
};
#endif
