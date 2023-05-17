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
#ifndef GPU_WORKUNIT_H
#define GPU_WORKUNIT_H

#ifndef __METAL_VERSION__
    #include "stdint.h"
#endif
namespace gpu {
template<int n> using fvector = float[n];
template <int n> struct ppBlk {
    int width() const {return n;}
    fvector<n> dx, dy, dz; // Offset from ilp->cx, cy, cz
    fvector<n> m;          // Mass
    fvector<n> fourh2;     // Softening: calculated
};

template <int n> struct pcBlk {
    int width() const {return n;}
    fvector<n> dx,dy,dz;
    fvector<n> xxxx,xxxy,xxxz,xxyz,xxyy,yyyz,xyyz,xyyy,yyyy;
    fvector<n> xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    fvector<n> xx,xy,xz,yy,yz;
#ifdef USE_DIAPOLE
    fvector<n> x,y,z;
#endif
    fvector<n> m,u;
};
typedef ppBlk<32> ppInteract;
typedef pcBlk<32> pcInteract;

// One of these entries for each interaction block
struct alignas(8) ppWorkUnit {
    uint32_t iP;   // Index of first particle
    uint16_t nP;   // Number of particles
    uint16_t nI;   // Number of interactions in the block
};
static_assert(sizeof(ppWorkUnit)==8,"check size of ppWorkUnit");

struct alignas(32) ppInput {
    float dx, dy, dz;
    float ax, ay, az;
    float fSoft2;
    float dImaga;
};
static_assert(sizeof(ppInput)==32,"Check isze of ppInput");

/* Each thread block outputs this for each particle */
struct alignas(32) ppResult {
    float ax;
    float ay;
    float az;
    float fPot;
    float dirsum;
    float normsum;
};
static_assert(sizeof(ppResult)==32,"Check size of ppResult");
}
#endif
