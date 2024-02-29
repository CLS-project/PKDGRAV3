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

template <int n> struct denBlk {
    int width() const {return n;}
    fvector<n> dx, dy, dz;
    fvector<n> m, iMat;
};

template <int n> struct denCorrBlk {
    int width() const {return n;}
    fvector<n> dx, dy, dz;
    fvector<n> T, P, expImb2;
};

template <int n> struct sphForceBlk {
    int width() const {return n;}
    fvector<n> dx, dy, dz;
    fvector<n> m, fBall, Omega;
    fvector<n> vx, vy, vz;
    fvector<n> rho, P, c, rung;
    fvector<n> Sxx, Syy, Sxy, Sxz, Syz;
};

typedef ppBlk<32> ppInteract;
typedef pcBlk<32> pcInteract;
typedef denBlk<32> denInteract;
typedef denCorrBlk<32> denCorrInteract;
typedef sphForceBlk<32> sphForceInteract;

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

struct alignas(32) denInput {
    float dx, dy, dz;
    float fBall, iMat;
};
static_assert(sizeof(denInput)==32,"Check size of denInput");

struct alignas(32) denResult {
    float rho;
    float drhodfball;
    float nden;
    float dndendfball;
    float nSmooth;
    float imbalanceX, imbalanceY, imbalanceZ;
};
static_assert(sizeof(denResult)==32,"Check size of denResult");

struct alignas(16) denCorrInput {
    float dx, dy, dz;
    float fBall;
};
static_assert(sizeof(denCorrInput)==16,"Check size of denCorrInput");

struct alignas(16) denCorrResult {
    float corrT;
    float corrP;
    float corr;
};
static_assert(sizeof(denCorrResult)==16,"Check size of denCorrResult");

struct alignas(64) sphForceInput {
    float dx, dy, dz;
    float fBall, Omega;
    float vx, vy, vz;
    float rho, P, c;
    float Sxx, Syy, Sxy, Sxz, Syz;
};
static_assert(sizeof(sphForceInput)==64,"Check size of sphForceInput");
static_assert(sizeof(sphForceInput)/sizeof(float)<=32,"sphForceInput cannot exceed 32 floats");

struct alignas(128) sphForceResult {
    float uDot;
    float ax, ay, az;
    float divv;
    float dtEst;
    float maxRung;
    float dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
    float Cinvxx, Cinvxy, Cinvxz, Cinvyx, Cinvyy, Cinvyz, Cinvzx, Cinvzy, Cinvzz;
};
static_assert(sizeof(sphForceResult)==128,"Check size of sphForceResult");
}
#endif
