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
#include "particle.h"
using blitz::TinyVector;
using blitz::dot;
using blitz::all;
#include <numeric>

/// @brief Calculate the bounding box for a range of particles
/// @param b starting iterator
/// @param e ending iterator
/// @return bounding box
Bound particleStore::bound(iterator b,iterator e) {
    return this->integerized() ? Bound(raw_bound<int32_t>(b,e),*this) : raw_bound<double>(b,e);
}
/// @brief Calculate the bounding box for all local particles
/// @return bounding box
Bound particleStore::bound() {
    return bound(begin(),end());
}

/*
** This function checks the predicate and returns a new value based on the flags.
** setIfTrue:    >0 -> return true if the predicate is true
**               =0 -> return the old value if the predicate is true
**               <0 -> return false if the predicate is true ("clear if true")
** clearIfFalse: >0 -> return false if the predicate is false
**               =0 -> return the old value if the predicate is false
**               <0 -> return true if the predicate is false ("set if false")
** A value of zero for either results in no action for the "IfTrue" or "IfFalse" flags.
** Conflicting options (e.g., setIfTrue and setIfFalse) result in a toggle.
*/
static inline bool isSelected( bool predicate, int setIfTrue, int clearIfFalse, bool value ) {
    int s = (predicate&(setIfTrue>0)) | (!predicate&(clearIfFalse<0));
    int c = (predicate&(setIfTrue<0)) | (!predicate&(clearIfFalse>0));
    return (~s&~c&value) | (s&~(c&value));
}

/// @brief Count the number of marked particles
/// @return total count of marked particles
int particleStore::CountSelected() {
    return std::accumulate(begin(),end(),0,
    [](int a,auto &p) {return a + p.marked(); });
}

/// @brief Mark active particles
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelActive(int setIfTrue, int clearIfFalse) {
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse](int a,auto &p) {
        return a + p.set_marked(isSelected(p.is_active(),setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark blackhole particles
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelBlackholes(int setIfTrue, int clearIfFalse) {
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse](int a,auto &p) {
        return a + p.set_marked(isSelected(p.is_star() && p.star().fTimer < 0,setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark particles of a certain species
/// @param mSpecies A bit mask indicating which species to mark
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelSpecies(uint64_t mSpecies, int setIfTrue, int clearIfFalse) {
    if (mSpecies&(1<<FIO_SPECIES_ALL)) mSpecies = 0xffffffffu;
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,mSpecies](int a,auto &p) {
        auto b = isSelected((1<<p.species()) & mSpecies,setIfTrue,clearIfFalse,p.marked());
        p.set_NN_flag(b); /* This is a bit clunky, but we only ever use this to reset the flags. */
        return a + p.set_marked(b);
    });
}

/// @brief Mark particles in a certain group
/// @param iGroup group number to mark
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelGroup(int iGroup, int setIfTrue, int clearIfFalse) {
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,iGroup](int a,auto &p) {
        return a + p.set_marked(isSelected(p.group()==iGroup,setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark particles in a certain mass range
/// @param dMinMass minimum mass
/// @param dMaxMass maximum mass
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelMass(double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse ) {
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,dMinMass,dMaxMass](int a,auto &p) {
        auto m = p.mass();
        return a + p.set_marked(isSelected(m >= dMinMass && m <=dMaxMass,setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark particles in a phase space density range
/// @param dMinDensity minimum phase space density
/// @param dMaxDensity maximum phase space density
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelPhaseDensity(double dMinDensity, double dMaxDensity, int setIfTrue, int clearIfFalse ) {
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,dMinDensity,dMaxDensity](int a,auto &p) {
        const auto &vel = p.VelSmooth();
        float density = p.density() * pow(vel.veldisp2,-1.5);
        return a + p.set_marked(isSelected(density >= dMinDensity && density <=dMaxDensity,setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark particles matching an ID range
/// @param idStart mimimum particle id
/// @param idEnd maximum particle id (inclusive)
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelById(uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse ) {
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,idStart,idEnd](int a,auto &p) {
        auto id = p.order();
        return a + p.set_marked(isSelected(id >= idStart && id <= idEnd,setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark particles inside a box
/// @param dCenter center coordinate of the box
/// @param dSize size of the box (half dimension)
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelBox(TinyVector<double,3> dCenter, TinyVector<double,3> dSize, int setIfTrue, int clearIfFalse ) {
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,dCenter,dSize](int a,auto &p) {
        TinyVector<double,3> dx = dCenter - p.position();
        return a + p.set_marked(isSelected(all(dx < dSize) && all(dx >= -dSize),setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark particles inside a sphere
/// @param r coordinate of the center of the sphere
/// @param dRadius radius of the sphere
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelSphere(TinyVector<double,3> r, double dRadius, int setIfTrue, int clearIfFalse ) {
    auto dRadius2 = dRadius * dRadius;
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,r,dRadius2](int a,auto &p) {
        TinyVector<double,3> dx = r - p.position();
        return a + p.set_marked(isSelected(dot(dx,dx)<=dRadius2,setIfTrue,clearIfFalse,p.marked()));
    });
}

/// @brief Mark particles inside a cylinder
/// @param dP1 first endpoint of the cylinder
/// @param dP2 second endpoint of the cylinder
/// @param dRadius radius of the cylinder
/// @param setIfTrue
/// @param clearIfFalse
/// @return total count of marked particles
int particleStore::SelCylinder(TinyVector<double,3> dP1, TinyVector<double,3> dP2, double dRadius, int setIfTrue, int clearIfFalse ) {
    TinyVector<double,3> dCyl = dP2 - dP1;
    auto dLength2 = dot(dCyl,dCyl);
    auto dRadius2 = dRadius*dRadius;
    return std::accumulate(begin(),end(),0,
    [setIfTrue,clearIfFalse,dRadius2,dLength2,dCyl,dP1](int a,auto &p) {
        TinyVector<double,3> dPart = p.position() - dP1;
        auto pdotr = dot(dPart,dCyl);
        bool predicate = pdotr >= 0.0 && pdotr <= dLength2 && (dot(dPart,dPart) - pdotr*pdotr/dLength2 <= dRadius2);
        return a + p.set_marked(isSelected(predicate,setIfTrue,clearIfFalse,p.marked()));
    });
}
