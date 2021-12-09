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

#ifndef BOUND_H
#define BOUND_H
#include "stdint.h"

typedef struct bndBound {
    double fCenter[3];
    double fMax[3];
    } BND;
typedef struct {
    int32_t fCenter[3];
    int32_t fMax[3];
    } IBND;

#ifdef __cplusplus
#include <ostream>
#include <utility>
#include <algorithm>
#include <type_traits>
#include "blitz/array.h"

class Bound : public BND {
public:
    typedef blitz::TinyVector<double,3> coord;
public:
    Bound() = default;
    Bound(const Bound &) = default;
    Bound(const BND &bnd) : BND(bnd) {}
    Bound(coord lower, coord upper)
	    : BND {{0.5 * (upper[0]+lower[0]), 0.5 * (upper[1]+lower[1]), 0.5 * (upper[2]+lower[2])},
	           {0.5 * (upper[0]-lower[0]), 0.5 * (upper[1]-lower[1]), 0.5 * (upper[2]-lower[2])}}
	  {}
    Bound(const double *lower, const double *upper)
	    : BND {{0.5 * (upper[0]+lower[0]), 0.5 * (upper[1]+lower[1]), 0.5 * (upper[2]+lower[2])},
	           {0.5 * (upper[0]-lower[0]), 0.5 * (upper[1]-lower[1]), 0.5 * (upper[2]-lower[2])}}
	  {}
    Bound(double lx, double ly, double lz, double ux, double uy, double uz)
	    : BND {{0.5 * (ux+lx), 0.5 * (uy+ly), 0.5 * (uz+lz)},
	           {0.5 * (ux-lx), 0.5 * (uy-ly), 0.5 * (uz-lz)}}
	  {}
public:
    void set(int d,double lower,double upper) {
    	fCenter[d] = 0.5*(upper+lower);
    	fMax[d] = 0.5*(upper-lower);
	}
    void set(coord lower,coord upper) {
	for(auto d=0; d<coord::length(); ++d)
	    set(d,lower[d],upper[d]);
	}
public:
    coord  width()       const {return coord(2.0*fMax[0],2.0*fMax[1],2.0*fMax[2]);}
    double width(int d)  const {return 2.0*fMax[d];}
    coord  lower()       const {return coord(fCenter[0]-fMax[0],fCenter[1]-fMax[1],fCenter[2]-fMax[2]);}
    double lower(int d)  const {return fCenter[d]-fMax[d];}
    coord  upper()       const {return coord(fCenter[0]+fMax[0],fCenter[1]+fMax[1],fCenter[2]+fMax[2]);}
    double upper(int d)  const {return fCenter[d]+fMax[d];}
    coord  center()      const {return coord(fCenter[0],fCenter[1],fCenter[2]);}
    double center(int d) const {return fCenter[d];}
    double minside()     const {return 2.0*std::min({fMax[2],fMax[1],fMax[0]});}
    double maxside()     const {return 2.0*std::max({fMax[2],fMax[1],fMax[0]});}
    int    maxdim()      const {return (fMax[0]>fMax[2]) ? (fMax[0]>fMax[1]?0:1) : (fMax[2]>fMax[1]?2:1);}
    int    mindim()      const {return (fMax[0]<fMax[2]) ? (fMax[0]<fMax[1]?0:1) : (fMax[2]<fMax[1]?2:1);}
    double maxdist(coord r) const {
    	coord x = blitz::abs(coord(fCenter)-r) + coord(fMax);
    	return blitz::dot(x,x);
	}
    double mindist(coord r) const {
        coord x = blitz::abs(coord(fCenter)-r) - coord(fMax);
        x = blitz::where(x>0.0,x,0.0);
        return blitz::dot(x,x);
        }
    // Return a new bound that includes both bounds
    Bound  combine(const Bound &rhs) const {
        Bound b;
        for(auto d=0; d<coord::length(); ++d) {
            auto l = std::min(lower(d),rhs.lower(d));
            auto u = std::max(upper(d),rhs.upper(d));
            b.fCenter[d] = 0.5*(u+l);
            b.fMax[d] = 0.5*(u-l);
            }
        return b;
        }
    // Return two new bounds split at "split" along dimension d
    std::pair<Bound,Bound> split(int d,double split) const {
	Bound l, r;
	for (auto j=0;j<3;++j) {
	    if (j == d) {
		l.fMax[j]    = 0.5 * (split - lower(d));
		l.fCenter[j] = 0.5 * (split + lower(d));
		r.fMax[j]    = 0.5 * (upper(d) - split);
		r.fCenter[j] = 0.5 * (upper(d) + split);
		}
	    else {
		l.fMax[j] = r.fMax[j] = fMax[j];
		l.fCenter[j] = r.fCenter[j] = fCenter[j];
		}
	    }
	return std::make_pair(l,r);
	}
    // Return two new bounds split in half along dimension d
    std::pair<Bound,Bound> split(int d) const {
	Bound l, r;
	for (auto j=0;j<3;++j) {
	    if (j == d) {
		r.fMax[j] = l.fMax[j] = 0.5*fMax[j];
		l.fCenter[j] = fCenter[j] - l.fMax[j];
		r.fCenter[j] = fCenter[j] + r.fMax[j];
		}
	    else {
		l.fMax[j] = r.fMax[j] = fMax[j];
		l.fCenter[j] = r.fCenter[j] = fCenter[j];
		}
	    }
	return std::make_pair(l,r);
	}
    friend std::ostream& operator<<(std::ostream& os, const Bound& b);
    };

inline std::ostream& operator<<(std::ostream& os, const Bound& b)
{
    os << "[" << b.lower() << "," << b.upper() << "]";
    return os;
}

static_assert(std::is_trivial<Bound>());
static_assert(std::is_standard_layout<Bound>());
#endif // __cplusplus
#endif
