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
#ifndef AWEIGHTS_HPP
#define AWEIGHTS_HPP
#include <cmath>
#include <vector>

template<int Order,typename F>
class AssignmentWeights {
    template <int A, typename B> struct identity {};
    void weights(identity<1,F> d, F r) {		// NGP
    	i = floorf(r);
    	H[0] = 1.0;
	}
    void weights(identity<2,F> d, F r) {		// CIC
	F rr = r - 0.5;
	i = floorf(rr);
	F h = rr - i;
	H[0] = 1.0 - h;
	H[1] = h;
	}
    void weights(identity<3,F> d, F r) {		// TSC
	auto K0 = [](F h) { return 0.75 - h*h; };
	auto K1 = [](F h) { return 0.50 * h*h; };
	i = floorf(r) - 1;
	F h = r - i - 1.5;
	H[0] = K1(0.5 - h);
	H[1] = K0(h);
	H[2] = K1(0.5 + h);
	}
    void weights(identity<4,F> d, F r) {		// PCS
	auto pow3 = [](F x) { return x*x*x; };
	auto K0   = [](F h) { return 1.0/6.0 * ( 4.0 - 6.0*h*h + 3.0*h*h*h); };
	auto K1   = [&pow3](F h) { return 1.0/6.0 * pow3(2.0 - h); };
	i = floorf(r-1.5);
	F h = r - (i+0.5);
	H[0] = K1(h);
	H[1] = K0(h-1);
	H[2] = K0(2-h);
	H[3] = K1(3-h);
        }
public:
    F H[Order];
    int i;
    AssignmentWeights(F r) { weights(identity<Order,F>(),r); }
    };

class AssignmentWindow : public std::vector<float> {
public:
    AssignmentWindow(int nGrid,int iAssignment) {
	reserve(nGrid);
	for( auto i=0; i<nGrid; ++i) {
	    float win = M_PI * i / nGrid;
	    if(win>0.1) win = win / sinf(win);
	    else win=1.0 / (1.0-win*win/6.0*(1.0-win*win/20.0*(1.0-win*win/76.0)));
	    push_back(powf(win,iAssignment));
	    }
	}
    };
#endif
