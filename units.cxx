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
#include "pkd_parameters.h"
#include "units.h"

UNITS::UNITS(pkd_parameters &parameters,double h) {
    // Input units (from the parameter file)
    dMsolUnit = parameters.get_dMsolUnit();
    dKpcUnit = parameters.get_dKpcUnit();
    dKBoltzUnit = parameters.get_dKBoltzUnit();
    dGasConst = parameters.get_dGasConst();

    // Convert kboltz/mhydrogen to system units, assuming that G == 1
    if (parameters.has_dMsolUnit() && parameters.has_dKpcUnit()) {
        /* code KBOLTZ/MHYDR */
        dGasConst = dKpcUnit*KPCCM*KBOLTZ/MHYDR/GCGS/dMsolUnit/MSOLG;
        /* code energy per unit mass --> erg per g */
        dErgPerGmUnit = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
        /* code energy --> erg */
        dErgUnit = GCGS*pow(dMsolUnit*MSOLG,2.0)/(dKpcUnit*KPCCM);
        /* code density --> g per cc */
        dGmPerCcUnit = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM,3.0);
        /* code time --> seconds */
        dSecUnit = sqrt(1/(dGmPerCcUnit*GCGS));
        /* code speed --> km/s */
        dKmPerSecUnit = sqrt(GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM))/1e5;
        /* code comove density -->g per cc = units.dGmPerCcUnit(1+z)^3*/
        dComovingGmPerCcUnit = dGmPerCcUnit;
    }
    else if (parameters.get_nGrid()) {
        // We need to properly set a unit system, we do so following the
        // convention: G=1, rho=Omega0 in code units
        dKpcUnit = parameters.get_dBoxSize()*1e3 / h;

        // The mass unit is set such that we recover a correct dHubble0
        // in code units and 100h in physical
        const double dHubbleCGS = 100.*h*1e5/(1e3*KPCCM); // 1/s
        dMsolUnit = pow(dKpcUnit * KPCCM, 3 ) / MSOLG * 3.0 * pow( dHubbleCGS, 2 ) * M_1_PI / 8.0 / GCGS;

        /* code KBOLTZ/MHYDR */
        dGasConst = dKpcUnit*KPCCM*KBOLTZ/MHYDR/GCGS/dMsolUnit/MSOLG;
        /* code energy per unit mass --> erg per g */
        dErgPerGmUnit = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
        /* code energy --> erg */
        dErgUnit = GCGS*pow(dMsolUnit*MSOLG,2.0)/(dKpcUnit*KPCCM);
        /* code density --> g per cc */
        dGmPerCcUnit = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM,3.0);
        /* code time --> seconds */
        dSecUnit = sqrt(1/(dGmPerCcUnit*GCGS));
        /* code speed --> km/s */
        dKmPerSecUnit = sqrt(GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM))/1e5;
        /* code comove density -->g per cc = units.dGmPerCcUnit(1+z)^3*/
        dComovingGmPerCcUnit = dGmPerCcUnit;

        // Some safety checks
        // double H0 = h * 100. / dKmPerSecUnit * dKpcUnit/1e3;
        // double rhoCrit = 3.*H0*H0/(8.*M_PI);
        // assert( fabs(H0-csm->val.dHubble0)/H0 < 0.01 );
        // assert( fabs(rhoCrit-1.0) < 0.01 );
    }
    else {
        dSecUnit = 1;
        dKmPerSecUnit = 1;
        dComovingGmPerCcUnit = 1;
        dGmPerCcUnit = 1;
        dErgPerGmUnit = 1;
        dErgUnit = 1;
    }

}
