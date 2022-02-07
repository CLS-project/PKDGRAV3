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

#ifndef ILLINOIS_H
#define ILLINOIS_H

typedef struct {
    double r, s;
    double fr, fs;
    double t;
} ILLINOIS;

#ifdef __cplusplus
extern "C" {
#endif
double illinoisInitialize(ILLINOIS *ctx,double r,double fr,double s,double fs);
double illinoisIterate(ILLINOIS *ctx,double ft);
double illinois(double (*func)(double,void *),void *ctx,double r,double s,double xacc,double yacc,int *pnIter);
#ifdef __cplusplus
}
#endif

#endif
