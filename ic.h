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

#ifndef IC_H
#define IC_H

typedef struct {
    float dr[3];
    float v[3];
    } basicParticle;

typedef struct {
    uint64_t ix : 21;
    uint64_t iy : 21;
    uint64_t iz : 21;
    float dr[3];
    float v[3];
    } expandParticle;

typedef union {
    basicParticle b;
    expandParticle e;
    } overlayedParticle;

#ifdef MDL_FFTW

typedef struct {
    FFTW3(real) x,y,z;
    } gridpos;

typedef struct {
    FFTW3(real) x,y,z;
    FFTW3(real) vx,vy,vz;
    } gridpsc;

#ifdef __cplusplus
extern "C"
#endif
int pkdGenerateIC(PKD pkd,MDLFFT fft,int iSeed,int bFixed,float fPhase,int nGrid,int b2LPT,double dBoxSize,
    struct csmVariables *cosmo,double a,int nTf, double *tk, double *tf,
    double *noiseMean, double *noiseCSQ);
#ifdef __cplusplus
extern "C"
#endif
int pkdGenerateClassICm(PKD pkd, MDLFFT fft, int iSeed, int bFixed, float fPhase, int nGrid,
    double dBoxSize, struct csmVariables *cosmo, double a, double *noiseMean, double *noiseCSQ);
#ifdef __cplusplus
extern "C"
#endif
void pkdGenerateNuGrid(PKD pkd, MDLFFT fft,double a,double LBox, int iSeed, int bFixed, float fPhase);
#endif
#endif
