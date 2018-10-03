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

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif
#include "cosmo.h"

//This version of the code requires GSL
#ifndef USE_GSL_COSMO
#error USE_GSL_COSMO must be defined!
#else
#define LIMIT 1000
#endif

/*
 * Cosmological module for PKDGRAV.
 * N.B.  This code is being shared with skid and the I.C. generator.
 */

void csmInitialize(CSM *pcsm) {
    CSM csm;

    csm = (CSM) malloc(sizeof(struct csmContext));
    assert(csm != NULL);

    csm->val.dHubble0 = 0.0;
    csm->val.dOmega0 = 0.0;
    csm->val.dLambda = 0.0;
    csm->val.dOmegaDE = 0.0;
    csm->val.w0 = 0.0;
    csm->val.wa = 0.0;
    csm->val.dOmegaRad = 0.0;
    csm->val.dOmegab = 0.0;
    csm->val.bComove = 0;
    csm->W = gsl_integration_workspace_alloc(LIMIT);
    csm->val.classData.bClass = 0;
    csm->val.classData.achFilename[0] = 0;
    csm->val.classData.achLinSpecies[0] = 0;
    csm->val.classData.background.size = 0;
    csm->val.classData.background.a[0] = 0.;
    csm->val.classData.background.t[0] = 0.;
    csm->val.classData.background.H[0] = 0.;
    csm->val.classData.perturbations.size_a = 0;
    csm->val.classData.perturbations.size_k = 0;
    csm->val.classData.perturbations.a[0] = 0.;
    csm->val.classData.perturbations.k[0] = 0.;
    csm->val.classData.perturbations.delta_m  [0] = 0.;
    csm->val.classData.perturbations.theta_m  [0] = 0.;
    csm->val.classData.perturbations.delta_lin[0] = 0.;
    csm->val.classData.perturbations.A_s = 0.;
    csm->val.classData.perturbations.n_s = 0.;
    csm->val.classData.perturbations.alpha_s = 0.;
    csm->val.classData.perturbations.k_pivot = 0.;
    csm->classGsl.initialized = 0;
    *pcsm = csm;
    }

void csmFinish(CSM csm) {
    if (csm->val.classData.bClass){
        gsl_interp_accel_free(csm->classGsl.background.logExp2logHub_acc);
        gsl_spline_free      (csm->classGsl.background.logExp2logHub_spline);
        gsl_interp_accel_free(csm->classGsl.background.logTime2logHub_acc);
        gsl_spline_free      (csm->classGsl.background.logTime2logHub_spline);
        gsl_interp_accel_free(csm->classGsl.background.logExp2logTime_acc);
        gsl_spline_free      (csm->classGsl.background.logExp2logTime_spline);
        gsl_interp_accel_free(csm->classGsl.background.logTime2logExp_acc);
        gsl_spline_free      (csm->classGsl.background.logTime2logExp_spline);
        gsl_interp_accel_free(csm->classGsl.background.logExp2logRho_m_acc);
        gsl_spline_free      (csm->classGsl.background.logExp2logRho_m_spline);
        gsl_interp_accel_free(csm->classGsl.background.logExp2logRho_lin_acc);
        gsl_spline_free      (csm->classGsl.background.logExp2logRho_lin_spline);
        gsl_interp_accel_free(csm->classGsl.perturbations.logk2delta_m_acc);
        gsl_interp_accel_free(csm->classGsl.perturbations.loga2delta_m_acc);
        gsl_spline2d_free    (csm->classGsl.perturbations.logkloga2delta_m_spline);
        gsl_interp_accel_free(csm->classGsl.perturbations.logk2theta_m_acc);
        gsl_interp_accel_free(csm->classGsl.perturbations.loga2theta_m_acc);
        gsl_spline2d_free    (csm->classGsl.perturbations.logkloga2theta_m_spline);
        gsl_interp_accel_free(csm->classGsl.perturbations.logk2delta_lin_acc);
        gsl_interp_accel_free(csm->classGsl.perturbations.loga2delta_lin_acc);
        gsl_spline2d_free    (csm->classGsl.perturbations.logkloga2delta_lin_spline);
    }
    gsl_integration_workspace_free(csm->W);
    free(csm);
    }

void csmClassRead(CSM csm, double dBoxSize){
    size_t i, j, l, index;
    int nLinSpecies;
    hid_t file, group, attr, string_type, rhocrit_dataset, rhocrit_dataspace, memspace;
    hsize_t size_bg, size_a, size_k, count[1], offset[1], offset_out[1];
    char *matter_name, hdf5_key[128], *unit_length,
        *linSpeciesNames[10], *linSpeciesName, *LinSpeciesParsing;
    double dOmegab, dOmega0, dOmegaDE, dOmegaRad, dOmegaNu, dSpectral,dNormalization, h;
    double a, k, rho_crit[1], unit_convertion_time, unit_convertion_density;
    double *loga, *logrho_lin, *deltarho_lin, *rho_lin;

    assert(csm->val.classData.bClass);
    if (strlen(csm->val.classData.achFilename) == 0){
        fprintf(stderr, "WARNING: No achClassFilename specified\n");
        abort();
    }

    /* Parse the achLinSpecies string */
    nLinSpecies = 0;
    LinSpeciesParsing = strdup(csm->val.classData.achLinSpecies);
    while ((linSpeciesName = strsep(&LinSpeciesParsing, "+")) != NULL){
        if (strlen(linSpeciesName)){
            linSpeciesNames[nLinSpecies] = strdup(linSpeciesName);
            nLinSpecies++;
        }
    }
    free(LinSpeciesParsing);
    if (strlen(csm->val.classData.achLinSpecies) && nLinSpecies == 0){
        fprintf(stderr,
            "WARNING: The non-empty achLinSpecies = \"%s\" was specified, "
            "which does not contain any names of species\n",
            csm->val.classData.achLinSpecies);
        abort();
    }
    printf("Reading background and linear perturbations");
    if (nLinSpecies){
        printf(", including the %d linear species ", nLinSpecies);
        for (l = 0; l < nLinSpecies; l++){
            if (l > 0){
                if (l == nLinSpecies - 1)
                    printf(" and ");
                else
                    printf(", ");
            }
            printf("\"%s\"", linSpeciesNames[l]);
        }
        printf(",");
    }
    printf(" from %s\n", csm->val.classData.achFilename);

    /* Create an HDF5 ASCII string type */
    string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, H5T_VARIABLE);

    /* Open the CLASS HDF5 file for reading */
    file = H5Fopen(csm->val.classData.achFilename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) abort();

    /* Check the unit system used by the HDF5 file.
    ** Though this system is specified by "unit length", "unit time" and
    ** "unit mass", we only rely on the length. Specifically, the
    ** "unit length" is assumed to be Mpc.
    */
    attr = H5Aopen_by_name(file, "/units", "unit length", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    unit_length = (char*)malloc(128*sizeof(char));
    if (H5Aread(attr, string_type, &unit_length) < 0) abort();
    H5Aclose(attr);
    if (strcmp(unit_length, "Mpc") != 0){
        fprintf(stderr,
            "WARNING: The HDF5 file uses a unit system in which %s = 1, "
            "but we need Mpc = 1.\n",
            unit_length);
        abort();
    }
    free(unit_length);

    /* Read in the Hubble parameter h = H0/(100*km/(s*Mpc)),
    ** which we will need to convert from 1/Mpc to h/Mpc.
    */
    attr = H5Aopen_by_name(file, "/background", "h", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &h) < 0) abort();
    H5Aclose(attr);

    /* Added by MK: read in full cosmology from HDF5 */
    /*
    ** Read in the baryonic matter density parameter dOmegab
    */
    attr = H5Aopen_by_name(file, "/background", "Omega_b", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &dOmegab) < 0) abort();
    H5Aclose(attr);
    // update csm->val.dOmegab
    csm->val.dOmegab=dOmegab;

    /*
    ** Read in the total (cdm+b) matter density parameter dOmega0
    */
    attr = H5Aopen_by_name(file, "/background", "Omega_cdm+b", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &dOmega0) < 0) abort();
    H5Aclose(attr);
    // update csm->val.dOmega0
    csm->val.dOmega0=dOmega0;
    
    /*
    ** Read in the dark energy density parameter dOmegaDE
    */
    attr = H5Aopen_by_name(file, "/background", "Omega_fld", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &dOmegaDE) < 0) abort();
    H5Aclose(attr);
    // update csm->val.dOmegaDE
    csm->val.dOmegaDE=dOmegaDE;
    
    /*
    ** Read in the radiation energy density (photons) parameter dOmegaRad
    */
    attr = H5Aopen_by_name(file, "/background", "Omega_g", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &dOmegaRad) < 0) abort();
    H5Aclose(attr);
    // update csm->val.dOmegaRad
    csm->val.dOmegaRad=dOmegaRad;

    /*
    ** Read in the neutrino density paremeter
    */
    attr = H5Aopen_by_name(file, "/background", "Omega_ncdm[0]", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &dOmegaNu) < 0) abort();
    H5Aclose(attr);
    // update csm->val.dOmegaRad
    csm->val.dOmegaRad += dOmegaNu;

    
    /*
    ** Read in the spectral index n_s
    */
    attr = H5Aopen_by_name(file, "/perturbations", "n_s", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &dSpectral) < 0) abort();
    H5Aclose(attr);
    // update csm->val.dSpectral
    csm->val.dSpectral = dSpectral;

    /*
    ** Read in the spectral amplitude
    */
    attr = H5Aopen_by_name(file, "/perturbations", "A_s", H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &dNormalization) < 0) abort();
    H5Aclose(attr);
    // update csm->val.dNormalization
    csm->val.dSpectral = dNormalization;

    //printf("Omega0 read in from HDF5 file (by MK): %.14f\n", dOmega0); 

    /* The matter species "m" is really the combination "cdm+b".
    ** Here we check whether this is written as "cdm+b" or "b+cdm"
    ** in the HDF5 file.
    */
    if (H5Lexists(file, "/background/rho_cdm+b", H5P_DEFAULT)){
        matter_name = "cdm+b";
    } else if (H5Lexists(file, "/background/rho_b+cdm", H5P_DEFAULT)){
        matter_name = "b+cdm";
    } else {
        fprintf(stderr,
            "WARNING: Could not find the matter species in %s\n",
            csm->val.classData.achFilename);
        abort();
    }

    /* Read in the background, excluding densities of linear species */
    if (H5LTget_dataset_info(file, "/background/a", &size_bg, NULL, NULL) < 0) abort();
    if (size_bg > CLASS_BACKGROUND_SIZE){
        fprintf(stderr,
            "WARNING: The CLASS background size (%zu) "
            "is larger than CLASS_BACKGROUND_SIZE (%d)\n",
            (size_t)size_bg, CLASS_BACKGROUND_SIZE);
        abort();
    }
    csm->val.classData.background.size = (size_t)size_bg;
    if (H5LTread_dataset_double(file, "/background/a",
        csm->val.classData.background.a) < 0) abort();
    if (H5LTread_dataset_double(file, "/background/t",
        csm->val.classData.background.t) < 0) abort();
    if (H5LTread_dataset_double(file, "/background/H",
        csm->val.classData.background.H) < 0) abort();
    snprintf(hdf5_key, sizeof(hdf5_key), "/background/rho_%s", matter_name);
    if (H5LTread_dataset_double(file, hdf5_key,
        csm->val.classData.background.rho_m) < 0) abort();

    /* Read in the perturbations, excluding density constrasts of linear species */
    /* a */
    if (H5LTget_dataset_info(file, "/perturbations/a", &size_a, NULL, NULL) < 0) abort();
    if (size_a > CLASS_PERTURBATIONS_A_SIZE){
        fprintf(stderr,
            "WARNING: The CLASS perturbations size_a (%zu) "
            "is larger than CLASS_PERTURBATIONS_A_SIZE (%d)\n",
            (size_t)size_a, CLASS_PERTURBATIONS_A_SIZE);
        abort();
    }
    csm->val.classData.perturbations.size_a = (size_t)size_a;
    if (H5LTread_dataset_double(file, "/perturbations/a",
        csm->val.classData.perturbations.a) < 0) abort();
    /* k
    ** This is given in 1/Mpc in the HDF5 file,
    ** but we want it in h/Mpc.
    */
    if (H5LTget_dataset_info(file, "/perturbations/k", &size_k, NULL, NULL) < 0) abort();
    if (size_k > CLASS_PERTURBATIONS_K_SIZE){
        fprintf(stderr,
            "WARNING: The CLASS perturbations size_k (%zu) "
            "is larger than CLASS_PERTURBATIONS_K_SIZE (%d)\n",
            (size_t)size_k, CLASS_PERTURBATIONS_K_SIZE);
        abort();
    }
    csm->val.classData.perturbations.size_k = (size_t)size_k;
    if (H5LTread_dataset_double(file, "/perturbations/k",
        csm->val.classData.perturbations.k) < 0) abort();
    for (i = 0; i < size_k; i++){
        csm->val.classData.perturbations.k[i] /= h;
    }
    /* delta_m[a, k] */
    snprintf(hdf5_key, sizeof(hdf5_key), "/perturbations/delta_%s", matter_name);
    if (H5LTread_dataset_double(file, hdf5_key,
        csm->val.classData.perturbations.delta_m) < 0) abort();
    /* theta_m[a, k] */
    snprintf(hdf5_key, sizeof(hdf5_key), "/perturbations/theta_%s", matter_name);
    if (H5LTread_dataset_double(file, hdf5_key,
        csm->val.classData.perturbations.theta_m) < 0) abort();

    /* If any linear species are requested, read in their background
    ** densities and perturbations and add them together to form
    ** the combined "lin" species.
    **/
    if (nLinSpecies){
        /* For each linear species, we read in its background density rho
        ** and its density contrast delta. We will keep running totals in
        ** the rho_lin and deltarho_lin arrays (deltarho since delta*rho
        ** is additive, unlike delta).
        **/
        rho_lin = (double*)calloc(size_bg, sizeof(double));
        for (i = 0; i < size_bg; i++){
            rho_lin[i] = 0.;
        }
        deltarho_lin = (double*)calloc(size_a*size_k, sizeof(double));
        for (i = 0; i < size_a*size_k; i++){
            deltarho_lin[i] = 0.;
        }
        loga = (double*)calloc(size_bg, sizeof(double));
        for (i = 0; i < size_bg; i++){
            loga[i] = log(csm->val.classData.background.a[i]);
        }
        logrho_lin = (double*)calloc(size_bg, sizeof(double));
        csm->classGsl.background.logExp2logRho_lin_acc = gsl_interp_accel_alloc();
        csm->classGsl.background.logExp2logRho_lin_spline = gsl_spline_alloc(
            gsl_interp_cspline, size_bg);
        for (l = 0; l < nLinSpecies; l++){
            /* Read in the background density of the l'th linear species,
            ** overwriting the previous data.
            */
            snprintf(hdf5_key, sizeof(hdf5_key), "/background/rho_%s", linSpeciesNames[l]);
            if (H5LTread_dataset_double(file, hdf5_key,
                csm->val.classData.background.rho_lin) < 0) abort();
            /* Add the background density of the l'th linear species to the running total */
            for (i = 0; i < size_bg; i++){
                rho_lin[i] += csm->val.classData.background.rho_lin[i];
            }
            /* Construct spline over the l'th linear species rho(a) */
            for (i = 0; i < size_bg; i++){
                logrho_lin[i] = log(csm->val.classData.background.rho_lin[i]);
            }
            gsl_spline_init(csm->classGsl.background.logExp2logRho_lin_spline,
                loga, logrho_lin, size_bg);
            /* Read in the density contrast of the l'th linear species,
            ** overwriting the previous data.
            */
            snprintf(hdf5_key, sizeof(hdf5_key), "/perturbations/delta_%s", linSpeciesNames[l]);
            if (H5LTread_dataset_double(file, hdf5_key,
                csm->val.classData.perturbations.delta_lin) < 0) abort();
            /* Add the density perturbation delta*rho of the l'th linear species
            ** to the running total.
            **/
            for (i = 0; i < size_a; i++){
                a = csm->val.classData.perturbations.a[i];
                for (j = 0; j < size_k; j++){
                    index = i*size_k + j;
                    deltarho_lin[index] +=
                        csm->val.classData.perturbations.delta_lin[index]*csmRhoBar_lin(csm, a);
                }
            }
        }
        /* Store the total background density in the classData struct */
        for (i = 0; i < size_bg; i++){
            csm->val.classData.background.rho_lin[i] = rho_lin[i];
        }
        /*  Construct spline over the total rho(a) */
        for (i = 0; i < size_bg; i++){
            logrho_lin[i] = log(rho_lin[i]);
        }
        gsl_spline_init(csm->classGsl.background.logExp2logRho_lin_spline,
            loga, logrho_lin, size_bg);
        /* Store the combined density contrast in the classData struct */
        for (i = 0; i < size_a; i++){
            a = csm->val.classData.perturbations.a[i];
            for (j = 0; j < size_k; j++){
                index = i*size_k + j;
                csm->val.classData.perturbations.delta_lin[index] =
                    deltarho_lin[index]/csmRhoBar_lin(csm, a);
            }
        }
        /* Cleanup */
        for (l = 0; l < nLinSpecies; l++){
            free(linSpeciesNames[l]);
        }
        free(rho_lin);
        free(deltarho_lin);
        free(loga);
        free(logrho_lin);
        gsl_interp_accel_free(csm->classGsl.background.logExp2logRho_lin_acc);
        gsl_spline_free      (csm->classGsl.background.logExp2logRho_lin_spline);
    }

    /* Convert background variables to PKDGRAV units.
    ** Here we make use of the fact that in PKDGRAV units,
    ** H0 = sqrt(8*pi/3) and rho_crit0 = 1.
    ** Since H0 and rho_crit0 are part of the CLASS data,
    ** we have enough information to be completely agnostic
    ** abouth the units actually used in the hdf5 file.
    */
    unit_convertion_time = sqrt(8*M_PI/3)/csm->val.classData.background.H[size_bg - 1];
    count[0] = 1;
    offset[0] = size_bg - 1;
    offset_out[0] = 0;
    rhocrit_dataset = H5Dopen(file, "/background/rho_crit", H5P_DEFAULT);
    if (rhocrit_dataset < 0) abort();
    rhocrit_dataspace = H5Dget_space(rhocrit_dataset);
    memspace = H5Screate_simple(1, count, NULL);
    H5Sselect_hyperslab(rhocrit_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count, NULL);
    if (H5Dread(rhocrit_dataset, H5T_NATIVE_DOUBLE, memspace, rhocrit_dataspace,
        H5P_DEFAULT, rho_crit) < 0) abort();
    H5Dclose(rhocrit_dataset);
    H5Sclose(rhocrit_dataspace);
    H5Sclose(memspace);
    unit_convertion_density = rho_crit[0];
    for (i = 0; i < size_bg; i++){
        csm->val.classData.background.t[i] /= unit_convertion_time;
    }
    for (i = 0; i < size_bg; i++){
        csm->val.classData.background.H[i] *= unit_convertion_time;
    }
    for (i = 0; i < size_bg; i++){
        csm->val.classData.background.rho_m[i] /= unit_convertion_density;
    }
    if (nLinSpecies){
        for (i = 0; i < size_bg; i++){
            csm->val.classData.background.rho_lin[i] /= unit_convertion_density;
        }
    }

    /* The delta and theta perturbations in the hdf5 file are given as
    ** transfer functions using the CLASS convention. To convert them
    ** to actual linear fields, the square of which equals the power
    ** spectrum, we use the CLASS convention
    **   delta(k) = T_delta(k)*zeta(k),
    ** (similar for theta) with
    **   zeta(k) = pi*sqrt(2*A_s)*k^(-3/2)*(k/k_pivot)^((n_s - 1)/2)
    **             *exp(alpha_s/4*(log(k/k_pivot))^2).
    ** the comoving curvature perturbation and T_delta(k) the CLASS
    ** transfer function. This comoving curvature perturbation
    ** is implemented in csmZeta(). What we store in
    ** csm->val.classData.perturbations is then just the
    ** transfer functions. When a perturbation is looked up by either
    ** the csmDelta_m(), csmTheta_m() or csmDelta_lin() function,
    ** the corresponding transfer function is looked up via spline
    ** interpolation and zeta(k) is constructed through a call
    ** to csmZeta(). We could also multiply  the zeta(k) values on now
    ** and be done with it, but this will make the interpolations
    ** slightly poorer.
    ** The transfer functions bears units. Specifically, the delta
    ** transfer functions are unitless while the theta transfer function
    ** has units of inverse time. From the above expression of zeta(k)
    ** we see that delta(k) has units of length^(3/2). After an inverse
    ** 3D FFT, a unitless delta(\vec{x}) is obtained by multiplying by
    ** the Fourier normalization boxsize^(-3/2). We include this
    ** normalization directly on the transfer functions.
    ** ACTUALLY, it turns out that we need a factor of boxsize^(-5/2),
    ** for some reason?
    */
    group = H5Gopen(file, "/perturbations", H5P_DEFAULT); if (group < 0) abort();
    attr = H5Aopen(group, "A_s", H5P_DEFAULT); if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &csm->val.classData.perturbations.A_s) < 0)
        abort();
    H5Aclose(attr);
    attr = H5Aopen(group, "n_s", H5P_DEFAULT); if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &csm->val.classData.perturbations.n_s) < 0)
        abort();
    H5Aclose(attr);
    attr = H5Aopen(group, "alpha_s", H5P_DEFAULT); if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &csm->val.classData.perturbations.alpha_s) < 0)
        abort();
    H5Aclose(attr);
    attr = H5Aopen(group, "k_pivot", H5P_DEFAULT); if (attr < 0) abort();
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, &csm->val.classData.perturbations.k_pivot) < 0)
        abort();
    H5Aclose(attr);
    csm->val.classData.perturbations.k_pivot /= h;  /* [1/Mpc] -> [h/Mpc] */
    H5Gclose(group);
    /* delta_m[a, k] */
    for (i = 0; i < size_a*size_k; i++){
        csm->val.classData.perturbations.delta_m[i] *= pow(dBoxSize, -2.5);
    }
    /* theta_m[a, k]
    ** Here we reuse unit_convertion_time to convert the unit of theta
    ** (inverse time) to PKDGRAV units.
    ** Also, theta is the comoving divergence of the peculiar velocity.
    ** To convert to the velocities used by PKDGRAV, we have to
    ** multiply this by the scale factor.
    */
    for (i = 0; i < size_a; i++){
        a = csm->val.classData.perturbations.a[i];
        for (j = 0; j < size_k; j++){
            csm->val.classData.perturbations.theta_m[i*size_k + j] *=
                unit_convertion_time*a*pow(dBoxSize, -2.5);
        }
    }
    /* delta_lin[a, k] */
    if (nLinSpecies){
        for (i = 0; i < size_a*size_k; i++){
            csm->val.classData.perturbations.delta_lin[i] *= pow(dBoxSize, -2.5);
        }
    }

    H5Tclose(string_type);
    H5Fclose(file);
}

void csmClassGslInitialize(CSM csm){
    if (csm->classGsl.initialized){
        return;
    }
    csm->classGsl.initialized = 1;

    size_t i, size, size_k, size_a;
    double *logx, *logy, *logk, *loga;

    /* Allocate and initialize background objects */
    size = csm->val.classData.background.size;
    logx = (double*)malloc(sizeof(double)*size);
    logy = (double*)malloc(sizeof(double)*size);
    /* Exp2Hub */
    csm->classGsl.background.logExp2logHub_acc = gsl_interp_accel_alloc();
    csm->classGsl.background.logExp2logHub_spline = gsl_spline_alloc(gsl_interp_cspline, size);
    for (i = 0; i < size; i++){
        logx[i] = log(csm->val.classData.background.a[i]);
        logy[i] = log(csm->val.classData.background.H[i]);
    }
    gsl_spline_init(csm->classGsl.background.logExp2logHub_spline, logx, logy, size);
    /* Time2Hub */
    csm->classGsl.background.logTime2logHub_acc = gsl_interp_accel_alloc();
    csm->classGsl.background.logTime2logHub_spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(csm->classGsl.background.logTime2logHub_spline, logx, logy, size);
    for (i = 0; i < size; i++){
        logx[i] = log(csm->val.classData.background.t[i]);
        logy[i] = log(csm->val.classData.background.H[i]);
    }
    gsl_spline_init(csm->classGsl.background.logTime2logHub_spline, logx, logy, size);
    /* Exp2Time */
    csm->classGsl.background.logExp2logTime_acc = gsl_interp_accel_alloc();
    csm->classGsl.background.logExp2logTime_spline = gsl_spline_alloc(gsl_interp_cspline, size);
    for (i = 0; i < size; i++){
        logx[i] = log(csm->val.classData.background.a[i]);
        logy[i] = log(csm->val.classData.background.t[i]);
    }
    gsl_spline_init(csm->classGsl.background.logExp2logTime_spline, logx, logy, size);
    /* Time2Exp */
    csm->classGsl.background.logTime2logExp_acc = gsl_interp_accel_alloc();
    csm->classGsl.background.logTime2logExp_spline = gsl_spline_alloc(gsl_interp_cspline, size);
    for (i = 0; i < size; i++){
        logx[i] = log(csm->val.classData.background.t[i]);
        logy[i] = log(csm->val.classData.background.a[i]);
    }
    gsl_spline_init(csm->classGsl.background.logTime2logExp_spline, logx, logy, size);
    /* Exp2Rho_m */
    csm->classGsl.background.logExp2logRho_m_acc = gsl_interp_accel_alloc();
    csm->classGsl.background.logExp2logRho_m_spline = gsl_spline_alloc(gsl_interp_cspline, size);
    for (i = 0; i < size; i++){
        logx[i] = log(csm->val.classData.background.a[i]);
        logy[i] = log(csm->val.classData.background.rho_m[i]);
    }
    gsl_spline_init(csm->classGsl.background.logExp2logRho_m_spline, logx, logy, size);
    /* Exp2Rho_lin */
    if (strlen(csm->val.classData.achLinSpecies)){
        csm->classGsl.background.logExp2logRho_lin_acc = gsl_interp_accel_alloc();
        csm->classGsl.background.logExp2logRho_lin_spline = gsl_spline_alloc(gsl_interp_cspline, size);
        for (i = 0; i < size; i++){
            logx[i] = log(csm->val.classData.background.a[i]);
            logy[i] = log(csm->val.classData.background.rho_lin[i]);
        }
        gsl_spline_init(csm->classGsl.background.logExp2logRho_lin_spline, logx, logy, size);
    }

    free(logx);
    free(logy);

    /* Allocate and initialize perturbation objects */
    size_k = csm->val.classData.perturbations.size_k;
    size_a = csm->val.classData.perturbations.size_a;
    logk = (double*)malloc(sizeof(double)*size_k);
    loga = (double*)malloc(sizeof(double)*size_a);
    for (i = 0; i < size_k; i++){
        logk[i] = log(csm->val.classData.perturbations.k[i]);
    }
    for (i = 0; i < size_a; i++){
        loga[i] = log(csm->val.classData.perturbations.a[i]);
    }
    loga[size_a - 1] = .0;  /* Ensure high accuracy at a = 1 boundary */
    /* delta_m */
    csm->classGsl.perturbations.logk2delta_m_acc = gsl_interp_accel_alloc();
    csm->classGsl.perturbations.loga2delta_m_acc = gsl_interp_accel_alloc();
    csm->classGsl.perturbations.logkloga2delta_m_spline = gsl_spline2d_alloc(
        gsl_interp2d_bicubic, size_k, size_a);
    gsl_spline2d_init(csm->classGsl.perturbations.logkloga2delta_m_spline, logk, loga,
        csm->val.classData.perturbations.delta_m, size_k, size_a);
    /* theta_m */
    csm->classGsl.perturbations.logk2theta_m_acc = gsl_interp_accel_alloc();
    csm->classGsl.perturbations.loga2theta_m_acc = gsl_interp_accel_alloc();
    csm->classGsl.perturbations.logkloga2theta_m_spline = gsl_spline2d_alloc(
        gsl_interp2d_bicubic, size_k, size_a);
    gsl_spline2d_init(csm->classGsl.perturbations.logkloga2theta_m_spline, logk, loga,
        csm->val.classData.perturbations.theta_m, size_k, size_a);
    /* delta_lin */
    if (strlen(csm->val.classData.achLinSpecies)){
        csm->classGsl.perturbations.logk2delta_lin_acc = gsl_interp_accel_alloc();
        csm->classGsl.perturbations.loga2delta_lin_acc = gsl_interp_accel_alloc();
        csm->classGsl.perturbations.logkloga2delta_lin_spline = gsl_spline2d_alloc(
            gsl_interp2d_bicubic, size_k, size_a);
        gsl_spline2d_init(csm->classGsl.perturbations.logkloga2delta_lin_spline, logk, loga,
            csm->val.classData.perturbations.delta_lin, size_k, size_a);
    }

    free(loga);
    free(logk);

    /* Testing */
    int do_background_test = 0;
    int do_D1_test = 0;
    double a, a_PKDGRAV, a_CLASS;
    double t, t_PKDGRAV, t_CLASS;
    double H, H_PKDGRAV, H_CLASS;
    double D1_PKDGRAV, f1_PKDGRAV, D1_CLASS, f1_CLASS;
    double D2_PKDGRAV, f2_PKDGRAV, D2_CLASS, f2_CLASS;
    if (do_background_test){
        /* H(a) */
        for (i = 1; i < csm->val.classData.background.size; i++){
            /* At tabulated point */
            a = csm->val.classData.background.a[i];
            csm->val.classData.bClass = 0; H_PKDGRAV = csmExp2Hub(csm, a);
            csm->val.classData.bClass = 1; H_CLASS   = csmExp2Hub(csm, a);
            printf("TEST H(a): a = %.17e H_PKDGRAV = %.17e H_CLASS = %.17e\n",
                a, H_PKDGRAV, H_CLASS);
            /* In between tabulated points */
            if (i + 1 < csm->val.classData.background.size){
                a = 0.5*(a + csm->val.classData.background.a[i + 1]);
                csm->val.classData.bClass = 0; H_PKDGRAV = csmExp2Hub(csm, a);
                csm->val.classData.bClass = 1; H_CLASS   = csmExp2Hub(csm, a);
                printf("TEST H(a): a = %.17e H_PKDGRAV = %.17e H_CLASS = %.17e\n",
                    a, H_PKDGRAV, H_CLASS);
            }
        }
        /* a(t) */
        for (i = 1; i < csm->val.classData.background.size; i++){
            /* At tabulated point */
            t = csm->val.classData.background.t[i];
            csm->val.classData.bClass = 0; a_PKDGRAV = csmTime2Exp(csm, t);
            csm->val.classData.bClass = 1; a_CLASS   = csmTime2Exp(csm, t);
            printf("TEST a(t): t = %.17e a_PKDGRAV = %.17e a_CLASS = %.17e\n",
                t, a_PKDGRAV, a_CLASS);
            /* In between tabulated points */
            if (i + 1 < csm->val.classData.background.size){
                t = 0.5*(t + csm->val.classData.background.t[i + 1]);
                csm->val.classData.bClass = 0; a_PKDGRAV = csmTime2Exp(csm, t);
                csm->val.classData.bClass = 1; a_CLASS   = csmTime2Exp(csm, t);
                printf("TEST a(t): t = %.17e a_PKDGRAV = %.17e a_CLASS = %.17e\n",
                    t, a_PKDGRAV, a_CLASS);
            }
        }
        /* H(t) */
        for (i = 1; i < csm->val.classData.background.size; i++){
           /* At tabulated point */
           t = csm->val.classData.background.t[i];
           csm->val.classData.bClass = 0; H_PKDGRAV = csmTime2Hub(csm, t);
           csm->val.classData.bClass = 1; H_CLASS   = csmTime2Hub(csm, t);
           printf("TEST H(t): t = %.17e H_PKDGRAV = %.17e H_CLASS = %.17e\n",
               t, H_PKDGRAV, H_CLASS);
           /* In between tabulated points */
           if (i + 1 < csm->val.classData.background.size){
               t = 0.5*(t + csm->val.classData.background.t[i + 1]);
               csm->val.classData.bClass = 0; H_PKDGRAV = csmTime2Hub(csm, t);
               csm->val.classData.bClass = 1; H_CLASS   = csmTime2Hub(csm, t);
               printf("TEST H(t): t = %.17e H_PKDGRAV = %.17e H_CLASS = %.17e\n",
                   t, H_PKDGRAV, H_CLASS);
           }
        }
        /* t(a) */
        for (i = 1; i < csm->val.classData.background.size; i++){
            /* At tabulated point */
            a = csm->val.classData.background.a[i];
            csm->val.classData.bClass = 0; t_PKDGRAV = csmExp2Time(csm, a);
            csm->val.classData.bClass = 1; t_CLASS   = csmExp2Time(csm, a);
            printf("TEST t(a): a = %.17e t_PKDGRAV = %.17e t_CLASS = %.17e\n",
                a, t_PKDGRAV, t_CLASS);
            /* In between tabulated points */
            if (i + 1 < csm->val.classData.background.size){
                a = 0.5*(a + csm->val.classData.background.a[i + 1]);
                csm->val.classData.bClass = 0; t_PKDGRAV = csmExp2Time(csm, a);
                csm->val.classData.bClass = 1; t_CLASS   = csmExp2Time(csm, a);
                printf("TEST t(a): a = %.17e t_PKDGRAV = %.17e t_CLASS = %.17e\n",
                    a, t_PKDGRAV, t_CLASS);
            }
        }
    }
    if (do_background_test || do_D1_test){
        /* Growth functions */
        for (i = 1; i < csm->val.classData.background.size; i++){
            /* At tabulated point */
            a = csm->val.classData.background.a[i];
            csm->val.classData.bClass = 0; csmComoveGrowth(csm, a,
                &D1_PKDGRAV, &D2_PKDGRAV, &f1_PKDGRAV, &f2_PKDGRAV);
            csm->val.classData.bClass = 1; csmComoveGrowth(csm, a,
                &D1_CLASS, &D2_CLASS, &f1_CLASS, &f2_CLASS);
            printf("TEST D_1(a): a = %.17e D1_PKDGRAV = %.17e D1_CLASS = %.17e\n",
                a, D1_PKDGRAV, D1_CLASS);
            printf("TEST f_1(a): a = %.17e f1_PKDGRAV = %.17e f1_CLASS = %.17e\n",
                a, f1_PKDGRAV, f1_CLASS);
            printf("TEST D_2(a): a = %.17e D2_PKDGRAV = %.17e D2_CLASS = %.17e\n",
                a, D2_PKDGRAV, D2_CLASS);
            printf("TEST f_2(a): a = %.17e f2_PKDGRAV = %.17e f2_CLASS = %.17e\n",
                a, f2_PKDGRAV, f2_CLASS);
            /* In between tabulated points */
            if (i + 1 < csm->val.classData.background.size){
                a = 0.5*(a + csm->val.classData.background.a[i + 1]);
                csm->val.classData.bClass = 0; csmComoveGrowth(csm, a,
                    &D1_PKDGRAV, &D2_PKDGRAV, &f1_PKDGRAV, &f2_PKDGRAV);
                csm->val.classData.bClass = 1; csmComoveGrowth(csm, a,
                    &D1_CLASS, &D2_CLASS, &f1_CLASS, &f2_CLASS);
                printf("TEST D_1(a): a = %.17e D1_PKDGRAV = %.17e D1_CLASS = %.17e\n",
                    a, D1_PKDGRAV, D1_CLASS);
                printf("TEST f_1(a): a = %.17e f1_PKDGRAV = %.17e f1_CLASS = %.17e\n",
                    a, f1_PKDGRAV, f1_CLASS);
                printf("TEST D_2(a): a = %.17e D2_PKDGRAV = %.17e D2_CLASS = %.17e\n",
                    a, D2_PKDGRAV, D2_CLASS);
                printf("TEST f_2(a): a = %.17e f2_PKDGRAV = %.17e f2_CLASS = %.17e\n",
                    a, f2_PKDGRAV, f2_CLASS);
            }
        }
        /* Done testing */
        abort();
    }
}

#define EPSCOSMO_FUTURE 1e-6
double csmRhoBar_m(CSM csm, double a){
    assert(csm->val.classData.bClass);
    double loga = log(a);
    if (loga > .0 && loga < EPSCOSMO_FUTURE){
        /* log(a) slightly in the future. Move back to the present. */
        loga = .0;
    }
    return exp(gsl_spline_eval(csm->classGsl.background.logExp2logRho_m_spline,
        loga, csm->classGsl.background.logExp2logRho_m_acc));
}
double csmRhoBar_lin(CSM csm, double a){
    assert(csm->val.classData.bClass);
    double loga = log(a);
    if (loga > .0 && loga < EPSCOSMO_FUTURE){
        /* log(a) slightly in the future. Move back to the present. */
        loga = .0;
    }
    return exp(gsl_spline_eval(csm->classGsl.background.logExp2logRho_lin_spline,
        loga, csm->classGsl.background.logExp2logRho_lin_acc));
}
double csmDelta_m(CSM csm, double a, double k){
    assert(csm->val.classData.bClass);
    double loga = log(a);
    if (loga > .0 && loga < EPSCOSMO_FUTURE){
        /* log(a) slightly in the future. Move back to the present. */
        loga = .0;
    }
    return csmZeta(csm, k)*gsl_spline2d_eval(
        csm->classGsl.perturbations.logkloga2delta_m_spline,
        log(k), loga,
        csm->classGsl.perturbations.logk2delta_m_acc,
        csm->classGsl.perturbations.loga2delta_m_acc);
}
double csmTheta_m(CSM csm, double a, double k){
    assert(csm->val.classData.bClass);
    double loga = log(a);
    if (loga > .0 && loga < EPSCOSMO_FUTURE){
        /* log(a) slightly in the future. Move back to the present. */
        loga = .0;
    }
    return csmZeta(csm, k)*gsl_spline2d_eval(
        csm->classGsl.perturbations.logkloga2theta_m_spline,
        log(k), loga,
        csm->classGsl.perturbations.logk2theta_m_acc,
        csm->classGsl.perturbations.loga2theta_m_acc);
}
double csmDelta_lin(CSM csm, double a, double k){
    assert(csm->val.classData.bClass);
    double loga = log(a);
    if (loga > .0 && loga < EPSCOSMO_FUTURE){
        /* log(a) slightly in the future. Move back to the present. */
        loga = .0;
    }
    return csmZeta(csm, k)*gsl_spline2d_eval(
        csm->classGsl.perturbations.logkloga2delta_lin_spline,
        log(k), loga,
        csm->classGsl.perturbations.logk2delta_lin_acc,
        csm->classGsl.perturbations.loga2delta_lin_acc);
}
double csmZeta(CSM csm, double k){
    double zeta;
    double A_s     = csm->val.classData.perturbations.A_s;
    double n_s     = csm->val.classData.perturbations.n_s;
    double alpha_s = csm->val.classData.perturbations.alpha_s;
    double k_pivot = csm->val.classData.perturbations.k_pivot;
    zeta = M_PI*sqrt(2*A_s)*pow(k, -1.5)*pow(k/k_pivot, 0.5*(n_s - 1));
    if (alpha_s != 0.0)
        zeta *= exp(0.25*alpha_s*pow(log(k/k_pivot), 2));
    return zeta;
}


#define EPSCOSMO 1e-7
#define EPSCOSMO_Time2Exp 1e-20
#define EPSCOSMO_Exp2TimeIntegrate 1e-13

/*
 * ** by MK: Computes the scale factor a at radiation-matter equivalence.
 * */
double csmRadMatEquivalence(CSM csm){
    return csm->val.dOmegaRad/csm->val.dOmega0;
}

double csmTime2Hub(CSM csm,double dTime) {
    if (csm->val.classData.bClass){
        return exp(gsl_spline_eval(
            csm->classGsl.background.logTime2logHub_spline,
            log(dTime),
            csm->classGsl.background.logTime2logHub_acc));
    }

    double a = csmTime2Exp(csm,dTime);

    assert(a > 0.0);
    return csmExp2Hub(csm, a);
    }

static double Exp2Time_integrand(double ak, void * params) {
    CSM csm = (CSM)params;

    double dExp = pow(ak,2.0/3.0);
    assert (dExp > 0.0);

    return 2.0/(3.0*ak*csmExp2Hub(csm,dExp));
    }

static double Exp2TimeIntegrate(CSM csm,double dExp) {
    gsl_function F;
    F.function = &Exp2Time_integrand;
    F.params = csm;
    double result,error;
    gsl_integration_qag(&F, 0.0, pow(dExp, 1.5),
	0.0, EPSCOSMO_Exp2TimeIntegrate, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
    //printf("a=%g,\t result of Exp2TimeIntegrate = %g\n", dExp, result);
    return result;
    }

double csmExp2Time(CSM csm,double dExp) {
    if (csm->val.classData.bClass){
        return exp(gsl_spline_eval(
            csm->classGsl.background.logExp2logTime_spline,
            log(dExp),
            csm->classGsl.background.logExp2logTime_acc));
    }

    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,eta;

    if (!csm->val.bComove) {
	/*
	 ** Invalid call!
	 */
	assert(0);
	}
    
    if (csm->val.dLambda == 0.0 && csm->val.dOmegaDE == 0.0 && csm->val.dOmegaRad == 0.0) {
	if (dOmega0 == 1.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    return(2.0/(3.0*dHubble0)*pow(dExp,1.5));
	    }
	else if (dOmega0 > 1.0) {
	    assert(dHubble0 >= 0.0);
	    if (dHubble0 == 0.0) {
		B = 1.0/sqrt(dOmega0);
		eta = acos(1.0-dExp);
		return(B*(eta-sin(eta)));
		}
	    if (dExp == 0.0) return(0.0);
	    a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
	    A = 0.5*dOmega0/(dOmega0-1.0);
	    B = A*a0;
	    eta = acos(1.0-dExp/A);
	    return(B*(eta-sin(eta)));
	    }
	else if (dOmega0 > 0.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
	    A = 0.5*dOmega0/(1.0-dOmega0);
	    B = A*a0;
	    eta = acosh(dExp/A+1.0);
	    return(B*(sinh(eta)-eta));
	    }
	else if (dOmega0 == 0.0) {
	    assert(dHubble0 > 0.0);
	    if (dExp == 0.0) return(0.0);
	    return(dExp/dHubble0);
	    }
	else {
	    /*
	     ** Bad value.
	     */
	    assert(0);
	    return(0.0);
	    }
	}
    else {
	return Exp2TimeIntegrate(csm,dExp);
	}
    }

#define MAX_ITER 100

double csmTime2Exp(CSM csm,double dTime) {
    if (csm->val.classData.bClass){
        if (dTime > csm->val.classData.background.t[csm->val.classData.background.size - 1]){
            /* dTime is in the future; do linear extrapolation */
            return csm->val.classData.background.a[csm->val.classData.background.size - 1]
                + (
                    csm->val.classData.background.a[csm->val.classData.background.size - 1]
                  - csm->val.classData.background.a[csm->val.classData.background.size - 2]
                )/(
                    csm->val.classData.background.t[csm->val.classData.background.size - 1]
                  - csm->val.classData.background.t[csm->val.classData.background.size - 2]
                )*(dTime - csm->val.classData.background.t[csm->val.classData.background.size - 1]);
        }
        return exp(gsl_spline_eval(
            csm->classGsl.background.logTime2logExp_spline,
            log(dTime),
            csm->classGsl.background.logTime2logExp_acc));
    }

    double al=0,ah=1,a0,a1=1,at,a;
    double th,f,f1,h,ho;
    int j;

    if (!csm->val.bComove) return(1.0);
    else {
	assert(dTime > 0);
	th = csmExp2Time(csm,ah);
	/*
	** Search for upper bracket if needed.
	*/
	while (dTime > th) {
	    a0 = a1;
	    a1 = ah;
	    ah = a1+a0;
	    th = csmExp2Time(csm,ah);
	    }
	a = 0.5*(al+ah);
	ho = ah-al;
	h = ho;

	f = dTime - Exp2TimeIntegrate(csm,a);
	f1 = 1/(a*csmExp2Hub(csm,a));
	for (j=0;j<MAX_ITER;++j) {
	    if (a+f/f1 < al || a+f/f1 > ah || fabs(2*f) > fabs(ho*f1)) {
		/*
		** Bisection Step.
		*/
		ho = h;
		h = 0.5*(ah-al);
		a = al+h;
		/*
				printf("bisect al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
		*/
		if (a == al) return a;
		}
	    else {
		/*
		** Newton Step.
		*/
		ho = h;
		h = f/f1;
		at = a;
		a += h;
		/*
				printf("newton al:%.14g ah:%.14g a:%.14g\n",al,ah,a);
		*/
		if (a == at) return a;
		}
	    if (fabs(h) < EPSCOSMO_Time2Exp) {
		/*
				printf("converged al:%.14g ah:%.14g a:%.14g t:%.14g == %.14g\n",
				       al,ah,a,dRombergO(csm, (double (*)(void *, double)) csmCosmoTint,0.0,pow(a,1.5),EPSCOSMO*1e-1),
				       dTime);
		*/
		return a;
		}

	    if (h == ho){
	        return a;
	    }
	    f = dTime - Exp2TimeIntegrate(csm,a);
	    f1 = 1/(a*csmExp2Hub(csm,a));
	    if (f < 0) ah = a;
	    else al = a;
	    }
	assert(0);
	}
    return 0.0; /* We never reach here, but this keeps the compiler happy */
    }


double csmComoveDriftInt(CSM csm, double dIExp) {
    return -dIExp/(csmExp2Hub(csm, 1.0/dIExp));
    }
static double ComoveDrift_integrand(double diExp, void * params) {
    return csmComoveDriftInt(params,diExp);
    }

/*
 ** Make the substitution y = 1/a to integrate da/(a^2*H(a))
 */
double csmComoveKickInt(CSM csm, double dIExp) {
    return -1.0/(csmExp2Hub(csm, 1.0/dIExp));
    }

static double ComoveKick_integrand(double diExp, void * params) {
    return csmComoveKickInt(params,diExp);
    }

/*
 ** This function integrates the time dependence of the "drift"-Hamiltonian.
 */
double csmComoveDriftFac(CSM csm,double dTime,double dDelta) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,a1,a2,eta1,eta2;

    if (!csm->val.bComove) return(dDelta);
 
    else if (bClass == 0) {
        if (csm->val.dLambda == 0.0 && csm->val.dOmegaDE == 0.0 && csm->val.dOmegaRad == 0.0) {
            a1 = csmTime2Exp(csm,dTime);
            a2 = csmTime2Exp(csm,dTime+dDelta);
            if (dOmega0 == 1.0) {
                return((2.0/dHubble0)*(1.0/sqrt(a1) - 1.0/sqrt(a2)));
                }
            else if (dOmega0 > 1.0) {
                assert(dHubble0 >= 0.0);
                if (dHubble0 == 0.0) {
                    A = 1.0;
                    B = 1.0/sqrt(dOmega0);
                    }
                else {
                    a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
                    A = 0.5*dOmega0/(dOmega0-1.0);
                    B = A*a0;
                    }
                eta1 = acos(1.0-a1/A);
                eta2 = acos(1.0-a2/A);
                return(B/A/A*(1.0/tan(0.5*eta1) - 1.0/tan(0.5*eta2)));
                }
            else if (dOmega0 > 0.0) {
                assert(dHubble0 > 0.0);
                a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
                A = 0.5*dOmega0/(1.0-dOmega0);
                B = A*a0;
                eta1 = acosh(a1/A+1.0);
                eta2 = acosh(a2/A+1.0);
                return(B/A/A*(1.0/tanh(0.5*eta1) - 1.0/tanh(0.5*eta2)));
                }
            else if (dOmega0 == 0.0) {
                /*
                 ** YOU figure this one out!
                 */
                assert(0);
                return(0.0);
                }
            else {
                /*
                 ** Bad value?
                 */
                assert(0);
                return(0.0);
                }
            }
        }
    else {
	gsl_function F;
	F.function = &ComoveDrift_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
	    1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
	    0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
	return result;
	}
    }


/*
 ** This function integrates the time dependence of the "kick"-Hamiltonian.
 */
double csmComoveKickFac(CSM csm,double dTime,double dDelta) {
    double dOmega0 = csm->val.dOmega0;
    double dHubble0 = csm->val.dHubble0;
    double a0,A,B,a1,a2,eta1,eta2;

    if (!csm->val.bComove) return(dDelta);
    else if (bClass == 0) {
        if (csm->val.dLambda == 0.0 && csm->val.dOmegaDE == 0.0 && csm->val.dOmegaRad == 0.0) {
            a1 = csmTime2Exp(csm,dTime);
            a2 = csmTime2Exp(csm,dTime+dDelta);
            if (dOmega0 == 1.0) {
                return((2.0/dHubble0)*(sqrt(a2) - sqrt(a1)));
                }
            else if (dOmega0 > 1.0) {
                assert(dHubble0 >= 0.0);
                if (dHubble0 == 0.0) {
                    A = 1.0;
                    B = 1.0/sqrt(dOmega0);
                    }
                else {
                    a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
                    A = 0.5*dOmega0/(dOmega0-1.0);
                    B = A*a0;
                    }
                eta1 = acos(1.0-a1/A);
                eta2 = acos(1.0-a2/A);
                return(B/A*(eta2 - eta1));
                }
            else if (dOmega0 > 0.0) {
                assert(dHubble0 > 0.0);
                a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
                A = 0.5*dOmega0/(1.0-dOmega0);
                B = A*a0;
                eta1 = acosh(a1/A+1.0);
                eta2 = acosh(a2/A+1.0);
                return(B/A*(eta2 - eta1));
                }
            else if (dOmega0 == 0.0) {
                /*
                 ** YOU figure this one out!
                 */
                assert(0);
                return(0.0);
                }
            else {
                /*
                 ** Bad value?
                 */
                assert(0);
                return(0.0);
                }
	}
    }
    else {
	gsl_function F;
	F.function = &ComoveKick_integrand;
	F.params = csm;
	double result,error;
	gsl_integration_qag(&F, 
	    1.0/csmTime2Exp(csm, dTime), 1.0/csmTime2Exp(csm, dTime + dDelta),
	    0.0, EPSCOSMO, LIMIT, GSL_INTEG_GAUSS61, csm->W, &result, &error);
	return result;
	}
    }

double csmComoveLookbackTime2Exp(CSM csm,double dComoveTime) {
    if (!csm->val.bComove) return(1.0);
    else {
	double dExpOld = 0.0;
	double dT0 = csmExp2Time(csm, 1.0);
	double dTime = dT0 - dComoveTime;
	double dExpNew;
	int it = 0;

	if (dTime < EPSCOSMO) dTime = EPSCOSMO;
	dExpNew = csmTime2Exp(csm, dTime);
	/*
	 * Root find with Newton's method.
	 */
	do {
	    double dTimeNew = csmExp2Time(csm, dExpNew);
	    double f = dComoveTime
		       - csmComoveKickFac(csm, dTimeNew, dT0 - dTimeNew);
	    double fprime = -1.0/(dExpNew*dExpNew*csmExp2Hub(csm, dExpNew));
	    dExpOld = dExpNew;
	    dExpNew += f/fprime;
	    it++;
	    assert(it < 20);
	    }
	while (fabs(dExpNew - dExpOld)/dExpNew > EPSCOSMO);
	return dExpNew;
	}
    }

static double RK4_f1(CSM csm, double lna, double G){
    double a = exp(lna);
    return G/csmExp2Hub(csm, a);
}

static double RK4_g1(CSM csm, double lna, double D, double G){
    double a = exp(lna);
    double inva = 1./a;
    return -2.0 * G + 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * D/csmExp2Hub(csm, a);
}

// This function is in principle redundant as it is exactly the same as RK4_f1
static double RK4_f2(CSM csm, double lna, double G){
    double a = exp(lna);
    return G/csmExp2Hub(csm, a);
}

static double RK4_g2(CSM csm, double lna, double D1, double D2, double G){
    double a = exp(lna);
    double inva = 1./a;
    return -2.0 * G + 1.5 * csm->val.dOmega0 * csm->val.dHubble0*csm->val.dHubble0 * inva*inva*inva * (D2 - D1*D1)/csmExp2Hub(csm, a);
}


#define NSTEPS 1000
void csmComoveGrowth(CSM csm, double a, double *D1LPT, double *D2LPT, double *f1LPT, double *f2LPT){
    /*
    ** Variable declarations & initializations
    */
    double a_init, lna_init = log(1e-12); // ln(a)=-12 ==> a = e^(-12) ~ 0 
    double stepwidth = (log(a)- lna_init)/NSTEPS;

    // NOTICE: Storing the following quantities into data structures is by far not optimal (we actually never need the old values after the update).
    double ln_timesteps[NSTEPS+1];
    // -- 1LPT 
    double D1[NSTEPS+1]; // 1LPT Growth factor D1(a)
    double G1[NSTEPS+1]; // G1(a) = dD1(a)/dln(a) *H  ==>  Growth rate: f1(a) = G1/(H*D1) 
    // -- 2LPT
    double D2[NSTEPS+1]; // 2LPT Growth factor D2(a)
    double G2[NSTEPS+1]; // G2(a) = dD2(a)/dln(a) *H  ==>  Growth rate: f2(a) = G1/(H*D1) 

    /* 
    ** Set boundary conditions
    */
    a_init = exp(lna_init);
    D1[0] = a_init + 2.0/3.0*csmRadMatEquivalence(csm); 
    G1[0] = csmExp2Hub(csm, a_init)*a_init;

    // This is the analytical approximation
    //double AnApprox = -3. * D1[0]*D1[0]/(7. * pow(csm->val.dOmega0,1./143.));

    D2[0] = -2.0/3.0*csmRadMatEquivalence(csm)*a_init;
    G2[0] = -2.0/3.0*csmRadMatEquivalence(csm)*a_init*csmExp2Hub(csm, a_init);


    //Classical RK4 Solver
    double k0, k1, k2, k3;
    double l0, l1, l2, l3;
    double m0, m1, m2, m3;
    double n0, n1, n2, n3;

    //FILE *fp;
    //fp = fopen("GrowthFactorTable.NewBC.dat","a");

    int i; // running loop variable
    for(i=0;i<NSTEPS;i++){
        ln_timesteps[i] = lna_init + i*stepwidth;
        //fprintf(file, "%.15f, %.5f,%.20f\n", exp(ln_timesteps[i]),1.0/exp(ln_timesteps[i])-1.0, D[i]+ 0.0001977011);
 
        //RK4 step 1
        k0 = stepwidth * RK4_f1(csm, ln_timesteps[i], G1[i]);
        l0 = stepwidth * RK4_g1(csm, ln_timesteps[i], D1[i], G1[i]);

	m0 = stepwidth * RK4_f2(csm, ln_timesteps[i], G2[i]);
        n0 = stepwidth * RK4_g2(csm, ln_timesteps[i], D1[i], D2[i], G2[i]);

	//RK4 step 2
	k1 = stepwidth * RK4_f1(csm, ln_timesteps[i] + stepwidth/2.0, G1[i] + l0/2.0); 
	l1 = stepwidth * RK4_g1(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k0/2.0, G1[i] + l0/2.0);
    
	m1 = stepwidth * RK4_f2(csm, ln_timesteps[i] + stepwidth/2.0, G2[i] + n0/2.0);
        n1 = stepwidth * RK4_g2(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k0/2.0, D2[i] + m0/2.0, G2[i] + n0/2.0);

      	//RK4 step 3
	k2 = stepwidth * RK4_f1(csm, ln_timesteps[i] + stepwidth/2.0, G1[i] + l1/2.0);
  	l2 = stepwidth * RK4_g1(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k1/2.0, G1[i] + l1/2.0);

	m2 = stepwidth * RK4_f2(csm, ln_timesteps[i] + stepwidth/2.0, G2[i] + n1/2.0);
        n2 = stepwidth * RK4_g2(csm, ln_timesteps[i] + stepwidth/2.0, D1[i] + k1/2.0, D2[i] + m1/2.0, G2[i] + n1/2.0);

	//RK4 step 4
	k3 = stepwidth * RK4_f1(csm, ln_timesteps[i] + stepwidth, G1[i] + l2);
	l3 = stepwidth * RK4_g1(csm, ln_timesteps[i] + stepwidth, D1[i] + k2, G1[i] + l2);

	m3 = stepwidth * RK4_f2(csm, ln_timesteps[i] + stepwidth, G2[i] + n2);
        n3 = stepwidth * RK4_g2(csm, ln_timesteps[i] + stepwidth, D1[i] + k2, D2[i] + m2, G2[i] + n2);

	//Update
	D1[i+1] = D1[i] + (k0 + 2*k1 + 2*k2 + k3)/6.0;
	G1[i+1] = G1[i] + (l0 + 2*l1 + 2*l2 + l3)/6.0; 

        D2[i+1] = D2[i] + (m0 + 2*m1 + 2*m2 + m3)/6.0;
        G2[i+1] = G2[i] + (n0 + 2*n1 + 2*n2 + n3)/6.0; 

	//fprintf(fp, "%.20g, %.20g, %.20g\n", exp(lna_init + i*stepwidth), D1[i], D2[i]);
    }
       
    //fclose(fp);

    *D1LPT = D1[NSTEPS];
    *f1LPT = G1[NSTEPS]/(csmExp2Hub(csm,a) * *D1LPT);

    *D2LPT = D2[NSTEPS];
    *f2LPT = G2[NSTEPS]/(csmExp2Hub(csm,a) * *D2LPT);
    return;
}

