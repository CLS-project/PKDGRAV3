#include "hydro/hydro.h"
#include "master.h"

void MSR::ComputeSmoothing(double dTime, double dDelta) {
    int nSmoothed = 1, it=0, maxit = 100;

    printf("Computing density... \n");
    TimerStart(TIMER_DENSITY);
    if (param.bIterativeSmoothingLength) {
#ifdef OPTIM_AVOID_IS_ACTIVE
        SelActives();
#else
        SelAll(); // We set all particles as "not converged"
#endif

        // Currently implemented but not used/tested
        //pstPredictSmoothing(msr->pst,&in,sizeof(in),NULL,NULL);

        while (nSmoothed>0 && it <= maxit) {
            // I think this is no longer used
            //msrSetFirstHydroLoop(msr, 1);
            // 1-> we care if the particle is marked
            // 0-> we dont
#ifdef OPTIM_SMOOTH_NODE
            nSmoothed = ReSmoothNode(dTime, dDelta, SMX_HYDRO_DENSITY,0);
#else
            nSmoothed = ReSmooth(dTime, dDelta, SMX_HYDRO_DENSITY,0);
#endif
            //msrSetFirstHydroLoop(msr, 0);
            it++;
        }
        if (nSmoothed >0) {
            // If we can not converge, at least make sure that we update the
            // density and other variables with the latest fBall
#ifdef OPTIM_SMOOTH_NODE
            nSmoothed = ReSmoothNode(dTime, dDelta, SMX_HYDRO_DENSITY_FINAL,0);
#else
            nSmoothed = ReSmooth(dTime, dDelta, SMX_HYDRO_DENSITY_FINAL,0);
#endif

            printf("Smoothing length did not converge for %d particles\n", nSmoothed);
        }

        TimerStop(TIMER_DENSITY);
        double dsec = TimerGet(TIMER_DENSITY);
        printf("Computing h took %d iterations and %.5f seconds \n", it, dsec);
    }
    else {
        bUpdateBall = 1;
        Smooth(dTime, dDelta, SMX_HYDRO_DENSITY,0,param.nSmooth);
        bUpdateBall = 0;
    }
}



static inline void densNodeOmegaE(NN *nnList, float fBall,
                                  float dx_node, float dy_node, float dz_node, int nCnt,
                                  double *omega, double *E, int *nSmooth) {
    float fBall2_p = 4.*fBall*fBall;
    *omega = 0.0;
    *nSmooth = 0.0;
    for (int j=0; j<6; ++j)
        E[j] = 0.;
    for (int pk=0; pk<nCnt; pk++) {
        // As both dr vector are relative to the cell, we can do:
        const float dx = dx_node - nnList[pk].dx;
        const float dy = dy_node - nnList[pk].dy;
        const float dz = dz_node - nnList[pk].dz;

        const float fDist2 = dx*dx + dy*dy + dz*dz;
        if (fDist2 <= fBall2_p) {
            const double rpq = sqrt(fDist2);
            const double Wpq = cubicSplineKernel(rpq, fBall);

            *omega += Wpq;


            E[XX] += dx*dx*Wpq;
            E[YY] += dy*dy*Wpq;
            E[ZZ] += dz*dz*Wpq;

            E[XY] += dy*dx*Wpq;
            E[XZ] += dz*dx*Wpq;
            E[YZ] += dy*dz*Wpq;

            *nSmooth += 1;

        }
    }
}


static inline double densNodeNcondB(PKD pkd, PARTICLE *p,
                                    double *E, double omega) {
    double B[6];

    // Normalize the matrix
    for (int j=0; j<6; ++j) {
        E[j] /= omega;
    }

    inverseMatrix(E, B);
    double Ncond = conditionNumber(E, B);
    assert(Ncond==Ncond);

    if (pkdIsGas(pkd,p)) {
        // We can already set this here, so it can be skipped in
        // hydroGradients
        SPHFIELDS *psph = pkdSph(pkd,p);
        psph->Ncond = Ncond;
        psph->B[XX] = B[XX];
        psph->B[YY] = B[YY];
        psph->B[ZZ] = B[ZZ];
        psph->B[XY] = B[XY];
        psph->B[XZ] = B[XZ];
        psph->B[YZ] = B[YZ];
    }

    return Ncond;
}

// Compute the density and derived variables simply given the fBall,
// without trying to converge to the correct value
void hydroDensityFinal(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
#ifdef OPTIM_UNION_EXTRAFIELDS
    double *omega = NULL;
    omega = pkdIsGas(pkd,p)  ? &(pkdSph(pkd,p)->omega) : omega;
    omega = pkdIsStar(pkd,p) ? &(pkdStar(pkd,p)->omega) : omega;
    omega = pkdIsBH(pkd,p)   ? &(pkdBH(pkd,p)->omega) : omega;
#else
    // Assuming *only* stars and gas
    double *omega = &(pkdSph(pkd,p)->omega);
#endif

    const double ph = pkdBall(pkd,p);

    *omega = 0.0;
    for (int i=0; i<nSmooth; ++i) {
        const float dx = -nnList[i].dx;
        const float dy = -nnList[i].dy;
        const float dz = -nnList[i].dz;
        const float fDist2 = dx*dx + dy*dy + dz*dz;
        const double rpq = sqrt(fDist2);


        *omega += cubicSplineKernel(rpq, ph);
    }


    if (pkdIsGas(pkd,p)) {
        double E[6];
        for (int i=0; i<6; ++i) {
            E[i] = 0.0;
        }

        for (int i=0; i<nSmooth; ++i) {
            const float dx = -nnList[i].dx;
            const float dy = -nnList[i].dy;
            const float dz = -nnList[i].dz;
            const float fDist2 = dx*dx + dy*dy + dz*dz;
            const double rpq = sqrt(fDist2);

            const double Wpq = cubicSplineKernel(rpq, ph);

            E[XX] += dx*dx*Wpq;
            E[YY] += dy*dy*Wpq;
            E[ZZ] += dz*dz*Wpq;

            E[XY] += dy*dx*Wpq;
            E[XZ] += dz*dx*Wpq;
            E[YZ] += dy*dz*Wpq;
        }

        /* Normalize the matrix */
        for (int i=0; i<6; ++i) {
            E[i] /= *omega;
        }
        inverseMatrix(E, pkdSph(pkd,p)->B);
    }

    pkdSetDensity(pkd, p, pkdMass(pkd,p)*(*omega));
    //printf("%" PRIu64 " final iteration omega %e\t density %e\t fBall %e\t nn %d \n", *pkdParticleID(pkd,p), *omega, pkdDensity(pkd,p), ph, nSmooth);
}


void hydroDensity_node(PKD pkd, SMF *smf, BND bnd_node, PARTICLE **sinks, NN *nnList,
                       int nCnt_own, int nCnt) {
    for (int pj=0; pj<nCnt_own; pj++) {
        PARTICLE *partj = sinks[pj];
#ifdef OPTIM_UNION_EXTRAFIELDS
        double *omega = NULL;
        omega = pkdIsGas(pkd,partj)  ? &(pkdSph(pkd,partj)->omega) : omega;
        omega = pkdIsStar(pkd,partj) ? &(pkdStar(pkd,partj)->omega) : omega;
        omega = pkdIsBH(pkd,partj)   ? &(pkdBH(pkd,partj)->omega) : omega;
#else
        // Assuming *only* stars and gas
        double *omega = &(pkdSph(pkd,partj)->omega);
#endif

        float dx_node = -pkdPos(pkd,partj,0)+bnd_node.fCenter[0];
        float dy_node = -pkdPos(pkd,partj,1)+bnd_node.fCenter[1];
        float dz_node = -pkdPos(pkd,partj,2)+bnd_node.fCenter[2];

        int niter = 0;
        int nSmooth;

        float Neff = smf->nSmooth;
        if (*omega < 0)
            Neff = - *omega;

        do {
            float ph = pkdBall(pkd,partj);
            double E[6];

            densNodeOmegaE(nnList, ph, dx_node, dy_node, dz_node,
                           nCnt,omega, E, &nSmooth);

            // Check if it has converged
            double c = 4.*M_PI/3. * (*omega) *ph*ph*ph*8.;
            //printf("%" PRIu64 " %d %d %e %e \n", *pkdParticleID(pkd,partj), nSmooth, nCnt, ph, Neff);
            if ((fabs(c-Neff) < smf->dNeighborsStd) ) {
                // Check if the converged density has a low enough condition number

                double Ncond = densNodeNcondB(pkd, partj, E, *omega);


                if (Ncond > 100) {
                    // In some configurations (outer particles in a isolated galaxy)
                    // the density could have converged but with a very high Ncond
                    // due to anisotropy.
                    // For those cases, we impose a maximum effective ngb number
                    if (Neff>200.) {
                        partj->bMarked=0;
                        printf("WARNING %d Maximum Neff reached: %e ; Ncond %e \n", nSmooth, Neff, Ncond);
                    }
                    else {
                        Neff *= 1.2;
                        niter = 0;
                        continue;
                    }
                }
                if (nSmooth < 20) {
                    // If we converge to a fBall that encloses too few neighbours,
                    // the effective number of nn is increased
                    Neff *= 1.2;
                    niter = 0;
                    continue;
                }

                partj->bMarked = 0;
                pkdSetDensity(pkd, partj, pkdMass(pkd,partj)*(*omega));
            }
            else {
                float newBall;

                newBall = (c!=0.0) ? ph * pow(  Neff/c,0.3333333333) : ph*4.0;
                if (newBall > 4.0*ph) newBall = ph*4.0;

                pkdSetBall(pkd,partj, 0.5*(newBall+ph));
                //printf("Setting new fBall %e %e %e \n", c, ph, pkdBall(pkd,partj));


                if (newBall>ph) {
                    float ph2 = 2.*pkdBall(pkd,partj);
                    // We check that the proposed ball is enclosed within the
                    // node search region
                    if (    (fabs(dx_node) + ph2 > bnd_node.fMax[0])||
                            (fabs(dy_node) + ph2 > bnd_node.fMax[1])||
                            (fabs(dz_node) + ph2 > bnd_node.fMax[2])) {
                        //printf("%" PRIu64 " \t nn %d \t omega %e \t Neff %e \t ph %e \t newBall %e \t fBall %e \n", *pkdParticleID(pkd,partj), nSmooth, *omega, Neff, ph, newBall, pkdBall(pkd,partj));
                        *omega = -Neff;
                        break;
                    }
                }

            }
            niter++;


            // At this point, we probably have a particle with plenty of
            //  neighbours, and a small increase/decrease in the radius causes
            //  omega to fluctuate around the desired value.
            //
            // In this cases, we stop iterating, making sure that we have the
            //  upper value, which should have more than the expected number of
            //  neighbours
            if (niter>1000 && partj->bMarked) {
                if (c > Neff) {
                    partj->bMarked = 0;
                    pkdSetBall(pkd,partj, ph);
                    densNodeNcondB(pkd, partj, E, *omega);
                    pkdSetDensity(pkd, partj, pkdMass(pkd,partj)*(*omega));
                    printf("WARNING %d Maximum iterations reached Neff %e c %e \n", pkdSpecies(pkd,partj), Neff, c);
                }
            }

        } while (partj->bMarked);

        // After a converged fBall is obtained, we limit fBall if needed
        if (smf->dhMinOverSoft > 0.) {
            if (pkdBall(pkd,partj) < smf->dhMinOverSoft*pkdSoft(pkd,partj)) {
                float newBall = smf->dhMinOverSoft*pkdSoft(pkd,partj);
                pkdSetBall(pkd, partj, newBall);

                double E[6];
                densNodeOmegaE(nnList, newBall, dx_node, dy_node, dz_node,
                               nCnt,omega, E, &nSmooth);

                densNodeNcondB(pkd, partj, E, *omega);

                pkdSetDensity(pkd, partj, pkdMass(pkd,partj)*(*omega));
            }
        }

    }
}



/* This function is now deprecated and may be removed in future versions of
 * the code.
 *
 * The supported density computation is that performed by hydroDensity_node,
 * which is automatically selected when using the OPTIM_SMOOTH_NODE flag
 */
void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    SPHFIELDS *psph;
    double ph, rpq, hpq, c;
    int i;

    /* Particle p data */
    psph = pkdSph(pkd,p);
    ph = fBall;

    /* Compute the \omega(x_i) normalization factor */
    psph->omega = 0.0;
    ph = pkdBall(pkd,p);
#ifdef OPTIM_DENSITY_REITER
    float maxr = 0.0;
#endif //OPTIM_DENSITY_REITER
    for (i=0; i<nSmooth; ++i) {
#ifdef OPTIM_DENSITY_REITER
        if (nnList[i].fDist2> maxr) maxr =nnList[i].fDist2;
#else

        rpq = sqrt(nnList[i].fDist2);
        hpq = ph;

        psph->omega += cubicSplineKernel(rpq, hpq);
#endif //OPTIM_DENSITY_REITER
    }

    /* If we are using a iterative procedure for computing the smoothing length:
     *    - Particles marked are those which still needs iterations
     *    - Particles not marked are those with a correct smoothing length
     */
#ifdef OPTIM_DENSITY_REITER
    maxr = sqrt(maxr);
    while (p->bMarked) {
        psph->omega=0.0;
        ph = pkdBall(pkd,p);
        for (i=0; i<nSmooth; ++i) {

            rpq = sqrt(nnList[i].fDist2);
            hpq = ph;

            psph->omega += cubicSplineKernel(rpq, hpq);
        }
#else
    if (smf->bIterativeSmoothingLength && p->bMarked) {
#endif //OPTIM_DENSITY_REITER

        c = 4.*M_PI/3. * psph->omega *ph*ph*ph*8.;
        if (fabs(c-smf->nSmooth) < smf->dNeighborsStd) {
            p->bMarked = 0;
        }
        else {
            float newBall;
            newBall = ph * pow(  smf->nSmooth/c,0.3333333333);
            //   if (nSmooth <= 1) newBall *= 2.*fBall;

            pkdSetBall(pkd,p, 0.5*(newBall+ph));
            //psph->fLastBall = ph;

#ifdef OPTIM_DENSITY_REITER
            // If the suggested new radius does not enclose all our neighbors,
            // we need to reiterate
            if (pkdBall(pkd,p)>maxr) break;
#endif

        }
    }


    /* We compute the density making use of Eq. 27 Hopkins 2015 */
    pkdSetDensity(pkd,p, pkdMass(pkd,p)*psph->omega);
}
