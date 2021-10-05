#include "hydro/hydro.h"

void msrComputeSmoothing(MSR msr,double dTime)
{
    int nSmoothed = 1, it=0, maxit = 100;

    printf("Computing density... \n");
    msrTimerStart(msr, TIMER_DENSITY);
    if (msr->param.bIterativeSmoothingLength) {
#ifdef OPTIM_AVOID_IS_ACTIVE
        msrSelActive(msr);
#else
        msrSelAll(msr); // We set all particles as "not converged"
#endif

        // Currently implemented but not used/tested
        //pstPredictSmoothing(msr->pst,&in,sizeof(in),NULL,NULL);

        while (nSmoothed>0 && it <= maxit) {
            msrSetFirstHydroLoop(msr, 1); // 1-> we care if the particle is marked
            // 0-> we dont
#ifdef OPTIM_SMOOTH_NODE
            nSmoothed = msrReSmoothNode(msr,dTime,SMX_FIRSTHYDROLOOP,0,0);
#else
            nSmoothed = msrReSmooth(msr,dTime,SMX_FIRSTHYDROLOOP,0,0);
#endif
            msrSetFirstHydroLoop(msr, 0);
            it++;
        }
        if (nSmoothed >0) {
            /* If after all this there are particles without a proper density...
             * we just hope for the best and print a warning message
             */

            printf("Smoothing length did not converge for %d particles\n", nSmoothed);
        }

        msrTimerStop(msr, TIMER_DENSITY);
        double dsec = msrTimerGet(msr, TIMER_DENSITY);
        printf("Computing h took %d iterations and %.5f seconds \n", it, dsec);
    } else {
        msrSetFirstHydroLoop(msr, 1);
        msrSmooth(msr,dTime,SMX_FIRSTHYDROLOOP,0,msr->param.nSmooth);
        msrSetFirstHydroLoop(msr, 0);
    }
}

void hydroDensity_node(PKD pkd, BND bnd_node, PARTICLE **sinks, NN *nnList,
                       int nCnt_own, int nCnt)
{
    for (int pj=0; pj<nCnt_own; pj++) {
        PARTICLE * partj = sinks[pj];
#ifdef OPTIM_UNION_EXTRAFIELDS
        double *omega = NULL;
        omega = pkdIsGas(pkd,partj)  ? &(pkdSph(pkd,partj)->omega) : omega;
        omega = pkdIsStar(pkd,partj) ? &(pkdStar(pkd,partj)->omega) : omega;
        omega = pkdIsBH(pkd,partj)   ? &(pkdBH(pkd,partj)->omega) : omega;
#else
        // Assuming *only* stars and gas
        double* omega = &(pkdSph(pkd,partj)->omega);
#endif

        float dx_node = -pkdPos(pkd,partj,0)+bnd_node.fCenter[0];
        float dy_node = -pkdPos(pkd,partj,1)+bnd_node.fCenter[1];
        float dz_node = -pkdPos(pkd,partj,2)+bnd_node.fCenter[2];

        int niter = 0;
        float Neff = pkd->param.nSmooth;
        do {
            float ph = pkdBall(pkd,partj);
            float fBall2_p = 4.*ph*ph;
            int nCnt_p = 0;
            *omega = 0.0;
            double E[6], B[6];
            for (int j=0; j<6; ++j) E[j] = 0.;
            for (int pk=0; pk<nCnt; pk++) {
                // As both dr vector are relative to the cell, we can do:
                float dx = dx_node - nnList[pk].dx;
                float dy = dy_node - nnList[pk].dy;
                float dz = dz_node - nnList[pk].dz;

                float fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2_p) {
                    double rpq = sqrt(fDist2);
                    double Wpq = cubicSplineKernel(rpq, ph);

                    *omega += Wpq;


                    E[XX] += dx*dx*Wpq;
                    E[YY] += dy*dy*Wpq;
                    E[ZZ] += dz*dz*Wpq;

                    E[XY] += dy*dx*Wpq;
                    E[XZ] += dz*dx*Wpq;
                    E[YZ] += dy*dz*Wpq;

                    nCnt_p++;
                }
            }



            // Check if it has converged
            double c = 4.*M_PI/3. * (*omega) *ph*ph*ph*8.;
            if ((fabs(c-Neff) < pkd->param.dNeighborsStd0) ) {
                // Check if the converged density has a low enough condition number

                // Normalize the matrix
                for (int j=0; j<6; ++j) {
                    E[j] /= *omega;
                }

                inverseMatrix(E, B);
                double Ncond = conditionNumber(E, B);
                assert(Ncond==Ncond);
                // TODO: Assign this here and void computing it in hydroGradients,
                // same for B
                //if (pkdIsSph(pkd,p)) psph->Ncond = Ncond;


                if (Ncond > 100) {
                    // In some configurations (outer particles in a isolated galaxy)
                    // the density could have converged but with a very high Ncond
                    // due to anisotropy.
                    // For those cases, we impose a maximum effective ngb number
                    if (Neff>200.) {
                        partj->bMarked=0;
                        printf("WARNING Neff %e Ncond %e \n", Neff, Ncond);
                    } else {
                        Neff *= 1.2;
                        niter = 0;
                        continue;
                    }
                }

                partj->bMarked = 0;
                pkdSetDensity(pkd, partj, pkdMass(pkd,partj)*(*omega));
            } else {
                float newBall;

                newBall = (c!=0.0) ? ph * pow(  Neff/c,0.3333333333) : ph*4.0;
                if (newBall > 4.0*ph) newBall = ph*4.0;

                pkdSetBall(pkd,partj, 0.5*(newBall+ph));
                //printf("Setting new fBall %e %e %e \n", c, ph, pkdBall(pkd,partj));


                if (newBall>ph) {
                    float ph = pkdBall(pkd,partj);
                    // We check that the proposed ball is enclosed within the
                    // node search region
                    if((fabs(dx_node) + ph > bnd_node.fMax[0])||
                       (fabs(dy_node) + ph > bnd_node.fMax[1])||
                       (fabs(dz_node) + ph > bnd_node.fMax[2])) {
                        // Removed this, see notes 29/10/20
                        //nSmoothed-=1; // We explicitly say that this particle was
                        //not succesfully smoothed
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
            if (niter>1000) {
                if (c > Neff) {
                    partj->bMarked = 0;
                    pkdSetBall(pkd,partj, ph);
                    printf("WARNING Neff %e c %e \n", Neff, c);
                }
            }

        } while(partj->bMarked);

        // After a converged fBall is obtained, we limit fBall if needed
        if (pkd->param.dhMinOverSoft > 0.) {
            float newBall = MAX(pkdBall(pkd,partj),
                                pkd->param.dhMinOverSoft*pkdSoft(pkd,partj));
            pkdSetBall(pkd, partj, newBall);
            float fBall2_p = 4.*newBall*newBall;
            *omega = 0.0;
            for (int pk=0; pk<nCnt; pk++) {
                float dx = dx_node - nnList[pk].dx;
                float dy = dy_node - nnList[pk].dy;
                float dz = dz_node - nnList[pk].dz;

                float fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2_p) {
                    double rpq = sqrt(fDist2);
                    double Wpq = cubicSplineKernel(rpq, newBall);

                    *omega += Wpq;
                }
            }
            pkdSetDensity(pkd, partj, pkdMass(pkd,partj)*(*omega));
        }

    }
}



/* This function is now deprecated and may be removed in future versions of
 * the code.
 *
 * The supported density computation is that performed by hydroDensity_node,
 * which is automatically selected when using the OPTIM_SMOOTH_NODE flag
 */
void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf)
{
    PKD pkd = smf->pkd;
    SPHFIELDS *psph;
    double ph, rpq, hpq, c;
    int i;

    /* Particle p data */
    psph = pkdSph(pkd,p);
    ph = fBall;

#ifndef FIXED_NSMOOTH_STRICT
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
    if (pkd->param.bIterativeSmoothingLength && p->bMarked) {
#endif //OPTIM_DENSITY_REITER

#ifndef FIXED_NSMOOTH_RELAXED
        c = 4.*M_PI/3. * psph->omega *ph*ph*ph*8.;
#else
        c = nSmooth;
#endif
        if (fabs(c-pkd->param.nSmooth) < pkd->param.dNeighborsStd0) {
            p->bMarked = 0;
        } else {
            float newBall;
            newBall = ph * pow(  pkd->param.nSmooth/c,0.3333333333);
            //   if (nSmooth <= 1) newBall *= 2.*fBall;

            pkdSetBall(pkd,p, 0.5*(newBall+ph));
            //psph->fLastBall = ph;

#ifdef OPTIM_DENSITY_REITER
            // If the suggested new radius does not enclose all our neighbors,
            // we need to reiterate
            if (pkdBall(pkd,p)>maxr) break;
#endif

#else // FIXED_NSMOOTH_STRICT
    double minR2, lastMin;
    if (pkd->param.bIterativeSmoothingLength && p->bMarked) {
        c = nSmooth;
        // Check if we have enough neighbors, otherwise increse fBall
        if (c <  pkd->param.nSmooth) {
            pkdSetBall(pkd,p,fBall * pow(  1.2*pkd->param.nSmooth/c,0.3333333333));
        } else if (c >= pkd->param.nSmooth) {
            // Now we look for the distance to the nSmooth-th neighbor
            minR2 = HUGE_VAL;
            lastMin = 0.;

            for (int n=0; n<pkd->param.nSmooth-2; n++) {
                minR2 = HUGE_VAL;
                for (i=0; i<nSmooth; i++) {
                    if (nnList[i].fDist2 < minR2)
                        if (nnList[i].fDist2 > lastMin)
                            minR2 = nnList[i].fDist2;
                }
                lastMin = minR2;

            }

            pkdSetBall(pkd,p, 0.5*sqrt(lastMin));

            p->bMarked=0;
#endif
        }
    }

#ifdef FIXED_NSMOOTH_STRICT
    psph->omega = 0.0;
    ph = pkdBall(pkd,p);
    for (i=0; i<nSmooth; ++i) {

        rpq = sqrt(nnList[i].fDist2);
        hpq = ph;

        psph->omega += cubicSplineKernel(rpq, hpq);
    }
#endif

    /* We compute the density making use of Eq. 27 Hopkins 2015 */
    pkdSetDensity(pkd,p, pkdMass(pkd,p)*psph->omega);
}
