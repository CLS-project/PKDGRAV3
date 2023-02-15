#include "hydro/hydro.h"
#include "master.h"
#include "blitz/array.h"
using blitz::TinyVector;
using blitz::abs;
using blitz::any;

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
            /* If after all this there are particles without a proper density...
             * we just hope for the best and print a warning message
             */

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



static inline void densNodeOmegaE(NN *nnList, double *rpqs, float fBall,
                                  TinyVector<float,3> dr_node, int nCnt,
                                  double *omega, blitz::TinyVector<double,6> &E) {
    float fBall2_p = 4.*fBall*fBall;
    *omega = 0.0;
    for (int j=0; j<6; ++j)
        E[j] = 0.;
    for (int pk=0; pk<nCnt; pk++) {
        // As both dr vector are relative to the cell, we can do:
        blitz::TinyVector<float,3> dr = dr_node[0] - nnList[pk].dr;

        float fDist2 = blitz::dot(dr,dr);
        if (fDist2 <= fBall2_p) {
            double rpq = rpqs[pk];
            double Wpq = cubicSplineKernel(rpq, fBall);

            *omega += Wpq;


            E[XX] += dr[0]*dr[0]*Wpq;
            E[YY] += dr[1]*dr[1]*Wpq;
            E[ZZ] += dr[2]*dr[2]*Wpq;

            E[XY] += dr[1]*dr[0]*Wpq;
            E[XZ] += dr[2]*dr[0]*Wpq;
            E[YZ] += dr[1]*dr[2]*Wpq;

        }
    }
}


static inline double densNodeNcondB(PKD pkd, particleStore::ParticleReference &p,
                                    blitz::TinyVector<double,6> &E, double omega) {
    blitz::TinyVector<double,6> B;

    // Normalize the matrix
    for (int j=0; j<6; ++j) {
        E[j] /= omega;
    }

    inverseMatrix(E.data(), B.data());
    double Ncond = conditionNumber(E.data(), B.data());
    assert(Ncond==Ncond);

    if (p.is_gas()) {
        // We can already set this here, so it can be skipped in
        // hydroGradients
        auto &sph = p.sph();
        sph.Ncond = Ncond;
        sph.B[XX] = B[XX];
        sph.B[YY] = B[YY];
        sph.B[ZZ] = B[ZZ];
        sph.B[XY] = B[XY];
        sph.B[XZ] = B[XZ];
        sph.B[YZ] = B[YZ];
    }

    return Ncond;
}



void hydroDensity_node(PKD pkd, SMF *smf, Bound bnd_node, const std::vector<PARTICLE *> &sinks, NN *nnList,int nCnt) {
    for (auto &P : sinks) {
        auto partj = pkd->particles[P];
#ifdef OPTIM_UNION_EXTRAFIELDS
        double *omega = NULL;
        omega = partj.is_gas()  ? &(partj.sph().omega) : omega;
        omega = partj.is_star() ? &(partj.star().omega) : omega;
        omega = partj.is_bh()   ? &(partj.BH().omega) : omega;
#else
        // Assuming *only* stars and gas
        double *omega = &(partj.sph().omega);
#endif
        auto r = partj.position();
        TinyVector<float,3> dr_node = bnd_node.center() - r;

        // The sqrt can be computed just once here, with higher probability
        // of being vectorized
        double rpqs[nCnt];
        for (int pk=0; pk<nCnt; pk++) {
            blitz::TinyVector<double,3> dr = dr_node - nnList[pk].dr;
            rpqs[pk] = sqrt(blitz::dot(dr,dr));
        }

        int niter = 0;
        float Neff = smf->nSmooth;
        do {
            float ph = partj.ball();
            blitz::TinyVector<double,6> E;

            densNodeOmegaE(nnList, rpqs, ph, dr_node, nCnt,omega, E);

            // Check if it has converged
            double c = 4.*M_PI/3. * (*omega) *ph*ph*ph*8.;
            if ((fabs(c-Neff) < smf->dNeighborsStd) ) {
                // Check if the converged density has a low enough condition number

                double Ncond = densNodeNcondB(pkd, partj, E, *omega);


                if (Ncond > 100) {
                    // In some configurations (outer particles in a isolated galaxy)
                    // the density could have converged but with a very high Ncond
                    // due to anisotropy.
                    // For those cases, we impose a maximum effective ngb number
                    if (Neff>200.) {
                        partj.set_marked(false);
                        printf("WARNING Neff %e Ncond %e \n", Neff, Ncond);
                    }
                    else {
                        Neff *= 1.2;
                        niter = 0;
                        continue;
                    }
                }

                partj.set_marked(false);
                partj.set_density(partj.mass()*(*omega));
            }
            else {
                float newBall;

                newBall = (c!=0.0) ? ph * pow(  Neff/c,0.3333333333) : ph*4.0;
                if (newBall > 4.0*ph) newBall = ph*4.0;

                partj.ball() = 0.5*(newBall+ph);
                //printf("Setting new fBall %e %e %e \n", c, ph, partj.ball());


                if (newBall>ph) {
                    float ph = partj.ball();
                    // We check that the proposed ball is enclosed within the node search region
                    if (any(abs(dr_node) + ph > bnd_node.apothem())) {
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
            if (niter>1000 && partj.marked()) {
                if (c > Neff) {
                    partj.set_marked(false);
                    partj.ball() = ph;
                    densNodeNcondB(pkd, partj, E, *omega);
                    partj.set_density(partj.mass()*(*omega));
                    printf("WARNING Neff %e c %e \n", Neff, c);
                }
            }

        } while (partj.marked());

        // After a converged fBall is obtained, we limit fBall if needed
        if (smf->dhMinOverSoft > 0.) {
            if (partj.ball() < smf->dhMinOverSoft*partj.soft()) {
                float newBall = smf->dhMinOverSoft*partj.soft();
                partj.ball() = newBall;

                blitz::TinyVector<double,6> E;
                densNodeOmegaE(nnList, rpqs, newBall, dr_node, nCnt,omega, E);

                densNodeNcondB(pkd, partj, E, *omega);

                partj.set_density(partj.mass()*(*omega));
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
void hydroDensity(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    double ph, rpq, hpq, c;
    int i;

    /* Particle p data */
    auto &sph = p.sph();
    ph = fBall;

    /* Compute the \omega(x_i) normalization factor */
    sph.omega = 0.0;
    ph = p.ball();
#ifdef OPTIM_DENSITY_REITER
    float maxr = 0.0;
#endif //OPTIM_DENSITY_REITER
    for (i=0; i<nSmooth; ++i) {
#ifdef OPTIM_DENSITY_REITER
        if (nnList[i].fDist2> maxr) maxr =nnList[i].fDist2;
#else

        rpq = sqrt(nnList[i].fDist2);
        hpq = ph;

        sph.omega += cubicSplineKernel(rpq, hpq);
#endif //OPTIM_DENSITY_REITER
    }

    /* If we are using a iterative procedure for computing the smoothing length:
     *    - Particles marked are those which still needs iterations
     *    - Particles not marked are those with a correct smoothing length
     */
#ifdef OPTIM_DENSITY_REITER
    maxr = sqrt(maxr);
    while (p.is_marked()) {
        sph.omega=0.0;
        ph = p->ball();
        for (i=0; i<nSmooth; ++i) {

            rpq = sqrt(nnList[i].fDist2);
            hpq = ph;

            sph.omega += cubicSplineKernel(rpq, hpq);
        }
#else
    if (smf->bIterativeSmoothingLength && p.marked()) {
#endif //OPTIM_DENSITY_REITER

        c = 4.*M_PI/3. * sph.omega *ph*ph*ph*8.;
        if (fabs(c-smf->nSmooth) < smf->dNeighborsStd) {
            p.set_marked(false);
        }
        else {
            float newBall;
            newBall = ph * pow(  smf->nSmooth/c,0.3333333333);
            //   if (nSmooth <= 1) newBall *= 2.*fBall;

            p.ball() = 0.5*(newBall+ph);
            //sph.fLastBall = ph;

#ifdef OPTIM_DENSITY_REITER
            // If the suggested new radius does not enclose all our neighbors,
            // we need to reiterate
            if (p.ball()>maxr) break;
#endif

        }
    }


    /* We compute the density making use of Eq. 27 Hopkins 2015 */
    p.set_density(p.mass()*sph.omega);
}
