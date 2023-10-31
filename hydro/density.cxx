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
    if (parameters.get_bIterativeSmoothingLength()) {
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
        Smooth(dTime, dDelta, SMX_HYDRO_DENSITY,0,parameters.get_nSmooth());
        bUpdateBall = 0;
    }
}

static inline void densNodeOmegaE(NN *nnList, float pH, TinyVector<double,3> dr_node, int nCnt,
                                  double *omega, TinyVector<double,6> &E, int *nSmooth) {
    const float pH2 = pH * pH;
    *omega = 0.0;
    *nSmooth = 0.0;
    E = 0.0;
    for (auto pk = 0; pk < nCnt; ++pk) {
        // As both dr vector are relative to the cell, we can do:
        const TinyVector<double,3> dr{dr_node - nnList[pk].dr};

        const auto fDist2 = blitz::dot(dr,dr);
        if (fDist2 <= pH2) {
            const double rpq = sqrt(fDist2);
            const double Wpq = cubicSplineKernel(rpq, pH);

            *omega += Wpq;

            E[XX] += dr[0]*dr[0]*Wpq;
            E[YY] += dr[1]*dr[1]*Wpq;
            E[ZZ] += dr[2]*dr[2]*Wpq;

            E[XY] += dr[1]*dr[0]*Wpq;
            E[XZ] += dr[2]*dr[0]*Wpq;
            E[YZ] += dr[1]*dr[2]*Wpq;

            *nSmooth += 1;
        }
    }
}

static inline double densNodeNcondB(PKD pkd, particleStore::ParticleReference &p,
                                    TinyVector<double,6> &E, double omega) {
    TinyVector<double,6> B;

    // Normalize the matrix
    E /= omega;

    inverseMatrix(E.data(), B.data());
    double Ncond = conditionNumber(E.data(), B.data());
    assert(Ncond==Ncond);

    if (p.is_gas()) {
        // We can already set this here, so it can be skipped in
        // hydroGradients
        auto &sph = p.sph();
        sph.Ncond = Ncond;
        sph.B = B;
    }

    return Ncond;
}

void packHydroDensity(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<hydroDensityPack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->fBall = p2.ball();
        p1->bMarked = p2.marked();
    }
}

void unpackHydroDensity(void *vpkd,void *dst,const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const hydroDensityPack *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_position(p2->position);
        p1.set_ball(p2->fBall);
        p1.set_marked(p2->bMarked);
    }
}

// Compute the density and derived variables simply given the fBall,
// without trying to converge to the correct value
void hydroDensityFinal(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
#ifdef OPTIM_UNION_EXTRAFIELDS
    double *omega = NULL;
    omega = p.is_gas()  ? &(p.sph().omega)  : omega;
    omega = p.is_star() ? &(p.star().omega) : omega;
    omega = p.is_bh()   ? &(p.BH().omega)   : omega;
#else
    // Assuming *only* stars and gas
    double *omega = &(p.sph().omega);
#endif

    const auto pH = p.ball();
    *omega = 0.0;

    if (p.is_gas()) {
        TinyVector<double,6> E{0.0};

        for (auto i = 0; i < nSmooth; ++i) {
            const double rpq = sqrt(nnList[i].fDist2);
            const double Wpq = cubicSplineKernel(rpq, pH);

            *omega += Wpq;

            const auto &dr = nnList[i].dr;

            E[XX] += dr[0]*dr[0]*Wpq;
            E[YY] += dr[1]*dr[1]*Wpq;
            E[ZZ] += dr[2]*dr[2]*Wpq;

            E[XY] += dr[1]*dr[0]*Wpq;
            E[XZ] += dr[2]*dr[0]*Wpq;
            E[YZ] += dr[1]*dr[2]*Wpq;
        }

        /* Normalize the matrix */
        E /= *omega;

        inverseMatrix(E.data(), p.sph().B.data());
    }
    else {
        for (auto i = 0; i < nSmooth; ++i) {
            const double rpq = sqrt(nnList[i].fDist2);
            *omega += cubicSplineKernel(rpq, pH);
        }
    }

    p.set_density(p.mass() * (*omega));
}

void hydroDensity_node(PKD pkd, SMF *smf, Bound bnd_node, const std::vector<PARTICLE *> &sinks,
                       NN *nnList,int nCnt) {
    for (auto &P : sinks) {
        auto partj = pkd->particles[P];
#ifdef OPTIM_UNION_EXTRAFIELDS
        double *omega = NULL;
        omega = partj.is_gas()  ? &(partj.sph().omega)  : omega;
        omega = partj.is_star() ? &(partj.star().omega) : omega;
        omega = partj.is_bh()   ? &(partj.BH().omega)   : omega;
#else
        // Assuming *only* stars and gas
        double *omega = &(partj.sph().omega);
#endif
        const TinyVector<double,3> dr_node {bnd_node.center() - partj.position()};

        int niter = 0;
        int nSmooth;
        int onLowerLimit = 0;

        float Neff = smf->nSmooth;
        if (*omega < 0)
            Neff = - *omega;

        double dConvFac = .5;
        do {
            const auto pH = partj.ball();
            TinyVector<double,6> E;

            densNodeOmegaE(nnList, pH, dr_node, nCnt, omega, E, &nSmooth);

            // Check if it has converged
            const double c = 4.*M_PI/3. * (*omega) * pH*pH*pH;
            if ((fabs(c-Neff) < smf->dNeighborsStd) ) {

                // Check if the converged density has a low enough condition number
                double Ncond = densNodeNcondB(pkd, partj, E, *omega);

                if (Ncond > 100) {
                    // In some configurations (outer particles in a isolated galaxy)
                    // the density could have converged but with a very high Ncond
                    // due to anisotropy.
                    if (Neff<200.) {
                        // To decrease anisotropy we increase the effective number
                        // of neighbours and repeat the iterative procedure
                        Neff *= 1.2;
                        niter = 0;
                        continue;
                    }
                    else {
                        // But we impose a maximun number of neighbours to avoid
                        // too long interaction lists
                        printf("WARNING %d Maximum Neff reached: %e ; Ncond %e \n", nSmooth, Neff, Ncond);
                    }
                }
                if (nSmooth < 20) {
                    // If we converge to a fBall that encloses too few neighbours,
                    // the effective number of ngb is increased
                    Neff *= 1.2;
                    niter = 0;
                    continue;
                }

                partj.set_marked(false);
                partj.set_density(partj.mass() * (*omega));
            }
            else {
                // We need to keep iterating
                constexpr float HMaxFac = 4.0f;
                const bool isMaxFac = Neff > c * pow(HMaxFac, 3);
                const float HFac = isMaxFac ? HMaxFac : pow(Neff/c, 1./3.);
                const float pNewH = pH * (1. + dConvFac*(HFac - 1.));
                partj.set_ball(pNewH);

                if (pNewH > pH) {
                    // We check that the proposed ball is enclosed within the
                    // node search region
                    if (any(abs(dr_node) + partj.ball() > bnd_node.apothem())) {
                        // A negative omega indicates the Neff for the next
                        // nbg search
                        *omega = -Neff;
                        break;
                    }
                }
                else {
                    // Put a hard lower limit based on the softening, even if we
                    // have not yet fully converged due to, e.g., an anisotropic
                    // particle distribution
                    if (smf->dhMinOverSoft > 0) {
                        const float pMinH = 2 * smf->dhMinOverSoft * partj.soft();
                        if (pNewH < pMinH) {
                            partj.set_ball(pMinH);
                            if (!onLowerLimit) {
                                // If in the next iteration still a lower smooth is
                                // preferred, we will skip this particle
                                onLowerLimit = 1;
                            }
                            else {
                                densNodeNcondB(pkd, partj, E, *omega);
                                partj.set_density(partj.mass() * (*omega));
                                partj.set_marked(false);
                            }
                        }
                    }
                    if (nSmooth < 20) {
                        if (!onLowerLimit) {
                            onLowerLimit = 1;
                        }
                        else {
                            // If we converge to a fBall that encloses too few neighbours,
                            // the effective number of nn is increased
                            Neff *= 1.2;
                            niter = 0;
                            onLowerLimit = 0;
                            continue;
                        }
                    }
                }
            }
            niter++;
            // Decrease the weight of the newly proposed smoothing length once
            // in a while to avoid getting into infinite loops
            if (niter%10==0 && niter<90)
                dConvFac *= 0.9;

            if (niter>100 && partj.marked()) {
                // At this point, we probably have a particle with plenty of
                // neighbours, and a small increase/decrease in the radius causes
                // omega to fluctuate around the desired value.
                //
                // In these cases, we stop iterating, making sure that we have the
                // upper value, which should have more than the expected number of
                // neighbours
                if (c > Neff) {
                    partj.set_marked(false);
                    partj.set_ball(pH);
                    densNodeNcondB(pkd, partj, E, *omega);
                    partj.set_density(partj.mass() * (*omega));
                    printf("WARNING %d Maximum iterations reached "
                           "Neff %e c %e \t %e %e %e %e \n",
                           partj.species(), Neff, c, partj.density(), partj.position(0),
                           partj.position(1),partj.position(2));
                }
            }
        } while (partj.marked());
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
    double pH, rpq, Hpq, c;

    /* Particle p data */
    auto &sph = p.sph();

    /* Compute the \omega(x_i) normalization factor */
    sph.omega = 0.0;
    pH = p.ball();
#ifdef OPTIM_DENSITY_REITER
    float maxr = 0.0;
#endif //OPTIM_DENSITY_REITER
    for (auto i = 0; i < nSmooth; ++i) {
#ifdef OPTIM_DENSITY_REITER
        if (nnList[i].fDist2> maxr) maxr =nnList[i].fDist2;
#else
        rpq = sqrt(nnList[i].fDist2);
        Hpq = pH;

        sph.omega += cubicSplineKernel(rpq, Hpq);
#endif //OPTIM_DENSITY_REITER
    }

    /* If we are using a iterative procedure for computing the smoothing length:
     *    - Particles marked are those which still needs iterations
     *    - Particles not marked are those with a correct smoothing length
     */
#ifdef OPTIM_DENSITY_REITER
    maxr = sqrt(maxr);
    while (p.marked()) {
        sph.omega=0.0;
        pH = p.ball();
        for (auto i = 0; i < nSmooth; ++i) {
            rpq = sqrt(nnList[i].fDist2);
            Hpq = pH;

            sph.omega += cubicSplineKernel(rpq, Hpq);
        }
#else
    if (smf->bIterativeSmoothingLength && p.marked()) {
#endif //OPTIM_DENSITY_REITER

        c = 4.*M_PI/3. * sph.omega *pH*pH*pH;
        if (fabs(c-smf->nSmooth) < smf->dNeighborsStd) {
            p.set_marked(false);
        }
        else {
            float pNewH;
            pNewH = pH * pow(smf->nSmooth/c,0.3333333333);
            //   if (nSmooth <= 1) pNewH *= 2.*fBall;

            p.set_ball(0.5*(pNewH+pH));
            //sph.fLastBall = pH;

#ifdef OPTIM_DENSITY_REITER
            // If the suggested new radius does not enclose all our neighbors,
            // we need to reiterate
            if (p.ball()>maxr) break;
#endif
        }
    }

    /* We compute the density making use of Eq. 27 Hopkins 2015 */
    p.set_density(p.mass() * sph.omega);
}
