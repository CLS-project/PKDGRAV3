#include <algorithm>
#include "hydro/hydro.h"
#include "master.h"
using blitz::TinyVector;
using blitz::dot;

void MSR::HydroStep(double dTime, double dDelta) {
    double dsec;

    printf("Computing hydro time step... ");

    TimerStart(TIMER_TIMESTEP);
#ifdef OPTIM_SMOOTH_NODE
#ifdef OPTIM_AVOID_IS_ACTIVE
    SelActives();
#endif
    ReSmoothNode(dTime, dDelta, SMX_HYDRO_STEP, 1);
#else
    ReSmooth(dTime, dDelta, SMX_HYDRO_STEP, 1);
#endif

    if (param.bGlobalDt) {
        uint8_t minDt;
        minDt = GetMinDt();
        SetGlobalDt(minDt);
    }

    TimerStop(TIMER_TIMESTEP);
    dsec = TimerGet(TIMER_TIMESTEP);
    printf("took %.5f seconds\n", dsec);

    if (param.bWakeUpParticles) {
        struct inDrift in;
        in.iRoot = 0; // Not used
        in.dTime = dTime;
        in.dDelta = 0; // Not used
        pstWakeParticles(pst,&in,sizeof(in),NULL,0);
    }
}

void initHydroStep(void *vpkd, void *vp) {
}

/* Compute the hydrodynamical time step of this particle, based:
 *    - Signal velocity
 *    - Acceleration
 *    - Timestep of the neighouring particles (in a scatter approach)
 */
void hydroStep(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    double dtEst;

    auto &psph = p.sph();
    const double ph = 0.5*p.ball();

    dtEst = HUGE_VAL;

    /*  Signal velocity criteria */
    for (auto i = 0; i < nSmooth; ++i) {
        if (pIn == nnList[i].pPart) continue;
        auto q = pkd->particles[nnList[i].pPart];
        const auto &qsph = q.sph();

        const auto &dr = nnList[i].dr;


        // From Eqs 24,25 Hopkins 2015, to limit deltaT
        const double dvDotdr = dot(dr,qsph.vPred - psph.vPred);
        double vsig_pq = psph.c + qsph.c;
        if (dvDotdr < 0.) vsig_pq -= dvDotdr/sqrt(nnList[i].fDist2);

        const double dt2 = 2.*smf->dEtaCourant * ph * smf->a / vsig_pq;
        if (dt2 < dtEst) dtEst = dt2;
    }

#ifdef ENTROPY_SWITCH
    // Gather the maximum relative kinetic energy
    psph.maxEkin = 0.0;
    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
        auto &qsph = q.sph();

        double dv2 = dot(qsph.vPred - psph.vPred, qsph.vPred - psph.vPred);
        double Ekin = 0.5*p.mass()*dv2;

        psph.maxEkin = (psph.maxEkin < Ekin) ? Ekin : psph.maxEkin;

    }
#endif

#ifdef HERNQUIST_POTENTIAL
    // Timestep criteria based on the Hernsquist potential
    const double const_reduced_hubble_cgs = 3.2407789e-18;
    //const double H0 = 0.704 * const_reduced_hubble_cgs * pkd->param.dSecUnit;
    const double H0 = 70.4/ smf->units.dKmPerSecUnit * ( smf->units.dKpcUnit / 1e3);

    const double concentration = 9.0;
    const double M200 = 135.28423603962767; //137.0 ; // / pkd->param.dMsolUnit;
    const double V200 = cbrt(M200*H0);
    //const double R200 = V200/(H0);
    const double R200 = cbrt(M200/(100.*H0*H0));
    const double RS = R200 / concentration;

    const double al = RS * sqrt(2. * (log(1. + concentration) -
                                      concentration / (1. + concentration)));

    const double mass = M200;
    (1.-0.041);

    /* Calculate the relative potential with respect to the centre of the
     * potential */
    auto dr = p.position(); //- potential->x;

    /* calculate the radius  */
    const double epsilon =  0.2/smf->units.dKpcUnit;
    const double epsilon2 = epsilon*epsilon;
    const float r = sqrtf(dot(dr,dr) + epsilon2);
    const float sqrtgm_inv = 1.f / sqrtf(mass);

    /* Calculate the circular orbital period */
    const float period = 2.f * M_PI * sqrtf(r) * al *
                         (1 + r / al) * sqrtgm_inv;

    /* Time-step as a fraction of the circular orbital time */
    const double time_step = 0.01 * period;

    if (time_step < dtEst) dtEst = time_step;
#endif





    // Timestep criteria based on the hydro+grav accelerations

    const TinyVector<double,3>
    a {p.acceleration() - smf->a *(-p.velocity()*psph.Frho + psph.Fmom)/p.mass()};

    // 1/a^3 from grav2.c:217
    // Explanation:
    // 1/a from different acceleration definition
    // 1/a in grav source term (ie. kick)
    // 1/a from the drift
    // But for the hydro, it is 1/a^2:
    // 1/a from hydro eqs (included normally in minDt)
    // 1/a from pkdVel, which is normally incorporated in the drift
    const double aFac = 1./(smf->a * smf->a * smf->a);
    const double acc = sqrt(dot(a,a)) * aFac;

    const float h = smf->bDoGravity ? std::min(ph, static_cast<double>(p.soft())) : ph;
    const double dtAcc = smf->dCFLacc * sqrt(2.*h/acc);


    if (dtAcc < dtEst) dtEst = dtAcc;
    const uint8_t uNewRung = pkdDtToRung(dtEst,smf->dDelta,MAX_RUNG);
    if (uNewRung > p.new_rung()) p.set_new_rung(uNewRung);


    /* Timestep limiter that imposes that the particle must have a dt
     * which is at most four times (i.e., 2 rungs) the smallest dt of the neighbours
     *
     * This is implemented as a scatter approach, in which the maximum dt of the
     * neighbours is set given the dt computed from this particle
     */
    if (smf->dDelta >= 0) {
        for (auto i = 0; i < nSmooth; ++i) {
            auto q = pkd->particles[nnList[i].pPart];
            if (p.new_rung() > q.new_rung() + 2) {
                q.set_new_rung(p.new_rung() - 1);
                if (!q.is_active()) {
                    // DEBUG: The number of wake up request and IDs must match!
                    //printf("Need to wake up! %"PRIu64" \n", *pkdParticleID(pkd,q));
                    q.sph().uWake = p.new_rung();
                }
            }
        }
    }
}

void combHydroStep(void *vpkd, void *v1,const void *v2) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(v1)];
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(v2)];
    assert(!pkd->bNoParticleOrder);
    if (p1.is_gas() && p2.is_gas()) {
        if (p2.new_rung() > p1.new_rung()) p1.set_new_rung(p2.new_rung());
    }
}

/* Wake the particles whose uWake is greater than zero (see hydroStep)
 *
 * TODO: clean unused function arguments
 */
void pkdWakeParticles(PKD pkd,int iRoot, double dTime, double dDelta) {
    for (auto &p : pkd->particles) {
        if (p.is_gas()) {
            auto &sph = p.sph();
            uint8_t uWake = sph.uWake;
            if (uWake) {
                p.set_rung(uWake);
                p.set_new_rung(uWake);
                sph.uWake = 0;
                sph.lastUpdateTime = dTime;

                // We revert to the state at the end of the previous timestep
                sph.E    = sph.lastE;
                sph.Uint = sph.lastUint;
                p.set_mass(sph.lastMass);
                sph.mom = sph.lastMom;

                /* NOTE: What do we do with the variables?
                 *  We could integrate the source terms without major problems.
                 *  However, with the hydrodynamics is not that easy.
                 *
                 *
                 *  Fene/Fmom are unusable at this time because the particle is
                 *  not synced. At most, we could use the previous hydro derivatives
                 *  to extrapolate up to this time... but they are also unusable
                 *
                 *  Storing them would be adding more variables that will be
                 *  rarely used... And even doing that... this will not be
                 *  conservative!!
                 */
            }
        }
    }
}

