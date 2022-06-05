#include "hydro/hydro.h"
#include "master.h"


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
void hydroStep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    PARTICLE *q;
    SPHFIELDS *psph, *qsph;
    double dt2, dtEst, vsig_pq, dvDotdr, dx, dy, dz;
    uint8_t uNewRung;
    int i,j;

    psph = pkdSph(pkd, p);

    dtEst = HUGE_VAL;

    /*  Signal velocity criteria */
    for (i=0; i<nSmooth; ++i) {

        dx = nnList[i].dx;
        dy = nnList[i].dy;
        dz = nnList[i].dz;

        if (dx==0 && dy==0 && dz==0) continue;

        q = nnList[i].pPart;
        qsph = pkdSph(pkd, q);

        // From Eqs 24,25 Hopkins 2015, to limit deltaT
        dvDotdr = (dx*(qsph->vPred[0] - psph->vPred[0]) +
                   dy*(qsph->vPred[1] - psph->vPred[1]) +
                   dz*(qsph->vPred[2] - psph->vPred[2]));

        if (dvDotdr < 0) {
            vsig_pq = psph->c + qsph->c - dvDotdr/sqrt(nnList[i].fDist2);
        }
        else {
            vsig_pq = psph->c + qsph->c;
        }

        dt2 = 2.*smf->dEtaCourant * fBall * smf->a /vsig_pq;
        if (dt2 < dtEst) dtEst=dt2;

    }

#ifdef ENTROPY_SWITCH
    // Gather the maximum relative kinetic energy
    psph->maxEkin = 0.0;
    for (i=0; i<nSmooth; i++) {
        q = nnList[i].pPart;
        qsph = pkdSph(pkd, q);

        double dv2 =(qsph->vPred[0] - psph->vPred[0])*(qsph->vPred[0] - psph->vPred[0]) +
                    (qsph->vPred[1] - psph->vPred[1])*(qsph->vPred[1] - psph->vPred[1]) +
                    (qsph->vPred[2] - psph->vPred[2])*(qsph->vPred[2] - psph->vPred[2]);
        double Ekin = 0.5*pkdMass(pkd,p)*dv2;

        psph->maxEkin = (psph->maxEkin < Ekin) ? Ekin : psph->maxEkin;

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
    dx = pkdPos(pkd,p,0); //- potential->x[0];
    dy = pkdPos(pkd,p,1); //- potential->x[1];
    dz = pkdPos(pkd,p,2); //- potential->x[2];

    /* calculate the radius  */
    const double epsilon =  0.2/smf->units.dKpcUnit;
    const double epsilon2 = epsilon*epsilon;
    const float r = sqrtf(dx * dx + dy * dy + dz * dz + epsilon2);
    const float sqrtgm_inv = 1.f / sqrtf(mass);

    /* Calculate the circular orbital period */
    const float period = 2.f * M_PI * sqrtf(r) * al *
                         (1 + r / al) * sqrtgm_inv;

    /* Time-step as a fraction of the circular orbital time */
    double time_step = 0.01 * period;

    if (time_step < dtEst) dtEst = time_step;
#endif





    // Timestep criteria based on the hydro+grav accelerations

    double a[3], acc;
    float *pa = pkdAccel(pkd,p);
    double cfl = smf->dCFLacc, dtAcc;


    for (j=0; j<3; j++) {
        a[j] = pa[j] -
               smf->a*(-pkdVel(pkd,p)[j]*psph->Frho + psph->Fmom[j])/pkdMass(pkd,p);
    }

    // 1/a^3 from grav2.c:217
    // Explanation:
    // 1/a from different acceleration definition
    // 1/a in grav source term (ie. kick)
    // 1/a from the drift
    // But for the hydro, it is 1/a^2:
    // 1/a from hydro eqs (included normally in minDt)
    // 1/a from pkdVel, which is normally incorporated in the drift
    double aFac = 1./(smf->a * smf->a * smf->a);

    acc = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) * aFac;


    float h = smf->bDoGravity ?
              (fBall < pkdSoft(pkd,p) ? fBall : pkdSoft(pkd,p) )
              : fBall;
    dtAcc = cfl*sqrt(2*h/acc);


    if (dtAcc < dtEst) dtEst = dtAcc;
    uNewRung = pkdDtToRung(dtEst,smf->dDelta,MAX_RUNG);
    if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung;


    /* Timestep limiter that imposes that the particle must have a dt
     * which is at most four times (i.e., 2 rungs) the smallest dt of the neighbours
     *
     * This is implemented as a scatter approach, in which the maximum dt of the
     * neighbours is set given the dt computed from this particle
     */
    if (smf->dDelta >= 0) {
        for (i=0; i<nSmooth; ++i) {
            q = nnList[i].pPart;

            //if ( (q->uNewRung - p->uNewRung) > 2) uNewRung = q->uNewRung-1;
            //if (uNewRung > p->uNewRung ) p->uNewRung = uNewRung;


            if ( (p->uNewRung - q->uNewRung) > 2) {
                q->uNewRung = p->uNewRung-1;
                if (!pkdIsActive(pkd,q)) {
                    // DEBUG: The number of wake up request and IDs must match!
                    //printf("Need to wake up! %"PRIu64" \n", *pkdParticleID(pkd,q));
                    pkdSph(pkd,q)->uWake = p->uNewRung;
                }
            }

        }

    }


}

void combHydroStep(void *vpkd, void *v1,const void *v2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = (PARTICLE *) v1;
    PARTICLE *p2 = (PARTICLE *) v2;
    assert(!pkd->bNoParticleOrder);
    if (pkdIsGas(pkd,p1) && pkdIsGas(pkd,p2)) {

        if (((PARTICLE *) p2)->uNewRung > ((PARTICLE *) p1)->uNewRung)
            ((PARTICLE *) p1)->uNewRung = ((PARTICLE *) p2)->uNewRung;
    }
}

/* Wake the particles whose uWake is greater than zero (see hydroStep)
 *
 * TODO: clean unused function arguments
 */
void pkdWakeParticles(PKD pkd,int iRoot, double dTime, double dDelta) {
    for (int i=0; i<pkd->Local(); ++i) {
        PARTICLE *p = pkd->Particle(i);
        if (pkdIsGas(pkd,p)) {
            uint8_t uWake = pkdSph(pkd,p)->uWake;
            if (uWake) {
                SPHFIELDS *psph = pkdSph(pkd,p);

                p->uRung = uWake;
                p->uNewRung = uWake;
                psph->uWake = 0;

                psph->lastUpdateTime = dTime;

                // We revert to the state at the end of the previous timestep
                psph->E    = psph->lastE;
                psph->Uint = psph->lastUint;
                float *mass = (float *) pkdField(p, pkd->oFieldOffset[oMass]);
                *mass = psph->lastMass;

                for (int j=0; j<3; j++) {
                    psph->mom[j] = psph->lastMom[j];
                }

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

                // DEBUG: The number of wake up request and IDs must match!
                //printf("Waking up %"PRIu64" \n", *pkdParticleID(pkd,p));
            }
        }
    }

}

