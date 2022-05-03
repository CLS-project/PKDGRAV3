
#include "blackhole/seed.h"
#include "smooth/smooth.h"
#include "master.h"

void MSR::PlaceBHSeed(double dTime, uint8_t uRungMax) {
    struct inPlaceBHSeed in;
    struct outPlaceBHSeed out;

    out.nBHs = 0;
    in.dTime = dTime;
    in.uRungMax = uRungMax;
    in.dScaleFactor = csmTime2Exp(csm,dTime);
    in.dBHMhaloMin = param.dBHMhaloMin;
    in.dTau = param.dTau;
    in.dInitialH = param.dInitialH;
    in.dBHSeedMass = param.dBHSeedMass;
#ifdef STAR_FORMATION
    if (csm->val.bComove)
        in.dDenMin = param.dSFThresholdDen*pow(in.dScaleFactor,3);
#else
    in.dDenMin = 0.0;
#endif

    printf("Planting seeds...\n");
    TimerStart(TIMER_BHS);
    assert(param.bFindGroups);
    pstPlaceBHSeed(pst, &in, sizeof(in), &out, sizeof(out));
    TimerStop(TIMER_BHS);

    nGas -= out.nBHs;
    nBH += out.nBHs;
    printf("Planted %d BH seeds \n", out.nBHs);


#ifdef OPTIM_REORDER_IN_NODES
    if (out.nBHs > 0)
        ReorderWithinNodes();
#endif

}

#ifdef __cplusplus
extern "C" {
#endif

int pkdPlaceBHSeed(PKD pkd, double dTime, double dScaleFactor,
                   uint8_t uRungMax, double dDenMin, double dBHMhaloMin,
                   double dTau, double dInitialH, double dBHSeedMass) {
    int newBHs = 0;
    SMX smx;
    smInitialize(&smx,pkd,NULL,32,1,0,SMX_NULL);

    // Look for any FoF group that do not contain any BH particle in it
    for (int gid=1; gid<=pkd->nLocalGroups; ++gid) {
        //printf("group %d mass %e \n", gid, pkd->veryTinyGroupTable[gid].fMass);
        //assert (pkd->ga[gid].id.iPid == pkd->idSelf);

        if (pkd->veryTinyGroupTable[gid].nBH==0 &&
                pkd->veryTinyGroupTable[gid].fMass > dBHMhaloMin) {
            printf("Group %d (mass %e) do not have any BH, adding one!\n",
                   gid, pkd->veryTinyGroupTable[gid].fMass);
            pkd->veryTinyGroupTable[gid].nBH++;



            // Fake particle for the neighbour search
            PARTICLE *p = NULL;

            smx->nnListSize = 0;
            smGather(smx,dTau,pkd->veryTinyGroupTable[gid].rPot, p);


            float minPot = HUGE_VAL;
            int index = 0;
            PARTICLE *pLowPot = NULL;
            for (int i=0; i<smx->nnListSize; i++) {
                if (pkd->idSelf != smx->nnList[i].iPid) continue;
                PARTICLE *q = smx->nnList[i].pPart;
                assert(pkdIsGas(pkd,q));
                pLowPot = (*pkdPot(pkd,q)<minPot)  ? q : pLowPot;
                index += 1;
                minPot = (pLowPot!=NULL) ? *pkdPot(pkd,pLowPot) : minPot;
            }

            // IA: I do not like this *at all*
            //  But maybe we are reading completely remote fof group??? TODO Check
            if (index==0) continue;

            // IA: We require the density to be above the SF threshold
#ifdef COOLING
            const double rho_H = pkdDensity(pkd,pLowPot) *
                                 pkdSph(pkd,pLowPot)->afElemMass[ELEMENT_H] / pkdMass(pkd,pLowPot);
#else
            // If no information, assume primoridal abundance
            const double rho_H = pkdDensity(pkd,pLowPot) * dInitialH;
#endif
            if (rho_H < dDenMin) continue;

            assert(pLowPot!=NULL);
            assert(pkdIsGas(pkd,pLowPot));

            // Now convert this particle into a BH
            // We just change the class of the particle
            double omega = pkdSph(pkd,pLowPot)->omega;
            pkdSetClass(pkd, pkdMass(pkd,pLowPot), pkdSoft0(pkd,pLowPot),
                        FIO_SPECIES_BH, pLowPot);

            BHFIELDS *pBH = pkdBH(pkd,pLowPot);
            // When changing the class, we have to take into account tht
            // the code velocity has different scale factor dependencies for
            // dm/star/bh particles and gas particles
            float *pv = pkdVel(pkd,pLowPot);
            for (int j=0; j<3; j++) {
                pv[j] *= dScaleFactor;
            }

            // We initialize the particle. Take into account that we need to
            // set EVERY variable, because as we use unions, there may be
            // dirty stuff
            pBH->omega = omega;
            pBH->dInternalMass = dBHSeedMass;
            pBH->newPos[0] = -1; // Ask for a reposition
            pBH->lastUpdateTime = dTime;
            pBH->dAccretionRate = 0.0;
            pBH->dEddingtonRatio = 0.0;
            pBH->dFeedbackRate = 0.0;
            pBH->dAccEnergy = 0.0;
            pBH->fTimer = dTime;
            pBH->doReposition = 2;

            // As the particle that was converted to a BH lies in a very
            // dense environment it will probably have a high rung, so
            // this is not required
            pLowPot->uNewRung = uRungMax;


            pkd->nBH++;
            pkd->nGas--;
            newBHs++;

            /* Old seeding process, which creates a BH rather than converting
             * a gas particle.
             * This is problematic because we would need to assing a
             * unique particleID to this newly
             * created BH, which is to completely trivial to do.
             *
             * Furthermore, Booth & Schaye also converted gas particles
             * into BH anyway
             */
            /*
            PARTICLE* p = (PARTICLE *) malloc(pkdParticleSize(pkd));
            vel_t *vel = pkdVel(pkd,p);
            float* acc = pkdAccel(pkd,p);
            for (int j=0; j<3; j++){
               pkdSetPos(pkd, p, j, pkd->veryTinyGroupTable[gid].rPot[j]);
               vel[j] = 0.0;
               acc[j] = 0.0;
            }



            p->uRung = uRungMax;
            p->uNewRung = uRungMax;

            // This class should have been created while reading the IC (at pkdReadFIO),
            // this won't be problematic as long as we use oMass and dSoft is set in the parameters file
            pkdSetClass(pkd,0,pkdSoft(pkd,p),FIO_SPECIES_BH,p);

            float *pmass = pkdField(p,pkd->oMass);
            *pmass = pkd->param.dBHSeedMass;
            pkdSetBall(pkd,p,pkd->param.dSoft);

            BHFIELDS* pBH = pkdBH(pkd,p);

            pBH->dAccretionRate = 0.0;
            pBH->dFeedbackRate = 0.0;
            pBH->dAccEnergy = 0.0;
            pBH->lastUpdateTime = dTime;
            pBH->fTimer = dTime;
            pBH->dInternalMass = *pmass;
            pBH->newPos[0] = -1;

            // TODO: Check in output that this is coherent with the rest of particles in the group
            pkdSetGroup(pkd,p,gid);

            pkdNewParticle(pkd,p);
            free(p);

            p = pkdParticle(pkd,pkd->nLocal-1);
            */
        }


    }
    smFinish(smx,NULL);

    return newBHs;

}
#ifdef __cplusplus
}
#endif
