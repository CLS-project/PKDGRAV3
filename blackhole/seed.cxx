#include <algorithm>
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

            smx->nnListSize = 0;
            smGather(smx,dTau,pkd->veryTinyGroupTable[gid].rPot);

            // Find the first local particle
            auto ii = std::find_if(smx->nnList,smx->nnList+smx->nnListSize,[pkd](const auto &nn) {return nn.iPid==pkd->Self();});
            // IA: I do not like this *at all*
            //  But maybe we are reading completely remote fof group??? TODO Check
            if (ii == smx->nnList+smx->nnListSize) continue;

            // Now find the minumum
            ii = std::min_element(ii,smx->nnList+smx->nnListSize,[pkd](const auto &a,const auto &b) {
                if (b.iPid != pkd->Self()) return true; // Well, "a" is valid
                auto p = pkd->particles[a.pPart];
                auto q = pkd->particles[b.pPart];
                return p.potential() < q.potential();
            });
            assert(ii < smx->nnList+smx->nnListSize);
            auto pLowPot = pkd->particles[ii->pPart];

            // IA: We require the density to be above the SF threshold
#ifdef COOLING
            const double rho_H = pLowPot.density() * pLowPot.sph().afElemMass[ELEMENT_H] / pLowPot.mass();
#else
            // If no information, assume primoridal abundance
            const double rho_H = pLowPot.density() * dInitialH;
#endif
            if (rho_H < dDenMin) continue;

            assert(pLowPot.is_gas());

            // Now convert this particle into a BH
            // We just change the class of the particle
            double omega = pLowPot.sph().omega;
            pLowPot.set_class(pLowPot.mass(),pLowPot.soft0(),0,FIO_SPECIES_BH);

            auto BH = pLowPot.BH();
            // When changing the class, we have to take into account tht
            // the code velocity has different scale factor dependencies for
            // dm/star/bh particles and gas particles
            pLowPot.velocity() *= dScaleFactor;

            // We initialize the particle. Take into account that we need to
            // set EVERY variable, because as we use unions, there may be
            // dirty stuff
            BH.omega = omega;
            BH.dInternalMass = dBHSeedMass;
            BH.newPos[0] = -1; // Ask for a reposition
            BH.lastUpdateTime = dTime;
            BH.dAccretionRate = 0.0;
            BH.dEddingtonRatio = 0.0;
            BH.dFeedbackRate = 0.0;
            BH.dAccEnergy = 0.0;
            BH.fTimer = dTime;

            // As the particle that was converted to a BH lies in a very
            // dense environment it will probably have a high rung, so
            // this is not required
            //pLowPot.uNewRung = uRungMax;


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
            auto p = pkd->particles.NewParticle();
            p.velocity() = 0.0;
            p.acceleration() = 0.0;
            p.set_position(pkd->veryTinyGroupTable[gid].rPot);
            p.set_rung(uRungMax);
            p.set_new_rung(uRungMax);

            // This class should have been created while reading the IC (at pkdReadFIO),
            // this won't be problematic as long as we use oMass and dSoft is set in the parameters file
            p.set_class(0,p.soft(),0,FIO_SPECIES_BH);

            p.mass() = pkd->param.dBHSeedMass;
            p.set_ball(pkd->param.dSoft);

            auto &pBH = p.BH();

            pBH.dAccretionRate = 0.0;
            pBH.dFeedbackRate = 0.0;
            pBH.dAccEnergy = 0.0;
            pBH.lastUpdateTime = dTime;
            pBH.fTimer = dTime;
            pBH.dInternalMass = *pmass;
            pBH.newPos[0] = -1;

            // TODO: Check in output that this is coherent with the rest of particles in the group
            p.set_group(gid);
            */
        }


    }
    smFinish(smx,NULL);

    return newBHs;

}
#ifdef __cplusplus
}
#endif
