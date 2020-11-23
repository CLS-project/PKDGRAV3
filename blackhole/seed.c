
#include "blackhole/seed.h"
#include "smooth.h"

void msrPlaceBHSeed(MSR msr, double dTime, uint8_t uRungMax) {
   struct inPlaceBHSeed in;
   struct outPlaceBHSeed out;

   out.nBHs = 0;
   in.dTime = dTime;
   in.uRungMax = uRungMax;
   in.dScaleFactor = csmTime2Exp(msr->csm,dTime);
#ifdef STAR_FORMATION
   if (msr->csm->val.bComove)
      in.dDenMin = msr->param.dSFThresholdDen*pow(in.dScaleFactor,3);
#else
   in.dDenMin = 0.0;
#endif

   printf("Planting seeds...\n");
   assert(msr->param.bFindGroups);
   pstPlaceBHSeed(msr->pst, &in, sizeof(in), &out, sizeof(out));

   msr->N += out.nBHs;
   msr->nBH += out.nBHs;
   printf("Planted %d BH seeds \n", out.nBHs);

}


int pkdPlaceBHSeed(PKD pkd, double dTime, double dScaleFactor, uint8_t uRungMax, double dDenMin){
   int newBHs = 0;
   SMX smx;
   smInitialize(&smx,pkd,NULL,32,1,0,SMX_NULL);

   // Look for any FoF group that do not contain any BH particle in it
   for (int gid=1;gid<=pkd->nLocalGroups;++gid) {
      //printf("group %d mass %e \n", gid, pkd->veryTinyGroupTable[gid].fMass);
      //assert (pkd->ga[gid].id.iPid == pkd->idSelf);

      if (pkd->veryTinyGroupTable[gid].nBH==0 && pkd->veryTinyGroupTable[gid].fMass > pkd->param.dMhaloMin){
         printf("Group %d (mass %e) do not have any BH, adding one!\n", gid, pkd->veryTinyGroupTable[gid].fMass);
         pkd->veryTinyGroupTable[gid].nBH++;



         // We look for gas particles close to the center of potential of the halo
         float fBall = pkdSoft(pkd, pkdParticle(pkd,0)  );

         // Fake particle for the neighbour search
         PARTICLE* p = NULL;// (PARTICLE *) malloc(pkdParticleSize(pkd));
         //for (int j=0; j<3; j++){
         //   pkdSetPos(pkd, p, j, pkd->veryTinyGroupTable[gid].rPot[j]);
         //}

         printf("Looking for ngbs.. %e \n", fBall);
         smx->nnListSize = 0;
	   smGather(smx,pkd->param.dTau,pkd->veryTinyGroupTable[gid].rPot, p);
         printf("End ngb search nSmooth %d %d \n", smx->nnListSize, smx->nSmooth);


         float minPot = HUGE_VAL;
         int index = 0;
         PARTICLE *pLowPot = NULL;
         for (int i=0; i<smx->nnListSize; i++){
            if (pkd->idSelf != smx->nnList[i].iPid) continue;
            PARTICLE* q = smx->nnList[i].pPart;
            assert(pkdIsGas(pkd,q));
            pLowPot = (*pkdPot(pkd,q)<minPot)  ? q : pLowPot;
            index += 1; 
            minPot = (pLowPot!=NULL) ? *pkdPot(pkd,pLowPot) : minPot;
         }
         printf("End pLowPot search %d \n", index);

         // IA: I do not like this *at all*
         //  But maybe we are reading completely remote fof group??? TODO Check
         if (index==0) continue;

         // IA: We require the density to be above the SF threshold
#ifdef COOLING
         const double rho_H = pkdDensity(pkd,pLowPot) * pkdSph(pkd,pLowPot)->chemistry[chemistry_element_H];
#else
         const double rho_H = pkdDensity(pkd,pLowPot) * 0.75; // If no information, assume primoridal abundance
#endif
         if (rho_H < dDenMin) continue;

         assert(pLowPot!=NULL);
         assert(pkdIsGas(pkd,pLowPot));

         // If, by chance (this is highly inprob), the most bound particle lies in another domain, then
         //  the other proc should take care of it
         double omega = pkdSph(pkd,pLowPot)->omega;
         float fmass = pkdMass(pkd,pLowPot);
         // Now convert this particle into a BH
         // We just change the class of the particle to stellar one
         
         pkdSetClass(pkd, 0., 0., FIO_SPECIES_BH, pLowPot);
         
         float *pfmass = pkdField(pLowPot, pkd->oMass);
         *pfmass = fmass;

         BHFIELDS* pBH = pkdBH(pkd,pLowPot);
         // When changing the class, we have to take into account that the code velocity
         // has different scale factor dependencies for dm/star/bh particles and gas particles
         float *pv = pkdVel(pkd,pLowPot);
         for (int j=0; j<3; j++){
            pv[j] *= dScaleFactor;
         }

         // We initialize the particle. Take into account that we need to set EVERY variable,
         // because as we use unions, there may be dirty stuff
         pBH->omega = omega;
         pBH->dInternalMass = pkd->param.dBHSeedMass;
         pBH->newPos[0] = -1; // Ask for a reposition
         pBH->lastUpdateTime = dTime;
         pBH->dAccretionRate = 0.0;
         pBH->dEddingtonRatio = 0.0;
         pBH->dFeedbackRate = 0.0;
         pBH->dAccEnergy = 0.0;
         pBH->fTimer = dTime;
         
         pLowPot->uNewRung = uRungMax;
         
         printf("New BH pos %e %e %e \n", pkdPos(pkd,pLowPot, 0), pkdPos(pkd,pLowPot, 1), pkdPos(pkd,pLowPot, 2));

         pkd->nBH++;
         pkd->nGas--;
         newBHs++;
         

         /* Old seeding process, which creates a BH rather than converting a gas particle.
          * This is problematic because we would need to assing a unique particleID to this newly
          * created BH, which is to completely trivial to do.
          *
          * Furthermore, Booth & Schaye also converted gas particles into BH anyway
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
