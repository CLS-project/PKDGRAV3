


#ifdef  BLACKHOLES
#include "blackhole/merger.h"


void msrBHmerger(MSR msr, double dTime) {
    struct inSmooth in;
    struct outSmooth out;
    struct outGetNParts Nout;
    int nOut;
    
    in.nSmooth = msr->param.nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = 0;
    in.iSmoothType = SMX_BH_MERGER;
    msrSmoothSetSMF(msr, &(in.smf), dTime);

    Nout.n = 0;
    Nout.nDark = 0;
    Nout.nGas = 0;
    Nout.nStar = 0;
    Nout.nBH = 0;

    if (msr->param.bVStep) {
	double sec,dsec;
	sec = msrTime();
	pstReSmoothNode(msr->pst,&in,sizeof(in),&out,sizeof(struct outSmooth));
      pstMoveDeletedParticles(msr->pst, NULL, 0, &Nout, sizeof(struct outGetNParts));
      pstRepositionBH(msr->pst, NULL, 0, NULL, 0);
	dsec = msrTime() - sec;
	printf("Merging %d BH particle pairs took %f secs\n\n", out.nSmoothed, dsec);
	}
    else {
	pstReSmoothNode(msr->pst,&in,sizeof(in),&out,sizeof(struct outSmooth));
      pstMoveDeletedParticles(msr->pst, NULL, 0, &Nout, sizeof(struct outGetNParts));
      pstRepositionBH(msr->pst, NULL, 0, NULL, 0);
	}
    msr->N = Nout.n;
    msr->nDark = Nout.nDark;
    msr->nGas = Nout.nGas;
    msr->nStar = Nout.nStar;
    msr->nBH = Nout.nBH;
    }


void pkdRepositionBH(PKD pkd){

   for (int i=0;i<pkdLocal(pkd);++i) { 
      PARTICLE* p = pkdParticle(pkd,i);
      if (pkdIsBH(pkd,p)){
         BHFIELDS* pBH = pkdBH(pkd,p);
         if (pBH->newPos[0]!=-1){
            for (int j=0; j<3; j++)
               pkdSetPos(pkd, p, j, pBH->newPos[j]);
            pBH->newPos[0] = -1;
         }
      }
   }
}



/* IA: This should be called when we have a merger across different domains... 
 * think about it! 
 * TODO
 **/
void combBHmerger(void *vpkd, void *p1,void *p2) {
    PKD pkd = (PKD) vpkd;

    assert(0);


    }



void smBHmerger(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    vel_t *pv, *qv;

    for (int i=0; i<nSmooth; ++i){
       PARTICLE* q = nnList[i].pPart;
       
       //printf("%d %d fBall %e pkdSoft %e %d %d \n", i, nSmooth, fBall, pkdSoft(pkd,p), p->bMarked, q->bMarked);

       /* IA: Merging across domain boundaries is a tricky thing,
        * This would require careful atomic operations in order to maintain
        * conservation *and* a deterministic behaviour independent of the number of cores.
        *
        * The easiest solution is to relax the deterministic and invariant behaviour and just
        * merge BH that are in the same domain.
        * This should be conservative, and in standard cosmological simulations this situation
        * will be very rare, and eventually both BH will lie in the same computational domain, so
        * I think this won't be a critical problem.
        *
        * Anyway, the most important thing is *conservation*, so that is what I should check for
        * in the test cases
        */
       if (pkd->idSelf != nnList[i].iPid) continue;


       if (p!=q && pkdIsBH(pkd,p) && pkdIsBH(pkd,q) && p->bMarked && q->bMarked){
          // At this point, we are sure that both particles are black holes and
          //  the distance between them is, at most, pkdSoft(pkd,p)
          // 
          // So, the following should not give any problem (remove for performance)
          assert( nnList[i].fDist2 < pkdSoft(pkd,p)*pkdSoft(pkd,p)  );
          assert( pkdIsBH(pkd,p) );
          assert( pkdIsBH(pkd,q) );
          assert( pkd->oMass );

          float qmass = pkdMass(pkd,q);
          float pmass = pkdMass(pkd,p);

          // As a rule, the particle that handles (and survives) the merger is the most massive one
          //if ( pmass >= qmass  ){
          if ( pkdPos(pkd,p,0) >= pkdPos(pkd,q,0)  ){
             pv = pkdVel(pkd,p);
             qv = pkdVel(pkd,q);

             vel_t dv2 = 0.0;
             for (int j=0; j<3; j++){
                dv2 += (pv[j]-qv[j])*(pv[j]-qv[j]);
             }

               //printf("%e \t %e \n", dv2, pkdMass(pkd,p)/pkdSoft(pkd,p));
             if (dv2 < pmass/pkdSoft(pkd,p)){
               // We have a merger!!
               printf("Merger!!! %e %e \n", pmass, qmass);

               //assert(pkd->idSelf == nnList[i].iPid);
              
               float newmass = pmass + qmass;
               float inv_newmass = 1./newmass;

               //printf("pos 1 %e \t pos 2 %e \t mean %e \n", pkdPos(pkd,p,0), pkdPos(pkd,q,0), pkdPos(pkd,p,0) - qmass*inv_newmass*nnList[i].dx);
               /* We can not update the position while in this loop, because we are still using 
                * information of the tree, which may be compromised.
                *
                * Instead, we save the new position and then update the merged BHs in a separated loop.
                * In the case that they are following the minimum potential of the nearby particles,
                * this update may not be needed.
                */
               BHFIELDS* pbh = pkdBH(pkd,p);
               pbh->newPos[0] = pkdPos(pkd,p,0) - qmass*inv_newmass*nnList[i].dx ;
               pbh->newPos[1] = pkdPos(pkd,p,1) - qmass*inv_newmass*nnList[i].dy ;
               pbh->newPos[2] = pkdPos(pkd,p,2) - qmass*inv_newmass*nnList[i].dz ;
               pbh->dInternalMass += pkdBH(pkd,q)->dInternalMass;
               //printf("Old velocity \t %e \t %e \t %e \n", pv[0], pv[1], pv[2]);
               for (int j=0; j<3; j++){
                  //pkdSetPos(pkd,p,j, (pmass*pkdPos(pkd,p,j) + qmass*pkdPos(pkd,q,j))*inv_newmass);
                  pv[j] = (pmass*pv[j] + qmass*qv[j])*inv_newmass;
               }
               //printf("New velocity \t %e \t %e \t %e \n", pv[0], pv[1], pv[2]);
               //printf("New position \t %e \t %e \t %e \n", pkdPos(pkd,p,0), pkdPos(pkd,p,1), pkdPos(pkd,p,2));

               float *mass_field = (float *)pkdField(p, pkd->oMass);
               *mass_field = newmass;

               //if (pkd->idSelf == nnList[i].iPid){
                  mass_field = (float *)pkdField(q, pkd->oMass);
                  *mass_field = 0.0f; // We do this to check that the deleted particles does not enter the gravity computation, as this
                                      //   will activate an assert

                  pkdDeleteParticle(pkd, q);
               //}
               
               // We only allow for one merger per particle per timestep
               p->bMarked=0;
               return;
               
             }
             // IA: this is not a good idea because we can not be sure that we are merging exactly the same particles
//          }else if (pkd->idSelf != nnList[i].iPid){
//             pv = pkdVel(pkd,p);
//             qv = pkdVel(pkd,q);
//
//             vel_t dv2 = 0.0;
//             for (int j=0; j<3; j++){
//                dv2 += (pv[j]-qv[j])*(pv[j]-qv[j]);
//             }
//
//             if (dv2 < 1.0/*qmass*//pkdSoft(pkd,q)){
//               float *mass_field = (float *)pkdField(p, pkd->oMass);
//               *mass_field = 0.0f;
//
//               pkdDeleteParticle(pkd, p);
//               return;
//             }
          }
       }
    }
}




int smReSmoothBHNode(SMX smx,SMF *smf, int iSmoothType) {
    PKD pkd = smx->pkd;
    int pj, pk, nCnt;
    double dx, dy, dz;
    double fDist2;
    PARTICLE* p;
    int nSmoothed=0;

   smx->nnListSize = 0;
   int nnListMax_p = NNLIST_INCREMENT;
    KDN *node;
    BND bnd_node;

    NN* nnList_p;
    nnList_p = (NN*)malloc(sizeof(NN)*nnListMax_p);

    // Here we store the pointers to the particle whose interaction need to be computed
    PARTICLE** sinks;
    sinks = malloc(64*sizeof(PARTICLE*)); // At most, the size of the bucket



    for (int i=NRESERVED_NODES; i<pkd->nNodes-1; i++){
      node = pkdTreeNode(pkd,i);
      if (!node->iLower && pkdNodeNbh(pkd,node)){ // We are in a bucket which contains a BH

         // Prepare the interaction list
         bnd_node = pkdNodeGetBnd(pkd,node);


         //printf("fBall %e nodeBall %e \n", fBall, pkdNodeBall(pkd,node));
         // Size of the ball that contains all possible particle interacting with this bucket

         double r[3];
         r[0] = bnd_node.fCenter[0];
         r[1] = bnd_node.fCenter[1];
         r[2] = bnd_node.fCenter[2];

         // First, we add all the particles whose interactions need to be computed
         int nActive = 0;
         float nodeBall = 0.;
#ifdef OPTIM_REORDER_IN_NODES
         int pStart = node->pLower + pkdNodeNgas(pkd,node);
#ifdef STAR_FORMATION
         pStart += pkdNodeNstar(pkd, node);
#endif
         int pEnd = node->pLower + pkdNodeNbh(pkd,node);
         //printf("pkdNodeNbh %d \n", pkdNodeNbh(pkd,node));
#else
         int pStart = node->pLower;
         int pEnd = node->pUpper+1;
#endif

         double fMax_shrink_x = 0.; 
         double fMax_shrink_y = 0.; 
         double fMax_shrink_z = 0.; 
         for (pj=pStart;pj<pEnd;++pj) {
             p = pkdParticle(pkd,pj);

             double dSoft = pkdSoft(pkd,p);

#ifdef OPTIM_AVOID_IS_ACTIVE
             if (p->bMarked){
#else
             if (pkdIsActive(pkd,p)) {
#endif

                const double x_disp = fabs(pkdPos(pkd,p,0) - bnd_node.fCenter[0]) + dSoft*2.;
                fMax_shrink_x = (x_disp > fMax_shrink_x) ? x_disp : fMax_shrink_x;

                const double y_disp = fabs(pkdPos(pkd,p,1) - bnd_node.fCenter[1]) + dSoft*2.;
                fMax_shrink_y = (y_disp > fMax_shrink_y) ? y_disp : fMax_shrink_y;

                const double z_disp = fabs(pkdPos(pkd,p,2) - bnd_node.fCenter[2]) + dSoft*2.;
                fMax_shrink_z = (z_disp > fMax_shrink_z) ? z_disp : fMax_shrink_z;

                if (nodeBall<dSoft) nodeBall=dSoft;
                sinks[nActive] = p; 
                nActive++;
             }
          }
         if (nActive==0) continue; // There is no elligible particle in this bucket, go to the next

         //printf("nodeBall %e nActive %d \n", nodeBall, nActive);
         nCnt = 0;
         //printf("%e %e \n", 2.*nodeBall, pkdNodeBall(pkd,node));
         int nCnt_own = nActive;
         //printf("start node %d %d \n", pkd->idSelf, i);

         //printf("nCnt_own %d \n", nCnt_own);

         double const fBall_x = bnd_node.fMax[0]+nodeBall;
         double const fBall_y = bnd_node.fMax[1]+nodeBall;
         double const fBall_z = bnd_node.fMax[2]+nodeBall;
         double fBall2 = fBall_x*fBall_x + fBall_y*fBall_y + fBall_z*fBall_z;
         double fBall2_shrink = fMax_shrink_x*fMax_shrink_x + fMax_shrink_y*fMax_shrink_y + fMax_shrink_z*fMax_shrink_z;


#ifdef OPTIM_EXTRA
         fBall2 = fBall2_shrink;
         bnd_node.fMax[0] = fMax_shrink_x;
         bnd_node.fMax[1] = fMax_shrink_y;
         bnd_node.fMax[2] = fMax_shrink_z;
#endif
    
         if (smx->bPeriodic) {
            double iStart[3], iEnd[3];
            for (int j=0;j<3;++j) {
                iStart[j] = d2i(floor((r[j] - bnd_node.fMax[j])/pkd->fPeriod[j] + 0.5));
                iEnd[j] = d2i(floor((r[j] + bnd_node.fMax[j])/pkd->fPeriod[j] + 0.5));
            }
            for (int ix=iStart[0];ix<=iEnd[0];++ix) {
                r[0] = bnd_node.fCenter[0] - ix*pkd->fPeriod[0];
                for (int iy=iStart[1];iy<=iEnd[1];++iy) {
                  r[1] = bnd_node.fCenter[1] - iy*pkd->fPeriod[1];
                  for (int iz=iStart[2];iz<=iEnd[2];++iz) {
                      r[2] = bnd_node.fCenter[2] - iz*pkd->fPeriod[2];
                          buildCandidateMergerList(smx, smf, node, bnd_node, &nCnt, r, fBall2, ix, iy, iz);
                  }
                }
            }
         }else{
            buildCandidateMergerList(smx, smf, node, bnd_node, &nCnt, r, fBall2, 0, 0, 0);
         }



              //printf("interaction list completed nCnt %d nCnt_own %d nActive  %d \n", nCnt, nCnt_own, nActive);

          // IA: Now we should have inside nnList all the particles in the bucket (sinks) and those of which can
          //  interact with them from other buckets (smx->nnList)
          //


             for (pj=0; pj<nCnt_own; pj++){
                PARTICLE * partj = sinks[pj];
                if (partj->bMarked){ // We need to double check this, as the particle may have been marked for deletion just before this
                   float fBall2_p = pkdSoft(pkd,partj)*pkdSoft(pkd,partj);
                   float dx_node = -pkdPos(pkd,partj,0)+bnd_node.fCenter[0];
                   float dy_node = -pkdPos(pkd,partj,1)+bnd_node.fCenter[1];
                   float dz_node = -pkdPos(pkd,partj,2)+bnd_node.fCenter[2];



                   int nCnt_p = 0;
                   for (pk=0;pk<nCnt;pk++){
                      dx = -dx_node + smx->nnList[pk].dx;
                      dy = -dy_node + smx->nnList[pk].dy;
                      dz = -dz_node + smx->nnList[pk].dz;

                      fDist2 = dx*dx + dy*dy + dz*dz;
                      //printf("fDist %e \t fBall2 %e \n", fDist2, fBall2_p);
                      if (fDist2 <= fBall2_p){
                         if (fDist2==0.)continue;

                         if (nCnt_p >= nnListMax_p) {
                             nnListMax_p += NNLIST_INCREMENT;
                             nnList_p = realloc(nnList_p,nnListMax_p*sizeof(NN));
                             assert(nnList_p != NULL);
                             }

                         nnList_p[nCnt_p].fDist2 = fDist2;
                         nnList_p[nCnt_p].dx = dx;
                         nnList_p[nCnt_p].dy = dy;
                         nnList_p[nCnt_p].dz = dz;
                         nnList_p[nCnt_p].pPart = smx->nnList[pk].pPart;
                         nnList_p[nCnt_p].iIndex = smx->nnList[pk].iIndex;
                         nnList_p[nCnt_p].iPid = smx->nnList[pk].iPid;



                         nCnt_p++;
                      }

                   }

                   smx->fcnSmooth(partj,pkdSoft(pkd,partj),nCnt_p,nnList_p,smf);

                }
          }



          nSmoothed += nCnt_own;

          for (pk=0;pk<nCnt;++pk) {
            if (smx->nnList[pk].iPid != pkd->idSelf) {
               // TODO: Do not forget to release previously acquired particles!!
               //mdlRelease(pkd->mdl,CID_PARTICLE,smx->nnList[pk].pPart);
            }
          }

		//mdlCacheCheck(pkd->mdl);
      }
   }
   free(nnList_p);
   free(sinks);
   return nSmoothed;
}

inline static KDN *getCell(PKD pkd, int iCell, int id) {
    if (id==pkd->idSelf) return pkdTreeNode(pkd,iCell);
    return mdlFetch(pkd->mdl,CID_CELL,iCell,id);
    }

/* This is mostly a version of buildInteractionList, but that only looks for BH particles.
 *
 * Probably, using C++ templates this could be done avoiding code duplications or extra overheads...
 */
void buildCandidateMergerList(SMX smx, SMF *smf, KDN* node, BND bnd_node, int *nCnt_tot, double r[3], double fBall2, int ix, int iy, int iz){
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    int id, sp, iCell, pStart, pEnd, pj;
    double dx, dy, dz, p_r[3], fDist2;
    KDN* kdn;
    BND bnd;
    struct stStack *S = smx->ST;
    int nCnt = *nCnt_tot;

   // We look for the biggest node that encloses the needed domain
   id = pkd->idSelf;

 // We can only take advantage of this if we are are in the original cell
#ifdef OPTIM_INVERSE_WALK
   kdn = pkdTreeNode(pkd,pkdNodeParent(pkd,node));
   bnd = pkdNodeGetBnd(pkd,kdn);
   if (ix==0 && iy==0 && iz==0){
      while((( fabs(bnd.fCenter[0] - r[0]) - bnd.fMax[0] + bnd_node.fMax[0] > 0  )||
             ( fabs(bnd.fCenter[1] - r[1]) - bnd.fMax[1] + bnd_node.fMax[1] > 0  )||
             ( fabs(bnd.fCenter[2] - r[2]) - bnd.fMax[2] + bnd_node.fMax[2] > 0  ))&&
             (pkdNodeParent(pkd,kdn)!=0)){
          kdn = pkdTreeNode(pkd, pkdNodeParent(pkd,kdn));
          bnd = pkdNodeGetBnd(pkd,kdn);
      }
   }else{
       kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = pkd->idSelf);
   }
#else
   kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = pkd->idSelf);
#endif

   //  Now we start the walk as usual
   sp = 0;
    while (1) {
        bnd = pkdNodeGetBnd(pkd, kdn);
      //MINDIST(&bnd,r,min2);
      //if (min2 > fBall2) {
      //    goto NoIntersect;
      //}
      for (int bnd_j=0; bnd_j<3; bnd_j++){
         if (fabs(bnd.fCenter[bnd_j]-r[bnd_j]) - bnd.fMax[bnd_j] - bnd_node.fMax[bnd_j] > 0. ) goto NoIntersect;
      }


      /*
      ** We have an intersection to test.
      */
      if (kdn->iLower) {
          int idUpper,iUpper;
          pkdGetChildCells(kdn,id,id,iCell,idUpper,iUpper);
          kdn = getCell(pkd,iCell,id);
          S[sp].id = idUpper;
          S[sp].iCell = iUpper;
          S[sp].min = 0.0;
          ++sp;
          continue;
      }
      else {
          if (id == pkd->idSelf) {
#ifdef OPTIM_REORDER_IN_NODES
            pStart = kdn->pLower+pkdNodeNgas(pkd,kdn);
#ifdef STAR_FORMATION
            pStart += pkdNodeNstar(pkd,kdn);
#endif
            pEnd = pStart + pkdNodeNbh(pkd,kdn);;
#else
            pStart = kdn->pLower;
            pEnd = kdn->pUpper+1;
#endif
            //printf("pEnd %d \n", pEnd);
            for (pj=pStart;pj<pEnd;++pj) {
                p = pkdParticle(pkd,pj);
#ifndef OPTIM_REORDER_IN_NODES
            if (!pkdIsBH(pkd,p)) continue;
#endif
                pkdGetPos1(pkd,p,p_r);
                dx = r[0] - p_r[0];
                dy = r[1] - p_r[1];
                dz = r[2] - p_r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2) {
                  if (nCnt >= smx->nnListMax) {
                      smx->nnListMax += NNLIST_INCREMENT;
                      smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
                      //printf("realloc \n");
                      assert(smx->nnList != NULL);
                      }
                  smx->nnList[nCnt].fDist2 = fDist2;
                  smx->nnList[nCnt].dx = dx;
                  smx->nnList[nCnt].dy = dy;
                  smx->nnList[nCnt].dz = dz;
                  smx->nnList[nCnt].pPart = p;
                  smx->nnList[nCnt].iIndex = pj;
                  smx->nnList[nCnt].iPid = pkd->idSelf;
                  ++nCnt;
                  }
                }
            }
          else {
#ifdef OPTIM_REORDER_IN_NODES
            pStart = kdn->pLower+pkdNodeNgas(pkd,kdn);
#ifdef STAR_FORMATION
            pStart += pkdNodeNstar(pkd,kdn);
#endif
            pEnd = pStart + pkdNodeNbh(pkd,kdn);;
#else
            pStart = kdn->pLower;
            pEnd = kdn->pUpper+1;
#endif
            for (pj=pStart;pj<pEnd;++pj) {
                p = mdlFetch(mdl,CID_PARTICLE,pj,id);
#ifndef OPTIM_REORDER_IN_NODES
            if (!pkdIsBH(pkd,p)) continue;
#endif
                pkdGetPos1(pkd,p,p_r);
                dx = r[0] - p_r[0];
                dy = r[1] - p_r[1];
                dz = r[2] - p_r[2];
                fDist2 = dx*dx + dy*dy + dz*dz;
                if (fDist2 <= fBall2) {
                  if (nCnt >= smx->nnListMax) {
                      smx->nnListMax += NNLIST_INCREMENT;
                      smx->nnList = realloc(smx->nnList,smx->nnListMax*sizeof(NN));
                      //printf("realloc \n");
                      assert(smx->nnList != NULL);
                      }

                  smx->nnList[nCnt].fDist2 = fDist2;
                  smx->nnList[nCnt].dx = dx;
                  smx->nnList[nCnt].dy = dy;
                  smx->nnList[nCnt].dz = dz;

                  // IA: As we do not allow merger across boundaries, we do not need to acquire the remote particle
                  smx->nnList[nCnt].pPart = p ;

                  smx->nnList[nCnt].iIndex = pj;
                  smx->nnList[nCnt].iPid = id;
                  ++nCnt;
                  }
                }
            }
          }
    NoIntersect:
      if (sp) {
          --sp;
          id = S[sp].id;
          iCell = S[sp].iCell;
          kdn = getCell(pkd,iCell,id);
          }
      else break;
    }

    *nCnt_tot = nCnt;

}

#endif
