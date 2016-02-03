#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include "fof.h"
#include "pkd.h"
#include "group.h"


static inline int getCell(PKD pkd,int iCell,int id,float *pcOpen,KDN **pc) {
    KDN *c;
    int nc;
    assert(iCell > 0);
    assert(id >= 0);
    if (id == pkd->idSelf) c = pkdTreeNode(pkd,iCell);
    else c = CAST(KDN *,mdlFetch(pkd->mdl,CID_CELL,iCell,id));
    *pc = c;
    if (c->bRemote|c->bTopTree) nc = 1000000000; /* we never allow pp with this cell */
    else nc = c->pUpper - c->pLower + 1;
    *pcOpen = c->bMax * pkd->fiCritTheta;
    return nc;
    }


void pkdFofGatherLocal(PKD pkd,int *S,FLOAT fBall2,FLOAT r[3],int32_t iGroup,
    int *piTail,uint32_t *Fifo,
    int *pbCurrFofContained,int *pnCurrFofParticles,
    FLOAT *fMinFofContained,FLOAT *fMaxFofContained) {
    KDN *kdn;
    PARTICLE *p;
    double p_r[3];
    FLOAT min2,dx,dy,dz,fDist2;
    int sp = 0;
    int iCell,pj,pEnd,j;
    const BND *bnd;
    int32_t iPartGroup;

    kdn = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
        bnd = pkdNodeBnd(pkd, kdn);
	MINDIST(bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	    }
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iLower);
	    S[sp++] = iCell+1;
	    continue;
	    }
	else {
	    pEnd = kdn->pUpper;
	    for (pj=kdn->pLower;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		iPartGroup = pkdGetGroup(pkd,p);
		if (iPartGroup) continue;		    
		pkdGetPos1(pkd,p,p_r);
		dx = r[0] - p_r[0];
		dy = r[1] - p_r[1];
		dz = r[2] - p_r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= fBall2) {
		    /*
		    **  Mark particle and add it to the do-fifo
		    */
		    pkdSetGroup(pkd,p,iGroup);
		    Fifo[(*piTail)++] = pj;
		    ++(*pnCurrFofParticles);
		    if (*pbCurrFofContained) {
			for (j=0;j<3;++j) {
			    if (p_r[j] < fMinFofContained[j]) {
				*pbCurrFofContained = 0;
				break;
				}
			    else if (p_r[j] > fMaxFofContained[j]) {
				*pbCurrFofContained = 0;
				break;
				}
			    }
			}
		    }
		}
	    }
    NoIntersect:
	if (sp) kdn = pkdTreeNode(pkd,iCell = S[--sp]);
	else break;
	}
    }



static void iOpenRemoteFof(PKD pkd,KDN *k,CL cl,CLTILE tile,float dTau2) {
    float dx,minbnd2,kOpen;
    CL_BLK *blk;
    int i,n,nLeft,iOpen;
    const BND *kbnd;

    kbnd = pkdNodeBnd(pkd,k);
    kOpen = kbnd->fMax[0] + kbnd->fMax[1] + kbnd->fMax[2]; /* Manhatten metric */
    blk = tile->blk;
    for(nLeft=tile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
	n = nLeft ? cl->lst.nPerBlock : tile->lstTile.nInLast;
	for(i=0; i<n; ++i) {
	    if (blk->idCell.i[i] > pkd->idSelf) iOpen = 10;  /* ignore this cell, but this never ignores the top tree */
	    else {
		minbnd2 = 0;
		dx = kbnd->fCenter[0] - kbnd->fMax[0] -  blk->xCenter.f[i] - blk->xOffset.f[i] - blk->xMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->xCenter.f[i] + blk->xOffset.f[i] - blk->xMax.f[i] - kbnd->fCenter[0] - kbnd->fMax[0];
		if (dx > 0) minbnd2 += dx*dx;
		dx = kbnd->fCenter[1] - kbnd->fMax[1] - blk->yCenter.f[i] - blk->yOffset.f[i] - blk->yMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->yCenter.f[i] + blk->yOffset.f[i] - blk->yMax.f[i] - kbnd->fCenter[1] - kbnd->fMax[1];
		if (dx > 0) minbnd2 += dx*dx;
		dx = kbnd->fCenter[2] - kbnd->fMax[2] - blk->zCenter.f[i] - blk->zOffset.f[i] - blk->zMax.f[i];
		if (dx > 0) minbnd2 += dx*dx;
		dx = blk->zCenter.f[i] + blk->zOffset.f[i] - blk->zMax.f[i] - kbnd->fCenter[2] - kbnd->fMax[2];
		if (dx > 0) minbnd2 += dx*dx;
       		if (minbnd2 > dTau2) iOpen = 10;  /* ignore this cell */
		else if (k->iLower == 0) {
		    if (blk->iLower.i[i] == 0) iOpen = 1;
		    else iOpen = 3;
		    }
		else if (kOpen > blk->cOpen.f[i] || blk->iLower.i[i] == 0) iOpen = 0;
		else iOpen = 3;
		}
	    blk->iOpen.i[i] = iOpen;
	    }
	}
    }


static void addChildFof(PKD pkd, CL cl, int iChild, int id, float *fOffset) {
    int idLower, iLower, idUpper, iUpper, iCache;
    float cOpen;
    KDN *c;
    int nc = getCell(pkd,iChild,id,&cOpen,&c);
    const BND *cbnd = pkdNodeBnd(pkd,c);
    iCache = 0;
    cOpen = cbnd->fMax[0] + cbnd->fMax[1] + cbnd->fMax[2]; /* Manhatten metric */
    pkdGetChildCells(c,id,idLower,iLower,idUpper,iUpper);
    clAppend(cl,iCache,id,iChild,idLower,iLower,idUpper,iUpper,nc,cOpen,
	pkdNodeMom(pkd,c)->m,4.0f*c->fSoft2,c->r,fOffset,cbnd->fCenter,cbnd->fMax);
    }


void pkdFofRemoteSearch(PKD pkd,double dTau2) {
    KDN *kdnSelf,*kdn,*c,*k;
    BND *bndSelf,*bnd;
    PARTICLE *p;
    CLTILE cltile;
    CL clTemp;
    double xj,yj,zj,d2;
    int npi;
    double *xi,*yi,*zi;
    uint32_t pjGroup;
    uint32_t *piGroup;
    uint32_t pi,pj;
    int iRemote;    
    int sp,i,j,ix,iy,iz,bRep;
    int idSelf,iTop,iCell,id,iCellLo,idLo,iCellUp,idUp,iSib,iCheckCell,iCheckLower;
    int jTile,M,iStack;
    float fOffset[3];


    /*
    ** Allocate the vectors to be large enough to handle all particles in a bucket.
    ** M should be aligned to 4*sizeof(double)
    */
    M = pkd->param.nBucket;
    xi = malloc(M*sizeof(double));
    mdlassert(pkd->mdl,xi != NULL);
    yi = malloc(M*sizeof(double));
    mdlassert(pkd->mdl,yi != NULL);
    zi = malloc(M*sizeof(double));
    mdlassert(pkd->mdl,zi != NULL);
    piGroup = malloc(M*sizeof(uint32_t));
    mdlassert(pkd->mdl,piGroup != NULL);

    iStack = 0;
    clClear(pkd->cl);

    kdnSelf = pkdTreeNode(pkd,ROOT);
    bndSelf = pkdNodeBnd(pkd, kdnSelf);
    idSelf = mdlSelf(pkd->mdl);
    iTop = pkd->iTopTree[ROOT];
    id = idSelf;
    for (j=0;j<3;++j) fOffset[j] = 0.0;
    /*
    ** Add all siblings of the top tree down to local root (but not including it) to 
    ** the checklist.
    */
    iCell = iTop;
    c = pkdTreeNode(pkd,iCell);
    while (c->bTopTree) {
	pkdGetChildCells(c,id,idLo,iCellLo,idUp,iCellUp);
	if (idLo == idSelf) { 
	    c = pkdTreeNode(pkd,iCellLo);
	    bnd = pkdNodeBnd(pkd,c);
	    for (j=0;j<3;++j) {
		if (fabs(bndSelf->fCenter[j] - bnd->fCenter[j]) > bnd->fMax[j]) {
		    addChildFof(pkd,pkd->cl,iCellLo,idLo,fOffset);
		    id = idUp;
		    assert(id == idSelf);
		    iCell = iCellUp;
		    c = pkdTreeNode(pkd,iCell);
		    goto NextCell;
		    }
		}
	    }
	addChildFof(pkd,pkd->cl,iCellUp,idUp,fOffset);
	iCell = iCellLo;
	id = idLo;
	assert(id == idSelf);
    NextCell:
	;
	}
    /*
    ** Add all replica global roots to the checklist for periodic BCs.
    */
    if (pkd->param.bPeriodic) {
	int nReps = pkd->param.nReplicas;
	for (ix=-nReps;ix<=nReps;++ix) {
	    fOffset[0] = ix*pkd->fPeriod[0];
	    for (iy=-nReps;iy<=nReps;++iy) {
		fOffset[1] = iy*pkd->fPeriod[1];
		for (iz=-nReps;iz<=nReps;++iz) {
		    fOffset[2] = iz*pkd->fPeriod[2];
		    bRep = ix || iy || iz;
		    if (bRep) addChildFof(pkd,pkd->cl,iTop,idSelf,fOffset);
		    }
		}
	    }
	}

    iCell = ROOT;
    id = idSelf;
    k = pkdTreeNode(pkd,iCell);
    /*
    ** The checklist is ready for a bound-bound walk of the remote particles.
    */
    while (1) {
	while (1) {
	    /*
	    ** Process the Checklist for the cell pointed to by k.
	    */
	    clClear(pkd->S[iStack+1].cl);
	    do {
		CL_LOOP(pkd->cl,cltile) {
#ifdef USE_SIMD_OPEN
		    iOpenRemoteFof(pkd,k,pkd->cl,cltile,dTau2);
#else
		    iOpenRemoteFof(pkd,k,pkd->cl,cltile,dTau2);
#endif
		    }
		clClear(pkd->clNew);
		CL_LOOP(pkd->cl,cltile) {
		    CL_BLK *blk = cltile->blk;
		    int nLeft;
		    for(nLeft=cltile->lstTile.nBlocks; nLeft>=0; --nLeft,blk++) {
			int n = nLeft ? pkd->cl->lst.nPerBlock : cltile->lstTile.nInLast;
			for (jTile=0;jTile<n;++jTile) {
			    switch (blk->iOpen.i[jTile]) {
			    case 0:
				/*
				** This checkcell stays on the checklist.
				*/
				clAppendItem(pkd->S[iStack+1].cl,blk,jTile);
				break;
			    case 1:
				/*
				** We check individual particles against each other here.
				*/
				iCheckCell = blk->iCell.i[jTile];
				id = blk->idCell.i[jTile];
				if (id == pkd->idSelf) {
				    printf("local checkcell:%d\n",iCheckCell);
				    c = pkdTreeNode(pkd,iCheckCell);
				    }
				else c = CAST(KDN *,mdlFetch(pkd->mdl,CID_CELL,iCheckCell,id));
				/*
				** Convert all the coordinates in the k-cell and store them in vectors.
				*/
				for (pi=k->pLower,npi=0;pi<=k->pUpper;++pi,++npi) {
				    p = pkdParticle(pkd,pi);
				    pkdGetPos3(pkd,p,xi[npi],yi[npi],zi[npi]);
				    piGroup[npi] = pkdGetGroup(pkd,p);
				    }
				for (pj=c->pLower;pj<=c->pUpper;++pj) {
				    if (id == pkd->idSelf) p = pkdParticle(pkd,pj);
				    else p = CAST(PARTICLE *,mdlFetch(pkd->mdl,CID_PARTICLE,pj,id));
				    pjGroup = pkdGetGroup(pkd,p);
				    pkdGetPos3(pkd,p,xj,yj,zj);
				    xj += blk->xOffset.f[jTile];
				    yj += blk->yOffset.f[jTile];
				    zj += blk->zOffset.f[jTile];
				    /*
				    ** The following could be vectorized over the vectors xi,yi and zi!
				    */
				    for (i=0;i<npi;++i) {
					d2 = (xj-xi[i])*(xj-xi[i]) + (yj-yi[i])*(yj-yi[i]) + (zj-zi[i])*(zj-zi[i]);
					if (d2 < dTau2) {
					    /*
					    ** We have found a remote group link!
					    ** Check if it is already linked to this group and if not add the link.
					    */
					    if (pjGroup == 0) printf("UNGROUPED PARTICLE FOUND at %.15g < %.15g\n",d2,dTau2);
					    assert(pjGroup > 0);
					    iRemote = pkd->ga[piGroup[i]].iLink;
					    while (iRemote) {
						if (pkd->tmpFofRemote[iRemote].key.iIndex == pjGroup && 
						    pkd->tmpFofRemote[iRemote].key.iPid == id) break;
						iRemote = pkd->tmpFofRemote[iRemote].iLink;
						}
					    if (!iRemote) {
						/*
						** Add this remote group to the list.
						*/
						iRemote = pkd->iRemoteGroup++;
						/*
						** Make sure we don't run out of ephemeral storage!
						** (we really should be ok as long as our nMinMembers is greater than 1 or 2)
						*/
						assert(pkd->iRemoteGroup < pkd->nMaxRemoteGroups);
						pkd->tmpFofRemote[iRemote].key.iIndex = pjGroup;
						assert(pjGroup > 0);
						pkd->tmpFofRemote[iRemote].key.iPid = id;
						pkd->tmpFofRemote[iRemote].iLink = pkd->ga[piGroup[i]].iLink;
						pkd->ga[piGroup[i]].iLink = iRemote;
						}
					    }
					}
				    }
				break;
			    case 3:
				/*
				** Open the cell.
				** We could do a prefetch here for non-local
				§				** cells.
				*/
				iCheckCell = blk->iCell.i[jTile];    assert(iCheckCell >= 0);
				iCheckLower = blk->iLower.i[jTile];  assert(iCheckLower > 0);

				fOffset[0] = blk->xOffset.f[jTile];
				fOffset[1] = blk->yOffset.f[jTile];
				fOffset[2] = blk->zOffset.f[jTile];

				addChildFof(pkd,pkd->clNew,blk->iLower.i[jTile],blk->idLower.i[jTile],fOffset);
				addChildFof(pkd,pkd->clNew,blk->iUpper.i[jTile],blk->idUpper.i[jTile],fOffset);
				break;
			    case 10:
				/*
				** This checkcell is removed from the checklist since it has no overlap with the current cell.
				*/
				break;
			    default:
				assert(0);
				} /* end of switch */
			    } /* end of for (jTile) */
			} /* end of for (nLeft) */
		    } /* end of CL_LOOP */
		clTemp = pkd->cl;
		pkd->cl = pkd->clNew;
		assert(pkd->cl!=NULL);
		pkd->clNew = clTemp;
		} while (clCount(pkd->cl));
	    /*
	    ** Done processing of the Checklist.
	    ** Now prepare to proceed to the next deeper
	    ** level of the tree.
	    */
	    if (k->iLower == 0) break;
	    iCell = k->iLower;
	    k = pkdTreeNode(pkd,iCell);
	    /*
	    ** Push the sibling onto the stack.
	    */
	    iSib = iCell+1;
	    c = pkdTreeNode(pkd,iSib);
	    clClone(pkd->cl,pkd->S[iStack+1].cl);
	    ++iStack;
	    assert(iStack < pkd->nMaxStack);
	    pkd->S[iStack].iNodeIndex = iSib;
	    }
	/*
	** Now the checklist should be empty and we should have dealt with all 
	** links going across the processor domains for this bucket.
	*/
	assert(clCount(pkd->S[iStack+1].cl)==0);
	/*
	** Get the next cell to process from the stack. 
	*/
	if (iStack == -1) break;  /* we are done! */
	k = pkdTreeNode(pkd,iCell = pkd->S[iStack].iNodeIndex);
	/*
	** Grab the checklist from the stack.
	*/
	clTemp = pkd->cl;
	assert(clCount(pkd->cl) == 0);
	pkd->cl = pkd->S[iStack].cl;
	assert(pkd->cl!=NULL);
	pkd->S[iStack].cl = clTemp;
	--iStack;
	}
    }


void updateFofBound(FLOAT *fMinFofContained,FLOAT *fMaxFofContained,BND *bnd) {
    int j;
    for (j=0;j<3;++j) {
	FLOAT bmin = bnd->fCenter[j] - bnd->fMax[j];
	FLOAT bmax = bnd->fCenter[j] + bnd->fMax[j];
	assert(bmin <= bmax);
	if (bmin > fMaxFofContained[j] || bmax < fMinFofContained[j]) continue;
	else if (bmax < fMaxFofContained[j] && bmin > fMinFofContained[j]) {
	    if (bmin - fMinFofContained[j] > fMaxFofContained[j] - bmax)
		fMaxFofContained[j] = bmin;
	    else 
		fMinFofContained[j] = bmax;
	    }
	else {
	    if (bmin < fMaxFofContained[j]) 
		fMaxFofContained[j] = bmin;
	    if (bmax > fMinFofContained[j]) 
		fMinFofContained[j] = bmax;		
	    }
	}
    }


int pkdNewFof(PKD pkd,double dTau2,int nMinMembers) {
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    double p_r[3];
    int32_t iGroup,*pGroup;
    int pn,i,j;
    KDN *kdnSelf;
    BND *bndSelf,*bnd,*bndTop;
    int *S;
    uint32_t iHead;
    uint32_t iTail;
    uint32_t *Fifo;
    uint32_t *iFofMap;
    int bCurrFofContained;
    uint32_t nCurrFofParticles;
    FLOAT fMinFofContained[3];
    FLOAT fMaxFofContained[3];    

    assert(pkd->oGroup); /* Validate memory model */
    S = malloc(1024*sizeof(int));
    assert(S != NULL);
    /*
    ** Set up the bounds for the FOF groups that are certainly contained in the domain.
    ** This is a little trickier for domains which could potentially overlap a bit.
    ** For now I assume that the domains do NOT overlap, but the calculation for overlapping
    ** domains just involves a tree walk.
    */
    kdnSelf = pkdTreeNode(pkd,ROOT);
    bndSelf = pkdNodeBnd(pkd, kdnSelf);
    for (j=0;j<3;++j) {
	fMinFofContained[j] = bndSelf->fCenter[j] - bndSelf->fMax[j];
	fMaxFofContained[j] = bndSelf->fCenter[j] + bndSelf->fMax[j];
	}
#if (0)
    /*
    ** The following code would allow doing fof on a substep, which isn't 
    ** forseen in the near future. We can test it at a later stage.
    ** Check bounds against all siblings of the top tree down to local root.
    */
    iCell = iTop;
    c = pkdTreeNode(pkd,iCell);
    bndTop = pkdNodeBnd(pkd,c);
    while (c->bTopTree) {
	pkdGetChildCells(c,id,idLo,iCellLo,idUp,iCellUp);
	if (idLo == idSelf) { 
	    c = pkdTreeNode(pkd,iCellLo);
	    bnd = pkdNodeBnd(pkd,c);
	    for (j=0;j<3;++j) {
		if (fabs(bndSelf->fCenter[j]-bnd->fCenter[j]) > bnd->fMax[j]) {
		    /*
		    ** Check bounds against this sibling.
		    */
		    updateFofBound(fMinFofContained,fMaxFofContained,bnd);
		    id = idUp;
		    assert(id == idSelf);
		    iCell = iCellUp;
		    c = pkdTreeNode(pkd,iCell);
		    goto NextCell;
		    }
		}
	    }
	assert(idUp == idSelf);
	c = pkdTreeNode(pkd,iCellUp);
	bnd = pkdNodeBnd(pkd,c);
	/*
	** Check bounds against this sibling.
	*/
	updateFofBound(fMinFofContained,fMaxFofContained,bnd);
	iCell = iCellLo;
	id = idLo;
    NextCell:
	;
	}
    /*
    ** Check bounds against first replica global roots for periodic BCs.
    */
    if (pkd->param.bPeriodic) {
	BND rbnd;
	for (j=0;j<3;++j) rbnd.fMax[j] = bndTop->fMax[j];
	for (ix=-1;ix<=1;++ix) {
	    fOffset[0] = ix*pkd->fPeriod[0];
	    for (iy=-1;iy<=1;++iy) {
		fOffset[1] = iy*pkd->fPeriod[1];
		for (iz=-1;iz<=1;++iz) {
		    fOffset[2] = iz*pkd->fPeriod[2];
		    bRep = ix || iy || iz;
		    if (bRep) {
			/*
			** Check bounds against this replica.
			*/
			for (j=0;j<3;++j) 
			    rbnd.fCenter[j] = bndTop->fCenter[j] + fOffset[j];
			updateFofBound(fMinFofContained,fMaxFofContained,&rbnd);
			}
		    }
		}
	    }
	}
#endif
    /*
    ** Finally make the contained region be dTau smaller on each side.
    */
    for (j=0;j<3;++j) {
	fMinFofContained[j] += sqrt(dTau2);
	fMaxFofContained[j] -= sqrt(dTau2);
	}
    /*
    ** Clear the group numbers!
    */
    for (pn=0;pn<pkd->nLocal;++pn) {
	p = pkdParticle(pkd,pn);
	pkdSetGroup(pkd,p,0);
	}
    /*
    ** The following *just* fits into ephemeral storage of 8bytes/particle.
    */
    assert(EPHEMERAL_BYTES >= 8);
    Fifo = (uint32_t *)(pkd->pLite);
    iFofMap = &Fifo[pkd->nLocal];
    pkd->nGroups = 0;    
    iGroup = 0;
    iFofMap[iGroup] = 0;
    pkd->nLocalGroups = 0;
    for (pn=0;pn<pkd->nLocal;++pn) {
	p = pkdParticle(pkd,pn);
	if (pkdGetGroup(pkd,p)) continue;
	++iGroup;
	/*
	** Mark particle and add it to the do-fifo
	*/
	iHead = iTail = 0;
	Fifo[iTail++] = pn;
	pkdSetGroup(pkd,p,iGroup);
	nCurrFofParticles = 1;
	bCurrFofContained = 1;
	pkdGetPos1(pkd,p,p_r);
	for (j=0;j<3;++j) {
	    if (p_r[j] < fMinFofContained[j]) {
		bCurrFofContained = 0;
		break;
		}
	    else if (p_r[j] > fMaxFofContained[j]) {
		bCurrFofContained = 0;
		break;
		}
	    }
	while (iHead != iTail) {
	    int pi = Fifo[iHead++];
	    p = pkdParticle(pkd,pi);
	    pkdGetPos1(pkd,p,p_r);
	    pkdFofGatherLocal(pkd,S,dTau2,p_r,iGroup,&iTail,Fifo,
		&bCurrFofContained,&nCurrFofParticles,fMinFofContained,fMaxFofContained);
	    }
	assert(iTail < pkd->nLocal);
	/*
	** Now check if this fof group is contained and has fewer than nMinFof particles.
	*/
	if (bCurrFofContained && nCurrFofParticles < nMinMembers) {
	    iFofMap[iGroup] = 0;
	    }
	else {
	    iFofMap[iGroup] = ++pkd->nLocalGroups;
	    }
	}
    /*
    ** Renumber the group assignments for the particles (having removed some small groups).
    */
    for (pn=0;pn<pkd->nLocal;++pn) {
	p = pkdParticle(pkd,pn);
	pkdSetGroup(pkd,p,iFofMap[pkdGetGroup(pkd,p)]);
	}
    printf("%3d:cull initial small groups from %d to nGroups=%d\n",pkd->idSelf,iGroup,pkd->nLocalGroups);
    pkd->nGroups = pkd->nLocalGroups + 1;
    free(S);  /* this stack is no longer needed */
    iFofMap = NULL; /* done with the temporary map of group numbers */ 
    Fifo = NULL;  /* done with the Fifo, can use the storage for other stuff now */
    /*
    ** Create initial group table. The assert below is a very minimal requirement as it doesn't account for remote
    ** links (tmpFofRemote). However, we check this again everytime we add a new remote link.
    */
    assert(sizeof(*pkd->ga)*pkd->nGroups+sizeof(*pkd->tmpFofRemote) <= EPHEMERAL_BYTES*pkd->nStore);
    pkd->nMaxRemoteGroups = (EPHEMERAL_BYTES*pkd->nStore - sizeof(*pkd->ga)*pkd->nGroups) / sizeof(*pkd->tmpFofRemote);
    pkd->ga = (struct smGroupArray *)(pkd->pLite);
    pkd->tmpFofRemote = (FOFRemote *)&pkd->ga[pkd->nGroups];
    for(i=0;i<pkd->nGroups;++i) {
	pkd->ga[i].id.iIndex = i;
	pkd->ga[i].id.iPid = pkd->idSelf;
	pkd->ga[i].iGid = i;
	pkd->ga[i].iLink = 0;   /* this is a linked list of remote groups linked to this local group */
	}
    pkd->iRemoteGroup = 1;  /* The first entry is a dummy one for a null index */
    /*
    ** Now lets go looking for local particles which have a remote neighbor that is part of 
    ** a group.
    */
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
    pkdFofRemoteSearch(pkd,dTau2);
    mdlFinishCache(mdl,CID_PARTICLE);
    }


/*
** When we virtual fetch a name of one of the groups we may already have fetched the
** same key on the same processor. For this reason we need to initialize the 
** virtual fetch to something that will for sure be updated on the first fetch.
** the second fetch will only update it if the new name is actually smaller than the
** one set by the first fetch.
*/
static void initNames(void *vctx, void *v) {
    struct smGroupArray *g = (struct smGroupArray *)v;
    g->id.iPid = INT32_MAX;
    g->id.iIndex = INT32_MAX;
    }
static void combNames(void *vctx, void *v1, void *v2) {
    struct smGroupArray *g1 = (struct smGroupArray *)v1;
    struct smGroupArray *g2 = (struct smGroupArray *)v2;
    if ( g1->id.iPid>g2->id.iPid || (g1->id.iPid==g2->id.iPid && g1->id.iIndex>g2->id.iIndex) ) {
	g1->id.iPid = g2->id.iPid;
	g1->id.iIndex = g2->id.iIndex;
	}
    }


int pkdFofPhases(PKD pkd) {
    MDL mdl = pkd->mdl;
    int bMadeProgress=0;
    int iIndex,iPid,iRemote,iLink,i;
    remoteID name;
    struct smGroupArray *pRemote;

    /*
    ** Phase 1: fetch remote names.
    */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->ga,sizeof(struct smGroupArray),pkd->nGroups);
    for (iRemote=1;iRemote<pkd->iRemoteGroup;++iRemote) {
	iIndex = pkd->tmpFofRemote[iRemote].key.iIndex;
	assert(iIndex > 0);
	iPid = pkd->tmpFofRemote[iRemote].key.iPid;
	pRemote = mdlFetch(mdl,CID_GROUP,iIndex,iPid);
	/*
	** Now update the name in the local table of remote group links.
	*/
	assert(pRemote->id.iIndex > 0);
	pkd->tmpFofRemote[iRemote].name.iIndex = pRemote->id.iIndex;
	pkd->tmpFofRemote[iRemote].name.iPid = pRemote->id.iPid;
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Phase 2: update to unique names.
    */
    for(i=1;i<pkd->nGroups;++i) {
	iLink = pkd->ga[i].iLink;
	if (iLink) {
	    name.iIndex = pkd->ga[i].id.iIndex;
	    name.iPid = pkd->ga[i].id.iPid;
	    /*
	    ** Find current master name (this is the lowest iPid,iIndex pair found).
	    */
	    while (iLink) {
		iPid = pkd->tmpFofRemote[iLink].name.iPid;
		iIndex = pkd->tmpFofRemote[iLink].name.iIndex;
		iLink = pkd->tmpFofRemote[iLink].iLink;
		if (iPid < name.iPid) {
		    name.iPid = iPid;
		    name.iIndex = iIndex;
		    bMadeProgress = 1;
		    }
		else if (iPid == name.iPid && iIndex < name.iIndex) {
		    name.iIndex = iIndex;
		    bMadeProgress = 1;
		    }
		}
	    assert(name.iIndex > 0);
	    pkd->ga[i].id.iIndex = name.iIndex;
	    pkd->ga[i].id.iPid = name.iPid;
	    iLink = pkd->ga[i].iLink;
	    while (iLink) {
		if (pkd->tmpFofRemote[iLink].name.iPid == name.iPid &&
		    pkd->tmpFofRemote[iLink].name.iIndex == name.iIndex) {
		    /*
		    ** There is no update to be made. We mark this with a 
		    ** a dummy iPid! It is ok to destroy it here since it 
		    ** will be refetched in phase 1 (all remote links names
		    ** are refetched).
		    */
		    pkd->tmpFofRemote[iLink].name.iPid = -1;
		    }
		else {
		    pkd->tmpFofRemote[iLink].name.iPid = name.iPid;
		    pkd->tmpFofRemote[iLink].name.iIndex = name.iIndex;
		    }
		iLink = pkd->tmpFofRemote[iLink].iLink;
		}
	    }
	}
    /*
    ** Phase 3: propagate.
    */
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->ga,sizeof(struct smGroupArray),pkd->nGroups,
	NULL, initNames, combNames );
    for(i=1;i<pkd->nGroups;++i) {
	iLink = pkd->ga[i].iLink;
	if (iLink) {
	    while (iLink) {
		name.iPid = pkd->tmpFofRemote[iLink].name.iPid;
		if (name.iPid >= 0) {
		    name.iIndex = pkd->tmpFofRemote[iLink].name.iIndex;
		    iPid = pkd->tmpFofRemote[iLink].key.iPid;
		    iIndex = pkd->tmpFofRemote[iLink].key.iIndex;
		    pRemote = mdlVirtualFetch(mdl,CID_GROUP,iIndex,iPid);
		    assert(pRemote);
		    if (name.iPid < pRemote->id.iPid) {
			pRemote->id.iPid = name.iPid;
			pRemote->id.iIndex = name.iIndex;
			}
		    else if (name.iPid == pRemote->id.iPid && name.iIndex < pRemote->id.iIndex) {
			pRemote->id.iIndex = name.iIndex;
			}
		    }
		iLink = pkd->tmpFofRemote[iLink].iLink;   
		}
	    }
	}
    mdlFinishCache(mdl,CID_GROUP);
    return(bMadeProgress);
    }


uint64_t pkdFofFinishUp(PKD pkd,int nMinGroupSize) {
    int i;

    /*
    ** Merge local groups in order to have only one group per process informing the master of its existence.
    ** For a given name in the table, we need to know if some other entry has the same name and remap that 
    ** index.
    */
    for(i=1; i<pkd->nGroups; ++i) {
	assert(pkd->ga[i].id.iIndex > 0);
	}
    pkd->nGroups = pkdGroupCombineDuplicateIds(pkd,pkd->nGroups,pkd->ga,1);
    pkd->nGroups = pkdPurgeSmallGroups(pkd,pkd->nGroups,pkd->ga,nMinGroupSize);
    /*
    ** This last call to pkdGroupCounts just gets the final group counts in pkd->ga restored
    ** which is optional if group stats are calculated later (which does the same thing again).
    ** pkd->nLocalGroups is kept consistent within pkdPurgeSmallGroups as well.
    */
    pkdGroupCounts(pkd,pkd->nGroups,pkd->ga);
    return pkd->nLocalGroups;
    }
