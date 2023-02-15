#ifdef FEEDBACK
#include <numeric>
#include "starformation/feedback.h"
#include "master.h"
using blitz::all;

void MSR::SetFeedbackParam() {
    const double dHydFrac = param.dInitialH;
    const double dnHToRho = MHYDR / dHydFrac / param.units.dGmPerCcUnit;
    param.dSNFBDu = param.dSNFBDT*dTuFac;
    if (!param.bRestart) {
        param.dSNFBDelay *=  SECONDSPERYEAR/param.units.dSecUnit ;
        param.dSNFBNumberSNperMass *= 8.73e15 / param.units.dErgPerGmUnit / 1.736e-2;
        param.dSNFBEffnH0 *= dnHToRho;
    }
}

/* Function that will be called with the information of all the neighbors.
 * Here we compute the probability of explosion, and we add the energy to the
 * surroinding gas particles
 */
void smSNFeedback(PARTICLE *pIn,float fBall,int nSmooth,NN *nnList,SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    int i;

    auto totMass = std::accumulate(nnList,nnList+nSmooth,0.0f,[pkd](auto &a,auto &b) { return a + pkd->particles[b.pPart].mass();}) - p.mass();

    // We could unify this factors in order to avoid more operations
    // (although I think wont make a huge change anyway)
    // There is no need to explicitly convert to solar masses
    // becuse we are doing a ratio
    const double prob = p.star().fSNEfficiency *
                        smf->dSNFBNumberSNperMass *
                        p.mass() / (smf->dSNFBDu*totMass);


    assert(prob<1.0);
    for (i=0; i<nSmooth; ++i) {
        if (all(nnList[i].dr==0)) continue;

        if (rand()<RAND_MAX*prob) { // We have a supernova explosion!
            auto q = pkd->particles[nnList[i].pPart];
            auto &qsph = q.sph();
            //printf("Uint %e extra %e \n",
            //   qsph->Uint, pkd->param.dFeedbackDu * pkdMass(pkd,q));

            const double feed_energy = smf->dSNFBDu * q.mass();
#ifdef OLD_FB_SCHEME
            qsph.Uint += feed_energy;
            qsph.E += feed_energy;
#ifdef ENTROPY_SWITCH
            qsph.S += feed_energy*(smf->dConstGamma-1.) *
                      pow(q.density(), -smf->dConstGamma+1);
#endif
#else // OLD_BH_SCHEME
            qsph.fAccFBEnergy += feed_energy;
#endif

            //printf("Adding SN energy! \n");
        }

    }
}



void initSNFeedback(void *vpkd, void *vp) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(vp)];

    if (p.is_gas()) {
        auto &sph = p.sph();

#ifdef OLD_FB_SCHEME
        sph.Uint = 0.;
        sph.E = 0.;
#ifdef ENTROPY_SWITCH
        sph.S = 0.;
#endif
#else // OLD_FB_SCHEME
        sph.fAccFBEnergy = 0.;
#endif
    }

}


void combSNFeedback(void *vpkd, void *v1, const void *v2) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(v1)];
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(v2)];

    if (p1.is_gas() && p2.is_gas()) {
        auto &sph1 = p1.sph();
        const auto &sph2 = p2.sph();

#ifdef OLD_FB_SCHEME
        sph1.Uint += sph2.Uint;
        sph1.E += sph2.E;
#ifdef ENTROPY_SWITCH
        sph1.S += sph2.S;
#endif
#else //OLD_FB_SCHEME
        sph1.fAccFBEnergy += sph2.fAccFBEnergy;
#endif

    }

}

#endif // FEEDBACK
