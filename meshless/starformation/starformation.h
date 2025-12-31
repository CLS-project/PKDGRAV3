#ifndef STARFORMATION_H
#define STARFORMATION_H

#include "TraversePST.h"
#include "feedback.h"
#include "eEOS/eEOS.h"

/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * ---------------------
 */
struct outStarForm {
    int nFormed = 0;
    int nDeleted = 0;
    double dMassFormed = 0.0;

    outStarForm& operator+=(outStarForm b) {
        nFormed += b.nFormed;
        nDeleted += b.nDeleted;
        dMassFormed += b.dMassFormed;
        return *this;
    }
};

struct inStarForm {
    int bSFThresholdDenSchaye2004;

    double dHubble;
    double dDelta;
    double dTime;
    double dScaleFactor;

    double dSFindexKS;
    double dSFnormalizationKS;
    double dConstGamma;
    double dSFGasFraction;
    double dSFThresholdDen;
    double dSFThresholdOD;
    double dSFThresholdu;
    double dSFEfficiency;

    double dEOSPolyFloorIndex;
    double dEOSPolyFloorDen;
    double dEOSPolyFlooru;

    int bCCSNFeedback;
    int bSNIaFeedback;
    double dSNFBEfficiency;
    double dSNFBMaxEff;
    double dSNFBEffIndex;
    double dSNFBEffnH0;

#ifdef STELLAR_EVOLUTION
    double dSNIaMaxMass;
    double dCCSNMinMass;
    double dCCSNMaxMass;
#endif
    struct eEOSparam eEOS;
};

class ServiceStarForm : public TraverseCount<struct outStarForm> {
  public:
    static constexpr std::string_view name = "StarForm";
    typedef struct inStarForm input;
    typedef struct outStarForm output;
    explicit ServiceStarForm(PST pst, std::string_view service_name = name, int nInBytes = 0) : TraverseCount(pst, service_name, nInBytes) {}

  protected:
    virtual int Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) override {
        auto pkd = pst->lcl()->pkd;
        auto& out = *reinterpret_cast<output*>(vout.data());
        auto& in = *reinterpret_cast<input*>(vin.data());
        pkdStarForm(pkd, in, &(out.nFormed), &(out.dMassFormed), &(out.nFormed));

        return sizeof(output);
    }

  private:
    void pkdStarForm(PKD pkd, struct inStarForm in, int* nFormed, /* number of stars formed */
                     double* dMassFormed,                         /* mass of stars formed */
                     int* nDeleted) {                             /* gas particles deleted */

        assert(pkd->decorators.present(DECORATOR_FIELD::oStar));
        assert(pkd->decorators.present(DECORATOR_FIELD::oMFM));
        assert(pkd->particles.present(PKD_FIELD::oMass));

        *nFormed = 0;
        *nDeleted = 0;
        *dMassFormed = 0.0;

        const double a_m1 = 1.0 / in.dScaleFactor;
        const double a_m2 = a_m1 * a_m1;
        const double a_m3 = a_m2 * a_m1;

        for (auto& p : pkd->particles) {

            if (p.is_active() && p.is_gas()) {
                const auto& mass = p.mass();
                auto& sph = p.sph();

#ifdef HAVE_METALLICITY
                double dThreshDen;
                if (in.bSFThresholdDenSchaye2004) {
                    const double dMetAbun = (double)sph.fMetalMass / (double)mass;
                    const double dThreshFac = dMetAbun < 1.5e-6 ? 100.0 : pow(2e-3 / dMetAbun, 0.64);
                    dThreshDen = in.dSFThresholdDen * dThreshFac;
                }
                else {
                    dThreshDen = in.dSFThresholdDen;
                }
#else
                double dThreshDen = in.dSFThresholdDen;
#endif

                dThreshDen = dThreshDen > in.dSFThresholdOD ? dThreshDen : in.dSFThresholdOD;

                double dmstar;
                if (in.dSFEfficiency > 0.0) {
                    dmstar = density_SFR(mass, p.density(), in.dHubble, a_m1, a_m3, dThreshDen, in.dSFThresholdu, in.dSFEfficiency, sph);
                }
                else {
                    dmstar = pressure_SFR(mass, p.density(), p.ball(), a_m3, dThreshDen, in.dSFnormalizationKS, in.dSFindexKS, in.dSFGasFraction,
                                          in.dConstGamma, in.eEOS, sph);
                }

                sph.SFR = dmstar;
                if (dmstar == 0.0) continue;

                const double dt = in.dTime - sph.lastUpdateTime;
                const double prob = 1.0 - exp(-dmstar * dt / mass);

                // Star formation event?
                if (rand() < RAND_MAX * prob) {
#ifdef STELLAR_EVOLUTION
                    blitz::TinyVector<float, ELEMENT_COUNT> ElemMass;
                    float fMetalMass;
                    ElemMass = sph.ElemMass;
                    fMetalMass = sph.fMetalMass;
#endif

                    // We just change the class of the particle to stellar one
                    p.set_class(mass, p.soft0(), 0, 1.0, FIO_SPECIES_STAR);
                    auto& star = p.star();

                    // When changing the class, we have to take into account that
                    // the code velocity has different scale factor dependencies for
                    // dm/star particles and gas particles
                    p.velocity() *= in.dScaleFactor;

                    // We log statistics about the formation time
                    star.fTimer = in.dTime;

#ifdef STELLAR_EVOLUTION
                    const float fMassInv = 1.0f / mass;
                    star.ElemAbun = ElemMass * fMassInv;
                    star.fMetalAbun = fMetalMass * fMassInv;

                    star.fInitialMass = mass;
                    star.fLastEnrichTime = 0.0f;

                    stevStarParticleInit(pkd, star, in.dSNIaMaxMass, in.dCCSNMinMass, in.dCCSNMaxMass);
#endif

                    star.bCCSNFBDone = in.bCCSNFeedback ? 0 : 1;
                    star.bSNIaFBDone = in.bSNIaFeedback ? 0 : 1;

                    // Compute the feedback efficiency for this particle based
                    // on the birth information (eq. 7 Schaye 2015)
#ifdef HAVE_METALLICITY
                    const double fMetalAbun = star.fMetalAbun;
#else
                    const double fMetalAbun = 0.00127; // 0.1 Zsolar
#endif
                    star.fSNEfficiency =
                        SNFeedbackEfficiency(fMetalAbun, p.density() * a_m3, in.dSNFBEfficiency, in.dSNFBMaxEff, in.dSNFBEffIndex, in.dSNFBEffnH0);

                    // Safety check
                    assert(p.is_star());
                    assert(!p.is_gas());

                    (*nFormed)++;
                    *dMassFormed += mass;
                }
            }
        }
    }

    inline double pressure_SFR(const float fMass, const float fDens, const double fBall, const double a_m3, const double dThreshDen,
                               const double dSFnormalizationKS, const double dSFindexKS, const double dSFGasFraction, const double dConstGamma, eEOSparam& eEOS,
                               decorators::MFM& sph) {
        // Two SF thresholds are applied:
        //      a) Minimum density, computed at the master level
        //      b) Maximum temperature of a factor 0.5 dex (i.e., 3.1622)
        //         above the polytropic eEOS
        double maxUint;
        if (eEOS.bDoEOS) {
            maxUint = 3.16228 * fMass * eEOSEnergyFloor<vec<double, double>, mmask<bool>>(a_m3, fDens, fBall, dConstGamma, eEOS);
        }
        else {
            // This configuration is allowed, but can easily produce numerical
            // fragmentation!!!
            maxUint = INFINITY;
        }

        const double dUint = sph.Uint + sph.fAccFBEnergy;

        if (std::max(dUint, sph.lastUint) > maxUint || fDens < dThreshDen) {
            return 0.0;
        }

        const double SFexp = 0.5 * (dSFindexKS - 1.);
        const double dmstar = dSFnormalizationKS * fMass * pow(dConstGamma * dSFGasFraction * sph.P * a_m3, SFexp);
        return dmstar;
    }

    inline double density_SFR(const float fMass, const float fDens, const double dHubble, const double a_m1, const double a_m3, const double dThreshDen,
                              const double dSFThresholdu, const double dSFEfficiency, decorators::MFM& sph) {
        // Three SF criteria are enforced:
        //      a) A minimum density, computed at the master level
        //      b) A maximum temperature, set by dSFThresholdTemp
        //      c) Converging flow (negative physical velocity divergence)
        const double dVelDiv = a_m1 * (sph.gradVx[0] + sph.gradVy[1] + sph.gradVz[2]) + 3. * dHubble;
        if (fDens < dThreshDen || sph.Uint > dSFThresholdu || dVelDiv > 0.0) {
            return 0.0;
        }

        const double tff = sqrt(3. * M_PI / (32. * fDens * a_m3));
        const double dmstar = dSFEfficiency * fMass / tff;
        return dmstar;
    }
};
#endif
