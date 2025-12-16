#ifndef GRAVITY_TRAITS_H
#define GRAVITY_TRAITS_H
#include "basetype.h"
#ifdef __CUDACC__
    #define CUDA_BOTH __host__ __device__
#else
    #define CUDA_BOTH
#endif

class GravityInteractionTraits {
public:
    struct name {
        static constexpr std::string_view value = "GravityInteraction";
    };

    struct Options {
        struct pkdKickParameters kick;
        struct pkdLightconeParameters lc;
        struct pkdTimestepParameters ts;
        double dTime;
        double dEwCut;
        double dEwhCut;
        double dTheta;
        int nReps;
        int bPeriodic;
        int bEwald;
        int bGPU;
        int iRoot1;
        int iRoot2;
        // SPHOptions SPHoptions;
        bool doExtensiveTest;
        uint8_t uRungLo;
        uint8_t uRungHi;
    };

    typedef struct PINFOIN {
        blitz::TinyVector<float, 3> r;
        blitz::TinyVector<float, 3> a;
        float fSmooth2;
        float fDensity;
        void init(particleStore::ParticlePointer p) {};
    } PINFOIN;

    typedef struct PINFOOUT {
        blitz::TinyVector<float, 3> a;
        float fPot;
        float dirsum, normsum;
        float rhopmax;
        float dtEst;
        float maxRung;

        void init(particleStore::ParticlePointer p, float dirLsum,
                  float normLsum) {
            a = 0.0f;
            fPot = 0.0f;
            dirsum = dirLsum;
            normsum = normLsum;
            rhopmax = 0.0f;
            dtEst = 1e14f;
            maxRung = 0.0f;
        }
    } PINFOOUT;

    using IL = ilpList;
    using ILC = ilcList;

    struct WORKPARTICLE : public workParticleCommon {
        // Use a constructor for convenience
        WORKPARTICLE(void *pkd, int knP, IL *il) {
            nRefs = 1; // I am using it currently
            ctx = pkd;
            dFlop = 0.0;
            dFlopSingleCPU = dFlopSingleGPU = 0.0;
            dFlopDoubleCPU = dFlopDoubleGPU = 0.0;
            nTilesTotal = nTilesCPU = 0;
            nP = 0;

            pPart = new PARTICLE*[knP];
            iPart = new uint32_t[knP];
            pInfoIn = new PINFOIN[knP];
            pInfoOut = new PINFOOUT[knP];
            c[0] = il->getReference(0);
            // assert(wp->c[0] == ilc.getReference(0));
            c[1] = il->getReference(1);
            // assert(wp->c[1] == ilc.getReference(1));
            c[2] = il->getReference(2);
            // assert(wp->c[2] == ilc.getReference(2));
        };
        ~WORKPARTICLE() {
            delete[] pPart;
            delete[] iPart;
            delete[] pInfoIn;
            delete[] pInfoOut;
        }
        int nTilesTotal, nTilesCPU;
        PINFOIN *pInfoIn;
        PINFOOUT *pInfoOut;
    };


    template <class F> struct ResultInteraction {};

    template <class F = float> struct ResultPC {
        F ax, ay, az, pot;
        F ir, norm;
        CUDA_BOTH void zero() { ax = ay = az = pot = ir = norm = 0; }
        CUDA_BOTH ResultPC<F> operator+=(const ResultPC<F> rhs) {
            ax += rhs.ax;
            ay += rhs.ay;
            az += rhs.az;
            pot += rhs.pot;
            ir += rhs.ir;
            norm += rhs.norm;
            return *this;
        }
    };

    template <class F = float> struct ResultPP {
        typedef F type;
        F ax, ay, az, pot;
        F ir, norm;
        CUDA_BOTH void zero() { ax = ay = az = pot = ir = norm = 0; }
        CUDA_BOTH ResultPP<F> operator+=(const ResultPP<F> rhs) {
            ax += rhs.ax;
            ay += rhs.ay;
            az += rhs.az;
            pot += rhs.pot;
            ir += rhs.ir;
            norm += rhs.norm;
            return *this;
        }
    };
};
#endif
