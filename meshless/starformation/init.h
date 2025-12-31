#ifndef SERVICE_STARFORMINIT_H
#define SERVICE_STARFORMINIT_H
#include "TraversePST.h"
#include "starformation.h"

struct inStarFormInit {
    double dTime;
    double dCCSNFBDelay;
    double dSNIaFBDelay;
    int bCCSNFeedback;
    int bSNIaFeedback;
    bool bFeedback;
};

class ServiceStarFormInit : public TraverseCount<struct outStarForm> {
  public:
    static constexpr std::string_view name = "StarFormInit";
    typedef struct inStarFormInit input;
    typedef struct outStarForm output;
    explicit ServiceStarFormInit(PST pst, std::string_view service_name = name, int nInBytes = 0) : TraverseCount(pst, service_name, nInBytes) {}

  protected:
    virtual int Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) override {
        auto pkd = pst->lcl()->pkd;
        auto& out = *reinterpret_cast<output*>(vout.data());
        auto& in = *reinterpret_cast<input*>(vin.data());
        pkdStarFormInit(pkd, in, &(out.nFormed));

        return sizeof(output);
    }

    void pkdStarFormInit(PKD pkd, struct inStarFormInit& in, int* nFormed) {
        if (in.bFeedback) {
            *nFormed = 0;
            for (auto& p : pkd->particles) {
                if (p.is_star()) {
                    auto& star = p.star();
                    if (star.fTimer >= 0) { // fTimer < 0 can be used for stars
                        // in the IC that are not supossed to explode
                        if ((in.dTime - star.fTimer) < in.dCCSNFBDelay) star.bCCSNFBDone = in.bCCSNFeedback ? 0 : 1;

                        if ((in.dTime - star.fTimer) < in.dSNIaFBDelay) star.bSNIaFBDone = in.bSNIaFeedback ? 0 : 1;

                        if (!(star.bCCSNFBDone && star.bSNIaFBDone)) *nFormed += 1;
                    }
                }
            }
        }
    }
};
#endif
