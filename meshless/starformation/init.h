#ifndef SERVICE_STARFORMINIT_H
#define SERVICE_STARFORMINIT_H
#include "TraversePST.h"

class ServiceStarFormInit : public TraverseCount<ServiceStarFormInit> {
  public:
    static constexpr std::string_view name = "StarFormInit";
    struct input {
        double dTime;
        double dCCSNFBDelay;
        double dSNIaFBDelay;
        bool bCCSNFeedback;
        bool bSNIaFeedback;
        bool bFeedback;
    };
    using output = uint64_t;
    explicit ServiceStarFormInit(PST pst, std::string_view service_name = name, int nInBytes = 0) : TraverseCount(pst, service_name, nInBytes) {}

  protected:
    auto Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int override {
        auto pkd = pst->pkd();
        auto& out = *reinterpret_cast<output*>(vout.data());
        auto& in = *reinterpret_cast<const input*>(vin.data());
        out = pkd->StarFormInit(in.dTime, in.dCCSNFBDelay, in.dSNIaFBDelay, in.bCCSNFeedback, in.bSNIaFeedback, in.bFeedback);
        return sizeof(output);
    }
};
#endif
