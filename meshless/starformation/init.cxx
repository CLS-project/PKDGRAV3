#include "init.h"

auto pkdContext::StarFormInit(double dTime, double dCCSNFBDelay, double dSNIaFBDelay, bool bCCSNFeedback, bool bSNIaFeedback, bool bFeedback) -> int {
    int nFormed = 0;
    if (bFeedback) {
        for (auto& p : particles) {
            if (p.is_star()) {
                auto& star = p.star();
                if (star.fTimer >= 0) { // fTimer < 0 can be used for stars
                    // in the IC that are not supossed to explode
                    if ((dTime - star.fTimer) < dCCSNFBDelay) star.bCCSNFBDone = bCCSNFeedback ? 0 : 1;
                    if ((dTime - star.fTimer) < dSNIaFBDelay) star.bSNIaFBDone = bSNIaFeedback ? 0 : 1;
                    if (!(star.bCCSNFBDone && star.bSNIaFBDone)) ++nFormed;
                }
            }
        }
    }
    return nFormed;
}
