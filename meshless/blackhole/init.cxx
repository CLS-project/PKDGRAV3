#include "blackhole/init.h"

void pkdContext::BHInit(uint8_t uRungMax) {
    for (auto& p : particles) {
        if (p.is_bh()) {
            p.set_rung(uRungMax);
            p.set_new_rung(uRungMax);
        }
    }
}
