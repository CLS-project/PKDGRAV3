#include "blackhole/init.h"

void pkdBHInit(PKD pkd, uint8_t uRungMax) {
    for (auto& p : pkd->particles) {
        if (p.is_bh()) {
            p.set_rung(uRungMax);
            p.set_new_rung(uRungMax);
        }
    }
}
