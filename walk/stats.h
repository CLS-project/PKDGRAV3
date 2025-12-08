#ifndef WALK_STATS_H
#define WALK_STATS_H

struct WalkStats {
    double dFlopSingleCPU = 0.0;
    double dFlopSingleGPU = 0.0;
    double dFlopDoubleCPU = 0.0;
    double dFlopDoubleGPU = 0.0;
    double dFlop = 0.0;
    uint64_t nTilesTotal = 0;
    uint64_t nTilesCPU = 0;
    int nActive = 0;
};

#endif
