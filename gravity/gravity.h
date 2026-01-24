// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "TraversePST.h"
#include <cstdint>

struct StatsCollector {
    double dSum = 0.0;
    double dSum2 = 0.0;
    double dMax = 0.0;
    int idMax = -1;
    int n = 0;

    StatsCollector() = default;
    StatsCollector(int id, double dSum) : dSum(dSum), dSum2(dSum * dSum), dMax(dSum), idMax(id), n(1) {}

    auto operator+=(const StatsCollector& other) -> StatsCollector& {
        dSum += other.dSum;
        dSum2 += other.dSum2;
        if (other.dMax > dMax) {
            dMax = other.dMax;
            idMax = other.idMax;
        }
        n += other.n;
        return *this;
    }
};
using STAT = StatsCollector;
static_assert(std::is_trivially_copyable_v<StatsCollector>);

class ServiceGravity : public TraverseCount<ServiceGravity> {
  public:
    static constexpr std::string_view name = "Gravity";
    struct input {
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
        struct pkdKickParameters kick;
        struct pkdLightconeParameters lc;
        struct pkdTimestepParameters ts;
        SPHOptions SPHoptions;
    };
    struct output {
        STAT sLocal;
        STAT sActive;
        STAT sPart;
        STAT sPartNumAccess;
        STAT sPartMissRatio;
        STAT sCell;
        STAT sCellNumAccess;
        STAT sCellMissRatio;
        STAT sFlop;
#ifdef INSTRUMENT
        STAT sComputing;
        STAT sWaiting;
        STAT sSynchronizing;
#endif
#ifdef __linux__
        STAT sRSS;
        STAT sFreeMemory;
#endif
        double dFlopSingleCPU;
        double dFlopDoubleCPU;
        double dFlopSingleGPU;
        double dFlopDoubleGPU;
        uint64_t nActive;
        blitz::TinyVector<uint64_t, IRUNGMAX + 1> nRung;
        uint64_t nTilesTotal;
        uint64_t nTilesCPU;

        auto operator+=(const output& other) -> output& {
            sLocal += other.sLocal;
            sActive += other.sActive;
            sPart += other.sPart;
            sPartNumAccess += other.sPartNumAccess;
            sPartMissRatio += other.sPartMissRatio;
            sCell += other.sCell;
            sCellNumAccess += other.sCellNumAccess;
            sCellMissRatio += other.sCellMissRatio;
            sFlop += other.sFlop;
#ifdef INSTRUMENT
            sComputing += other.sComputing;
            sWaiting += other.sWaiting;
            sSynchronizing += other.sSynchronizing;
#endif
#ifdef __linux__
            sRSS += other.sRSS;
            sFreeMemory += other.sFreeMemory;
#endif
            dFlopSingleCPU += other.dFlopSingleCPU;
            dFlopDoubleCPU += other.dFlopDoubleCPU;
            dFlopSingleGPU += other.dFlopSingleGPU;
            dFlopDoubleGPU += other.dFlopDoubleGPU;
            nActive += other.nActive;
            nRung += other.nRung;
            nTilesTotal += other.nTilesTotal;
            nTilesCPU += other.nTilesCPU;
            return *this;
        }
    };
    explicit ServiceGravity(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int final;
};
