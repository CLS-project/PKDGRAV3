// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "pkd_config.h"
#include "gravity.h"

auto ServiceGravity::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->pkd();
    auto in = reinterpret_cast<input*>(vin.data());
    auto outr = reinterpret_cast<output*>(vout.data());

    double dPart, dPartNumAccess, dPartMissRatio, dCell, dCellNumAccess, dCellMissRatio, dFlop;
    pkdGravAll(pkd, &in->kick, &in->lc, &in->ts, in->dTime, in->nReps, in->bPeriodic, in->bGPU, in->bEwald, in->iRoot1, in->iRoot2, in->dEwCut, in->dEwhCut,
               in->dTheta, &in->SPHoptions, &outr->nActive, &dPart, &dPartNumAccess, &dPartMissRatio, &dCell, &dCellNumAccess, &dCellMissRatio, &dFlop,
               outr->nRung);
    outr->dFlopSingleCPU = 1e-9 * pkd->dFlopSingleCPU;
    outr->dFlopDoubleCPU = 1e-9 * pkd->dFlopDoubleCPU;
    outr->dFlopSingleGPU = 1e-9 * pkd->dFlopSingleGPU;
    outr->dFlopDoubleGPU = 1e-9 * pkd->dFlopDoubleGPU;
    outr->nTilesTotal = pkd->nTilesTotal;
    outr->nTilesCPU = pkd->nTilesCPU;
    outr->sLocal = STAT(pst->Self(), pkd->Local());
    outr->sActive = STAT(pst->Self(), outr->nActive);
#ifdef __linux__
    FILE* fp;
    char buffer[512], *save, *f, *v;
    double dRSS = 0.0;
    fp = fopen("/proc/self/stat", "r");
    if (fp != NULL) {
        if (fgets(buffer, sizeof(buffer), fp) == NULL) buffer[0] = '\0';
        fclose(fp);
        f = strtok_r(buffer, " ", &save);
        for (auto i = 0; i <= 36 && f; i++) {
            switch (i) {
            case 23:
    #ifdef HAVE_GETPAGESIZE
                dRSS = atol(f) * 1.0 * getpagesize() / 1024 / 1024 / 1024;
    #else
                dRSS = atol(f) * 1.0 / 1024 / 1024;
    #endif
                break;
            default:
                break;
            }
            f = strtok_r(NULL, " ", &save);
        }
    }
    outr->sRSS = STAT(pst->Self(), dRSS);
    fp = fopen("/proc/meminfo", "r");
    double dFreeMemory = 0;
    if (fp != NULL) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            f = strtok_r(buffer, ":", &save);
            v = strtok_r(NULL, " ", &save);
            if (strcmp(f, "MemFree") == 0 || strcmp(f, "Buffers") == 0 || strcmp(f, "Cached") == 0) {
                dFreeMemory -= atol(v) * 1.0 / 1024 / 1024;
            }
        }
        fclose(fp);
    }
    outr->sFreeMemory = STAT(pst->Self(), dFreeMemory);
#endif
    if (outr->nActive) { // Some stats are valid only if there are active particles
        outr->sPart = STAT(pst->Self(), dPart);
        outr->sPartNumAccess = STAT(pst->Self(), dPartNumAccess);
        outr->sPartMissRatio = STAT(pst->Self(), dPartMissRatio);
        outr->sCell = STAT(pst->Self(), dCell);
        outr->sCellNumAccess = STAT(pst->Self(), dCellNumAccess);
        outr->sCellMissRatio = STAT(pst->Self(), dCellMissRatio);
    }
    outr->sFlop = STAT(pst->Self(), dFlop);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    outr->sComputing = STAT(pst->Self(), pst->mdl()->TimeComputing());
    outr->sWaiting = STAT(pst->Self(), pst->mdl()->TimeWaiting());
    outr->sSynchronizing = STAT(pst->Self(), pst->mdl()->TimeSynchronizing());
#endif
    return sizeof(output);
}