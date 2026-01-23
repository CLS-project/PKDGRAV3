// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "pkd_config.h"
#include "memstatus.h"

auto ServiceMemStatus::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->lcl()->pkd;
    auto out = reinterpret_cast<output*>(vout.data());
#ifdef __linux__
    FILE* fp;
    char buffer[512], *save, *f, *v;
    int i;

    fp = fopen("/proc/self/stat", "r");
    if (fp != NULL) {
        if (fgets(buffer, sizeof(buffer), fp) == NULL) buffer[0] = '\0';
        fclose(fp);
        f = strtok_r(buffer, " ", &save);
        for (i = 0; i <= 36 && f; i++) {
            switch (i) {
            case 23:
    #ifdef HAVE_GETPAGESIZE
                out->rss = atol(f) * getpagesize() / 1024 / 1024;
    #else
                out->rss = atol(f) / 1024;
    #endif
                break;
            default:
                break;
            }
            f = strtok_r(NULL, " ", &save);
        }
    }
    out->freeMemory = 0;
    fp = fopen("/proc/meminfo", "r");
    if (fp != NULL) {
        while (fgets(buffer, sizeof(buffer), fp)) {
            f = strtok_r(buffer, ":", &save);
            v = strtok_r(NULL, " ", &save);
            if (strcmp(f, "MemFree") == 0 || strcmp(f, "Buffers") == 0 || strcmp(f, "Cached") == 0) {
                out->freeMemory += atol(v) / 1024;
            }
        }
        fclose(fp);
    }
#endif
    out->nBytesTree = pkd->tree.memory();
    out->nBytesCl = pkdClMemory(pkd);
    out->nBytesIlp = pkdIlpMemory(pkd);
    out->nBytesIlc = pkdIlcMemory(pkd);
    return sizeof(output);
}