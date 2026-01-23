// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#include "services/checkpoint.h"
#include "fmt/core.h"

void ServiceCheckpoint::IO(PST pst, std::span<std::byte> vin, int iGroup, int iSegment, int nSegment) {
    auto in = reinterpret_cast<input*>(vin.data());
    PKD pkd = pst->lcl()->pkd;
    auto achOutFile = fmt::format(in->filename, fmt::arg("i", pkd->mdl->Self()));
    auto achOutDecoratorFile = fmt::format(in->filename, fmt::arg("i", pkd->mdl->Self())) + ".dec";

    // Need to add prefix as with previous "makeName" implementation
    pkd->Checkpoint(achOutFile.c_str(), achOutDecoratorFile.c_str());
}
