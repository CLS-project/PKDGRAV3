/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "core/partclass.h"
#include "pkd_config.h"
#include "classes.h"
#include <algorithm>

void pkdContext::swap_classes(int id) {
    auto classes = particles.getClasses();
    auto rID = mdl->ReqService<ServiceSwapClasses>(id, classes);
    std::vector<PARTCLASS> in_classes(PKD_MAX_CLASSES);
    size_t n = mdl->GetReply(rID, in_classes);
    in_classes.resize(n / sizeof(PARTCLASS));
    particles.setClasses(in_classes, false);
}

// This is a special case: we simply desecend the PST to the leaf and swap the classes there.
auto ServiceSwapClasses::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    return Traverse(pst->lower(), vin, vout);
}

auto ServiceSwapClasses::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = span_cast<const PARTCLASS>(vin);
    auto out = span_cast<PARTCLASS>(vout);
    auto pkd = pst->pkd();
    auto classes = pkd->particles.getClasses();
    std::ranges::copy(classes, out.begin());
    pkd->particles.setClasses(in, false);
    return classes.size_bytes();
}

auto ServiceGetClasses::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto out = span_cast<PARTCLASS>(vout);
    auto pkd = pst->pkd();
    // Get the list of particle classes, copy to output and sort
    auto classes = pkd->particles.getClasses();
    assert(out.size() >= classes.size());
    out = out.subspan(0, classes.size());
    std::ranges::copy(classes, out.begin());
    std::ranges::sort(out);
    return out.size_bytes();
}

auto ServiceGetClasses::Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST) -> int {
    auto out1 = span_cast<PARTCLASS>(vout1);
    auto out2 = span_cast<const PARTCLASS>(vout2);
    // Combine the two sorted lists into one sorted list without duplicates
    std::vector<PARTCLASS> combined;
    combined.reserve(PKD_MAX_CLASSES);
    std::ranges::set_union(out1, out2, std::back_inserter(combined));
    assert(combined.size() <= PKD_MAX_CLASSES);
    out1 = out1.subspan(0, combined.size());
    std::ranges::copy(combined, out1.begin());
    return out1.size_bytes();
}

auto ServiceSetClasses::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = span_cast<const PARTCLASS>(vin);
    auto pkd = pst->pkd();
    pkd->particles.setClasses(in, true); // update existing particles
    return 0;
}
