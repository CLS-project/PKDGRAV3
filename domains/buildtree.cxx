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

#include "pkd_config.h"
#include "buildtree.h"

auto ServiceBuildTree::Recurse(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    LCL* plcl = pst->lcl();
    PKD pkd = plcl->pkd;
    auto pTop = pkd->tree[reinterpret_cast<KDN*>(vout.data())];

    auto out = vout.subspan(pkd->tree.ElementSize());                  // The two subtrees
    auto loBytes = (pst->n_lower() * 2 - 1) * pkd->tree.ElementSize(); // left subtree size
    auto hiBytes = (pst->n_upper() * 2 - 1) * pkd->tree.ElementSize(); // right subtree size
    auto lo = out.subspan(0, loBytes);                                 // left subtree
    auto hi = out.subspan(loBytes, hiBytes);                           // right subtree

    /* We will accumulate the top tree here */
    int rID = ReqService(pst, vin);
    auto nLower = Traverse(pst->lower(), vin, lo);
    assert(nLower == loBytes);
    auto nUpper = GetReply(rID, hi);
    assert(nUpper == hiBytes);

    /*
    ** Combine Cell1 and Cell2 into pCell
    ** to find cell CoM, bounds and multipoles.
    ** This also computes the opening radius for gravity.
    */
    auto pCell1 = pTop + 1;
    auto pCell2 = pTop + pst->n_lower() * 2;
    pTop->bMax() = HUGE_VAL; /* initialize bMax for CombineCells */
    auto minside = pst->bnd.minside();
    pkd->CombineCells1(pTop, pCell1, pCell2);
    pTop->calc_open(minside);
    pkd->CombineCells2(pTop, pCell1, pCell2);

    /* Get our cell ready */
    pTop->set_top_tree(pst->n_lower() * 2);

    return nLower + nUpper + pkd->tree.ElementSize();
}

auto ServiceBuildTree::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    LCL* plcl = pst->lcl();
    PKD pkd = plcl->pkd;
    auto in = reinterpret_cast<const input*>(vin.data());
    auto pTop = pkd->tree[reinterpret_cast<KDN*>(vout.data())];

    auto pRoot = pkd->tree[in->uRoot];
    pkd->TreeAlignNode();
    build_tree(pkd, vin);
    *pTop = *pRoot;
    /* Get our cell ready */
    pTop->set_remote(pst->Self());
    return pkd->tree.ElementSize();
}

auto ServiceBuildTree::build_tree(PKD pkd, std::span<std::byte> in) -> void {
    auto input_data = reinterpret_cast<const input*>(in.data());
    pkd->TreeBuild(input_data->nBucket, input_data->nGroup, input_data->uRoot, input_data->utRoot, input_data->ddHonHLimit);
}

auto ServiceTreeUpdateFlagBounds::build_tree(PKD pkd, std::span<std::byte> in) -> void {
    auto input_data = reinterpret_cast<const input*>(in.data());
    pkd->TreeUpdateFlagBounds(input_data->base.uRoot, &input_data->SPHoptions);
}

auto ServiceTreeInitMarked::Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int {
    auto pkd = pst->lcl()->pkd;
    pkd->TreeInitMarked();
    return 0;
}