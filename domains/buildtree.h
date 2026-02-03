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
#ifndef BUILDTREE_H
#define BUILDTREE_H
#include "TraversePST.h"
#include <string_view>

class ServiceBuildTree : public TraversePST {
  public:
    static constexpr std::string_view name = "BuildTree";
    struct input {
        int nBucket;     /* Bucket Size */
        int nGroup;      /* Group Size */
        uint32_t uRoot;  /* Which root node to use */
        uint32_t utRoot; /* Template tree */
        double ddHonHLimit;
    };
    explicit ServiceBuildTree(PST pst, std::string_view service_name = name, size_t input_size = sizeof(input))
        : TraversePST(pst, service_name, input_size, (pst->mdl()->Threads() == 1 ? 1 : 2 * pst->mdl()->Threads() - 1) * pkdContext::MaxNodeSize()) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    virtual auto build_tree(PKD pkd, std::span<std::byte> in) -> void;
};

class ServiceTreeUpdateFlagBounds : public ServiceBuildTree {
  public:
    static constexpr std::string_view name = "TreeUpdateFlagBounds";
    struct input {
        ServiceBuildTree::input base;
        SPHOptions SPHoptions;
    };
    explicit ServiceTreeUpdateFlagBounds(PST pst) : ServiceBuildTree(pst, name, sizeof(input)) {}

  protected:
    // auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto build_tree(PKD pkd, std::span<std::byte> in) -> void override;
};

class ServiceTreeInitMarked : public TraversePST {
  public:
    static constexpr std::string_view name = "TreeInitMarked";
    using input = void;
    using output = void;
    explicit ServiceTreeInitMarked(PST pst) : TraversePST(pst, name) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

#endif // BUILDTREE_H