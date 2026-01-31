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
#pragma once

#include "pkd.h"
#include "mdl.h"
#include "core/dataregistry.h"

struct LCL : private DataRegistryOwner {
    using DataRegistryOwner::data;
    using DataRegistryOwner::data_emplace;
    using DataRegistryOwner::data_ensure;
    using DataRegistryOwner::data_erase;
    using DataRegistryOwner::data_present;
    PKD pkd;
};

class pstNode : private DataRegistryOwner {
  public:
    pstNode() = delete;
    explicit pstNode(mdl::mdlClass* mdl, LCL* lcl) : pstLower_(nullptr), mdl_(mdl), lcl_(lcl), idUpper_(-1), nLeaves_(1), nLower_(0), nUpper_(0), iLvl_(0) {}
    using DataRegistryOwner::data;
    using DataRegistryOwner::data_emplace;
    using DataRegistryOwner::data_ensure;
    using DataRegistryOwner::data_erase;
    using DataRegistryOwner::data_present;

  private:
    pstNode* pstLower_;  // Left child node: found locally
    mdl::mdlClass* mdl_; // MDL class
    LCL* lcl_;           // Variables that are identical on all PST nodes
    int idUpper_;        // ID of the thread with the right child node
    int nLeaves_;        // Number of leaves in this subtree (nLower_ + nUpper_, or 1 for a leaf)
    int nLower_;         // Number of leaves in the lower subtree (0 if leaf)
    int nUpper_;         // Number of leaves in the upper subtree (0 if leaf)
    int iLvl_;           // Depth in the PST tree

  public:
    [[nodiscard]] auto OffNode() const { return nLeaves_ > mdl_->Cores(); }
    [[nodiscard]] auto OnNode() const { return nLeaves_ <= mdl_->Cores(); }
    [[nodiscard]] auto AmNode() const { return nLeaves_ == mdl_->Cores(); }
    [[nodiscard]] auto NotNode() const { return nLeaves_ != mdl_->Cores(); }
    [[nodiscard]] auto AmCore() const { return nLeaves_ == 1; }
    [[nodiscard]] auto NotCore() const { return nLeaves_ > 1; }

    [[nodiscard]] auto Self() const { return mdl_->Self(); }
    [[nodiscard]] auto n_lower() const { return nLower_; }
    [[nodiscard]] auto n_upper() const { return nUpper_; }
    [[nodiscard]] auto n_total() const { return nLeaves_; }
    [[nodiscard]] auto lower() const { return pstLower_; }
    [[nodiscard]] auto upper() const { return idUpper_; }
    [[nodiscard]] auto depth() const { return iLvl_; }

    [[nodiscard]] auto mdl() const { return mdl_; }
    [[nodiscard]] auto lcl() const { return lcl_; }

  private:
    friend class ServiceSetAdd;
    auto create_lower(int idMiddle, int idUp) {
        assert(pstLower_ == nullptr);
        idUpper_ = idMiddle;
        nLower_ = idMiddle - Self();
        nUpper_ = idUp - idMiddle;
        nLeaves_ = nLower_ + nUpper_;
        pstLower_ = new pstNode(mdl_, lcl_);
        pstLower_->iLvl_ = iLvl_ + 1;
        return pstLower_;
    }
};
using PST = pstNode*;
