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
#ifndef TRAVERSEPST_H
#define TRAVERSEPST_H
#include "pst.h"
#include <cstddef>
#include <string>
#include <type_traits>

class TraversePST; // forward declaration is sufficient

template <typename T>
concept TraverseService = std::derived_from<T, TraversePST>;

// This class will traverse the PST calling Recurse() at each level.
// The Recurse method is responsible for making a request,
// calling and calling Traverse() recursively. When a leaf is reach
// the Service() function is called.
class TraversePST : public mdl::BasicService {
    PST node_pst;

  protected:
    template <typename T>
    std::span<T> span_cast(std::span<std::byte> bytes, size_t nElements) {
        assert(bytes.size() >= sizeof(T) * nElements);
        assert(reinterpret_cast<std::uintptr_t>(bytes.data()) % alignof(T) == 0);
        return {reinterpret_cast<T*>(bytes.data()), nElements};
    }
    template <typename T>
    std::span<T> span_cast(std::span<std::byte> bytes) {
        return span_cast<T>(bytes, bytes.size() / sizeof(T));
    }

    template <typename T>
    std::span<const T> span_cast(std::span<const std::byte> bytes, size_t nElements) {
        assert(bytes.size() >= sizeof(T) * nElements);
        assert(reinterpret_cast<std::uintptr_t>(bytes.data()) % alignof(T) == 0);
        return {reinterpret_cast<const T*>(bytes.data()), nElements};
    }
    template <typename T>
    std::span<const T> span_cast(std::span<const std::byte> bytes) {
        return span_cast<T>(bytes, bytes.size() / sizeof(T));
    }

  public:
    explicit TraversePST(PST node_pst, std::string_view service_name, int nInBytes = 0, int nOutBytes = 0)
        : BasicService(service_name, nInBytes, nOutBytes), node_pst(node_pst) {}
    virtual ~TraversePST() = default;

    template <typename I = std::byte>
    int ReqService(PST pst, std::span<I> in = {}) {
        return pst->mdl()->ReqService(pst->upper(), getServiceID(), mdl::as_bytes(in));
    }
    template <typename I>
    int ReqService(PST pst, const I& in) {
        return pst->mdl()->ReqService(pst->upper(), getServiceID(), mdl::as_bytes(in));
    }

    template <typename O = std::byte>
    int GetReply(int rID, std::span<O> out = {}) {
        return node_pst->mdl()->GetReply(rID, mdl::as_writable_bytes(out));
    }

    template <typename O>
    int GetReply(int rID, O& out) {
        return node_pst->mdl()->GetReply(rID, mdl::as_writable_bytes(out));
    }

  protected:
    virtual int operator()(std::span<std::byte> In, std::span<std::byte> Out) final;
    virtual int Traverse(PST pst, std::span<std::byte> in, std::span<std::byte> out);

    virtual int OffNode(PST pst, std::span<std::byte> in, std::span<std::byte> out) { return Recurse(pst, in, out); }
    virtual int AtNode(PST pst, std::span<std::byte> in, std::span<std::byte> out) { return Recurse(pst, in, out); }
    virtual int Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out);
    virtual int Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) = 0;

  protected:
    template <typename I = std::byte, typename O = std::byte>
        requires(!TraverseService<I>)
    int Traverse(PST pst, std::span<I> in = {}, std::span<O> out = {}) {
        return Traverse(pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }
    template <typename I, typename O>
        requires(!TraverseService<I>)
    int Traverse(PST pst, I& in, O& out) {
        return Traverse(pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }
    template <typename I, typename O = std::byte>
        requires(!TraverseService<I>)
    int Traverse(PST pst, I& in, std::span<O> out = {}) {
        return Traverse(pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }
    template <typename I = std::byte, typename O>
        requires(!TraverseService<I>)
    int Traverse(PST pst, std::span<I> in, O& out) {
        return Traverse(pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }

  protected:
    static int Traverse(unsigned sid, PST pst, std::span<std::byte> in, std::span<std::byte> out);

    static int Traverse(std::string_view service_name, PST pst, std::span<std::byte> in, std::span<std::byte> out);
    template <TraverseService T, typename I = std::byte, typename O = std::byte>
    static int Traverse(PST pst, std::span<I> in = {}, std::span<O> out = {}) {
        return Traverse(T::name, pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }
    template <TraverseService T, typename I, typename O>
    static int Traverse(PST pst, I& in, O& out) {
        return Traverse(T::name, pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }
    template <TraverseService T, typename I, typename O = std::byte>
    static int Traverse(PST pst, I& in, std::span<O> out = {}) {
        return Traverse(T::name, pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }
    template <TraverseService T, typename I = std::byte, typename O>
    static int Traverse(PST pst, std::span<I> in, O& out) {
        return Traverse(T::name, pst, mdl::as_writable_bytes(in), mdl::as_writable_bytes(out));
    }
};

// This class is for services where the input does not change as the PST is walked,
// but where a customized Combine for a fixed size output is needed.
class TraverseCombinePST : public TraversePST {
  public:
    explicit TraverseCombinePST(PST node_pst, std::string_view service_name, int nInBytes, int nOutBytes)
        : TraversePST(node_pst, service_name, nInBytes, nOutBytes) {}
    virtual ~TraverseCombinePST() = default;

  protected:
    virtual int Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) final;
    virtual int Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, [[maybe_unused]] std::span<const std::byte> in,
                        [[maybe_unused]] PST pst) = 0;
};

// This class is for services where the input does not change as the PST is walked,
// but where a field is gathered from each thread.
template <typename Derived>
class TraverseGatherPST : public TraversePST {
  public:
    explicit TraverseGatherPST(PST node_pst, std::string_view service_name, int nInBytes, int nOutBytes)
        : TraversePST(node_pst, service_name, nInBytes, nOutBytes) {}
    explicit TraverseGatherPST(PST node_pst, std::string_view service_name, int nInBytes = 0)
        : TraversePST(node_pst, service_name, nInBytes, node_pst->mdl()->Threads() * sizeof(typename Derived::output)) {}
    virtual ~TraverseGatherPST() = default;

  protected:
    virtual int Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> vout) final {
        using output = typename Derived::output;
        int nBytesLower = pst->n_lower() * sizeof(output);
        int nBytesUpper = pst->n_upper() * sizeof(output);
        assert(nBytesLower + nBytesUpper <= vout.size());
        auto out = vout.subspan(0, nBytesLower + nBytesUpper);
        auto out1 = out.subspan(0, nBytesLower);
        auto out2 = out.subspan(nBytesLower, nBytesUpper);
        auto nOut = nBytesLower + nBytesUpper;
        auto rID = ReqService(pst, in);
        Traverse(pst->lower(), in, out1);
        nBytesUpper = GetReply(rID, out2);
        return nOut;
    }
};

// This function more generally handles the common case where the return
// needs to simply add the lower and upper values. This works for any
// type that provides a "+=" function.
template <class Derived>
class TraverseCount : public TraverseCombinePST {
  public:
    explicit TraverseCount(PST pst, std::string_view service_name, int nInBytes = 0)
        : TraverseCombinePST(pst, service_name, nInBytes, sizeof(typename Derived::output)) {}

  protected:
    int Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST) final {
        using output = typename Derived::output;
        auto out1 = reinterpret_cast<output*>(vout1.data());
        auto out2 = reinterpret_cast<const output*>(vout2.data());
        *out1 += *out2;
        return sizeof(output);
    }
    virtual int Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) override = 0;
};

#endif
