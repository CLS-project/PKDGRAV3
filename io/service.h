#ifndef SERVICE_H
#define SERVICE_H
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
#include "TraversePST.h"
#include <string>

class ServiceIO : public TraversePST {
  public:
    // This should be initialized with the number of parallel readers/writers
    // in nSimultaneous. The number of files is nFiles. Immediately after this
    // header is an array of file_info[nFiles] with information about each file.
    // The nElements in the header is the sum of the nElements in all files.
    // Finally iBeg is set to zero and iEnd is set to the total nElements.
    // As the service descends the PST, nSimultaneous, iBeg and iEnd are updated.
    static constexpr std::string_view name = "InputOutput";
    struct input {
        uint32_t nSimultaneous; // Number at this level of the PST
        uint32_t iReaderWriter; // Index of the reader/write [0,nTotalActive)
        uint32_t nSegment;      // Number of threads that act serially together
        uint32_t iThread;       // Index of the first thread in this group
        char filename[PATH_MAX];
    };
    typedef void output;
    explicit ServiceIO(PST pst, std::string_view service_name = name, int nInBytes = 0) : TraversePST(pst, service_name, nInBytes + sizeof(input)) {}

  private:
    virtual int Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) final;
    virtual int Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) final;

  protected:
    virtual int do_lower(PST pst, std::span<std::byte> in, std::span<std::byte> out);
    virtual int do_upper(PST pst, std::span<std::byte> in);
    /// @brief Perform the I/O operation
    /// @param pst The left PST
    /// @param in Input data
    /// @param iGroup The group index between 0 and nSimultaneous-1 (on service start)
    /// @param iSegment The segment within this group (0 to nSegment-1)
    /// @param nSegment The total number of segments in this group (at least one)
    /// This function is called on every thread, but only one segment at a time per group.
    /// Groups are independent and will run in parallel.
    virtual void IO(PST pst, std::span<std::byte> in, int iGroup, int iSegment, int nSegment) = 0;
};

class ServiceFileSizes : public TraversePST {
  public:
    static constexpr std::string_view name = "FileSizes";
    struct input {
        uint32_t nTotalActive;  // Number of simultaneous readers/writers
        uint32_t nSimultaneous; // Number at this level of the PST
        uint32_t iReaderWriter; // Index of the reader/write [0,nTotalActive)
        uint32_t nElementSize;  // Size of each element (1 for file size)
        char filename[PATH_MAX];
    };
    static constexpr int max_files = 100'000;
    struct output {
        uint64_t nFileBytes : 40;
        uint64_t iFileIndex : 24;
    };
    explicit ServiceFileSizes(PST pst, std::string_view service_name = name) : TraversePST(pst, service_name, sizeof(input), max_files * (sizeof(output))) {}

  private:
    virtual int Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) final;
    virtual int Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) final;

  protected:
    virtual uint64_t GetSize(const std::string& filename, uint64_t file_size);
};

class ServiceInput : public ServiceIO {
  public:
    // This should be initialized with the number of parallel readers/writers
    // in nSimultaneous. The number of files is nFiles. Immediately after this
    // header is an array of file_info[nFiles] with information about each file.
    // The nElements in the header is the sum of the nElements in all files.
    // Finally iBeg is set to zero and iEnd is set to the total nElements.
    // As the service descends the PST, nSimultaneous, iBeg and iEnd are updated.
    static constexpr std::string_view name = "Input";
    struct input {
        ServiceIO::input io;
        int nFiles;         // Number of different files
        uint64_t nElements; // Total number of elements
        uint64_t iBeg;      // Index of first element to read
        uint64_t iEnd;      // Index of one past last element to read
    };
    static constexpr int max_files = 100'000;
    using io_elements = uint64_t; // Read this many elements
    typedef void output;
    explicit ServiceInput(PST pst, std::string_view service_name = name, int nInBytes = 0)
        : ServiceIO(pst, service_name, nInBytes + sizeof(input) + max_files * (sizeof(io_elements))) {}

  protected:
    virtual int do_lower(PST pst, std::span<std::byte> in, std::span<std::byte> out) override;
    virtual int do_upper(PST pst, std::span<std::byte> in) override;
    virtual void IO(PST pst, std::span<std::byte> in, int iGroup, int iSegment, int nSegment) override;

  protected:
    virtual void start(PST pst, uint64_t nElements, std::span<std::byte> in);
    virtual void finish(PST pst, uint64_t nElements, std::span<std::byte> in);
    virtual void Read(PST pst, uint64_t iElement, const std::string& filename, uint64_t iBeg, uint64_t iEnd) = 0;
};

class ServiceOutput : public ServiceIO {
  public:
    // This should be initialized with the number of parallel readers/writers
    // in nSimultaneous. The number of files is nFiles. Immediately after this
    // header is an array of file_info[nFiles] with information about each file.
    // The nElements in the header is the sum of the nElements in all files.
    // Finally iBeg is set to zero and iEnd is set to the total nElements.
    // As the service descends the PST, nSimultaneous, iBeg and iEnd are updated.
    static constexpr std::string_view name = "Output";
    struct input {
        ServiceIO::input io;
    };
    typedef void output;
    explicit ServiceOutput(PST pst, std::string_view service_name = name, int nInBytes = 0) : ServiceIO(pst, service_name, nInBytes + sizeof(input)) {}

  protected:
    virtual void IO(PST pst, std::span<std::byte> in, int iGroup, int iSegment, int nSegment) override;

  protected:
    virtual void Write(PST pst, std::span<std::byte> in, int iGroup, const std::string& filename, int iSegment, int nSegment) = 0;
};

class ServiceSetTotal : public TraverseCombinePST {
  public:
    static constexpr std::string_view name = "SetTotal";
    using input = void;
    using output = uint64_t;
    explicit ServiceSetTotal(PST pst) : TraverseCombinePST(pst, name, 0, sizeof(output)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Combine(std::span<std::byte> vout1, std::span<const std::byte> vout2, std::span<const std::byte>, PST) -> int override;
};

class ServiceSetWriteStart : public TraversePST {
  public:
    static constexpr std::string_view name = "SetWriteStart";
    using input = uint64_t;
    using output = void;
    explicit ServiceSetWriteStart(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto Recurse(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};
#endif /* SERVICE_H */