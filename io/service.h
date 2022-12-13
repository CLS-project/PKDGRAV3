#ifndef B52D0DB0_06A9_416B_A559_F5010B417974
#define B52D0DB0_06A9_416B_A559_F5010B417974
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
    struct input {
        uint32_t nSimultaneous; // Number at this level of the PST
        uint32_t iReaderWriter; // Index of the reader/write [0,nTotalActive)
        uint32_t nSegment;      // Number of threads that act serially together
        uint32_t iThread;       // Index of the first thread in this group
        char filename[256];
    };
    typedef void output;
    explicit ServiceIO(PST pst,int service_id,int nInBytes=0,const char *name="InputOutput")
        : TraversePST(pst,service_id,nInBytes+sizeof(input),name) {}
private:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) final;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) final;
protected:
    virtual int do_lower(PST pst,void *vin,int nIn,void *vout,int nOut);
    virtual int do_upper(PST pst,void *vin,int nIn);
    virtual void IO(PST pst,void *vin,int nIn,int iGroup,int iSegment,int nSegment) = 0;
};

class ServiceFileSizes : public TraversePST {
public:
    struct input {
        uint32_t nTotalActive;  // Number of simultaneous readers/writers
        uint32_t nSimultaneous; // Number at this level of the PST
        uint32_t iReaderWriter; // Index of the reader/write [0,nTotalActive)
        uint32_t nElementSize;  // Size of each element (1 for file size)
        char filename[256];
    };
    static constexpr int max_files = 100'000;
    struct output {
        uint64_t nFileBytes : 40;
        uint64_t iFileIndex : 24;
    };
    explicit ServiceFileSizes(PST pst,int service_id=PST_FILE_SIZES,const char *name="FileSizes")
        : TraversePST(pst,service_id,sizeof(input),max_files*(sizeof(output)),name) {}
private:
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) final;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) final;
protected:
    virtual uint64_t GetSize(const std::string &filename,uint64_t file_size);
};

class ServiceInput : public ServiceIO {
public:
    // This should be initialized with the number of parallel readers/writers
    // in nSimultaneous. The number of files is nFiles. Immediately after this
    // header is an array of file_info[nFiles] with information about each file.
    // The nElements in the header is the sum of the nElements in all files.
    // Finally iBeg is set to zero and iEnd is set to the total nElements.
    // As the service descends the PST, nSimultaneous, iBeg and iEnd are updated.
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
    explicit ServiceInput(PST pst,int service_id,int nInBytes=0,const char *name="Input")
        : ServiceIO(pst,service_id,nInBytes+sizeof(input)+max_files*(sizeof(io_elements)),name) {}
protected:
    virtual int do_lower(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    virtual int do_upper(PST pst,void *vin,int nIn) override;
    virtual void IO(PST pst,void *vin,int nIn,int iGroup,int iSegment,int nSegment) override;
protected:
    virtual void start(PST pst,uint64_t nElements,void *vin,int nIn);
    virtual void finish(PST pst,uint64_t nElements,void *vin,int nIn);
    virtual void Read(PST pst,uint64_t iElement,const std::string &filename,uint64_t iBeg,uint64_t iEnd) = 0;
};

class ServiceOutput : public ServiceIO {
public:
    // This should be initialized with the number of parallel readers/writers
    // in nSimultaneous. The number of files is nFiles. Immediately after this
    // header is an array of file_info[nFiles] with information about each file.
    // The nElements in the header is the sum of the nElements in all files.
    // Finally iBeg is set to zero and iEnd is set to the total nElements.
    // As the service descends the PST, nSimultaneous, iBeg and iEnd are updated.
    struct input {
        ServiceIO::input io;
    };
    typedef void output;
    explicit ServiceOutput(PST pst,int service_id,int nInBytes=0,const char *name="Output")
        : ServiceIO(pst,service_id,nInBytes+sizeof(input),name) {}
protected:
    virtual void IO(PST pst,void *vin,int nIn,int iGroup,int iSegment,int nSegment) override;
protected:
    virtual void Write(PST pst,void *vin,int nIn,int iGroup,const std::string &filename,int iSegment,int nSegment) = 0;
};

#endif /* B52D0DB0_06A9_416B_A559_F5010B417974 */
