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
#include "service.h"
#include "fmt/format.h"
using namespace fmt::literals;
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

static_assert(std::is_void<ServiceFileSizes::input>()  || std::is_standard_layout<ServiceFileSizes::input>());

int ServiceFileSizes::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = pst->mdl;
    auto hdr  = static_cast<ServiceFileSizes::input *>(vin);
    auto out  = static_cast<ServiceFileSizes::output *>(vout);
    assert(hdr->nSimultaneous>0 && hdr->nSimultaneous<=pst->nLeaves);

    // We are still trying to launch parallel I/O so split the work
    if (hdr->nSimultaneous>1) {
        int nLower = hdr->nSimultaneous * pst->nLower / pst->nLeaves;
        if (nLower==0) nLower=1;
        hdr->nSimultaneous -= nLower;
        hdr->iReaderWriter += nLower;
        auto rID = ReqService(pst,vin,nIn);;
        hdr->nSimultaneous = nLower;
        hdr->iReaderWriter -= nLower;
        nOut = Traverse(pst->pstLower,vin,nIn,vout,nOut);
        nOut += mdl->GetReply(rID,out + nOut/sizeof(output));
    }
    // Now calculate the file sizes
    else {
        nOut = Service(pst,vin,nIn,vout,nOut);
    }
    return nOut;
}

int ServiceFileSizes::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto hdr  = static_cast<ServiceFileSizes::input *>(vin);
    auto out  = static_cast<ServiceFileSizes::output *>(vout);
    nOut = 0;
    struct stat s;
    std::string filename_template(hdr->filename);
    std::string filename;
    for (auto i = hdr->iReaderWriter; stat((filename = fmt::format(filename_template,"i"_a=i)).c_str(),&s)==0; i+=hdr->nTotalActive) {
        out->iFileIndex = i;
        out->nFileBytes = GetSize(filename,s.st_size) / hdr->nElementSize;
        nOut += sizeof(output);
        ++out;
    }
    return nOut;
}

uint64_t ServiceFileSizes::GetSize(const std::string &filename,uint64_t file_size) {
    return file_size;
}


static_assert(std::is_void<ServiceIO::input>()  || std::is_standard_layout<ServiceIO::input>());

int ServiceIO::Recurse(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = pst->mdl;
    auto hdr  = static_cast<ServiceIO::input *>(vin);
    assert(hdr->nSimultaneous>0 && hdr->nSimultaneous<=pst->nLeaves);

    // We are still trying to launch parallel I/O so split the work
    if (hdr->nSimultaneous>1) {
        int nLower = hdr->nSimultaneous * pst->nLower / pst->nLeaves;
        if (nLower==0) nLower=1;
        hdr->nSimultaneous -= nLower;
        hdr->iReaderWriter += nLower;
        hdr->iThread = pst->idUpper;
        auto rID = do_upper(pst,vin,nIn);
        hdr->nSimultaneous = nLower;
        hdr->iReaderWriter -= nLower;
        hdr->iThread = mdl->Self();
        do_lower(pst,vin,nIn,vout,nOut);
        return mdl->GetReply(rID,vout);
    }
    // Now we go sequentially
    else {
        if (hdr->nSegment==0) hdr->nSegment = pst->nLeaves;
        do_lower(pst,vin,nIn,vout,nOut);
        auto rID = do_upper(pst,vin,nIn);
        return mdl->GetReply(rID,vout);
    }
}

int ServiceIO::do_lower(PST pst,void *vin,int nIn,void *vout,int nOut) {
    return Traverse(pst->pstLower,vin,nIn,vout,nOut);
}

int ServiceIO::do_upper(PST pst,void *vin,int nIn) {
    return ReqService(pst,vin,nIn);
}

int ServiceIO::Service(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto mdl = pst->mdl;
    auto hdr  = static_cast<ServiceIO::input *>(vin);
    if (hdr->nSegment==0) hdr->nSegment = pst->nLeaves;
    auto iSegment = mdl->Self() - hdr->iThread;
    IO(pst,vin,nIn,hdr->iReaderWriter,iSegment,hdr->nSegment);
    return 0;
}


// Make sure that the communication structure is "standard" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void<ServiceInput::input>()  || std::is_standard_layout<ServiceInput::input>());
static_assert(std::is_void<ServiceInput::io_elements>()  || std::is_standard_layout<ServiceInput::io_elements>());
static_assert(std::is_void<ServiceInput::output>() || std::is_standard_layout<ServiceInput::output>());

auto get_elements = [](auto pst,auto hdr) {
    auto iBeg = hdr->iBeg, iEnd = hdr->iEnd, nElements = hdr->iEnd-hdr->iBeg;
    auto nElementsLower = (nElements + pst->nLeaves - 1) / pst->nLeaves * pst->nLower;
    auto iMid = iBeg + nElementsLower;
    return std::make_tuple(iBeg,iMid,iEnd);
};

int ServiceInput::do_lower(PST pst,void *vin,int nIn,void *vout,int nOut) {
    auto hdr  = static_cast<ServiceInput::input *>(vin);
    auto [iBeg,iMid,iEnd] = get_elements(pst,hdr);
    hdr->iEnd = iMid;
    auto n = Traverse(pst->pstLower,vin,nIn,vout,nOut);
    hdr->iEnd = iEnd;
    return n;
}

int ServiceInput::do_upper(PST pst,void *vin,int nIn) {
    auto hdr  = static_cast<ServiceInput::input *>(vin);
    auto [iBeg,iMid,iEnd] = get_elements(pst,hdr);
    hdr->iBeg = iMid;
    auto rID = pst->mdl->ReqService(pst->idUpper,getServiceID(),vin,nIn);
    hdr->iBeg = iBeg;
    return rID;
}

void ServiceInput::IO(PST pst,void *vin,int nIn,int iGroup,int iSegment,int nSegment) {
    auto hdr  = static_cast<ServiceInput::input *>(vin);
    auto elements = reinterpret_cast<ServiceInput::io_elements *>(hdr+1);
    auto next = reinterpret_cast<void *>(elements+hdr->nFiles);
    auto cnt = nIn - (reinterpret_cast<char *>(next) - reinterpret_cast<char *>(vin));

    start(pst,hdr->iEnd - hdr->iBeg,next,cnt);
    uint64_t iBeg = 0, iElement=0;
    for (auto iFile = 0; iFile<hdr->nFiles; iBeg+=elements[iFile++]) {
        if (iBeg > hdr->iEnd) break;                            // We are past our part
        uint64_t iEnd = iBeg + elements[iFile];                 // for this file
        if (iEnd <= hdr->iBeg) continue;                        // We aren't there yet
        uint64_t iFileBeg = std::max(iBeg,hdr->iBeg) - iBeg;    // our part of the file
        uint64_t iFileEnd = std::min(iEnd,hdr->iEnd) - iBeg;
        if (iFileBeg == iFileEnd) continue;
        assert(iFileBeg<=iFileEnd);
        // Now perform the I/O for element [iFileBeg,iFileEnd) on iFile
        auto filename = fmt::format(hdr->io.filename,"i"_a=iFile);
        Read(pst,iElement,filename,iFileBeg,iFileEnd);
        iElement += (iFileEnd-iFileBeg);
    }
    finish(pst,hdr->iEnd - hdr->iBeg,next,cnt);
}

void ServiceInput::start(PST pst,uint64_t nElements,void *vin,int nIn) {}
void ServiceInput::finish(PST pst,uint64_t nElements,void *vin,int nIn) {}

void ServiceOutput::IO(PST pst,void *vin,int nIn,int iGroup,int iSegment,int nSegment) {
    auto hdr  = static_cast<ServiceOutput::input *>(vin);
    auto filename = fmt::format(hdr->io.filename,"i"_a=iGroup);
    Write(pst,vin,nIn,iGroup,filename,iSegment,nSegment);
}