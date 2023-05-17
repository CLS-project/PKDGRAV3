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
#include "rshalocount.h"
#include <fstream>
using std::ifstream;
using std::ios_base;
namespace rs {
#include "rockstar/io/io_internal.h"
}

uint64_t ServiceRsHaloCount::GetSize(const std::string &filename,uint64_t file_size) {
    std::ifstream fhalo(filename,std::ios_base::in | std::ios_base::binary);
    rs::binary_output_header hdr;
    fhalo.read(reinterpret_cast<std::ifstream::char_type *>(&hdr),sizeof(hdr));      assert(fhalo.good());
    return hdr.num_halos;
}
