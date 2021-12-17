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

#ifndef OUTPUT_H
#define OUTPUT_H
#include "pkd.h"

typedef enum {
    OUT_TINY_GROUP,
    OUT_KGRID,
    OUT_RGRID,
} outType;

struct inOutputSend {
    int iPartner;     /* Who to send the data to */
    int iGrid;
    outType eOutputType;  /* What kind of output */
};

void pkdOutputSend(PKD pkd, outType eOutputType, int iPartner);
void pkdOutput(PKD pkd, outType eOutputType, int iProcessor,int nProcessor,
               int iPartner,int nPartner, const char *fname );

#endif
