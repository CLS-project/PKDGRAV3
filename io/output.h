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
#include "io/service.h"

typedef enum {
    OUT_TINY_GROUP,
    OUT_KGRID,
    OUT_RGRID,
} outType;

struct inOutput {
    int iProcessor;   /* Output number: 0 to nParaWrite */
    int nProcessor;   /* Number of processors left in parallel */
    int iPartner;     /* Who to send the data to */
    int nPartner;     /* How many partners there are */
    outType eOutputType;  /* What kind of output */
    int iGrid;        /* Which grid (for grid output) */
    char achOutFile[PATH_MAX];
};

/* PST_OUTPUT */
int pstOutput(PST,void *,int,void *,int);
struct inOutputSend {
    int iPartner;     /* Who to send the data to */
    int iGrid;
    outType eOutputType;  /* What kind of output */
};

void pkdOutputSend(PKD pkd, outType eOutputType, int iPartner);
void pkdOutput(PKD pkd, outType eOutputType, int iProcessor,int nProcessor,
               int iPartner,int nPartner, const char *fname );

class ServiceSaveGroupStats : public ServiceOutput {
public:
    static constexpr std::string_view name = "SaveGroupStats";
    using input = ServiceOutput::input;
    using output = void;

    explicit ServiceSaveGroupStats(PST pst)
        : ServiceOutput(pst,PST_SAVE_GROUP_STATS,name,sizeof(input)) {}
protected:
    virtual void Write(PST pst,void *vin,int nIn,int iGroup,const std::string &filename,int iSegment,int nSegment) override;
};

#endif
