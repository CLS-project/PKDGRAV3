#ifndef IO_HINCLUDED
#define IO_HINCLUDED
#include "iohdf5.h"
#include "mdl.h"

typedef struct ioContext {
    MDL mdl;
    double dTime;

    total_t iMinOrder;
    total_t iMaxOrder;

    local_t nAllocated;  /* Total allocated on this processor */
    local_t N;           /* Total in use on this processor */
    local_t nExpected;   /* Number left to be received */

    ioV3 *r;       /* Position */
    ioV3 *v;       /* Velocity */
    //FLOAT *m;      /* Mass */
    //FLOAT *s;      /* Softening */
    float *d;      /* Density */
    float *p;      /* Potential */

} * IO;

void ioInitialize(IO *,MDL);
void ioAddServices(IO io,MDL mdl);

enum io_services {
    /* These are requested by the master */
    IO_SRV_STOP,
    IO_SETUP,
    IO_START_SAVE,

    /* These are requested by the I/O master */
    IO_ALLOCATE,
    IO_START_RECV,
    IO_MAKE_PNG
};

/* IO_SETUP */
struct inIOSetup {
    total_t N;
    };
void ioSetup(IO,void *,int,void *,int *);

/* IO_START_SAVE */
struct inStartSave {
    double dTime;
    double dEcosmo;
    double dTimeOld;
    double dUOld;
    total_t N;
    int    bCheckpoint;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioStartSave(IO,void *,int,void *,int *);

/* IO_START_RECV */
struct inStartRecv {
    double dTime;
    double dEcosmo;
    double dTimeOld;
    double dUOld;
    total_t iIndex;
    local_t nCount;
    int    bCheckpoint;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioStartRecv(IO,void *,int,void *,int *);

/* IO_ALLOCATE */
struct inIOAllocate {
    local_t nCount;
    };
void ioAllocate(IO,void *,int,void *,int *);

#ifdef USE_PNG
/* IO_MAKE_PNG */
struct inMakePNG {
    uint_fast32_t iResolution;  /* Image resolution RxR */
    float minValue;
    float maxValue;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioMakePNG(IO,void *,int,void *,int *);
#endif

#endif
