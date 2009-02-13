#ifndef IO_HINCLUDED
#define IO_HINCLUDED
#define IO_H_MODULE_ID "$Id$"
#ifdef USE_HDF5
#include "iohdf5.h"
#endif
#include "mdl.h"

#define MAX_IO_CLASSES 256
typedef struct {
    uint64_t iMinOrder;
    uint64_t iMaxOrder;
    FLOAT    dMass;
    FLOAT    dSoft;
    } ioClass;

typedef struct {
    FLOAT v[3];
    } ioV3;

typedef struct ioContext {
    MDL mdl;
    double dTime;

    total_t iMinOrder;
    total_t iMaxOrder;

    local_t nAllocated;  /* Total allocated on this processor */
    local_t N;           /* Total in use on this processor */

    local_t nExpected;   /* Number left to be received */

    total_t nTotal;      /* Total on all processors */
    total_t iOrder;

    ioV3 *r;             /* Position */
    ioV3 *v;             /* Velocity */
    float *d;            /* Density */
    float *p;            /* Potential */
    uint8_t *vClass;     /* Class */

    int      nClasses;   /* Number of particle classes */
    ioClass ioClasses[MAX_IO_CLASSES];
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
    IO_START_SEND,
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
    int    iStandard;
    int    bHDF5;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioStartSave(IO,void *,int,void *,int *);

/* IO_START_RECV */
struct inStartRecv {
    double dTime;
    double dEcosmo;
    double dTimeOld;
    double dUOld;
    total_t N;
    total_t iIndex;
    local_t nCount;
    int    bCheckpoint;
    int    iStandard;
    int    bHDF5;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioStartRecv(IO,void *,int,void *,int *);

/* IO_START_SEND */
struct inStartSend {
    char achInName[PST_FILENAME_SIZE];
    total_t N;            /* total number of particles */
    total_t iMinOrder;
    total_t iMaxOrder;
    int iFirstFile;       /* index of first file to read */
    int iLastFile;        /* index of last file */
    total_t iFirstOffset; /* index of first particle in the first file */
    total_t iLastOffset;  /* index of the last particle in the last file */
    };
void ioStartSend(IO,void *,int,void *,int *);


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
    float scale;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioMakePNG(IO,void *,int,void *,int *);
#endif

#endif
