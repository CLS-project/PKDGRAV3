#ifndef IO_HINCLUDED
#define IO_HINCLUDED
#include "iohdf5.h"
#include "mdl.h"

typedef struct ioContext {
    MDL mdl;
    int nExpected;
    int nReceived;
    double dTime;

    int N;
    int iMinOrder;
    int iMaxOrder;

    ioV3 *r;       /* Position */
    ioV3 *v;       /* Velocity */
    FLOAT *m;      /* Mass */
    FLOAT *s;      /* Softening */

} * IO;

void ioInitialize(IO *,MDL);
void ioAddServices(IO io,MDL mdl);

enum io_services {
    IO_SRV_STOP,
    IO_START_SAVE,
    IO_START_RECV,
    IO_RECV_DATA,
};

/* IO_START_SAVE */
struct inStartSave {
    double dTime;
    int    N;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioStartSave(IO,void *,int,void *,int *);

/* IO_START_RECV */
struct inStartRecv {
    double dTime;
    int iIndex;
    int nCount;
    char achOutName[PST_FILENAME_SIZE];
    };
void ioStartRecv(IO,void *,int,void *,int *);

/* IO_RECV_DATA */
struct inRecvData {
    };
void ioStartRecv(IO,void *,int,void *,int *);




#endif
