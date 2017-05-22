#ifndef IOMODULE_H
#define IOMODULE_H
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <sys/types.h>

#if defined(HAVE_LIBAIO)
#include <libaio.h>
#elif defined(HAVE_AIO_H)
#include <aio.h>
#endif
#define IO_MAX_ASYNC_COUNT 8
typedef struct {
#if defined(HAVE_LIBAIO)
    struct iocb cb[IO_MAX_ASYNC_COUNT];
    struct io_event events[IO_MAX_ASYNC_COUNT];
    io_context_t ctx;
#elif defined(HAVE_AIO_H)
    struct aiocb cb[IO_MAX_ASYNC_COUNT];
    struct aiocb const * pcb[IO_MAX_ASYNC_COUNT];
#endif
    char *pBuffer[IO_MAX_ASYNC_COUNT];
    off_t iFilePosition;   /* File position */
    size_t nBufferSize;    /* Total size of the buffer */
    size_t iByte;          /* Index into current buffer */
    int nBuffers;          /* Total number of buffers */
    int iBuffer;           /* Current buffer */
    int nPending;
    int nPageSize;
    int fd;
    } asyncFileInfo;


void io_init(asyncFileInfo *info, size_t nBuffers,size_t nBufferSize);
void io_free(asyncFileInfo *info);
int io_create(asyncFileInfo *info, const char *pathname);
void io_write(asyncFileInfo *info, void *buf, size_t count);
void io_close(asyncFileInfo *info);

#endif
