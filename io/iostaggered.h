#ifndef IO_STAGGERED
#define IO_STAGGERED
#include <stdint.h>

void io_staggered_read(const char *filename, void *buffer1, void *buffer2, uint64_t nElements, uint64_t nBytes1, uint64_t nBytes2, uint64_t iOffset=0);
void io_staggered_write(const char *filename, void *buffer1, void *buffer2, uint64_t nElements, uint64_t nBytes1, uint64_t nBytes2, bool bAppend=false);

#endif /* IO_STAGGERED */

