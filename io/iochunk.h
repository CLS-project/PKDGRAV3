#ifndef A5ADF993_BA19_4C3E_B6C9_7FE3E709953D
#define A5ADF993_BA19_4C3E_B6C9_7FE3E709953D
#include <stdint.h>

void io_chunk_read(const char *filename, void *buffer, uint64_t nBytes, uint64_t iOffset=0);
void io_chunk_write(const char *filename, const void *buffer, uint64_t nBytes, bool bAppend=false);

#endif /* A5ADF993_BA19_4C3E_B6C9_7FE3E709953D */
