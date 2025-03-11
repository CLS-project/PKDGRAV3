/*
* IA: Extracted from align.h and memuse.h
*/

#ifndef INLINE
    #if defined(__INTEL_COMPILER)
        #define INLINE extern inline
    #else
        #define INLINE extern inline
    #endif
#endif


/**
 * @brief allocate aligned memory. The use and results are the same as the
 *        posix_memalign function. This function should be used for any
 *        significant allocations and consistently labelled.
 *
 * @param label a symbolic label for the memory, i.e. "parts".
 * @param memptr pointer to the allocated memory.
 * @param alignment alignment boundary.
 * @param size the quantity of bytes to allocate.
 * @result zero on success, otherwise an error code.
 */
INLINE int swift_memalign(const char *label,
                          void **memptr,
                          size_t alignment,
                          size_t size) {
    int result = posix_memalign(memptr, alignment, size);
#ifdef SWIFT_MEMUSE_REPORTS
    if (result == 0) {
        memuse_log_allocation(label, *memptr, 1, size);
    }
    else {
        /* Failed allocations are interesting as well. */
        memuse_log_allocation(label, NULL, -1, size);
    }
#endif
    return result;
}



#define SWIFT_STRUCT_ALIGNMENT 32


/**
 * @brief Macro to tell the compiler that a given array has the specified
 * alignment.
 *
 * Note that this turns into a no-op but gives information to the compiler.
 * For GCC versions older than 4.6 this is ignored as the builtin does not
 * exist.
 *
 * @param type The type of the array.
 * @param array The array.
 * @param alignment The alignment in bytes of the array.
 */
#if defined(__ICC)
#define swift_align_information(type, array, alignment) \
  __assume_aligned(array, alignment);
#elif (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ > 6)
#define swift_align_information(type, array, alignment) \
  array = (type *)__builtin_assume_aligned(array, alignment);
#else
#define swift_align_information(type, array, alignment) ;
#endif





