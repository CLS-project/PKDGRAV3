from libc.stdint cimport uint64_t

cdef extern from "core/memory.h":
    cdef cppclass EphemeralMemory:
        uint64_t per_particle
        uint64_t per_process

        EphemeralMemory()

        EphemeralMemory(uint64_t per_particle, uint64_t per_process)

        # Overload | operator for union-like behavior
        EphemeralMemory operator|(const EphemeralMemory &other) const

        # Overload + operator for sum-like behavior
        EphemeralMemory operator+(const EphemeralMemory &other) const