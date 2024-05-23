cimport ephemeral

cdef class PyEphemeralMemory:
    cdef uint64_t per_particle
    cdef uint64_t per_process

    def __init__(self, uint64_t per_particle=0, uint64_t per_process=0):
        self.per_particle = per_particle
        self.per_process = per_process

    def __or__(self, other):
        cdef EphemeralMemory result = EphemeralMemory(self.per_particle, self.per_process) | EphemeralMemory(other.per_particle, other.per_process)
        return PyEphemeralMemory(result.per_particle, result.per_process)

    def __add__(self, other):
        cdef EphemeralMemory result = EphemeralMemory(self.per_particle, self.per_process) + EphemeralMemory(other.per_particle, other.per_process)
        return PyEphemeralMemory(result.per_particle, result.per_process)

    def __ior__(self, other):
        cdef EphemeralMemory result = EphemeralMemory(self.per_particle, self.per_process) | EphemeralMemory(other.per_particle, other.per_process)
        self.per_particle = result.per_particle
        self.per_process = result.per_process
        return self

    def __iadd__(self, other):
        cdef EphemeralMemory result = EphemeralMemory(self.per_particle, self.per_process) + EphemeralMemory(other.per_particle, other.per_process)
        self.per_particle = result.per_particle
        self.per_process = result.per_process
        return self

    @property
    def per_particle(self):
        return self.per_particle

    @property
    def per_process(self):
        return self.per_process
