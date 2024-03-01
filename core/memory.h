#ifndef PKD_CORE_MEMORY_H
#define PKD_CORE_MEMORY_H

#include <algorithm> // For std::max
#include <cstdint> // For std::uint64_t

class EphemeralMemory {
public:
    std::uint64_t per_particle = 0;
    std::uint64_t per_process = 0;

    // Default constructor
    EphemeralMemory(std::uint64_t per_particle = 0, std::uint64_t per_process = 0) : per_particle(per_particle), per_process(per_process) {}

    // Constructor with MSR, grid, and count
    EphemeralMemory(class MSR *msr, int grid, int count = 1);

    // Overload | operator for union-like behavior
    EphemeralMemory operator|(const EphemeralMemory &other) const {
        return EphemeralMemory(std::max(this->per_particle, other.per_particle), std::max(this->per_process, other.per_process));
    }

    // Overload + operator for sum-like behavior
    EphemeralMemory operator+(const EphemeralMemory &other) const {
        return EphemeralMemory(this->per_particle + other.per_particle, this->per_process + other.per_process);
    }

    // Overload |= operator for in-place union-like behavior
    EphemeralMemory &operator|=(const EphemeralMemory &other) {
        this->per_process = std::max(this->per_process, other.per_process);
        this->per_particle = std::max(this->per_particle, other.per_particle);
        return *this;
    }

    // Overload += operator for in-place sum-like behavior
    EphemeralMemory &operator+=(const EphemeralMemory &other) {
        this->per_process += other.per_process;
        this->per_particle += other.per_particle;
        return *this;
    }

};

#endif // PKD_CORE_MEMORY_H