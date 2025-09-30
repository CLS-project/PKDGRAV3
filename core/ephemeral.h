/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef EPHEMERAL_H
#define EPHEMERAL_H
#include <memory_resource>
#include <memory>

class ephemeral_buffer_resource : public std::pmr::memory_resource {
public:
    ephemeral_buffer_resource() = delete;
    explicit
    ephemeral_buffer_resource(void* __buffer, size_t __buffer_size,
                              memory_resource* __upstream) noexcept
        : M_current_buf(__buffer), M_avail(__buffer_size),
          M_orig_buf(__buffer), M_orig_size(__buffer_size) {
    }

    // Reset the buffer to its original state, releasing any allocated space
    void release() noexcept {
        M_current_buf = M_orig_buf;
        M_avail = M_orig_size;
    }

    // Returns the available number of elements of type T (for reserve())
    template<typename T>
    size_t elements() {
        return M_avail / sizeof(T);
    }

    // Reclaim buffer space after a raw pointer
    // Any data structures after this pointer cannot be used anymore.
    void reclaim(void* end_ptr) {
        char *new_current = static_cast<char *>(end_ptr);
        char *buffer_start = static_cast<char *>(M_orig_buf);
        char *buffer_end = static_cast<char *>(M_orig_buf) + M_orig_size;
        if (new_current < buffer_start || new_current > buffer_end) {
            throw std::runtime_error("ephemeral_buffer_resource::reclaim: pointer out of bounds");
        }
        if (new_current > static_cast<char*>(M_current_buf)) {
            throw std::runtime_error("ephemeral_buffer_resource::reclaim: reclaim pointer is after current buffer position");
        }
        M_current_buf = new_current;
        M_avail = buffer_end - new_current;
    }

    // Template helper for containers with .data() and .size()
    template<typename Container>
    void reclaim(const Container& container) {
        void *end_ptr = static_cast<void *>(
                            const_cast<typename Container::value_type *>(container.data() + container.size()));
        reclaim(end_ptr);
    }

protected:
    void *
    do_allocate(size_t bytes, size_t alignment) override {
        if (bytes==0) bytes=1; // Ensures we don't return the same pointer twice.
        void *p = std::align(alignment, bytes, M_current_buf, M_avail);
        if (p==nullptr) throw std::bad_alloc();
        M_current_buf = static_cast<char *>(M_current_buf) + bytes;
        M_avail -= bytes;
        return p;
    }

    void
    do_deallocate(void*, size_t, size_t) override { }

    bool
    do_is_equal(const memory_resource& __other) const noexcept override
    { return this == &__other; }

private:
    void *M_current_buf;
    size_t M_avail;
    void *M_orig_buf;
    size_t M_orig_size;


};

#endif