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

#ifndef ABORT_ERROR_H
#define ABORT_ERROR_H
#include <string_view>
#include <type_traits>
#include <utility>   // For std::forward
#include <stdexcept> // for std::runtime_error

class pkdgrav_error : public std::runtime_error {
  public:
    explicit pkdgrav_error(int error_code);
    explicit pkdgrav_error(int error_code, const char* func, const char* file, int line);
    explicit pkdgrav_error(int error_code, const char* func, const char* file, int line, const std::string_view message);
};

[[noreturn]] void abort_error_print(int error_code, const char* func, const char* file, int line);
#define abort_assert(code, condition)                                                                                                                          \
    if (!(condition)) {                                                                                                                                        \
        abort_error_print(code, __func__, __FILE__, __LINE__);                                                                                                 \
    }

[[noreturn]] void abort_error_print(int code, const char* func, const char* file, int line, const std::string_view message);
#define abort_false(code, condition, message)                                                                                                                  \
    if (!(condition)) {                                                                                                                                        \
        abort_error_print(code, __func__, __FILE__, __LINE__, message);                                                                                        \
    }

#endif
