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
#include "abort_error.h"
#include "mdlbt.h"
#include "fmt/format.h"  // This will be part of c++20
#include <fmt/color.h>
using namespace fmt::literals; // Gives us ""_a and ""_format literals
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <mutex>

static std::mutex error_mutex;

[[noreturn]] void abort_error_print(int code, const char *func, const char *file, int line) {
    error_mutex.lock();
    std::cerr << std::flush;
    fmt::print(stderr, fg(fmt::color::crimson) | fmt::emphasis::bold,
               "INTERNAL error {code}, function {func}, file {file}, line {line}\n"
               "see: https://pkdgrav3.readthedocs.io/en/latest/error_codes.html#error-{code}\n",
               "code"_a=code, "func"_a=func, "file"_a=file, "line"_a=line);
    mdlbt::show_backtrace(63,1);
    std::cerr << std::flush;
    error_mutex.unlock();
    std::exit(EXIT_FAILURE);
}

[[noreturn]] void abort_error_print_message(int code, const char *func, const char *file, int line, const std::string_view message) {
    error_mutex.lock();
    std::cerr << std::flush;
    fmt::print(stderr, fg(fmt::color::crimson) | fmt::emphasis::bold,
               "error {code}, function {func}, file {file}, line {line}\n"
               "see: https://pkdgrav3.readthedocs.io/en/latest/error_codes.html#error-{code}\n"
               "{message}: {error}\n"
               ,"code"_a=code, "func"_a=func, "file"_a=file, "line"_a=line, "message"_a=message,"error"_a=strerror(errno));
    mdlbt::show_backtrace(63,1);
    std::cerr << std::flush;
    error_mutex.unlock();
    std::exit(EXIT_FAILURE);
}