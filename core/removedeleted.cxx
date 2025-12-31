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
#include "removedeleted.h"

// Make sure that the communication structure is "trivial" so that it
// can be moved around with "memcpy" which is required for MDL.
static_assert(std::is_void_v<ServiceRemoveDeleted::input> || std::is_standard_layout_v<ServiceRemoveDeleted::input>);
static_assert(std::is_void_v<ServiceRemoveDeleted::output> || std::is_standard_layout_v<ServiceRemoveDeleted::output>);

int ServiceRemoveDeleted::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) {
    static_assert(std::is_void_v<input>);
    auto pkd = pst->lcl()->pkd;
    // Move the deleted particles to the end of the array
    auto i = std::partition(pkd->particles.begin(), pkd->particles.end(), [](auto& p) { return !p.is_deleted(); });

    // If any, put their decorators also at the end of the array
    int remove_decorator = 0;
    for (auto p = i; p < pkd->particles.end(); p++)
        if (p->have_decorator()) remove_decorator++;

    auto objective_decorators = pkd->decorators.size() - remove_decorator;
    if (remove_decorator) {
        auto ephem = pkd->Ephemeral<uint32_t>();

        // The ephemeral storage is populated as a map between decorators and
        // their particles. This way, there is no need to do complex searches
        for (uint32_t i = 0; i < pkd->Local(); ++i) {
            if (pkd->particles[i].have_decorator()) {
                ephem[pkd->particles[i].decorator_index()] = i;
            }
        }

        for (auto p = i; p < pkd->particles.end(); p++) {
            auto prev_decor = p->decorator_field();
            // Objective index for the decorator
            auto iDecorator = pkd->decorators.size() - remove_decorator;

            // Get the particle
            auto particle_to_swap = pkd->particles[ephem[iDecorator]];
            auto new_decor = particle_to_swap.decorator_field();

            // Swap the decorators
            auto a = pkd->decorators[prev_decor];
            auto b = pkd->decorators[new_decor];
            swap(a, b);

            // Update the mapping
            ephem[p->decorator_index()] = ephem[iDecorator];
            ephem[iDecorator] = ephem[p->decorator_index()];

            // Update the particles' decorators pointers
            particle_to_swap.set_decorator(prev_decor);
            p->set_decorator(new_decor);
            remove_decorator--;
        }
    }

    // ... and remove them
    pkd->SetLocal(i - pkd->particles.begin());
    // Count the particles species
    ServiceCountSpecies::Service(pst, vin, vout);
    auto out = reinterpret_cast<output*>(vout.data());
    auto nDecorator = out->counts[FIO_SPECIES_SPH] + out->counts[FIO_SPECIES_STAR] + out->counts[FIO_SPECIES_BH];

    pkd->decorators.SetLocal(nDecorator);
    assert(nDecorator == objective_decorators);
    return sizeof(output);
}
