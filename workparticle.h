#include "basetype.h"
#include "interaction/gravity_traits.h"

typedef std::variant<
workParticle *,
             GravityInteractionTraits::WORKPARTICLE *
             > workParticleVariants;
