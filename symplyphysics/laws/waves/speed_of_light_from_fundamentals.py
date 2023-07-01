from sympy import (Eq, sqrt)
from sympy.physics.units import speed_of_light, magnetic_constant, electric_constant
from symplyphysics import units, print_expression, convert_to

# Description
## Speed of light in vacuum is fundamental but still might be calculated from other fundamentals.

# Law: c = 1/sqrt(e0 * u0), where
## C is speed of light in vacuum,
## e0 is electric constant or vacuum permittivity,
## u0 is magnetic constant or vacuum permeability.

law = Eq(speed_of_light, 1 / sqrt(magnetic_constant * electric_constant))


def print_law() -> str:
    return print_expression(law)


assert convert_to(law.lhs, units.meter / units.second) == convert_to(law.rhs,
    units.meter / units.second)
