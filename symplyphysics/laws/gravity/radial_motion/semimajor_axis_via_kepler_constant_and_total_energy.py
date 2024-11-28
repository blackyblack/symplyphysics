r"""
Semimajor axis via Kepler's constant and total energy
=====================================================

The semi-major axis of an orbiting planet depends on the Kepler's constant of
the star—planet system and the total energy of the planet per unit of its mass.

**Notes:**

#. Works for both elliptical (:math:`\varepsilon < 0`) and hyperbolical
   (:math:`\varepsilon > 0`) orbits.
"""

from sympy import Eq, pi
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
)

semimajor_axis = symbols.semimajor_axis
"""
The :symbols:`semimajor_axis` of the planet's orbit.
"""

kepler_constant = SymbolNew("K", units.length**3 / units.time**2, display_latex="\\mathfrak{K}")
"""
The Kepler's constant, whose value is determined by the :symbols:`mass`
of the orbited star.
"""

specific_energy = symbols.specific_energy
"""
:symbols:`specific_energy` of the planet. Can be negative or positive depending on the sign of
the planet's energy.
"""

law = Eq(semimajor_axis, 2 * pi**2 * kepler_constant / abs(specific_energy))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    kepler_constant_=kepler_constant,
    total_energy_per_unit_mass_=specific_energy,
)
@validate_output(semimajor_axis)
def calculate_semimajor_axis(
    kepler_constant_: Quantity,
    total_energy_per_unit_mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        kepler_constant: kepler_constant_,
        specific_energy: total_energy_per_unit_mass_,
    })
    return Quantity(result)
