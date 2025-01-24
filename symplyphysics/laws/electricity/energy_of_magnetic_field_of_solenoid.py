"""
Energy of magnetic field of coil
================================

Energy of a coil depends on the material of the core, the mangetic field strength within
the coil and volume of the core.

**Links:**

#. `Physics LibreTexts, derivable from 14.4.1 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/14%3A_Inductance/14.04%3A_Energy_in_a_Magnetic_Field>`__.

..
    NOTE replace `H` with `B`?
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols

energy = symbols.energy
"""
:symbols:`energy` of the coil.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium.
"""

mangetic_field_strength = symbols.magnetic_field_strength
"""
:symbols:`magnetic_field_strength`.
"""

volume = symbols.volume
"""
:symbols:`volume` of the coil.
"""

law = Eq(energy, quantities.vacuum_permeability * relative_permeability * mangetic_field_strength**2 * volume / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permeability_=relative_permeability, intensity_=mangetic_field_strength, volume_=volume)
@validate_output(energy)
def calculate_energy(relative_permeability_: float, intensity_: Quantity,
    volume_: Quantity) -> Quantity:
    result_energy_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_energy_expr.subs({
        relative_permeability: relative_permeability_,
        mangetic_field_strength: intensity_,
        volume: volume_
    })
    return Quantity(result_expr)
