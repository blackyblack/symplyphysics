"""
Wave impedance from permeability and permittivity
=================================================

The impedance of a wave is a value determined by the ratio of the transverse component
of the electric field to the transverse component of the magnetic field of a traveling
wave.

**Notation:**

#. :quantity_notation:`vacuum_impedance`.

**Conditions:**

#. The dielectric medium is ideal, i.e. its conductivity :math:`sigma = 0`.

**Links:**

#. `Wikipedia, derivable from last formula in paragraph <https://en.wikipedia.org/wiki/Wave_impedance#Definition>`__.

#. `Physics LibreTexts — Plane Waves, formula 7.7 <https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Essential_Graduate_Physics_-_Classical_Electrodynamics_(Likharev)/07%3A_Electromagnetic_Wave_Propagation/7.01%3A_Plane_Waves>`__.

..
    TODO: rename file
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

wave_impedance = symbols.wave_impedance
"""
:symbols:`wave_impedance`.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the medium.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium.
"""

law = Eq(wave_impedance,
    quantities.vacuum_impedance * sqrt(relative_permeability / relative_permittivity))
"""
:laws:symbol::

:laws:latex::
"""

# For the derivation see the Physics LibreTexts link
# TODO: add the definition of wave impedance as `eta = E / H`


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability)
@validate_output(wave_impedance)
def calculate_resistance(relative_permittivity_: float, relative_permeability_: float) -> Quantity:
    result_expr = solve(law, wave_impedance, dict=True)[0][wave_impedance]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
    })
    return Quantity(result_expr)
