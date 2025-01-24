"""
Characteristic resistance of medium
===================================

The characteristic resistance of a wave is a value determined by the ratio of the
transverse component of the electric field strength to the transverse component of the
magnetic field strength of a traveling wave.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, sqrt
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

resistance = symbols.electrical_resistance
"""
:symbols:`electrical_resistance` of the medium.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the medium.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the medium.
"""

impedance_constant = Quantity(120 * pi * units.ohm, display_symbol="Z_0")
"""
Constant equal to :math:`120 \\Omega`.
"""

law = Eq(resistance, impedance_constant * sqrt(relative_permeability / relative_permittivity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability)
@validate_output(resistance)
def calculate_resistance(relative_permittivity_: float, relative_permeability_: float) -> Quantity:
    result_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
    })
    return Quantity(result_expr)
