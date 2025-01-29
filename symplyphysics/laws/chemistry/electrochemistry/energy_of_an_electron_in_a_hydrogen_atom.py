"""
Energy of electron in hydrogen atom per Bohr
============================================

There is an expression for the total energy of the hydrogen atom according to Bohr's
theory.

**Notation:**

#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`vacuum_permittivity`.

**Links:**

#. `Wikipedia, derivable from first formula and formula for kinetic energy <https://en.wikipedia.org/wiki/Classical_electron_radius#>`__.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (Quantity, validate_input, validate_output, quantities, symbols)

energy = symbols.energy
"""
Electron :symbols:`energy`.
"""

radius = symbols.radius
"""
Electron :symbols:`radius`.
"""

law = Eq(energy,
    quantities.elementary_charge**2 / (8 * pi * quantities.vacuum_permittivity * radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(radius_=radius)
@validate_output(energy)
def calculate_energy_of_electron(radius_: Quantity) -> Quantity:
    result_expr = solve(law, energy, dict=True)[0][energy]
    result = result_expr.subs(radius, radius_)
    return Quantity(result)
