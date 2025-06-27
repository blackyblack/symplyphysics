"""
Energy of electron in hydrogen atom per Bohr
============================================

In the Bohr's model of the Hydrogen atom, the electron is viewed as moving along a circular orbit
around the nucleus (which was further generalized onto elliptical orbits by Sommerfeld). In the
circular case, the energy of the atom can be expressed directly as a function of the electron's
orbit radius.

**Notation:**

#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`vacuum_permittivity`.

**Links:**

#. Formula 13.8 on p. 71 of "General Course of Physics" (Obschiy kurs fiziki), vol. 5, part 1 by
   Sivukhin D.V. (1979).

..
    TODO: find English link
"""

from sympy import Eq, solve, pi
from symplyphysics import Quantity, validate_input, validate_output, quantities, symbols

energy = symbols.energy
"""
:symbols:`energy` of the Hydrogen atom.
"""

radius = symbols.radius
"""
:symbols:`radius` of the electron's orbit.
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
