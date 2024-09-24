r"""
Avogadro constant is particle count over amount of substance
============================================================

The Avogadro constant is the constant of proportionality between particle count
and amount of substance.

**Notation:**

#. :quantity_notation:`avogadro_constant`.
"""

from sympy import Eq, solve
from symplyphysics import (
    quantities,
    Quantity,
    convert_to_float,
    validate_input,
    validate_output,
    symbols,
)

particle_count = symbols.particle_count
"""
Number of particles.
"""

amount_of_substance = symbols.amount_of_substance
"""
Amount of substance.
"""

law = Eq(quantities.avogadro_constant, particle_count / amount_of_substance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(mole_count_=amount_of_substance)
@validate_output(particle_count)
def calculate_particles_count(mole_count_: Quantity) -> int:
    solved = solve(law, particle_count, dict=True)[0][particle_count]
    result_expr = solved.subs(amount_of_substance, mole_count_)
    result = Quantity(result_expr)
    return int(convert_to_float(result))
