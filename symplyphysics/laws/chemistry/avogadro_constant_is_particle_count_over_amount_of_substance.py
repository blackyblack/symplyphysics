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
    units,
    quantities,
    Quantity,
    Symbol,
    dimensionless,
    convert_to_float,
    validate_input,
    validate_output,
)

particle_count = Symbol("particle_count", dimensionless)
"""
Number of particles.

Symbol:
    :code:`N`
"""

amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
"""
Amount of substance.

Symbol:
    :code:`n`
"""

law = Eq(quantities.avogadro_constant, particle_count / amount_of_substance)
r"""
:code:`N_A = N / n`

Latex:
    .. math::
        N_\text{A} = \frac{N}{n}
"""


@validate_input(mole_count_=amount_of_substance)
@validate_output(particle_count)
def calculate_particles_count(mole_count_: Quantity) -> int:
    solved = solve(law, particle_count, dict=True)[0][particle_count]
    result_expr = solved.subs(amount_of_substance, mole_count_)
    result = Quantity(result_expr)
    return int(convert_to_float(result))
