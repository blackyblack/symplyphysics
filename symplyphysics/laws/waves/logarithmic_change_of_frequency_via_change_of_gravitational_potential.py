"""
Logarithmic change of wave frequency via change of gravitational potential
==========================================================================

TODO
"""

from sympy import Eq, log, solve
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

frequency_after = Symbol("frequency_after", units.frequency)
r"""
TODO

Symbol:
    :math:`\nu_2`
"""

frequency_before = Symbol("frequency_before", units.frequency)
r"""
TODO

Symbol:
    :math:`\nu_1`
"""

gravitational_potential_change = Symbol("gravitational_potential_change", units.energy)
r"""
TODO

Symbol:
    :math:`\Delta \phi`
"""

law = Eq(
    log(frequency_after / frequency_before),
    gravitational_potential_change / units.speed_of_light**2,
)
r"""
:math:`\log{\frac{\nu_2}{\nu_1}} = \frac{\Delta \phi}{c^2}`
"""


@validate_input(
    frequency_before_=frequency_before,
    gravitational_potential_change_=gravitational_potential_change,
)
@validate_output(frequency_after)
def calculate_frequency_after(
    frequency_before_: Quantity,
    gravitational_potential_change_: Quantity,
) -> Quantity:
    expr = solve(law, frequency_after)[0]
    result = expr.subs({
        frequency_before: frequency_before_,
        gravitational_potential_change: gravitational_potential_change_,
    })
    return Quantity(result)
