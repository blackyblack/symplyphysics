r"""
Light frequency change is proportional to gravitational potential change
========================================================================

When light is propagating in a gravitational field, its frequency changes proportionally to the
change in the potential of the gravitational field.

Let us consider an infinitesimally small section :math:`d \vec r` of the light's path, such that
the frequency of light is constant within that section. In that case we can obtain a dependency
between the change in light's frequency and the change in the gravitational potential.

**Notes:**

    .. _gravitational potential:

#. The *gravitational potential* :math:`\varphi` is defined as a scalar quantity such that the equation
   :math:`\vec g = - \nabla \varphi` holds where :math:`\vec g` is the vector of acceleration due to gravity
   and :math:`\nabla` is the nabla operator.

#. :math:`d \varphi = - \left( \vec g, d \vec r \right)` where :math:`\left( \vec a_1, \vec a_2 \right)` is
   the dot product between :math:`\vec a_1` and :math:`\vec a_2`.

**Links:**

#. Formula 72.4 on p. 378 of "General Course of Physics" (Obschiy kurs fiziki), vol. 1 by Sivukhin D.V. (1979).
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    SymbolNew,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)

frequency_change = clone_as_symbol(
    symbols.temporal_frequency,
    display_symbol="d(f)",
    display_latex="df",
)
r"""
The infinitesimal change in :symbols:`temporal_frequency` after passing an infinitesimal section
:math:`d \vec r`.
"""

frequency = symbols.temporal_frequency
r"""
The :symbols:`temporal_frequency` of light within an infinitesimal section :math:`d \vec r`.
"""

gravitational_potential_change = SymbolNew("d(phi)", units.velocity**2, display_latex="d \\phi")
r"""
The infinitesimal change in :ref:`gravitational potential <gravitational potential>` after passing an
infinitesimal section :math:`d \vec r`.
"""

law = Eq(
    frequency_change / frequency,
    -1 * gravitational_potential_change / quantities.speed_of_light**2,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    frequency_=frequency,
    gravitational_potential_change_=gravitational_potential_change,
)
@validate_output(frequency_change)
def calculate_frequency_change(
    frequency_: Quantity,
    gravitational_potential_change_: Quantity,
) -> Quantity:
    expr = solve(law, frequency_change)[0]
    result = expr.subs({
        frequency: frequency_,
        gravitational_potential_change: gravitational_potential_change_,
    })
    return Quantity(result)
