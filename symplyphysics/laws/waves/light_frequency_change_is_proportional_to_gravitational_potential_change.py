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
   :math:`\vec g = - \nabla \varphi` is hold where :math:`\vec g` is the vector of acceleration due to gravity
   and :math:`\nabla` is the nabla operator.

#. :math:`d \varphi = - \left( \vec g, d \vec r \right)` where :math:`\left( \vec a_1, \vec a_2 \right)` is
   the dot product between :math:`\vec a_1` and :math:`\vec a_2`.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

frequency_change = Symbol("frequency_change", units.frequency)
r"""
The infinitesimal change in frequency after passing an infinitesimal section :math:`d \vec r`.

Symbol:
    :code:`dnu`

Latex:
    :math:`d \nu`
"""

frequency = Symbol("frequency", units.frequency)
r"""
The frequency of light within an infinitesimal section :math:`d \vec r`.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

gravitational_potential_change = Symbol("gravitational_potential_change", units.velocity**2)
r"""
The infinitesimal change in :ref:`gravitational potential <gravitational potential>` after passing an infinitesimal section :math:`d \vec r`.

Symbol:
    :code:`dphi`

Latex:
    :math:`d \varphi`
"""

law = Eq(
    frequency_change / frequency,
    -1 * gravitational_potential_change / units.speed_of_light**2,
)
r"""
:code:`dnu / nu = -1 * dphi / c^2`

Latex:
    .. math::
        \frac{d \nu}{\nu} = - \frac{d \varphi}{c^2}
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
