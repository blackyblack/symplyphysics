"""
Capacitance of spherical capacitor
==================================

A spherical capacitor is composed of two concentric spheres with the space between them
filled with a dielectric medium. See `Figure`_.

.. _Figure: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/capsph.html
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

capacitance = Symbol("capacitance", units.capacitance)
"""
Capacitance of the capacitor.

Symbol:
    :code:`C`
"""

absolute_permittivity = Symbol("absolute_permittivity", units.capacitance / units.length)
r"""
Absolute permittivity of the medium between the spheres.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

inner_radius = Symbol("inner_radius", units.length)
r"""
Radius of the inner sphere.

Symbol:
    :code:`r_1`

Latex:
    :math:`r_1`
"""

outer_radius = Symbol("outer_radius", units.length)
r"""
Radius of the outer sphere.

Symbol:
    :code:`r_2`

Latex:
    :math:`r_2`
"""

law = Eq(
    capacitance, 4 * pi * absolute_permittivity * inner_radius * outer_radius /
    (outer_radius - inner_radius))
r"""
:code:`C = 4 * pi * epsilon * r_1 * r_2 / (r_2 - r_1)`

Latex:
    .. math::
        C = \frac{4 \pi \varepsilon r_1 r_2}{r_2 - r_1}
"""


@validate_input(relative_permittivity_=absolute_permittivity,
    inner_radius_=inner_radius,
    outer_radius_=outer_radius)
@validate_output(capacitance)
def calculate_capacity(relative_permittivity_: float, inner_radius_: Quantity,
    outer_radius_: Quantity) -> Quantity:
    result_expr = solve(law, capacitance, dict=True)[0][capacitance]
    result_expr = result_expr.subs({
        absolute_permittivity: relative_permittivity_,
        inner_radius: min([inner_radius_, outer_radius_], key=lambda x: x.scale_factor),
        outer_radius: max([inner_radius_, outer_radius_], key=lambda x: x.scale_factor)
    })
    return Quantity(result_expr)
