"""
Capacitance of spherical capacitor
==================================

A spherical capacitor is composed of two concentric spheres with the space between them
filled with a dielectric medium. See `Figure`_.

.. _Figure: http://hyperphysics.phy-astr.gsu.edu/hbase/electric/capsph.html

**Links:**

#. `Spherical capacitor <http://hyperphysics.phy-astr.gsu.edu/hbase/electric/capsph.html>`__.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
)

capacitance = SymbolNew("C", units.capacitance)
"""
Capacitance of the capacitor.
"""

absolute_permittivity = SymbolNew("epsilon", units.capacitance / units.length, display_latex="\\varepsilon")
"""
Absolute permittivity of the medium between the spheres.
"""

inner_radius = SymbolNew("r_in", units.length, display_latex="r_\\text{in}")
"""
Radius of the inner sphere.
"""

outer_radius = SymbolNew("r_out", units.length, display_latex="r_\\text{out}")
"""
Radius of the outer sphere.
"""

law = Eq(
    capacitance, 4 * pi * absolute_permittivity * inner_radius * outer_radius /
    (outer_radius - inner_radius))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=absolute_permittivity,
    inner_radius_=inner_radius,
    outer_radius_=outer_radius)
@validate_output(capacitance)
def calculate_capacity(absolute_permittivity_: Quantity, inner_radius_: Quantity,
    outer_radius_: Quantity) -> Quantity:
    result_expr = solve(law, capacitance, dict=True)[0][capacitance]
    result_expr = result_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        inner_radius: min([inner_radius_, outer_radius_], key=lambda x: x.scale_factor),
        outer_radius: max([inner_radius_, outer_radius_], key=lambda x: x.scale_factor)
    })
    return Quantity(result_expr)
