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
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import (
    electric_field_outside_charged_sphere as _electric_field_law,
    voltage_is_line_integral_of_electric_field as _voltage_law,
    capacitance_from_charge_and_voltage as _capacitance_law,
)

capacitance = symbols.capacitance
"""
:symbols:`capacitance` of the capacitor.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium between the spheres.
"""

inner_radius = clone_as_symbol(symbols.distance_to_origin, display_symbol="r_in", display_latex="r_\\text{in}")
"""
Radius of the inner sphere. See :symbols:`distance_to_origin`.
"""

outer_radius = clone_as_symbol(symbols.distance_to_origin, display_symbol="r_out", display_latex="r_\\text{out}")
"""
Radius of the outer sphere. See :symbols:`distance_to_origin`.
"""

law = Eq(
    capacitance, 4 * pi * absolute_permittivity * inner_radius * outer_radius /
    (outer_radius - inner_radius))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law for spherical capacitor filled with vacuum

_distance = _voltage_law.distance

_electric_field_expr = _electric_field_law.law.rhs.subs({
    _electric_field_law.charge: _capacitance_law.charge,
    _electric_field_law.distance: _distance,
})

# TODO: make use of symbols with proper assumptions
_voltage_expr = _voltage_law.law.rhs.subs({
    _voltage_law.electric_field_component(_distance): _electric_field_expr,
    _voltage_law.initial_distance: outer_radius,
    _voltage_law.final_distance: inner_radius,
}).doit()

_capacitance_expr = _capacitance_law.definition.rhs.subs({
    _capacitance_law.voltage: _voltage_expr,
}).simplify()

_capacitance_from_law = law.rhs.subs({
    absolute_permittivity: quantities.vacuum_permittivity,
})

assert expr_equals(_capacitance_expr, _capacitance_from_law)


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
