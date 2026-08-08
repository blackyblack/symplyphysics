"""
Electrostatic potential due to point charge
===========================================

Electrostatic potential of electric field due to a point charge is inversely proportional
to the distance to the point charge. Also see :ref:`Electrostatic potential is work to bring from reference point over charge`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Electric_potential#Electric_potential_due_to_a_point_charge>`__.
"""

from sympy import (Eq, solve, pi, oo)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.electromagnetism.electrostatics.point_charge import electric_field_due_to_point_charge as field_law
from symplyphysics.electromagnetism.electrostatics.voltage import voltage_is_line_integral_of_electric_field as voltage_integral_law

electrostatic_potential = symbols.electric_potential
"""
Electrostatic potential at given point. See :symbols:`electric_potential`.
"""

absolute_permittivity = symbols.absolute_permittivity
"""
:symbols:`absolute_permittivity` of the medium.
"""

distance = symbols.euclidean_distance
"""
:symbols:`euclidean_distance` to the point charge.
"""

charge = symbols.charge
"""
Electric :symbols:`charge`.
"""

law = Eq(electrostatic_potential, charge / (4 * pi * absolute_permittivity * distance))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the electric field due to a point charge and the line integral of the
# electric field.

## The integration path is radial, so the tangent component of the electric field along the
## path is the field strength of the point charge at the current distance.
_field_expr = field_law.law.rhs.subs({
    field_law.charge: charge,
    field_law.distance: voltage_integral_law.distance,
})

## The electrostatic potential at a point is the voltage between the reference point and the
## given point. Assumption: the reference point is located at infinity, where the potential
## of the point charge vanishes.
_potential_derived = voltage_integral_law.law.rhs.subs(
    voltage_integral_law.electric_field_component(voltage_integral_law.distance),
    _field_expr,
).subs({
    voltage_integral_law.initial_distance: oo,
    voltage_integral_law.final_distance: distance,
}).doit()

## The electric field law is stated for the vacuum. In a linear dielectric medium the vacuum
## permittivity is replaced by the absolute permittivity of the medium.
## TODO: use a law for the electric field of a point charge in a dielectric medium when it
## appears in the repository.
assert expr_equals(_potential_derived,
    law.rhs.subs(absolute_permittivity, quantities.vacuum_permittivity))


@validate_input(absolute_permittivity_=absolute_permittivity, distance_=distance, charge_=charge)
@validate_output(electrostatic_potential)
def calculate_electrostatic_potential(absolute_permittivity_: Quantity, distance_: Quantity,
    charge_: Quantity) -> Quantity:
    result_expr = solve(law, electrostatic_potential, dict=True)[0][electrostatic_potential]
    result_expr = result_expr.subs({
        absolute_permittivity: absolute_permittivity_,
        distance: distance_,
        charge: charge_,
    })
    return Quantity(result_expr)
