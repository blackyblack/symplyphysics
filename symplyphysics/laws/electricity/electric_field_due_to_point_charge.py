r"""
Electric field due to point charge
==================================

The strength of the electric field set up by a point charge is linearly proportional
to the value of the charge and the square inverse of the distance to it.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import electrostatic_force_via_charges_and_distance as coulombs_law
from symplyphysics.laws.electricity import electric_field_is_force_over_test_charge as electric_field_def

electric_field_strength = symbols.electric_field_strength
r"""
:symbols:`electric_field_strength` due to point charge :math:`q`.
"""

charge = symbols.charge
"""
Value of the point :symbols:`charge`.
"""

distance = symbols.distance
"""
:symbols:`distance` to the charge.
"""

law = Eq(electric_field_strength, 1 / (4 * pi * quantities.vacuum_permittivity) * charge / distance**2)
"""
:laws:symbol::

:laws:latex::
"""

# Derive this law from Coulomb's law and the definition of electric field

_test_charge = clone_as_symbol(symbols.charge, display_symbol="test_charge")

_coulombs_law_sub = coulombs_law.law.subs({
    coulombs_law.first_charge: charge,
    coulombs_law.second_charge: _test_charge,
    coulombs_law.distance: distance,
})
_force = solve(_coulombs_law_sub, coulombs_law.electrostatic_force)[0]

_electric_field_def_sub = electric_field_def.law.subs({
    electric_field_def.electrostatic_force: _force,
    electric_field_def.test_charge: _test_charge,
})
_electric_field_derived = solve(_electric_field_def_sub, electric_field_def.electric_field_strength)[0]
_electric_field_from_law = law.rhs

assert expr_equals(_electric_field_from_law, _electric_field_derived)


@validate_input(point_charge_=charge, distance_=distance)
@validate_output(electric_field_strength)
def calculate_electric_field(point_charge_: Quantity, distance_: Quantity) -> Quantity:
    result = solve(law, electric_field_strength)[0]
    result_field = result.subs({
        charge: point_charge_,
        distance: distance_,
    })
    return Quantity(result_field)
