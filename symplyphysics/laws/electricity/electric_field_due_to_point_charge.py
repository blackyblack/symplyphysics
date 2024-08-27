r"""
Electric field due to point charge
==================================

The value of the electric field set up by a point charge is linearly proportional
to the value of the charge and the square inverse of the distance to it.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import electrostatic_force_via_charges_and_distance as coulombs_law
from symplyphysics.laws.electricity import electric_field_value_is_force_over_test_charge as electric_field_def

electric_field = Symbol("electric_field", units.force / units.charge)
r"""
Value of the electric field due to point charge :math:`q`.

Symbol:
    :code:`E`
"""

charge = Symbol("charge", units.charge)
"""
Value of the point charge.

Symbol:
    :code:`q`
"""

distance = Symbol("distance", units.length)
"""
Distance to the charge.

Symbol:
    :code:`r`
"""

law = Eq(electric_field, 1 / (4 * pi * units.vacuum_permittivity) * charge / distance**2)
r"""
:code:`E = 1 / (4 * pi * epsilon_0) * q / r^2`

Latex:
    .. math::
        E = \frac{1}{4 \pi \varepsilon_0} \frac{q}{r^2}
"""

# Derive this law from Coulomb's law and the definition of electric field

_test_charge = Symbol("_test_charge", units.charge)

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
_electric_field_derived = solve(_electric_field_def_sub, electric_field_def.electric_field)[0]
_electric_field_from_law = law.rhs

assert expr_equals(_electric_field_from_law, _electric_field_derived)


@validate_input(point_charge_=charge, distance_=distance)
@validate_output(electric_field)
def calculate_electric_field(point_charge_: Quantity, distance_: Quantity) -> Quantity:
    result = solve(law, electric_field)[0]
    result_field = result.subs({
        charge: point_charge_,
        distance: distance_,
    })
    return Quantity(result_field)
