from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import force_from_charge_and_distance as coulombs_law
from symplyphysics.laws.electricity import electric_field_value_is_force_over_test_charge as electric_field_def

# Description
## The value of the electric field set up by a point charge is linearly proportional
## to the value of the charge and the square inverse of the distance to it.

# Law: E = k_e * q / r^2
## E - value of the electric field, which has the same sign as the charge producing it
## k_e - Coulomb's constant
## q - point charge
## r - distance to point charge

electric_field = Symbol("electric_field", units.force / units.charge)
point_charge = Symbol("point_charge", units.charge)
distance = Symbol("distance", units.length)

law = Eq(electric_field, units.coulomb_constant * point_charge / distance**2)

# Derive this law from Coulomb's law and the definition of electric field

test_charge = Symbol("test_charge", units.charge)

coulombs_law_sub = coulombs_law.law.subs({
    coulombs_law.first_charge: point_charge,
    coulombs_law.second_charge: test_charge,
    coulombs_law.distance: distance,
})
force = solve(coulombs_law_sub, coulombs_law.electrostatic_force)[0]

electric_field_def_sub = electric_field_def.law.subs({
    electric_field_def.electrostatic_force: force,
    electric_field_def.test_charge: test_charge,
})
electric_field_derived = solve(electric_field_def_sub, electric_field_def.electric_field)[0]
electric_field_from_law = law.rhs

assert expr_equals(electric_field_from_law, electric_field_derived)


def print_law() -> str:
    return print_expression(law)


@validate_input(point_charge_=point_charge, distance_=distance)
@validate_output(electric_field)
def calculate_electric_field(point_charge_: Quantity, distance_: Quantity) -> Quantity:
    result = solve(law, electric_field)[0]
    result_field = result.subs({
        point_charge: point_charge_,
        distance: distance_,
    })
    return Quantity(result_field)
