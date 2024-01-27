from sympy import Eq, solve, series
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as point_field
from symplyphysics.laws.electricity import electric_dipole_moment

# Description
## The value of the electric field set up by a dipole at a distant point on the dipole axis
## (which runs through both particles) is proportional to the inverse cube of the distance
## to the dipole and the value of the dipole moment.

# Law:
## E = 2 * k_e * p / z^3
## E - value of electric field of dipole
## k_e - Coulomb's constant
## p = ql - value of the dipole moment
## z - distance to dipole on the dipole axis

# Condition
## z/l >> 1 - the point of measuring the electric field should be far enough from the dipole itself

electric_field = Symbol("electric_field", units.force / units.charge)
dipole_moment = Symbol("dipole_moment", units.charge * units.length)
distance_to_dipole = Symbol("distance_to_dipole", units.length)

law = Eq(electric_field, 2 * units.coulomb_constant * dipole_moment / distance_to_dipole**3)

# Derive the law from the expression for the electric field of point charges.
# Assuming the dipole is made up of two point charges q and -q (q > 0) with distance d in between.
# Assume the z-axis running through both charges, and let its origin be at their middle point,
# the negative charge located below and the positive one above the origin.

charge = Symbol("charge", units.charge)
distance_between_charges = Symbol("distance_between_charges", units.length)
distance_to_origin = Symbol("distance_to_origin", units.length)

positive_charge_field = point_field.law.rhs.subs({
    point_field.point_charge: charge,
    point_field.distance: distance_to_origin - distance_between_charges / 2
})
negative_charge_field = point_field.law.rhs.subs({
    point_field.point_charge: -charge,
    point_field.distance: distance_to_origin + distance_between_charges / 2
})

# The net electric field is the sum of electric fields due to both charges
net_field = positive_charge_field + negative_charge_field

# The condition that distance_to_origin/distance_between_charges >> 1 can be analyzed as such:
# Let distance_between_charges = factor * distance_to_origin, where factor -> 0.
factor = Symbol("factor", dimensionless)

# Use the above definition of factor to substitute distance_to_origin in the net_field formula
net_field_sub = net_field.subs(distance_between_charges, factor * distance_to_origin)

# Expand net_field_sub with respect to factor around 0 up to the first power
# in order to find the first approximation of the current formula.
net_field_approx = series(net_field_sub, factor, 0, 2).removeO()

# Substitute distance_to_origin back into the net_field formula
net_field_approx_sub = net_field_approx.subs(factor, distance_between_charges / distance_to_origin)

# Replace charge*distance_between_charges back with dipole_moment
net_field_derived = solve(
    [
    Eq(electric_field, net_field_approx_sub),
    electric_dipole_moment.law.subs({
    electric_dipole_moment.electric_moment: dipole_moment,
    electric_dipole_moment.charge: charge,
    electric_dipole_moment.distance: distance_between_charges,
    }),
    ],
    (charge, electric_field),
    dict=True,
)[0][electric_field]

net_field_from_law = law.rhs.subs(distance_to_dipole, distance_to_origin)

assert expr_equals(net_field_from_law, net_field_derived)


def print_law() -> str:
    return print_expression(law)


@validate_input(dipole_moment_=dipole_moment, distance_to_dipole_=distance_to_dipole)
@validate_output(electric_field)
def calculate_electric_field(dipole_moment_: Quantity, distance_to_dipole_: Quantity) -> Quantity:
    result = solve(law, electric_field)[0]
    result_field = result.subs({
        dipole_moment: dipole_moment_,
        distance_to_dipole: distance_to_dipole_,
    })
    return Quantity(result_field)
