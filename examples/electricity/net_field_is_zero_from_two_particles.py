from sympy import solve, Eq, symbols, sqrt
from symplyphysics import print_expression, units, convert_to, Quantity, prefixes
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as point_charge

# Description
## Two charged particles are fixed to an x-axis. Particle 1 of charge q1 = 2.1*10^-8 C is at position x = 20 cm;
## particle 2 of charge q2 = -4.0*q1 is at position x = 70 cm. At what coordinate on the axis other than the
## infinity is the net electric field produces by the two particles equal to zero?

first_charge, second_charge = symbols("first_charge second_charge")
first_position, second_position, zero_position = symbols("first_position second_position zero_position")

law_electric_field_first_charge = point_charge.law.subs({
    point_charge.point_charge: first_charge,
    point_charge.distance: zero_position - first_position,
})
electric_field_first_charge = solve(law_electric_field_first_charge, point_charge.electric_field)[0]

law_electric_field_second_charge = point_charge.law.subs({
    point_charge.point_charge: second_charge,
    point_charge.distance: zero_position - second_position,
})
electric_field_second_charge = solve(law_electric_field_second_charge, point_charge.electric_field)[0]

# By condition the total field is zero.
# Note that the point-charge law for electric field employs the magnitude of the field.
# But here we take a projection of the field to the axis, and since the two charges are
# oppositely charged, one of the two will be negative and the other positive. This explains
# the minus sign at the second term.
eqn = Eq(electric_field_first_charge - electric_field_second_charge, 0)
# This equation has 2 solutions
solved = solve(eqn, zero_position)

values = {
    first_charge: Quantity(2.1e-8 * units.coulomb),
    first_position: Quantity(20.0 * units.centimeter),
    second_position: Quantity(70.0 * units.centimeter)
}
values[second_charge] = -4.0 * values[first_charge]

for name, solution in zip(["first", "second"], solved):
    quantity = Quantity(solution.subs(values))
    value = convert_to(quantity, units.centimeter).evalf(3)
    print(f"The {name} coordinate is at position x = {value} {units.centimeter}")
