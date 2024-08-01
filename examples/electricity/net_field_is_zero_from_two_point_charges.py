from sympy import solve, Eq, symbols
from symplyphysics import units, convert_to, Quantity
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as point_charge

# Description
## Two charged particles are fixed to an x-axis. Particle 1 of charge q1 = 2.1*10^-8 C is at position x = 20 cm;
## particle 2 of charge q2 = -4.0*q1 is at position x = 70 cm. Find the coordinate on the axis other than the
## infinity at which the net electric field produced by the two particles equals zero?

first_charge, second_charge = symbols("first_charge second_charge")
first_position, second_position, zero_position = symbols(
    "first_position second_position zero_position")

law_electric_field_first_charge = point_charge.law.subs({
    point_charge.point_charge: first_charge,
    point_charge.distance: zero_position - first_position,
})
electric_field_first_charge = solve(law_electric_field_first_charge, point_charge.electric_field)[0]

law_electric_field_second_charge = point_charge.law.subs({
    point_charge.point_charge: second_charge,
    point_charge.distance: zero_position - second_position,
})
electric_field_second_charge = solve(law_electric_field_second_charge,
    point_charge.electric_field)[0]

# By condition the total field is zero.
eqn = Eq(electric_field_first_charge + electric_field_second_charge, 0)

# This equation has 2 solutions
solved = solve(eqn, zero_position)

FIRST_CHARGE_VALUE = 2.1e-8
SECOND_CHARGE_VALUE = -4.0 * FIRST_CHARGE_VALUE
values = {
    first_charge: Quantity(FIRST_CHARGE_VALUE * units.coulomb),
    second_charge: Quantity(SECOND_CHARGE_VALUE * units.coulomb),
    first_position: Quantity(20.0 * units.centimeter),
    second_position: Quantity(70.0 * units.centimeter)
}

for name, solution in zip(["first", "second"], solved):
    quantity = Quantity(solution.subs(values))
    value = convert_to(quantity, units.centimeter).evalf(3)
    print(f"The {name} coordinate is at position x = {value} {units.centimeter}")
