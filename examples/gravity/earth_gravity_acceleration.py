#!/usr/bin/env python3
r"""
In this example, we calculate the free fall acceleration on Earth's surface using the law of
gravitational force and Newton's second law of motion and show that the free fall acceleration is
independent of the probe mass.

Earth's radius is :math:`6371 \, \text{km}`; Earth's mass is :math:`5.9722 \cdot 10^{24} \, \text{kg}`.

**Notes**:

#. *Free fall acceleration* is also referred to as *acceleration due to gravity*.
"""

from sympy import solve
from symplyphysics import print_expression, units, convert_to_si, Quantity, clone_as_symbol, symbols
from symplyphysics.laws.gravity import gravity_force_from_mass_and_distance as gravity_force_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as newtons_second_law

planet_mass = clone_as_symbol(symbols.mass, display_symbol="M")
probe_mass = clone_as_symbol(symbols.mass)

gravity_force_expr = gravity_force_law.law.rhs.subs({
    gravity_force_law.first_mass: planet_mass,
    gravity_force_law.second_mass: probe_mass,
})

acceleration_expr = solve(
    newtons_second_law.law,
    newtons_second_law.acceleration,
    dict=True,
)[0][newtons_second_law.acceleration]

free_fall_acceleration_expr = acceleration_expr.subs({
    newtons_second_law.force: gravity_force_expr,
    newtons_second_law.mass: probe_mass,
})

assert probe_mass not in free_fall_acceleration_expr.free_symbols

print(
    "Free fall acceleration formula:",
    print_expression(free_fall_acceleration_expr),
    sep="\n",
    end="\n\n",
)

earth_mass = Quantity(5.9722e24 * units.kilogram)
earth_radius = Quantity(6371 * units.kilometer)

earth_acceleration_expr = free_fall_acceleration_expr.subs({
    planet_mass: earth_mass,
    gravity_force_law.distance_between_mass_centers: earth_radius
})
earth_acceleration_qty = Quantity(earth_acceleration_expr)

earth_acceleration_float = convert_to_si(earth_acceleration_qty).evalf(4)
print(f"Free fall acceleration on Earth is {earth_acceleration_float} m/s^2")
