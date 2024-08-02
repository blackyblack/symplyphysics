#!/usr/bin/env python3

from sympy import solve
from symplyphysics import units, convert_to, Quantity, print_expression
from symplyphysics.laws.conservation import momentum_after_collision_equals_to_momentum_before as momentum_conservation_law
from symplyphysics.definitions import momentum_is_mass_times_speed as momentum_def

# Inelastic collision transforms two objects with masses m1, m2 and velocities v1 and v2 to one object with mass m = m1 + m2 and unknown velocity to be resolved
# Assert 1-dimentional environment.
# Let the first object be the bullet with mass of 10gram and velocity of 700m/s hits the second object - steady 75kg piece of plasticine, stucks in it and makes this plasticine move with some unknown velocity.
# The space is 1-dimensional, all vectors are collinear with axis.

bullet_mass = Quantity(10 * units.gram)
bullet_velocity = Quantity(700 * units.meter / units.second)
body_mass = Quantity(75 * units.kilogram)
# Hint dimension so 'convert_to' is able to convert it to m/s
body_velocity = Quantity(0 * units.meter / units.second, dimension=units.velocity)

print(
    f"Formula for momentum conservation law is:\n{print_expression(momentum_conservation_law.law)}")
print(f"Formula for momentum is:\n{print_expression(momentum_def.definition)}")

# initial body velocity is 0 so the momentum is 0 as well. So the momentum of system before collision equals to momentum of bullet
momentum_before = momentum_def.calculate_momentum(bullet_mass, bullet_velocity)
momentum_after = momentum_conservation_law.calculate_momentum_after(momentum_before)

# Mass of the resulting object is a sum of masses of a bullet and plasticine
solved = solve(momentum_def.definition, momentum_def.speed, dict=True)[0][momentum_def.speed]
result_expr = solved.subs({
    momentum_def.mass: Quantity(bullet_mass + body_mass),
    momentum_def.momentum: momentum_after
})
result_velocity = Quantity(result_expr)

result_velocity_meter_per_second = convert_to(result_velocity, units.meter / units.second).evalf(2)
bullet_mass_gram = convert_to(bullet_mass, units.gram).evalf(2)
body_mass_kg = convert_to(body_mass, units.kilogram).evalf(2)
bullet_velocity_meter_per_second = convert_to(bullet_velocity, units.meter / units.second).evalf(3)
body_velocity_meter_per_second = convert_to(body_velocity, units.meter / units.second).evalf(2)

print(f"Velocity = {result_velocity_meter_per_second} {units.meter / units.second}\n")
print(
    f"for two inelasticly collided objects with masses = {bullet_mass_gram} {units.gram}, {body_mass_kg} {units.kilogram}\n"
)
print(
    f"velocities = {bullet_velocity_meter_per_second} {units.meter / units.second}, {body_velocity_meter_per_second} {units.meter / units.second}\n"
)
