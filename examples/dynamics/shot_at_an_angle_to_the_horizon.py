#!/usr/bin/env python3

import math
from sympy import solve, Symbol, Eq
from symplyphysics import print_expression
from symplyphysics.laws.geometry import planar_projection_is_cosine as projection_velocity
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy

# The initial velocity of the bullet is 600 m/s, its mass is 10 g.
# At what angle to the horizon did it fly out of the muzzle of the gun,
# if its kinetic energy at the highest point of the trajectory is 450 J?

start_velocity = Symbol("start_velocity")
kinetic_energy_in_peak = Symbol("kinetic_energy_in_peak")
mass_of_bullet = Symbol("mass_of_bullet")

angle_of_shot = Symbol("angle_of_shot")

# We take only the horizontal component of the velocity,
# since there is no vertical component of the velocity at the highest point of the trajectory.
velocity_projection_equation = projection_velocity.law.subs(({
    projection_velocity.vector_length: start_velocity,
    projection_velocity.vector_angle: angle_of_shot
})).rhs

kinetic_energy_equation = kinetic_energy.law.subs({
    kinetic_energy.kinetic_energy_of_body: kinetic_energy_in_peak,
    kinetic_energy.body_velocity: velocity_projection_equation,
    kinetic_energy.symbols.basic.mass: mass_of_bullet
})
print(f"Final equation: {print_expression(kinetic_energy_equation)}")

# The last solve is correct for this example
angle_of_shot_equation = solve(kinetic_energy_equation, angle_of_shot, dict=True)[-1][angle_of_shot]
answer = Eq(angle_of_shot, angle_of_shot_equation)
print(f"Total angle of shot equation:\n{print_expression(answer)}")

angle_of_shot_rad = angle_of_shot_equation.subs({
    start_velocity: 600,
    mass_of_bullet: 0.01,
    kinetic_energy_in_peak: 450
})
# 1 rad = 180 / pi degree
print(f"Angle of shot is: {angle_of_shot_rad * 180 / math.pi} degree")
