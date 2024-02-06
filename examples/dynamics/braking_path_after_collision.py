#!/usr/bin/env python3

from sympy import solve, Symbol, Eq
from symplyphysics import Quantity, convert_to, units, print_expression
from symplyphysics.definitions import (
    momentum_is_mass_times_velocity as momentum_def,
)
from symplyphysics.laws.conservation import (
    mechanical_energy_after_equals_to_mechanical_energy_before as energy_conservation_law,
    momentum_after_collision_equals_to_momentum_before as momentum_conservation_law
)
from symplyphysics.laws.dynamics import (
    kinetic_energy_from_mass_and_velocity as kinetic_energy_def,
    potential_energy_from_mass_and_height as potential_energy_def,
    friction_force_from_normal_force as friction_law,
    braking_path as braking_path_law,
)

# Description
## Block 1 of mass m1 slides from rest along a frictionless ramp from height h = 2.50 m and then
## collides with stationary block 2, which has mass m2 = 2.00*m1. After the collision, block 2
## slides into a region where the coefficient of kinetic friction is mu = 0.500 and comes to a
## stop in distance d within that region. What is the value of distance d if the collision was
## (a) elastic and (b) completely inelastic?

mass_1 = Symbol("mass_1")
mass_2 = 2 * mass_1
height = Symbol("height")
friction_coefficient = Symbol("friction_coefficient")

values = {
    height: Quantity(2.50 * units.meter),
    friction_coefficient: Quantity(0.500),
    units.acceleration_due_to_gravity: Quantity(9.81 * units.meter / units.second**2),
}

speed_before_1 = Symbol("speed_before_1")
speed_after_1 = Symbol("speed_after_1")
speed_after_2 = Symbol("speed_after_2")

# Find speed of block 1 just before the collision

energy_before_1 = potential_energy_def.law.rhs.subs({
    potential_energy_def.height: height,
    potential_energy_def.body_mass: mass_1,
    potential_energy_def.free_fall_acceleration: units.acceleration_due_to_gravity,
})

energy_after_1 = kinetic_energy_def.law.rhs.subs({
    kinetic_energy_def.body_mass: mass_1,
    kinetic_energy_def.body_velocity: speed_before_1,
})

energy_conservation_eqn_1 = energy_conservation_law.law.subs({
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_before): energy_before_1,
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_after): energy_after_1,
})

speed_before_1 = solve(energy_conservation_eqn_1, speed_before_1)[0]

momentum_before_1 = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass_1,
    momentum_def.velocity: speed_before_1,
})

# (a) Elastic collision

momentum_after_1_elastic = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass_1,
    momentum_def.velocity: speed_after_1,
})

momentum_after_2_elastic = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass_2,
    momentum_def.velocity: speed_after_2,
})

momentum_after_elastic = momentum_after_1_elastic + momentum_after_2_elastic

momentum_conservation_eqn_collision = momentum_conservation_law.law.subs({
    momentum_conservation_law.momentum(momentum_conservation_law.time_before): momentum_before_1,
    momentum_conservation_law.momentum(momentum_conservation_law.time_after): momentum_after_elastic,
})

energy_after_1_elastic = kinetic_energy_def.law.rhs.subs({
    kinetic_energy_def.body_mass: mass_1,
    kinetic_energy_def.body_velocity: speed_after_1,
})

energy_after_2_elastic = kinetic_energy_def.law.rhs.subs({
    kinetic_energy_def.body_mass: mass_2,
    kinetic_energy_def.body_velocity: speed_after_2,
})

energy_after_elastic = energy_after_1_elastic + energy_after_2_elastic

energy_conservation_eqn_collision = energy_conservation_law.law.subs({
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_before): energy_before_1,
    energy_conservation_law.mechanical_energy(energy_conservation_law.time_after): energy_after_elastic,
})

speed_after_2_elastic = solve(
    [momentum_conservation_eqn_collision, energy_conservation_eqn_collision],
    (speed_after_1, speed_after_2),
    dict=True
)[1][speed_after_2]

# (b) Inelastic collision

speed_after = Symbol("speed_after")

momentum_after_inelastic = momentum_def.definition.rhs.subs({
    momentum_def.mass: mass_1 + mass_2,
    momentum_def.velocity: speed_after,
})

speed_after_inelastic = solve(
    Eq(momentum_before_1, momentum_after_inelastic),
    speed_after
)[0]

# Find braking path from friction force in general case

mass_after = Symbol("mass_after")

normal_reaction = mass_after * units.acceleration_due_to_gravity

friction_force = friction_law.law.rhs.subs({
    friction_law.friction_factor: friction_coefficient,
    friction_law.normal_reaction: normal_reaction,
})

braking_path = braking_path_law.law.rhs.subs({
    braking_path_law.mass: mass_after,
    braking_path_law.velocity: speed_after,
    braking_path_law.friction_force: friction_force,
})

braking_path_elastic = braking_path.subs({
    mass_after: mass_2,
    speed_after: speed_after_2_elastic,
}).simplify()

braking_path_inelastic = braking_path.subs({
    mass_after: mass_1 + mass_2,
    speed_after: speed_after_inelastic,
}).simplify()

print(f"General expression for braking path:\n{print_expression(braking_path)}\n")
print(f"Braking path after elastic collision:\n{print_expression(braking_path_elastic)}\n")
print(f"Braking path after completely inelastic collision:\n{print_expression(braking_path_inelastic)}\n")

braking_path_elastic_value = convert_to(
    Quantity(braking_path_elastic.subs(values)),
    units.meter,
).evalf(3)

braking_path_inelastic_value = convert_to(
    Quantity(braking_path_inelastic.subs(values)),
    units.meter,
).evalf(3)

print(f"Braking path after elastic collision is {braking_path_elastic_value} m.")
print(f"Braking path after completely inelastic collision is {braking_path_inelastic_value} m.")
