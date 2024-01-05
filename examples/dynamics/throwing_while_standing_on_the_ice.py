from sympy import solve
from symplyphysics import print_expression, units, Symbol, dimensionless
from symplyphysics.laws.dynamics import friction_force_from_normal_force as friction_force
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as work_friction
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy
from symplyphysics.definitions import momentum_is_mass_times_velocity as impuls
from symplyphysics.laws.conservation import momentum_after_collision_equals_to_momentum_before as impuls_conservation
from symplyphysics.laws.conservation import mechanical_energy_after_equals_to_mechanical_energy_before as energy_conservation

# Example from https://uchitel.pro/%D0%B7%D0%B0%D0%B4%D0%B0%D1%87%D0%B8-%D0%BD%D0%B0-%D0%B7%D0%B0%D0%BA%D0%BE%D0%BD-%D1%81%D0%BE%D1%85%D1%80%D0%B0%D0%BD%D0%B5%D0%BD%D0%B8%D1%8F-%D0%B8%D0%BC%D0%BF%D1%83%D0%BB%D1%8C%D1%81%D0%B0/
# A skater with a mass of M = 70 kg, standing on the ice,
# throws a puck with a mass of m = 0.3 kg horizontally at a speed of v = 40 m/s.
# How far s will the skater roll back if the coefficient of friction of the skates
# on the ice is μ = 0.02?

mass_of_skater = Symbol("mass_of_skater", units.mass)
mass_of_puck = Symbol("mass_of_puck", units.mass)
friction_factor = Symbol("friction_factor", dimensionless)
velocity_of_puck = Symbol("velocity_of_puck", units.velocity)
distance = Symbol("distance", units.length)

gravity_acceleration = Symbol("gravity_acceleration", units.acceleration)

velocity_of_skater = Symbol("velocity_of_scater", units.velocity)

impuls_of_skater = impuls.definition.subs({
    impuls.mass: mass_of_skater,
    impuls.velocity: velocity_of_skater
}).rhs
impuls_of_puck = impuls.definition.subs({
    impuls.mass: mass_of_puck,
    impuls.velocity: velocity_of_puck
}).rhs

impuls_conservation_law = impuls_conservation.law.subs({
    impuls_conservation.momentum(impuls_conservation.time_before): 0,
    impuls_conservation.momentum(impuls_conservation.time_after): impuls_of_skater - impuls_of_puck
})
velocity_of_skater_law = solve(impuls_conservation_law, velocity_of_skater, dict=True)[0][velocity_of_skater]

friction_force_law = friction_force.law.subs({
    friction_force.friction_factor: friction_factor,
    friction_force.normal_reaction: mass_of_skater * gravity_acceleration
}).rhs
work_friction_law = work_friction.law.subs({
    work_friction.force: friction_force_law,
    work_friction.distance: distance
}).rhs

kinetic_energy_law = kinetic_energy.law.subs({
    kinetic_energy.body_mass: mass_of_skater,
    kinetic_energy.body_velocity: -velocity_of_skater_law
}).rhs

conservation_energy = energy_conservation.law.subs({
    energy_conservation.mechanical_energy(energy_conservation.time_before): 0,
    energy_conservation.mechanical_energy(energy_conservation.time_after): kinetic_energy_law - work_friction_law,
})
print(f"Final equation: {print_expression(conservation_energy)}")
distance_law = solve(conservation_energy, distance, dict=True)[0][distance]
print(f"Total distance equation: {print_expression(distance_law)}")
