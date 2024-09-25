from sympy import solve, Symbol, Eq
from symplyphysics import print_expression
from symplyphysics.laws.dynamics import friction_force_from_normal_force as friction_force
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_distance as work_friction
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy
from symplyphysics.definitions import momentum_is_mass_times_speed as momentum
from symplyphysics.laws.conservation import momentum_after_collision_equals_to_momentum_before as momentum_conservation
from symplyphysics.laws.conservation import mechanical_energy_after_equals_to_mechanical_energy_before as energy_conservation
from symplyphysics.laws.dynamics import reaction_force_from_action_force as third_newton_law
from symplyphysics.laws.dynamics import acceleration_is_force_over_mass as second_newton_law

# Example from https://uchitel.pro/%D0%B7%D0%B0%D0%B4%D0%B0%D1%87%D0%B8-%D0%BD%D0%B0-%D0%B7%D0%B0%D0%BA%D0%BE%D0%BD-%D1%81%D0%BE%D1%85%D1%80%D0%B0%D0%BD%D0%B5%D0%BD%D0%B8%D1%8F-%D0%B8%D0%BC%D0%BF%D1%83%D0%BB%D1%8C%D1%81%D0%B0/
# A skater with a mass of M = 70 kg, standing on the ice,
# throws a puck with a mass of m = 0.3 kg horizontally at a speed of v = 40 m/s.
# How far s will the skater roll back if the coefficient of friction of the skates
# on the ice is μ = 0.02?

mass_of_skater = Symbol("mass_of_skater")
mass_of_puck = Symbol("mass_of_puck")
coefficient_of_friction = Symbol("coefficient_of_friction")
velocity_of_puck = Symbol("velocity_of_puck")
distance = Symbol("distance")

gravity_acceleration = Symbol("gravity_acceleration")

momentum_of_skater = momentum.definition.subs({
    momentum.mass: mass_of_skater,
}).rhs
momentum_of_puck = momentum.definition.subs({
    momentum.mass: mass_of_puck,
    momentum.speed: velocity_of_puck
}).rhs
momentum_conservation_law = momentum_conservation.law.subs({
    momentum_conservation.momentum(momentum_conservation.time_before): momentum_of_skater,
    momentum_conservation.momentum(momentum_conservation.time_after): momentum_of_puck
})
velocity_of_skater_equation = solve(momentum_conservation_law, momentum.speed,
    dict=True)[0][momentum.speed]

acceleration_equation = second_newton_law.law.subs({
    second_newton_law.mass: mass_of_skater,
    second_newton_law.acceleration: gravity_acceleration
})
gravity_force_value = solve(acceleration_equation, second_newton_law.force,
    dict=True)[0][second_newton_law.force]
reaction_force_value = third_newton_law.law.subs({
    third_newton_law.action_force: gravity_force_value
}).rhs
# We take only the modulus of the reaction force vector,
# since we set the direction of the vertical axis coinciding with the direction of this vector
friction_force_value = friction_force.law.subs({
    friction_force.coefficient_of_friction: coefficient_of_friction,
    friction_force.normal_force: abs(reaction_force_value)
}).rhs

work_friction_value = work_friction.law.subs({
    work_friction.force: friction_force_value,
    work_friction.distance: distance
}).rhs
kinetic_energy_value = kinetic_energy.law.subs({
    kinetic_energy.mass: mass_of_skater,
    kinetic_energy.speed: velocity_of_skater_equation
}).rhs
conservation_energy = energy_conservation.law.subs({
    energy_conservation.mechanical_energy(energy_conservation.time_before): kinetic_energy_value,
    energy_conservation.mechanical_energy(energy_conservation.time_after): work_friction_value,
})
print(f"Final equation: {print_expression(conservation_energy)}")
distance_equation = solve(conservation_energy, distance, dict=True)[0][distance]
answer = Eq(distance, distance_equation)
print(f"Total distance equation: {print_expression(answer)}")

distance_m = distance_equation.subs({
    gravity_acceleration: 9.8,
    mass_of_skater: 70,
    mass_of_puck: 0.3,
    coefficient_of_friction: 0.02,
    velocity_of_puck: 40
})
print(f"Distance is: {distance_m} m")
