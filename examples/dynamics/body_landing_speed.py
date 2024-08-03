from sympy import solve, Eq, Symbol
from symplyphysics import print_expression, units
from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as kinetic_energy
from symplyphysics.laws.conservation import mechanical_energy_after_equals_to_mechanical_energy_before as conservation
from symplyphysics.definitions import mechanical_energy_is_kinetic_and_potential_energy as mechanical_energy

# From http://kornev-school.ru/f9_law_of_energy_conservation.html first example
# A soccer ball falls without initial velocity from a height.
# Find its speed before hitting the ground.

body_height = Symbol("body_height")
gravity_acceleration = Symbol("gravity_acceleration")
body_mass = Symbol("mass_of_the_body")
landing_speed = Symbol("landing_speed")

kinetic_energy_before = kinetic_energy.law.subs({
    kinetic_energy.speed: 0,
    kinetic_energy.mass: body_mass
}).rhs
potential_energy_before = potential_energy.law.subs({
    potential_energy.height: body_height,
    potential_energy.mass: body_mass
}).rhs
mechanical_energy_before = mechanical_energy.definition.subs({
    mechanical_energy.kinetic_energy: kinetic_energy_before,
    mechanical_energy.potential_energy: potential_energy_before,
}).rhs
kinetic_energy_after = kinetic_energy.law.subs({
    kinetic_energy.speed: landing_speed,
    kinetic_energy.mass: body_mass
}).rhs
potential_energy_after = potential_energy.law.subs({
    potential_energy.height: 0,
    potential_energy.mass: body_mass
}).rhs
mechanical_energy_after = mechanical_energy.definition.subs({
    mechanical_energy.kinetic_energy: kinetic_energy_after,
    mechanical_energy.potential_energy: potential_energy_after,
}).rhs

conservation_law = conservation.law.subs({
    conservation.mechanical_energy(conservation.time_before): mechanical_energy_before,
    conservation.mechanical_energy(conservation.time_after): mechanical_energy_after
})

# First solution is negative - ignore it
solved = solve(conservation_law, landing_speed, dict=True)[1][landing_speed]
answer = Eq(landing_speed, solved)

print(f"\nSolution:\nIF potential_energy = kinetic_energy\nTHEN {print_expression(answer)}")

landing_speed_ms = solved.subs({
    units.acceleration_due_to_gravity: 9.8,
    body_height: 11.25
}).evalf(3)

print(f"\nLanding speed when falling from 11.25 meters is: {landing_speed_ms} m/s")
