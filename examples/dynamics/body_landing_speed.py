from sympy import solve, Eq
from symplyphysics import print_expression, units, Symbol

from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy

# From http://kornev-school.ru/f9_law_of_energy_conservation.html first example
# A soccer ball falls without initial velocity from a height.
# Find its speed before hitting the ground.

height = Symbol("height", units.length)
gravity_acceleration = Symbol("gravity_acceleration", units.acceleration)
body_mass = Symbol("mass_of_the_body", units.mass)
landing_speed = Symbol("landing_speed", units.velocity)

E_k = kinetic_energy.law.subs({
    kinetic_energy.body_mass: body_mass,
    kinetic_energy.body_velocity: landing_speed
})
E_p = potential_energy.law.subs({
    potential_energy.body_mass: body_mass,
    potential_energy.height: height,
    potential_energy.free_fall_acceleration: gravity_acceleration
})

law = [E_p, E_k]

solved = solve(law, (landing_speed, kinetic_energy.kinetic_energy_of_body), dict=True)[0][kinetic_energy.body_velocity]
answer = Eq(landing_speed, solved)

print(f"\nSolution:\nIF potential_energy = kinetic_energy\nTHEN {print_expression(answer)}")
