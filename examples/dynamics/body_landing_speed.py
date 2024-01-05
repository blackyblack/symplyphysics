from sympy import solve, Eq
from symplyphysics import print_expression, units, Symbol, Function, Quantity

from symplyphysics.laws.dynamics import potential_energy_from_mass_and_height as potential_energy
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy
from symplyphysics.laws.conservation import mechanical_energy_after_equals_to_mechanical_energy_before as conservation
from symplyphysics.definitions import mechanical_energy_is_kinetic_and_potential as mechanical_energy

# From http://kornev-school.ru/f9_law_of_energy_conservation.html first example
# A soccer ball falls without initial velocity from a height.
# Find its speed before hitting the ground.

height = Function("height", units.length)
velocity = Function("velocity", units.velocity)

gravity_acceleration = Symbol("gravity_acceleration", units.acceleration)
body_mass = Symbol("mass_of_the_body", units.mass)

time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)

landing_speed = velocity(time_after)
velocity_start = velocity(time_before)

height_start = height(time_before)
height_land = height(time_after)


def get_mechanical_energy_at_the_moment(time: Symbol) -> Quantity:
    kinetic_energy_of_body_at_the_moment = kinetic_energy.law.subs({
        kinetic_energy.body_velocity: velocity(time),
        kinetic_energy.body_mass: body_mass
    })

    potential_energy_of_body_at_the_moment = potential_energy.law.subs({
        potential_energy.free_fall_acceleration: gravity_acceleration,
        potential_energy.height: height(time),
        potential_energy.body_mass: body_mass
    })

    mechanical_energy_of_body_at_the_moment = mechanical_energy.definition.subs({
        mechanical_energy.kinetic_energy: kinetic_energy_of_body_at_the_moment,
        mechanical_energy.potential_energy: potential_energy_of_body_at_the_moment
    })

    return mechanical_energy_of_body_at_the_moment


mechanical_energy_of_body_before = get_mechanical_energy_at_the_moment(time_before)
mechanical_energy_of_body_after = get_mechanical_energy_at_the_moment(time_after)

conservation_law = conservation.law.subs({
    conservation.mechanical_energy(conservation.time_before): mechanical_energy_of_body_before,
    conservation.mechanical_energy(conservation.time_after): mechanical_energy_of_body_after
})


solved = solve(conservation_law, landing_speed, dict=True)[0][landing_speed]
answer = Eq(landing_speed, solved)

print(f"\nSolution:\nIF potential_energy = kinetic_energy\nTHEN {print_expression(answer)}")
