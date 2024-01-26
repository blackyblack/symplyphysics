from sympy import Eq, solve, integrate
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import (
    acceleration_from_force as newtons_second_law,
    kinetic_energy_from_mass_and_velocity as kinetic_energy_def,
    mechanical_work_from_force_and_move as work_def,
)
from symplyphysics.laws.kinematic import (
    accelerated_velocity_from_time as velocity_law,
    constant_acceleration_movement_is_parabolic as distance_law,
    distance_from_constant_velocity,
)

# Description
## The work-energy principle states that the work done by all forces acting on a particle
## (the work of the resultant force) equals the change in the kinetic energy of the particle.

# Law: W = K_after - K_before
## W - total work on particle
## K_after, K_before - kinetic energy of particle before and after work is done, respectively

total_work = Symbol("total_work", units.energy)
kinetic_energy = Function("kinetic_energy", units.energy)
time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)

law = Eq(total_work, kinetic_energy(time_after) - kinetic_energy(time_before))


# Derive the law in case of rectilinear motion with constant total force acting on particle

total_force = Symbol("total_force", units.force)
particle_mass = Symbol("particle_mass", units.mass)

acceleration = solve(
    newtons_second_law.law,
    newtons_second_law.acceleration
)[0].subs({
    newtons_second_law.force: total_force,
    newtons_second_law.mass: particle_mass,
})

velocity = Function("velocity", units.velocity)

movement_time = solve(
    velocity_law.law,
    velocity_law.time
)[0].subs({
    velocity_law.velocity: velocity(time_after),
    velocity_law.initial_velocity: velocity(time_before),
    velocity_law.acceleration: acceleration,
})

displacement = distance_law.law.rhs.subs({
    distance_law.movement_time: movement_time,
    distance_law.constant_acceleration: acceleration,
    distance_law.initial_velocity: velocity(time_before),
})

work = work_def.law.rhs.subs({
    work_def.force: total_force,
    work_def.distance: displacement,
})

kinetic_energy_before = kinetic_energy_def.law.rhs.subs({
    kinetic_energy_def.body_mass: particle_mass,
    kinetic_energy_def.body_velocity: velocity(time_before),
})

kinetic_energy_after = kinetic_energy_def.law.rhs.subs({
    kinetic_energy_def.body_mass: particle_mass,
    kinetic_energy_def.body_velocity: velocity(time_after),
})

work_sub = solve(
    [
        Eq(total_work, work),
        Eq(kinetic_energy(time_before), kinetic_energy_before),
        Eq(kinetic_energy(time_after), kinetic_energy_after),
    ],
    (total_work, velocity(time_before), velocity(time_after)),
    dict=True,
)[0][total_work]

assert expr_equals(work_sub, law.rhs)


# Derive the law for the general case in one dimension

velocity_ = Symbol("velocity", units.velocity)
infinitesimal_time = Symbol("infinitesimal_time", units.time)
infinitesimal_velocity = Symbol("infinitesimal_velocity", units.velocity)

# dx = v*dt
infinitesimal_displacement = distance_from_constant_velocity.law.rhs.subs({
    distance_from_constant_velocity.initial_position: 0,
    distance_from_constant_velocity.constant_velocity: velocity_,
    distance_from_constant_velocity.movement_time: infinitesimal_time,
})

# a = dv/dt
acceleration_ = solve(
    velocity_law.law, velocity_law.acceleration
)[0].subs({
    velocity_law.initial_velocity: 0,
    velocity_law.velocity: infinitesimal_velocity,
    velocity_law.time: infinitesimal_time,
})

# F = m*a = m*(dv/dt)
force_ = solve(
    newtons_second_law.law, newtons_second_law.force
)[0].subs({
    newtons_second_law.mass: particle_mass,
    newtons_second_law.acceleration: acceleration_,
})

# dW = F*dx = m*(dv/dt)*v*dt = m*v*dv
infinitesimal_work = work_def.law.rhs.subs({
    work_def.force: force_,
    work_def.distance: infinitesimal_displacement,
})

# W = m*(v1**2)/2  - m*(v0**2)/2
finite_work = integrate(
    infinitesimal_work.subs(infinitesimal_velocity, 1),
    (velocity_, velocity(time_before), velocity(time_after))
)

finite_work_sub = solve(
    [
        Eq(total_work, finite_work),
        Eq(kinetic_energy(time_before), kinetic_energy_before),
        Eq(kinetic_energy(time_after), kinetic_energy_after),
    ],
    (total_work, velocity(time_before), velocity(time_after)),
    dict=True,
)[0][total_work]

assert expr_equals(finite_work_sub, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(kinetic_energy_before_=kinetic_energy, kinetic_energy_after_=kinetic_energy)
@validate_output(total_work)
def calculate_total_work(kinetic_energy_before_: Quantity, kinetic_energy_after_: Quantity) -> Quantity:
    result = law.rhs.subs({
        kinetic_energy(time_before): kinetic_energy_before_,
        kinetic_energy(time_after): kinetic_energy_after_,
    })
    return Quantity(result)
