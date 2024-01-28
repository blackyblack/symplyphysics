from sympy import Derivative, Eq, solve, integrate
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
from symplyphysics.definitions import velocity_is_movement_derivative as velocity_def
from symplyphysics.definitions import acceleration_is_velocity_derivative as acceleration_def

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

# Derive the law in case of rectilinear motion with infinitesimal constant total force acting on particle and
# integrating the resulting infinitesimal work over time

total_force = Symbol("total_force", units.force)
particle_mass = Symbol("particle_mass", units.mass)
infinitesimal_time = Symbol("infinitesimal_time", units.time)
infinitesimal_velocity = Function("infinitesimal_velocity", units.velocity)

infinitesimal_displacement = solve(
    velocity_def.definition,
    Derivative(velocity_def.movement(velocity_def.moving_time), velocity_def.moving_time))[0]
infinitesimal_displacement = infinitesimal_displacement.subs(velocity_def.moving_time,
    infinitesimal_time)
infinitesimal_displacement = infinitesimal_displacement.subs(
    velocity_def.velocity(infinitesimal_time), infinitesimal_velocity(infinitesimal_time))
infinitesimal_acceleration = solve(acceleration_def.definition,
    acceleration_def.acceleration(acceleration_def.time))[0]
infinitesimal_acceleration = infinitesimal_acceleration.subs(acceleration_def.time,
    infinitesimal_time)
infinitesimal_acceleration = infinitesimal_acceleration.subs(
    acceleration_def.velocity(infinitesimal_time), infinitesimal_velocity(infinitesimal_time))

# F = m*a = m*(dv/dt)
force_ = solve(newtons_second_law.law, newtons_second_law.force)[0].subs({
    newtons_second_law.mass: particle_mass,
    newtons_second_law.acceleration: infinitesimal_acceleration,
})

# dW = F*dx = m*(dv/dt)*v*dt = m*v*dv
infinitesimal_work = work_def.law.rhs.subs({
    work_def.force: force_,
    work_def.distance: infinitesimal_displacement,
})

# W = m*(v1**2)/2  - m*(v0**2)/2
finite_work = integrate(infinitesimal_work, (infinitesimal_time, time_before, time_after))

kinetic_energy_before_eq = kinetic_energy_def.law.subs({
    kinetic_energy_def.body_mass: particle_mass,
    kinetic_energy_def.body_velocity: infinitesimal_velocity(time_before),
    kinetic_energy_def.kinetic_energy_of_body: kinetic_energy(time_before)
})

kinetic_energy_after_eq = kinetic_energy_def.law.subs({
    kinetic_energy_def.body_mass: particle_mass,
    kinetic_energy_def.body_velocity: infinitesimal_velocity(time_after),
    kinetic_energy_def.kinetic_energy_of_body: kinetic_energy(time_after)
})

finite_work_sub = solve(
    [
    Eq(total_work, finite_work),
    kinetic_energy_before_eq,
    kinetic_energy_after_eq,
    ],
    (total_work, infinitesimal_velocity(time_before), infinitesimal_velocity(time_after)),
    dict=True,
)[0][total_work]

assert expr_equals(finite_work_sub, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(kinetic_energy_before_=kinetic_energy, kinetic_energy_after_=kinetic_energy)
@validate_output(total_work)
def calculate_total_work(kinetic_energy_before_: Quantity,
    kinetic_energy_after_: Quantity) -> Quantity:
    result = law.rhs.subs({
        kinetic_energy(time_before): kinetic_energy_before_,
        kinetic_energy(time_after): kinetic_energy_after_,
    })
    return Quantity(result)
