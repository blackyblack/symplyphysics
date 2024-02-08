from sympy import (
    Eq,
    pi,
    sqrt,
    solve,
    Symbol as SymSymbol,
    Function as SymFunction,
)
from sympy.physics.units import acceleration_due_to_gravity
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.definitions import (
    angular_velocity_is_angle_derivative as angular_velocity_def,
    angular_acceleration_is_angular_velocity_derivative as angular_acceleration_def,
)
from symplyphysics.laws.dynamics import (
    acceleration_from_force as newtons_second_law,
    torque_due_to_twisting_force as torque_def,
    moment_of_force_from_moment_of_inertia_and_angular_acceleration as torque_law,
)

# Description
## A physical pendulum is a pendulum with an arbitrary distribution of mass that
## oscillates about a given pivot point.

# Law: T = 2 * pi * sqrt(I / (m * g * h))
## T - period of physical pendulum
## I - rotational inertia of pendulum
## m - mass of pendulum
## g - acceleration due to gravity
## h - distance between pivot and pendulum's center of mass

oscillation_period = Symbol("oscillation_period", units.time)
rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
pendulum_mass = Symbol("pendulum_mass", units.mass)
distance_to_pivot = Symbol("distance_to_pivot", units.length)

law = Eq(
    oscillation_period, 
    2 * pi * sqrt(rotational_inertia / (pendulum_mass * acceleration_due_to_gravity * distance_to_pivot))
)


# Derive from torque definition

time = SymSymbol("time")
angle = SymFunction("angle")
torque = SymSymbol("torque")

gravitational_force = solve(newtons_second_law.law, newtons_second_law.force)[0].subs({
    newtons_second_law.mass: pendulum_mass,
    newtons_second_law.acceleration: acceleration_due_to_gravity,
})

angular_velocity = (
    angular_velocity_def.definition.rhs
    .subs(angular_velocity_def.time, time)
    .subs(angular_velocity_def.angle_function(time), angle(time))
)

angular_acceleration = (
    angular_acceleration_def.definition.rhs
    .subs(angular_acceleration_def.time, time)
    .subs(angular_acceleration_def.angular_velocity(time), angular_velocity)
)

angle_sym = SymSymbol("angle")

# The factor of -1 indicates that the torque acts to reduce the angle.
torque_from_def = -1 * torque_def.law.rhs.subs({
    torque_def.force: gravitational_force,
    torque_def.distance_to_axis: distance_to_pivot,
    torque_def.angle: angle(time),
})

torque_from_def = (
    torque_from_def
    .subs(angle(time), angle_sym)
    .series(angle_sym, 0, 2)
    .removeO()
    .subs(angle_sym, angle(time))
)

torque_from_law = torque_law.law.rhs.subs({
    torque_law.moment_of_inertia: rotational_inertia,
    torque_law.angular_acceleration: angular_acceleration,
})

diff_eqn = torque_from_def - torque_from_law


def print_law() -> str:
    return print_expression(law)


@validate_input(
    rotational_inertia_=rotational_inertia,
    pendulum_mass_=pendulum_mass,
    distance_to_pivot_=distance_to_pivot,
)
@validate_output(oscillation_period)
def calculate_period(
    rotational_inertia_: Quantity,
    pendulum_mass_: Quantity,
    distance_to_pivot_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        pendulum_mass: pendulum_mass_,
        distance_to_pivot: distance_to_pivot_,
    })
    return Quantity(result)
