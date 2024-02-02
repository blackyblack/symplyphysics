from sympy import Eq, dsolve, solve, Symbol as SymSymbol
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import (
    momentum_after_collision_equals_to_momentum_before as momentum_conservation_law,
)
from symplyphysics.definitions import (
    momentum_is_mass_times_velocity as momentum_def,
    mass_flow_rate as flow_rate_def
)
from symplyphysics.laws.kinematic import (
    accelerated_velocity_from_time as acceleration_def
)

# Description
## Assuming we are at rest relative to an inertial reference frame, we observe a rocket
## accelerate through space with no gravitational or atmospheric drag forces acting on it.
## The mass of the rocket changes as it burns fuel and releases the products of burning,
## but the total mass of the system does not change.

# Law: R * v_rel = M * a
## R = -dM/dt - rate of fuel consumption
## v_rel - velocity of rocket relative to products
## M - rocket mass
## a - rocket acceleration

# Note: (R * v_rel) is called thrust of rocket engine

# Conditions:
## - Fuel consumption rate is constant
## - Non-relativistic velocities

fuel_consumption_rate = Symbol("fuel_consumption_rate", units.mass / units.time)
relative_velocity = Symbol("relative_velocity", units.velocity)
rocket_mass = Symbol("rocket_mass", units.mass)
rocket_acceleration = Symbol("rocket_acceleration", units.acceleration)

law = Eq(fuel_consumption_rate * relative_velocity, rocket_mass * rocket_acceleration)


# Derive this law from the law of conservation of momentum.

rocket_speed = SymSymbol("rocket_speed")
rocket_speed_change = SymSymbol("rocket_speed_change")
fuel_mass_thrusted = SymSymbol("fuel_mass_thrusted")

rocket_momentum_before_release = momentum_def.definition.rhs.subs({
    momentum_def.mass: rocket_mass,
    momentum_def.velocity: rocket_speed,
})

rocket_momentum_after_release = momentum_def.definition.rhs.subs({
    momentum_def.mass: rocket_mass - fuel_mass_thrusted,
    momentum_def.velocity: rocket_speed + rocket_speed_change,
})

products_speed = SymSymbol("products_speed")

products_momentum = momentum_def.definition.rhs.subs({
    momentum_def.mass: fuel_mass_thrusted,
    momentum_def.velocity: products_speed,
})

# Summary momentum is conserved
final_momentum = rocket_momentum_after_release + products_momentum

rocket_speed_relative_to_frame = rocket_speed + rocket_speed_change
rocket_speed_relative_to_products = relative_velocity
products_speed_relative_to_frame = products_speed

# NOTE: probably needs a separate law - this is Galilean law of addition of velocities
relative_speed_eqn = Eq(
    rocket_speed_relative_to_frame,
    # velocities are directed in opposite sides
    rocket_speed_relative_to_products + products_speed_relative_to_frame,
)

momentum_conservation_eqn = momentum_conservation_law.law.subs({
    momentum_conservation_law.momentum(momentum_conservation_law.time_before): rocket_momentum_before_release,
    momentum_conservation_law.momentum(momentum_conservation_law.time_after): final_momentum,
})

relative_velocity_expr = solve(
    [relative_speed_eqn, momentum_conservation_eqn],
    (products_speed, relative_velocity),
    dict=True,
)[0][relative_velocity]

time_change = SymSymbol("time_change")

# solve differential equation with constant fuel_consumption_rate
dsolved_fuel_mass = dsolve(flow_rate_def.definition.subs(
    flow_rate_def.mass_flow_rate(flow_rate_def.time), fuel_consumption_rate),
    flow_rate_def.mass(flow_rate_def.time))
fuel_consumption_eqn = Eq(fuel_mass_thrusted, dsolved_fuel_mass.rhs)
# C1 is initial fuel mass thrusted
fuel_consumption_eqn = fuel_consumption_eqn.subs({"C1": 0, flow_rate_def.time: time_change})

rocket_acceleration_expr = solve(
    acceleration_def.law,
    acceleration_def.acceleration,
)[0].subs({
    acceleration_def.velocity: rocket_speed_change,
    acceleration_def.time: time_change,
    acceleration_def.initial_velocity: 0,
})

relative_velocity_eqn_system = [
    Eq(relative_velocity, relative_velocity_expr),
    fuel_consumption_eqn,
    Eq(rocket_acceleration, rocket_acceleration_expr),
]

relative_velocity_derived = solve(
    relative_velocity_eqn_system,
    (relative_velocity, fuel_mass_thrusted, rocket_speed_change),
    dict=True
)[0][relative_velocity]

relative_velocity_from_law = solve(law, relative_velocity)[0]

assert expr_equals(relative_velocity_derived, relative_velocity_from_law)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    fuel_consumption_rate_=fuel_consumption_rate,
    rocket_mass_=rocket_mass,
    rocket_acceleration_=rocket_acceleration,
)
@validate_output(relative_velocity)
def calculate_relative_velocity(
    fuel_consumption_rate_: Quantity,
    rocket_mass_: Quantity,
    rocket_acceleration_: Quantity,
) -> Quantity:
    result = solve(law, relative_velocity)[0].subs({
        fuel_consumption_rate: fuel_consumption_rate_,
        rocket_mass: rocket_mass_,
        rocket_acceleration: rocket_acceleration_,
    })
    return Quantity(result)
