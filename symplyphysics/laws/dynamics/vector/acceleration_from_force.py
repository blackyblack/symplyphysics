from sympy import solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    QuantityVector,
    Vector,
    scale_vector,
    validate_input,
    validate_output,
    list_of_quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    momentum_is_mass_times_velocity as momentum_def,
    acceleration_is_velocity_derivative as acceleration_def,
)
from symplyphysics.laws.dynamics.vector import force_is_derivative_of_momentum as force_momentum_law

# Description
## Newton's second law in vector form: a = 1/m * F
## Where:
## F - force vector,
## m - mass,
## a - acceleration vector,
## * - scalar multiplication (scale vector).

mass = Symbol("mass", units.mass)


def acceleration_law(force_: Vector) -> Vector:
    return scale_vector(1 / mass, force_)


def force_law(acceleration_: Vector) -> Vector:
    return scale_vector(mass, acceleration_)


# Derive this law from law of force and momentum
# Condition: mass is constant

time = force_momentum_law.time
momentum_x = Function("momentum_x", units.mass * units.velocity)
momentum_y = Function("momentum_y", units.mass * units.velocity)
momentum_z = Function("momentum_z", units.mass * units.velocity)

momentum_vec = Vector([momentum_x(time), momentum_y(time), momentum_z(time)])
force_derived = force_momentum_law.force_law(momentum_vec)

momentum_def_sub = momentum_def.definition.subs(momentum_def.mass, mass)
velocity_x = solve(momentum_def_sub, momentum_def.velocity)[0].subs(momentum_def.momentum, momentum_x(time))
velocity_y = solve(momentum_def_sub, momentum_def.velocity)[0].subs(momentum_def.momentum, momentum_y(time))
velocity_z = solve(momentum_def_sub, momentum_def.velocity)[0].subs(momentum_def.momentum, momentum_z(time))

acceleration_def_sub = acceleration_def.definition.rhs.subs(acceleration_def.time, time)
acceleration_x = acceleration_def_sub.subs(acceleration_def.velocity(time), velocity_x)
acceleration_y = acceleration_def_sub.subs(acceleration_def.velocity(time), velocity_y)
acceleration_z = acceleration_def_sub.subs(acceleration_def.velocity(time), velocity_z)

acceleration_vec = Vector([acceleration_x, acceleration_y, acceleration_z])
force_from_law = force_law(acceleration_vec)

for component_derived, component_from_law in zip(force_derived.components, force_from_law.components):
    assert expr_equals(component_derived, component_from_law)


@validate_input(mass_=mass, acceleration_=units.acceleration)
@validate_output(units.force)
def calculate_force(mass_: Quantity, acceleration_: QuantityVector) -> QuantityVector:
    result_force = force_law(acceleration_)
    force_components = list_of_quantities(result_force.components, {mass: mass_})
    return QuantityVector(force_components, acceleration_.coordinate_system)


@validate_input(mass_=mass, force_=units.force)
@validate_output(units.acceleration)
def calculate_acceleration(mass_: Quantity, force_: QuantityVector) -> QuantityVector:
    result_acceleration = acceleration_law(force_)
    acceleration_components = list_of_quantities(
        result_acceleration.components, {mass: mass_}
    )
    return QuantityVector(acceleration_components, force_.coordinate_system)
