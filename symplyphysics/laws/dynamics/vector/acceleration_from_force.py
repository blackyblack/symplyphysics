from sympy import solve
from symplyphysics import (
    units,
    Quantity,
    Function,
    QuantityVector,
    Vector,
    scale_vector,
    validate_input,
    validate_output,
    symbols,
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

mass = symbols.basic.mass


def acceleration_law(force_: Vector) -> Vector:
    return scale_vector(1 / mass, force_)


def force_law(acceleration_: Vector) -> Vector:
    return scale_vector(mass, acceleration_)


# Derive this law from law of force and momentum
# Condition: mass is constant

time = force_momentum_law.time

momentum_x = Function("momentum_x", units.momentum)
momentum_y = Function("momentum_y", units.momentum)
momentum_z = Function("momentum_z", units.momentum)
momentum_vec = Vector([momentum_x(time), momentum_y(time), momentum_z(time)])

force_derived = force_momentum_law.force_law(momentum_vec)

momentum_def_sub = momentum_def.definition.subs(momentum_def.mass, mass)
velocity_from_momentum = solve(momentum_def_sub, momentum_def.speed)[0]
velocity_vec = Vector([
    velocity_from_momentum.subs(momentum_def.momentum, momentum_component)
    for momentum_component in momentum_vec.components
])

acceleration_def_sub = acceleration_def.definition.rhs.subs(acceleration_def.time, time)
acceleration_vec = Vector([
    acceleration_def_sub.subs(acceleration_def.speed(time), velocity_component)
    for velocity_component in velocity_vec.components
])

force_from_law = force_law(acceleration_vec)

for component_derived, component_from_law in zip(force_derived.components,
    force_from_law.components):
    assert expr_equals(component_derived, component_from_law)


@validate_input(mass_=mass, acceleration_=units.acceleration)
@validate_output(units.force)
def calculate_force(mass_: Quantity, acceleration_: QuantityVector) -> QuantityVector:
    result_vector = force_law(acceleration_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={mass: mass_})


@validate_input(mass_=mass, force_=units.force)
@validate_output(units.acceleration)
def calculate_acceleration(mass_: Quantity, force_: QuantityVector) -> QuantityVector:
    result_vector = acceleration_law(force_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={mass: mass_})
