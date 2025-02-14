from sympy import solve
from symplyphysics import (
    units,
    Quantity,
    QuantityVector,
    Vector,
    scale_vector,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import (
    acceleration_is_speed_derivative as acceleration_def,
    momentum_is_mass_times_speed as momentum_def,
)
from symplyphysics.laws.dynamics.vector import force_is_derivative_of_momentum as force_momentum_law

# Description
## Newton's second law in vector form: a = 1/m * F
## Where:
## F - force vector,
## m - mass,
## a - acceleration vector,
## * - scalar multiplication (scale vector).

mass = symbols.mass


def acceleration_law(force_: Vector) -> Vector:
    return scale_vector(1 / mass, force_)


def force_law(acceleration_: Vector) -> Vector:
    return scale_vector(mass, acceleration_)


# Derive this law from law of force and momentum
# Condition: mass is constant

time = force_momentum_law.time

_momentum_x = clone_as_function(symbols.momentum, [time], subscript="x")
_momentum_y = clone_as_function(symbols.momentum, [time], subscript="y")
_momentum_z = clone_as_function(symbols.momentum, [time], subscript="z")
_momentum_vec = Vector([_momentum_x(time), _momentum_y(time), _momentum_z(time)])

_force_derived = force_momentum_law.force_law(_momentum_vec)

_momentum_def_sub = momentum_def.definition.subs(momentum_def.mass, mass)
_velocity_from_momentum = solve(_momentum_def_sub, momentum_def.speed)[0]
_velocity_vec = Vector([
    _velocity_from_momentum.subs(momentum_def.momentum, momentum_component)
    for momentum_component in _momentum_vec.components
])

_acceleration_def_sub = acceleration_def.definition.rhs.subs(acceleration_def.time, time)
_acceleration_vec = Vector([
    _acceleration_def_sub.subs(acceleration_def.speed(time), velocity_component)
    for velocity_component in _velocity_vec.components
])

_force_from_law = force_law(_acceleration_vec)

for _component_derived, _component_from_law in zip(_force_derived.components,
    _force_from_law.components):
    assert expr_equals(_component_derived, _component_from_law)


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
