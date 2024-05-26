from sympy import sqrt, symbols, sympify
from sympy.physics.units import speed_of_light
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    scale_vector,
    QuantityVector,
    add_cartesian_vectors,
    dot_vectors,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.relativistic.vector import force_via_acceleration

# Description
# In special relativity, the Newton's second law does not hold in the form `F = m * a`. There still exists a relation
## between force `a` and acceleration `F`, albeit more complex.

# Law: a = 1 / (m0 * gamma) * (F - dot(v, F) * v / c**2)
## a - vector of acceleration
## m0 - rest mass
## v - vector of velocity
## c - speed of light
## F - vector of force
## gamma = 1 / sqrt(1 - dot(v, v) / c**2) - Lorentz factor
## dot(a, b) - dot product between vectors `a` and `b`

# Notes
## - For force and rest mass expressions, see [force via acceleration law](./force_via_acceleration.py)

rest_mass = Symbol("rest_mass", units.mass)


def acceleration_law(force_: Vector, velocity_: Vector) -> Vector:
    force_parallel_to_velocity_ = scale_vector(
        -1 * dot_vectors(velocity_, force_) / speed_of_light**2,
        velocity_,
    )

    resulting_force = add_cartesian_vectors(force_, force_parallel_to_velocity_)

    mass_factor_ = sqrt(1 - dot_vectors(velocity_, velocity_) / speed_of_light**2) / rest_mass

    return scale_vector(mass_factor_, resulting_force)

# Show that the force-to-acceleration and acceleration-to-force laws are consistent with each other

_acceleration = Vector(symbols("acceleration_x:z", real=True))
_force = Vector(symbols("force_x:z", real=True))
_velocity = Vector(symbols("velocity_x:z", real=True))

_acceleration_via_force = acceleration_law(_force, _velocity).simplify()
_force_derived = force_via_acceleration.force_law(_acceleration_via_force, _velocity).simplify()

for _force_derived_component, _force_component in zip(_force_derived.components, _force.components):
    assert expr_equals(
        sympify(_force_derived_component).subs(force_via_acceleration.rest_mass, rest_mass),
        _force_component,
    )

_force_via_acceleration = force_via_acceleration.force_law(_acceleration, _velocity).simplify()
_acceleration_derived = acceleration_law(_force_via_acceleration, _velocity).simplify()

for _acceleration_derived_component, _acceleration_component in zip(
    _acceleration_derived.components, _acceleration.components
):
    assert expr_equals(
        sympify(_acceleration_derived_component).subs(force_via_acceleration.rest_mass, rest_mass),
        _acceleration_component,
    )

# TODO: do the same for the rest mass formula


@validate_input(
    rest_mass_=rest_mass,
    force_=units.force,
    velocity_=units.velocity,
)
@validate_output(units.acceleration)
def calculate_acceleration(
    rest_mass_: Quantity,
    force_: QuantityVector,
    velocity_: QuantityVector,
) -> QuantityVector:
    result = acceleration_law(force_.to_base_vector(), velocity_.to_base_vector())
    return QuantityVector.from_base_vector(
        result,
        subs={rest_mass: rest_mass_},
    )
