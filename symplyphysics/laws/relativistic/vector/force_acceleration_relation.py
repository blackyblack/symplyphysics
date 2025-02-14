from sympy import Expr, symbols as sym_symbols, sympify
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    Vector,
    vector_magnitude,
    add_cartesian_vectors,
    dot_vectors,
    scale_vector,
    QuantityVector,
    symbols,
)
from symplyphysics.quantities import speed_of_light
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.vectors.arithmetics import project_vector, reject_cartesian_vector
from symplyphysics.definitions import lorentz_factor as lorentz_factor_law

# Description
## In special relativity, the Newton's second law does not hold in the form `F = m * a`. There still exists a relation
## between force `F` and acceleration `a`, albeit more complex.

# Force law: F = gamma**3 * m0 * a_parallel + gamma * m0 * a_orthogonal
# Acceleration law: a = 1 / (m0 * gamma) * (F - dot(v, F) * v / c**2)
## F - force vector
## a - acceleration vector
## v - velocity vector
## a = a_parallel + a_orthogonal - acceleration vector
## a_parallel = dot(a, v) / dot(v, v) * v - acceleration component parallel to velocity
## a_orthogonal = a - a_parallel - acceleration component orthogonal to velocity
## m0 - rest mass
## gamma = 1 / sqrt(1 - dot(v, v) / c**2) - Lorentz factor
## dot(a, b) - dot product between vectors `a` and `b`

# Conditions
## - This law applies to special relativity.

# Links: Wikipedia, see paragraph <https://en.wikipedia.org/wiki/Acceleration_(special_relativity)#Acceleration_and_force>

rest_mass = symbols.rest_mass


def acceleration_law(force_: Vector, velocity_: Vector) -> Vector:
    force_parallel_to_velocity_ = scale_vector(
        -1 * dot_vectors(velocity_, force_) / speed_of_light**2,
        velocity_,
    )

    resulting_force = add_cartesian_vectors(force_, force_parallel_to_velocity_)

    mass_factor_ = rest_mass * lorentz_factor_law.definition.rhs.subs({
        lorentz_factor_law.speed: vector_magnitude(velocity_),
    })

    return scale_vector(1 / mass_factor_, resulting_force)


def force_law(acceleration_: Vector, velocity_: Vector) -> Vector:
    acceleration_parallel_ = project_vector(acceleration_, velocity_)
    acceleration_orthogonal_ = reject_cartesian_vector(acceleration_, velocity_)

    lorentz_factor_ = lorentz_factor_law.definition.rhs.subs({
        lorentz_factor_law.speed: vector_magnitude(velocity_),
    })

    force_parallel_ = scale_vector(
        lorentz_factor_**3 * rest_mass,
        acceleration_parallel_,
    )
    force_orthogonal_ = scale_vector(
        lorentz_factor_ * rest_mass,
        acceleration_orthogonal_,
    )

    return add_cartesian_vectors(
        force_parallel_,
        force_orthogonal_,
    )


def rest_mass_law(
    force_: Vector,
    acceleration_: Vector,
    velocity_: Vector,
) -> Expr:
    acceleration_parallel_ = project_vector(acceleration_, velocity_)
    acceleration_orthogonal_ = reject_cartesian_vector(acceleration_, velocity_)

    lorentz_factor_ = lorentz_factor_law.definition.rhs.subs({
        lorentz_factor_law.speed: vector_magnitude(velocity_),
    })

    lhs = vector_magnitude(force_)
    rhs = vector_magnitude(
        add_cartesian_vectors(
        scale_vector(lorentz_factor_**3, acceleration_parallel_),
        scale_vector(lorentz_factor_, acceleration_orthogonal_),
        ))

    return lhs / rhs


# Show that the force-to-acceleration and acceleration-to-force laws are consistent with each other

_acceleration = Vector(sym_symbols("acceleration_x:z", real=True))
_force = Vector(sym_symbols("force_x:z", real=True))
_velocity = Vector(sym_symbols("velocity_x:z", real=True))

_acceleration_via_force = acceleration_law(_force, _velocity).simplify()
_force_derived = force_law(_acceleration_via_force, _velocity).simplify()

for _force_derived_component, _force_component in zip(_force_derived.components, _force.components):
    assert expr_equals(
        sympify(_force_derived_component),
        _force_component,
    )

_force_via_acceleration = force_law(_acceleration, _velocity).simplify()
_acceleration_derived = acceleration_law(_force_via_acceleration, _velocity).simplify()

for _acceleration_derived_component, _acceleration_component in zip(
        _acceleration_derived.components, _acceleration.components):
    assert expr_equals(
        sympify(_acceleration_derived_component),
        _acceleration_component,
    )


@validate_input(
    rest_mass_=rest_mass,
    acceleration_=units.acceleration,
    velocity_=units.velocity,
)
@validate_output(units.force)
def calculate_force(
    rest_mass_: Quantity,
    acceleration_: QuantityVector,
    velocity_: QuantityVector,
) -> QuantityVector:
    result = force_law(
        acceleration_.to_base_vector(),
        velocity_.to_base_vector(),
    )
    return QuantityVector.from_base_vector(
        result,
        subs={rest_mass: rest_mass_},
    )
