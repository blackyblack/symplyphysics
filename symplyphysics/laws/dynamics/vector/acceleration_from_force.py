from sympy import sympify
from symplyphysics import (units, Quantity, Symbol, QuantityVector, Vector, scale_vector,
    validate_input, validate_output)

# Description
## Newton's second law in vector form: a = 1/m * F
## Where:
## F - force vector.
## m - mass.
## a - acceleration vector
## * - scalar multiplication (scale vector)

mass = Symbol("mass", units.mass)


def acceleration_law(force_: Vector) -> Vector:
    return scale_vector(1 / mass, force_)


def force_law(acceleration_: Vector) -> Vector:
    return scale_vector(mass, acceleration_)


@validate_input(mass_=units.mass, acceleration_=units.acceleration)
@validate_output(units.force)
def calculate_force(mass_: Quantity, acceleration_: QuantityVector) -> QuantityVector:
    quantities_vector = Vector(acceleration_.components, acceleration_.coordinate_system)
    result_force = force_law(quantities_vector)
    force_components = []
    for c in result_force.components:
        with_mass = sympify(c).subs(mass, mass_)
        force_components.append(Quantity(with_mass))
    return QuantityVector(force_components, acceleration_.coordinate_system)
