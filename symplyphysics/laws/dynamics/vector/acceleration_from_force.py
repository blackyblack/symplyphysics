from sympy import sympify
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)
from symplyphysics.core.vectors.vectors import Vector
from symplyphysics.core.vectors.vector_arithmetics import scale_vector

# Description
## Newton's second law in vector form: a = 1/m * F
## Where:
## F - force vector.
## m - mass.
## a - acceleration vector
## * - scalar multiplication (scale vector)

mass = Symbol("mass", units.mass)

def acceleration_law(force_: Vector) -> Vector:
    return scale_vector(1/mass, force_)

def force_law(acceleration_: Vector) -> Vector:
    return scale_vector(mass, acceleration_)


# TODO: use one dimension for entire vector, not for each component
# TODO: add helpers to generate Vector easily
# TODO: add Vector and arithmetics to __init__

@validate_input(mass_=units.mass, acceleration_=units.acceleration)
@validate_output(units.force)
def calculate_force(mass_: Quantity, acceleration_: Vector) -> Vector:
    result_force = force_law(acceleration_)
    force_components = []
    for c in result_force.components:
        with_mass = sympify(c).subs(mass, mass_)
        force_components.append(Quantity(with_mass))
    return Vector(force_components, result_force.coordinate_system)
