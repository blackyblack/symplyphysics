from sympy import sympify
from symplyphysics import (Vector, scale_vector, units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.vectors.vectors import QuantityVector

# Description
## Deformed spring is about to return back to it's undeformed state and responds with some force. Law is:
## F = -kx, where
## F is force of response vector,
## k is elastic coefficient,
## x is vector of deformation.

# Condition
## Deformation is elactic (reversible).

elastic_coefficient = Symbol("elastic_coefficient", units.force / units.length)


def force_law(deformation_: Vector) -> Vector:
    return scale_vector(-1 * elastic_coefficient, deformation_)


def deformation_law(force_: Vector) -> Vector:
    return scale_vector(-1 / elastic_coefficient, force_)


@validate_input(coefficient_=elastic_coefficient, deformation_=units.length)
@validate_output(units.force)
def calculate_force(coefficient_: Quantity, deformation_: QuantityVector) -> QuantityVector:
    quantities_vector = Vector(deformation_.components, deformation_.coordinate_system)
    result_force = force_law(quantities_vector)
    force_components = []
    for c in result_force.components:
        with_coefficient = sympify(c).subs(elastic_coefficient, coefficient_)
        force_components.append(Quantity(with_coefficient))
    return QuantityVector(force_components, deformation_.coordinate_system)


@validate_input(coefficient_=elastic_coefficient, force_=units.force)
@validate_output(units.length)
def calculate_deformation(coefficient_: Quantity, force_: QuantityVector) -> QuantityVector:
    quantities_vector = Vector(force_.components, force_.coordinate_system)
    result_deformation = deformation_law(quantities_vector)
    deformation_components = []
    for c in result_deformation.components:
        with_coefficient = sympify(c).subs(elastic_coefficient, coefficient_)
        deformation_components.append(Quantity(with_coefficient))
    return QuantityVector(deformation_components, force_.coordinate_system)
