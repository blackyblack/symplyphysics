from symplyphysics import (Vector, QuantityVector, scale_vector, units, Quantity, Symbol,
    validate_input, validate_output, list_of_quantities)

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
    result_force = force_law(deformation_)
    force_components = list_of_quantities(result_force.components,
        {elastic_coefficient: coefficient_})
    return QuantityVector(force_components, deformation_.coordinate_system)


@validate_input(coefficient_=elastic_coefficient, force_=units.force)
@validate_output(units.length)
def calculate_deformation(coefficient_: Quantity, force_: QuantityVector) -> QuantityVector:
    result_deformation = deformation_law(force_)
    deformation_components = list_of_quantities(result_deformation.components,
        {elastic_coefficient: coefficient_})
    return QuantityVector(deformation_components, force_.coordinate_system)
