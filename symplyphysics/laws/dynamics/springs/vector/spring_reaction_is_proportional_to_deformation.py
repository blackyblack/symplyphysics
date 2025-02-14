from symplyphysics import (Vector, QuantityVector, scale_vector, units, Quantity,
    validate_input, validate_output, symbols)

# Description
## Deformed spring is about to return back to it's undeformed state and responds with some force. Law is:
## F = -k * x, where
## F is force of response vector,
## k is spring stiffness,
## x is vector of deformation.

# Condition
## Deformation is elactic (reversible).

# Also see its [scalar counterpart](../spring_reaction_is_proportional_to_deformation.py)

stiffness = symbols.stiffness


def force_law(deformation_: Vector) -> Vector:
    return scale_vector(-1 * stiffness, deformation_)


def deformation_law(force_: Vector) -> Vector:
    return scale_vector(-1 / stiffness, force_)


@validate_input(coefficient_=stiffness, deformation_=units.length)
@validate_output(units.force)
def calculate_force(coefficient_: Quantity, deformation_: QuantityVector) -> QuantityVector:
    result_vector = force_law(deformation_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={stiffness: coefficient_})


@validate_input(coefficient_=stiffness, force_=units.force)
@validate_output(units.length)
def calculate_deformation(coefficient_: Quantity, force_: QuantityVector) -> QuantityVector:
    result_vector = deformation_law(force_.to_base_vector())
    return QuantityVector.from_base_vector(result_vector, subs={stiffness: coefficient_})
