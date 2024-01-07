from sympy import Expr
from symplyphysics import (units, Quantity, Vector,
    validate_input, validate_output)
from symplyphysics.core.vectors.arithmetics import dot_vectors
from symplyphysics.core.vectors.vectors import QuantityVector

# Description
## Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
## Work is scalar value equal to dot product of force and movement vectors.
## Law: A = F * S, where
## A is mechanical work
## F is force vector applied to object
## S is movement vector caused by this force
## * is a scalar multiplication of vectors (dot product)


def work_law(force_: Vector, distance_: Vector) -> Expr:
    return dot_vectors(force_, distance_)


@validate_input(force_=units.force, distance_=units.length)
@validate_output(units.energy)
def calculate_work(force_: QuantityVector, distance_: QuantityVector) -> Quantity:
    result_work = work_law(force_, distance_)
    return Quantity(result_work)
