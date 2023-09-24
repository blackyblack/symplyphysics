from typing import Sequence
from sympy import (Expr, Integral, Derivative, Symbol as SymSymbol, Eq, simplify, sympify)
from symplyphysics import (dot_vectors, print_expression, Quantity, Vector)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.vector_field import VectorField

# Description
## Field circulation along closed curve is curvilinear integral of field by this curve.
## A space curve is parametrized by a vector-valued function that depends upon a single parameter t that varies over some interval.
## We shall always assume that x(t) is continuously differentiable. The curve is smooth provided its tangent
## vector is continuous and everywhere nonzero.
## Physically, we can think of a curve as the trajectory described by a particle moving in
## space. At each time t, the tangent vector x'(t) represents the instantaneous velocity of the
## particle. Thus, as long as the particle moves with nonzero speed, its trajectory is necessarily a smooth curve.
## When curve is not smooth (eg trajectory is a square) one should represent this curve as a sum of its smooth parts
## (eg sum of sides of a square).

# Definition
## C = CurveIntegral(dot(F, dr), Curve)
# Where:
## C is circulation
## F is vector field
## dr is curve tangent vector
## dot is dot product

# Conditions
## - Curve is a function of a single parameter (eg y(x) = x**2), or parametrized with a single parameter
## (eg x(t) = cos(t), y(t) = sin(t))
## - Curve is smooth and continuous

# These are not physical symbols - SymPy 'Symbol' is good enough.

# trajectory is a function of the moving particle
circulation = SymSymbol("circulation")
trajectory = SymSymbol("trajectory")
# trajectory_element (dl) is trajectory derivative by parameter
# parameter is an argument of trajectory function, eg x coordinate
trajectory_element = SymSymbol("trajectory_element")
parameter = SymSymbol("parameter")
parameter_from = SymSymbol("parameter_from")
parameter_to = SymSymbol("parameter_to")
field_dot_tangent = SymSymbol("field_dot_tangent")

trajectory_element_definition = Eq(trajectory_element, Derivative(trajectory, parameter))
definition = Eq(
    circulation,
    Integral(field_dot_tangent,
    (parameter, parameter_from, parameter_to)))


def print_law() -> str:
    return print_expression(definition)


def _calculate_dot_tangent(field_: VectorField, trajectory_: Vector) -> Expr:
    trajectory_sympy_vector = trajectory_.to_sympy_vector()
    trajectory_element_sympy_vector = trajectory_element_definition.rhs.subs(
        trajectory, trajectory_sympy_vector).doit()
    trajectory_element_vector = Vector.from_sympy_vector(trajectory_element_sympy_vector,
        trajectory_.coordinate_system)
    trajectory_components = [sympify(c) for c in trajectory_.components]
    field_applied = field_.apply(trajectory_components)
    return dot_vectors(field_applied, trajectory_element_vector)


# trajectory_ should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def calculate_circulation(field_: VectorField, trajectory_: Sequence[Expr], parameter_limits: tuple[ScalarValue, ScalarValue]) -> Quantity:
    trajectory_vector = Vector(trajectory_, field_.coordinate_system)
    field_dot_tangent_value = _calculate_dot_tangent(field_, trajectory_vector)
    result_expr = definition.rhs.subs({
        field_dot_tangent: field_dot_tangent_value,
        parameter_from: parameter_limits[0],
        parameter_to: parameter_limits[1]
    }).doit()
    return Quantity(simplify(result_expr))
