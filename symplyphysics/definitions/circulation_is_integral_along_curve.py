from typing import List
from sympy import Integral, Derivative
from sympy.vector import Dot
from symplyphysics import (
    symbols, Eq, pretty, simplify, CoordSys3D, array_to_sympy_vector, VectorField, apply_field
)

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
## C = CurveIntegral(F * dl, Curve), where
## C is circulation
## F is vector field
## l is radius-vector of the curve
## * is dot product

# Conditions
## - Curve is a function of a single parameter (eg y(x) = x**2), or parametrized with a single parameter
## (eg x(t) = cos(t), y(t) = sin(t))
## - Curve is smooth and continuous

# trajectory is a function of the moving particle
circulation, field, trajectory, = symbols('circulation field trajectory')
# trajectory_element (dl) is trajectory derivative by parameter
# parameter is an argument of trajectory function, eg x coordinate
trajectory_element, parameter = symbols('trajectory_element parameter')
parameter_from, parameter_to = symbols('parameter_from parameter_to')

trajectory_element_definition = Eq(trajectory_element, Derivative(trajectory, parameter))
definition = Eq(circulation, Integral(Dot(field, trajectory_element), (parameter, parameter_from, parameter_to)))

def print():
    return pretty(definition, use_unicode=False)

# field_ should be VectorField
# trajectory_ should be array with projections to coordinates, eg [3 * cos(parameter), 3 * sin(parameter)]
def calculate_circulation(
    coord_system_: CoordSys3D,
    field_: VectorField,
    trajectory_: List,
    parameter_from_,
    parameter_to_):

    field_app = apply_field(field_, trajectory_)
    field_as_vector = array_to_sympy_vector(coord_system_, field_app)
    trajectory_as_vector = array_to_sympy_vector(coord_system_, trajectory_)
    trajectory_element_result = trajectory_element_definition.rhs.subs(trajectory, trajectory_as_vector).doit()
    result_expr = definition.rhs.subs({field: field_as_vector, trajectory_element: trajectory_element_result, parameter_from: parameter_from_, parameter_to: parameter_to_}).doit()
    # some expressions are invalid without simplifying them first
    return simplify(result_expr)
