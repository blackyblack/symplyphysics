from sympy import Integral, Derivative
from sympy.vector import Dot
from symplyphysics import (
    symbols, Eq, pretty, simplify
)

# Description
## Field circulation along closed curve is curvilinear integral of field by this curve.
## C = CurveIntegral(F * dl, Curve), where
## C is circulation
## F is field vector function
## l if radius-vector of the curve

## For the sake of simplicity, this module defines Circulation as an integral by one parameter. Fields and trajectories
## with more dimensions should be parametrized by one parameter or calculated separately, eg by dx and dy.
## Or use sympy.vector.vector_integrate().

# trajectory is a function of the moving particle. Should be a function of a single parameter.
circulation, field, trajectory, = symbols('circulation field trajectory')
# trajectory_element (dl) is trajectory derivative by parameter
# parameter is an argument of trajectory function, eg x coordinate
trajectory_element, parameter = symbols('trajectory_element parameter')
parameter_from, parameter_to = symbols('parameter_from parameter_to')

trajectory_element_definition = Eq(trajectory_element, Derivative(trajectory, parameter))
definition = Eq(circulation, Integral(Dot(field, trajectory_element), (parameter, parameter_from, parameter_to)))

def print():
    return pretty(definition, use_unicode=False)

# field_ and trajectory_ should be instances of CoordSys3D
def calculate_circulation(
    field_,
    trajectory_,
    parameter_from_,
    parameter_to_):

    trajectory_element_result = trajectory_element_definition.rhs.subs(trajectory, trajectory_).doit()
    result_expr = definition.rhs.subs({field: field_, trajectory_element: trajectory_element_result, parameter_from: parameter_from_, parameter_to: parameter_to_}).doit()
    # some expressions are invalid without simplifying them first
    return simplify(result_expr)
