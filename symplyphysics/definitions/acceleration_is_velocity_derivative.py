from sympy import symbols, Function, Eq, pretty, solve, dsolve
from sympy.utilities.lambdify import lambdify, implemented_function
from sympy.physics.units import Quantity
import sympy.physics.units as phy_units
from symplyphysics.quantity_decorator import validate_input, validate_output
import symplyphysics.expr_to_quantity

# Description
## Acceleration definition: A = dv/dt

time = symbols('time')
acceleration, velocity_function = symbols('acceleration velocity', cls = Function)
definition = Eq(acceleration(time), velocity_function(time).diff(time))
definition_dimension_SI = phy_units.meter / phy_units.second**2

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(definition_dimension_SI, use_unicode=False)

@validate_input(velocity_start_=phy_units.velocity, velocity_end_=phy_units.velocity, time_=phy_units.time)
@validate_output(phy_units.acceleration)
def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity, time_: Quantity) -> Quantity:
    velocity_function_ = implemented_function('velocity_function', lambda x: x * (velocity_end_ - velocity_start_) / time_)
    velocity_function_lambda = lambdify(time, velocity_function_(time))
    # solve differential equation with custom function
    dsolved = dsolve(definition.subs(velocity_function(time), velocity_function_lambda(time)), acceleration(time))
    # calculate acceleration at the given point of time
    ## Note: since acceleration is constant, any time argument can be passed here.
    ## For nonlinear functions 'time' should be passed for a correct answer.
    solved = solve(dsolved.subs(time, time_))
    result_expr = solved[0][acceleration(time_)]
    return symplyphysics.expr_to_quantity.convert(result_expr, 'acceleration')
