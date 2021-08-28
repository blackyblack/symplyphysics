from sympy import symbols, Function, Eq, pretty, solve, dsolve
from sympy.utilities.lambdify import lambdify, implemented_function
import sympy.physics.units as phy_units
from sympy.physics.units import Quantity
from sympy.physics.units.systems.si import SI

# Description
## Acceleration definition: A = dv/dt

time = symbols('time')
acceleration, velocity_function = symbols('acceleration velocity', cls = Function)
definition = Eq(acceleration(time), velocity_function(time).diff(time))
dim_definition_SI = phy_units.meter/phy_units.second**2

def print():
    return pretty(definition, use_unicode=False)

def print_dimension():
    return pretty(dim_definition_SI, use_unicode=False)

def calculate_linear_acceleration(velocity_start_: Quantity, velocity_end_: Quantity, time_: Quantity) -> Quantity:
    dimsys_SI = SI.get_dimension_system()
    assert dimsys_SI.equivalent_dims(velocity_start_.dimension, phy_units.velocity)
    assert dimsys_SI.equivalent_dims(velocity_end_.dimension, phy_units.velocity)
    assert dimsys_SI.equivalent_dims(time_.dimension, phy_units.time)
    
    velocity_function_ = implemented_function('velocity_function', lambda x: x * (velocity_end_ - velocity_start_) / time_)
    velocity_function_lambda = lambdify(time, velocity_function_(time))
    # solve differential equation with custom function
    dsolved = dsolve(definition.subs(velocity_function(time), velocity_function_lambda(time)), acceleration(time))
    # calculate acceleration at the given point of time
    ## Note: since acceleration is constant, any time argument can be passed here.
    ## For nonlinear functions 'time' should be passed for a correct answer.
    solved = solve(dsolved.subs(time, time_))
    result_expr = solved[0][acceleration(time_)]

    quantity_scale = SI._collect_factor_and_dimension(result_expr)
    result_quantity = Quantity('acceleration')
    dimsys_SI.set_quantity_dimension(result_quantity, quantity_scale[1])
    dimsys_SI.set_quantity_scale_factor(result_quantity, quantity_scale[0])
    assert dimsys_SI.equivalent_dims(result_quantity.dimension, phy_units.acceleration)
    return result_quantity