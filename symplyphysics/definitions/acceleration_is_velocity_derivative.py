from sympy import symbols, Function, Eq, pretty, solve, dsolve
from sympy.utilities.lambdify import lambdify, implemented_function

# Description
## Acceleration definition: A = dv/dt

time = symbols('time')
acceleration, velocity_function = symbols('acceleration velocity', cls = Function)
law = Eq(acceleration(time), velocity_function(time).diff(time))

def print():
    return pretty(law, use_unicode=False)

def calculate_linear_acceleration(velocity_start_, velocity_end_, time_):
    velocity_function_ = implemented_function('velocity_function', lambda x: x * (velocity_end_ - velocity_start_) / time_)
    velocity_function_lambda = lambdify(time, velocity_function_(time))
    # solve differential equation with custom function
    dsolved = dsolve(law.subs(velocity_function(time), velocity_function_lambda(time)), acceleration(time))
    # calculate acceleration at the given point of time
    ## Note: since acceleration is constant, any time argument can be passed here.
    ## For nonlinear functions 'time' should be passed for a correct answer.
    solved = solve(dsolved.subs(time, time_))
    return solved[0][acceleration(time_)].evalf()