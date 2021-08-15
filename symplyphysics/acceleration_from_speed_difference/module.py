from sympy import symbols, Function, Eq, pretty, solve, dsolve
from sympy.utilities.lambdify import lambdify, implemented_function

# Description
## Acceleration definition: A = dv/dt

time_ = symbols('time')
acceleration_, velocity_function_ = symbols('acceleration velocity', cls = Function)
law = Eq(acceleration_(time_), velocity_function_(time_).diff(time_))

def print():
    return pretty(law, use_unicode=False)

def calculate_linear_acceleration(velocity_start, velocity_end, time):
    velocity_function = implemented_function('velocity_function', lambda x: x * (velocity_end - velocity_start) / time)
    velocity_function_lambda = lambdify(time_, velocity_function(time_))
    # solve differential equation with custom function
    dsolved = dsolve(law.subs(velocity_function_(time_), velocity_function_lambda(time_)), acceleration_(time_))
    # calculate acceleration at the given point of time
    ## Note: since acceleration is constant, any time argument can be passed here.
    ## For nonlinear functions 'time' should be passed for a correct answer.
    solved = solve(dsolved.subs(time_, time))
    return solved[0][acceleration_(time)].evalf()