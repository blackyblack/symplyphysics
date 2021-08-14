from sympy import symbols, Function, Eq, pretty, solve, dsolve

# Description
## Acceleration definition: A = dv/dt

acceleration_, time_ = symbols('acceleration time')
velocity_function_ = Function('velocity')
law = Eq(acceleration_, velocity_function_(time_).diff(time_))

def print():
    return pretty(law, use_unicode=False)

def calculate_linear_acceleration(velocity_start, velocity_end, time):
    # solve differential equation with initial conditions
    dsolved = dsolve(law, velocity_function_(time_), ics = {velocity_function_(0): velocity_start})
    # calculate acceleration for given velocity change over time
    solved = solve(dsolved.subs(velocity_function_(time_), velocity_end).subs(time_, time))
    return solved[0].evalf()