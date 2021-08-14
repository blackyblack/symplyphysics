import sympy

# Description
## Acceleration definition: A = dv/dt

acceleration_, time_ = sympy.symbols('acceleration time')
velocity_function_ = sympy.Function('velocity')
law = sympy.Eq(acceleration_, sympy.Derivative(velocity_function_(time_), time_))

def print():
    return sympy.pretty(law, use_unicode=False)

def calculate_linear_acceleration(velocity_start, velocity_end, time):
    # solve differential equation with initial conditions
    dsolved = sympy.solvers.ode.dsolve(law, velocity_function_(time_), ics = {velocity_function_(0): velocity_start})
    # calculate acceleration for given velocity change over time
    solved = sympy.solvers.solve(dsolved.subs(velocity_function_(time_), velocity_end).subs(time_, time))
    return solved[0].evalf()