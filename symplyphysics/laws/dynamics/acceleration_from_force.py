from sympy import symbols, Eq, pretty, solve

# Description
## Newton's second law: F = m * a

force, mass, acceleration = symbols('force mass acceleration')
law = Eq(acceleration, force / mass)

def print():
    return pretty(law, use_unicode=False)

def calculate_force(mass_, acceleration_):
    return solve(law.subs(mass, mass_).subs(acceleration, acceleration_))[0].evalf()