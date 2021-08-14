from sympy import symbols, Eq, pretty, solve

# Description
## Newton's second law: F = m * a

force_, mass_, acceleration_ = symbols('force mass acceleration')
law = Eq(force_, mass_ * acceleration_)

def print():
    return pretty(law, use_unicode=False)

def calculate_force(mass, acceleration):
    return solve(law.subs(mass_, mass).subs(acceleration_, acceleration))[0].evalf()