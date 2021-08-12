import sympy

# Description
## Newton's second law: F = m * a

force_, mass_, acceleration_ = sympy.symbols('force mass acceleration')
law = sympy.Eq(force_, mass_ * acceleration_)

def print():
    return sympy.pretty(law, use_unicode=False)

def calculate_force(mass, acceleration):
    return sympy.solvers.solve(law.subs(mass_, mass).subs(acceleration_, acceleration))[0]