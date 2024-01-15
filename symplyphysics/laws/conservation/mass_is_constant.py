from sympy import (Eq, dsolve, Derivative)
from symplyphysics import (units, Quantity, Symbol, print_expression, Function, validate_input,
    validate_output)

# Description
## The mass is constant in a system that is closed, and mass is not transformed to energy

# Law: dm/dt = 0
## Where:
## m - summary mass of a system,
## t - time,
## dm/dt - mass derivative over time.

# Conditions:
## - System in a closed impenetrable volume, that is, molecules/atoms cannot leave it
##   and they are always inside;
## - Mass is not transformed to energy, for example due to annihilation.

# Note:
## SymPy does not have a proper way to represent constant mass. We use it's derivative over time instead. Derivative
## of the constant value is zero.

time = Symbol("time", units.time)
mass = Function("mechanical_energy", units.mass)

law = Eq(Derivative(mass(time), time), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_before_=mass)
@validate_output(mass)
def calculate_mass_after(mass_before_: Quantity) -> Quantity:
    solved = dsolve(law, mass(time))
    result_expr = solved.subs("C1", mass_before_).rhs
    return Quantity(result_expr)
