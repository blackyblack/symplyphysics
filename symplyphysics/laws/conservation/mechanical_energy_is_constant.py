from sympy import (Eq, dsolve, Derivative)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression, Function,
    validate_input, validate_output)

# Description
## Mechanical energy, sum of the kinetic energy, or energy of motion, and the potential energy, or energy stored in a system by
## reason of the position of its parts.
## Mechanical energy is constant in a system that has only gravitational forces or in an otherwise idealized system — that is,
## one lacking dissipative forces, such as friction and air resistance, or one in which such forces can be reasonably neglected.

# Law: dE/dt = 0
## Where:
## E - summary mechanical energy of a system,
## t - time,
## dE/dt - energy derivative over time.

# Note:
## SymPy does not have a proper way to represent constant energy. We use it's derivative over time instead. Derivative
## of the constant value is zero.

time = Symbol("time", units.time)
mechanical_energy = Function("mechanical_energy", units.energy)

law = Eq(Derivative(mechanical_energy(time), time), 0)


def print_law() -> str:
    return print_expression(law)


@validate_input(mechanical_energy_before_=mechanical_energy)
@validate_output(mechanical_energy)
def calculate_energy_after(mechanical_energy_before_: Quantity) -> Quantity:
    solved = dsolve(law, mechanical_energy(time))
    result_expr = solved.subs("C1", mechanical_energy_before_).rhs
    return expr_to_quantity(result_expr)
