from sympy import (Eq, solve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression, Function,
    validate_input_symbols, validate_output_symbol)

# Description
## The total mechanical energy of a system is conserved i.e., the energy can neither be created nor be destroyed.
## It can only be internally converted from one form to another if the forces doing work on the system are conservative in nature.

# Law: E(t1) = E(t0)
## Where:
## E - summary mechanical energy of a system,
## t1 - point of time, when interaction of a system with conservative forces occured,
## t0 - initial time.

time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)
mechanical_energy = Function("mechanical_energy", units.energy)

law = Eq(mechanical_energy(time_after), mechanical_energy(time_before))


def print() -> str:
    return print_expression(law)


@validate_input_symbols(mechanical_energy_before_=mechanical_energy)
@validate_output_symbol(mechanical_energy)
def calculate_energy_after(mechanical_energy_before_: Quantity) -> Quantity:
    solved = solve(law, mechanical_energy(time_after), dict=True)[0][mechanical_energy(time_after)]
    result_expr = solved.subs(mechanical_energy(time_before), mechanical_energy_before_)
    return expr_to_quantity(result_expr)
