from sympy import (Eq, solve, dsolve)
from symplyphysics import (units, expr_to_quantity, Quantity, Symbol, print_expression, Function,
    validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import mechanical_energy_is_constant

# Description
## The total mechanical energy of an isolated system is conserved i.e., the energy can neither be created nor be destroyed.
## It can only be internally converted from one form to another if the forces doing work on the system are conservative in nature.
## See [mechanical_energy_is_constant](./mechanical_energy_is_constant.py) for more information.

# Law: E(t1) = E(t0)
## Where:
## E - summary mechanical energy of a system,
## t1 - point of time, when interaction of a system with conservative forces occured,
## t0 - initial time.

time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)
mechanical_energy = Function("mechanical_energy", units.energy)

law = Eq(mechanical_energy(time_after), mechanical_energy(time_before))

# Derive the same law from constant mechanical energy

## dsolve() shows that solution is constant C1
dsolved = dsolve(mechanical_energy_is_constant.law, mechanical_energy_is_constant.mechanical_energy(mechanical_energy_is_constant.time))

energy_before_eq = dsolved.subs(mechanical_energy_is_constant.time, time_before)
energy_before_eq = energy_before_eq.subs(mechanical_energy_is_constant.mechanical_energy(time_before), mechanical_energy(time_before))
energy_after_eq = dsolved.subs(mechanical_energy_is_constant.time, time_after)
energy_after_eq = energy_after_eq.subs(mechanical_energy_is_constant.mechanical_energy(time_after), mechanical_energy(time_after))

## Show that when energy is constant, energy_before equals to energy_after
energy_after_solved = solve([energy_after_eq, energy_before_eq], (mechanical_energy(time_after), "C1"), dict=True)[0][mechanical_energy(time_after)]
assert expr_equals(energy_after_solved, law.rhs)


def print() -> str:
    return print_expression(law)


@validate_input(mechanical_energy_before_=mechanical_energy)
@validate_output(mechanical_energy)
def calculate_energy_after(mechanical_energy_before_: Quantity) -> Quantity:
    solved = solve(law, mechanical_energy(time_after), dict=True)[0][mechanical_energy(time_after)]
    result_expr = solved.subs(mechanical_energy(time_before), mechanical_energy_before_)
    return expr_to_quantity(result_expr)
