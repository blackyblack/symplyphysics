from sympy import (Eq, solve, dsolve)
from symplyphysics import (units, Quantity, Symbol, print_expression, Function, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import mass_is_constant

# Description
## The total mass of the closed system is preserved.
## For more information, see [mass_is_constant](./mass_is_constant.py ).

# Law: m(t1) = m(t0)
## Where:
## m - summary mass of a system,
## t1 - the moment of time after the action in the system,
## t0 - initial time.

time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)
mass = Function("mas", units.mass)

law = Eq(mass(time_after), mass(time_before))

# Derive the same law from constant mass

## dsolve() shows that solution is constant C1
dsolved = dsolve(
    mass_is_constant.law,
    mass_is_constant.mass(mass_is_constant.time))

energy_before_eq = dsolved.subs(mass_is_constant.time, time_before)
energy_before_eq = energy_before_eq.subs(
    mass_is_constant.mass(time_before), mass(time_before))
energy_after_eq = dsolved.subs(mass_is_constant.time, time_after)
energy_after_eq = energy_after_eq.subs(mass_is_constant.mass(time_after),
    mass(time_after))

## Show that when mass is constant, mass_before equals to mass_after
energy_after_solved = solve([energy_after_eq, energy_before_eq],
    (mass(time_after), "C1"),
    dict=True)[0][mass(time_after)]
assert expr_equals(energy_after_solved, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_before_=mass)
@validate_output(mass)
def calculate_mass_after(mass_before_: Quantity) -> Quantity:
    solved = solve(law, mass(time_after), dict=True)[0][mass(time_after)]
    result_expr = solved.subs(mass(time_before), mass_before_)
    return Quantity(result_expr)
