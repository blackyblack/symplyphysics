from sympy import (Eq, solve, dsolve)
from symplyphysics import (units, Quantity, Symbol, Function, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import mass_is_constant

# Description
## The total mass of the closed system is preserved.
## For more information, see [mass_is_constant](./mass_is_constant.py).

# Law: m(t1) = m(t0)
## Where:
## m - summary mass of a system,
## t1 - the moment of time after the action in the system,
## t0 - initial time.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Conservation_of_mass>

time_before = Symbol("time_before", units.time)
time_after = Symbol("time_after", units.time)
mass_function = Function("mass_function", units.mass)

law = Eq(mass_function(time_after), mass_function(time_before))

# Derive the same law from constant mass

## dsolve() shows that solution is constant C1
dsolved = dsolve(mass_is_constant.law, mass_is_constant.mass_function(mass_is_constant.time))

mass_before_eq = dsolved.subs(mass_is_constant.time, time_before)
mass_before_eq = mass_before_eq.subs(mass_is_constant.mass_function(time_before),
    mass_function(time_before))
mass_after_eq = dsolved.subs(mass_is_constant.time, time_after)
mass_after_eq = mass_after_eq.subs(mass_is_constant.mass_function(time_after),
    mass_function(time_after))

## Show that when mass is constant, mass_before equals to mass_after
mass_after_solved = solve([mass_after_eq, mass_before_eq], (mass_function(time_after), "C1"),
    dict=True)[0][mass_function(time_after)]
assert expr_equals(mass_after_solved, law.rhs)


@validate_input(mass_before_=mass_function)
@validate_output(mass_function)
def calculate_mass_after(mass_before_: Quantity) -> Quantity:
    solved = solve(law, mass_function(time_after), dict=True)[0][mass_function(time_after)]
    result_expr = solved.subs(mass_function(time_before), mass_before_)
    return Quantity(result_expr)
