"""
Initial mass equals final mass
==============================

The total mass of a closed system is conserved. For more information, see
:ref:`Mass is constant`.

**Conditions:**

#. The system is isolated, i.e. no particles can leave it.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Conservation_of_mass>`__.
"""

from sympy import Eq, solve, dsolve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
    clone_as_symbol,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.conservation import mass_is_constant

initial_time = clone_as_symbol(symbols.time, subscript="0")
"""
Initial :symbols:`time`.
"""

final_time = clone_as_symbol(symbols.time, subscript="1")
"""
Final :symbols:`time`.
"""

mass = clone_as_function(symbols.mass, [symbols.time])
"""
:symbols:`mass` as a function of :symbols:`time`.
"""

law = Eq(mass(final_time), mass(initial_time))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from constant mass

## dsolve() shows that solution is constant C1
_dsolved = dsolve(mass_is_constant.law, mass_is_constant.mass(mass_is_constant.time))

_mass_before_eq = _dsolved.subs(mass_is_constant.time, initial_time)
_mass_before_eq = _mass_before_eq.subs(mass_is_constant.mass(initial_time),
    mass(initial_time))
_mass_after_eq = _dsolved.subs(mass_is_constant.time, final_time)
_mass_after_eq = _mass_after_eq.subs(mass_is_constant.mass(final_time),
    mass(final_time))

## Show that when mass is constant, mass_before equals to mass_after
_mass_after_solved = solve([_mass_after_eq, _mass_before_eq], (mass(final_time), "C1"),
    dict=True)[0][mass(final_time)]
assert expr_equals(_mass_after_solved, law.rhs)


@validate_input(mass_before_=mass)
@validate_output(mass)
def calculate_mass_after(mass_before_: Quantity) -> Quantity:
    solved = solve(law, mass(final_time), dict=True)[0][mass(final_time)]
    result_expr = solved.subs(mass(initial_time), mass_before_)
    return Quantity(result_expr)
