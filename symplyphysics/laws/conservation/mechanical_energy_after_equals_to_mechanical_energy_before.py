"""
Initial mechanical energy equals final mechanical energy
========================================================

The total mechanical energy of an isolated system is conserved, i.e. the energy can
neither be created nor be destroyed. It can only be internally converted from one form
to another if the forces doing work on the system are conservative in nature.

**Notes:**

#. Also see :ref:`Mechanical energy is constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Conservation_of_energy>`__.

..
    TODO: rename file
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
from symplyphysics.laws.conservation import mechanical_energy_is_constant

initial_time = clone_as_symbol(symbols.time, subscript="0")
"""
Initial :symbols:`time`.
"""

final_time = clone_as_symbol(symbols.time, subscript="1")
"""
Final :symbols:`time`.
"""

mechanical_energy = clone_as_function(symbols.mechanical_energy, [symbols.time])
"""
:symbols:`mechanical_energy` of the system as a function of :symbols:`time`.
"""

law = Eq(mechanical_energy(final_time), mechanical_energy(initial_time))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the same law from constant mechanical energy

## dsolve() shows that solution is constant C1
dsolved = dsolve(
    mechanical_energy_is_constant.law,
    mechanical_energy_is_constant.mechanical_energy(mechanical_energy_is_constant.time))

energy_before_eq = dsolved.subs(mechanical_energy_is_constant.time, initial_time)
energy_before_eq = energy_before_eq.subs(
    mechanical_energy_is_constant.mechanical_energy(initial_time), mechanical_energy(initial_time))
energy_after_eq = dsolved.subs(mechanical_energy_is_constant.time, final_time)
energy_after_eq = energy_after_eq.subs(mechanical_energy_is_constant.mechanical_energy(final_time),
    mechanical_energy(final_time))

## Show that when energy is constant, energy_before equals to energy_after
energy_after_solved = solve([energy_after_eq, energy_before_eq],
    (mechanical_energy(final_time), "C1"),
    dict=True)[0][mechanical_energy(final_time)]
assert expr_equals(energy_after_solved, law.rhs)


@validate_input(mechanical_energy_before_=mechanical_energy)
@validate_output(mechanical_energy)
def calculate_energy_after(mechanical_energy_before_: Quantity) -> Quantity:
    solved = solve(law, mechanical_energy(final_time), dict=True)[0][mechanical_energy(final_time)]
    result_expr = solved.subs(mechanical_energy(initial_time), mechanical_energy_before_)
    return Quantity(result_expr)
