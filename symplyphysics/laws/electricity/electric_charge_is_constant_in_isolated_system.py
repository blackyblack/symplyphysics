"""
Electric charge is constant in isolated system
==============================================

The electric charge of an isolated system is conserved. As a result, its charge is constant
at all times.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

initial_charge = clone_as_symbol(symbols.charge, subscript="0")
"""
Initial :symbols:`charge` of the system.
"""

final_charge = clone_as_symbol(symbols.charge, subscript="1")
"""
Final :symbols:`charge` of the system.
"""

law = Eq(final_charge, initial_charge)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(charge_before_=initial_charge)
@validate_output(final_charge)
def calculate_charge_after(charge_before_: Quantity) -> Quantity:
    solved = solve(law, final_charge, dict=True)[0][final_charge]
    result_expr = solved.subs(initial_charge, charge_before_)
    return Quantity(result_expr)
