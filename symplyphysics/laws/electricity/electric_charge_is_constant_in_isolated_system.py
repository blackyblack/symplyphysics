"""
Electric charge is constant in isolated system
==============================================

The electric charge of an isolated system is conserved. As a result, its charge is constant
at all times.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

initial_charge = Symbol("initial_charge", units.charge)
r"""
Initial charge of the system.

Symbol:
    :code:`q_0`

Latex:
    :math:`q_1`
"""

final_charge = Symbol("final_charge", units.charge)
r"""
Final charge of the system.

Symbol:
    :code:`q_1`

Latex:
    :math:`q_1`
"""

law = Eq(final_charge, initial_charge)


@validate_input(charge_before_=initial_charge)
@validate_output(final_charge)
def calculate_charge_after(charge_before_: Quantity) -> Quantity:
    solved = solve(law, final_charge, dict=True)[0][final_charge]
    result_expr = solved.subs(initial_charge, charge_before_)
    return Quantity(result_expr)
