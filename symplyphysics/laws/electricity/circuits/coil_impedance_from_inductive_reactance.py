"""
Coil impedance from inductive reactance
=======================================

The impedance of ideal coil depends on its inductive reactance. Having zero resistivity,
the real part of coil impedance is zero.

..
    TODO: find link
"""

from sympy import (I, Eq, solve)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols
)

impedance = symbols.electrical_impedance
"""
:symbols:`electrical_impedance` of a coil.
"""

reactance = symbols.electrical_reactance
"""
:symbols:`electrical_reactance` of a coil.
"""

law = Eq(impedance, I * reactance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(reactance_=reactance)
@validate_output(impedance)
def calculate_impedance(reactance_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, impedance, dict=True)[0][impedance]
    result_expr = result_impedance_expr.subs({
        reactance: reactance_,
    })
    return Quantity(result_expr)
