"""
Input impedance of lossy transmission line
==========================================

Knowing the length of the transmission line, the loss factor of the transmission line
and its surge impedance, as well as the propagation constant and load resistance, it is
possible to determine the input impedance of the transmission line.

..
    TODO: find link
"""

from sympy import Eq, solve, sinh, cosh, evaluate
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

input_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_in", display_latex="Z_\\text{in}")
"""
Input :symbols:`electrical_impedance` of the transmission line.
"""

load_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_L", display_latex="Z_\\text{L}")
"""
Load :symbols:`electrical_impedance`.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the transmission line.
"""

length = symbols.length
"""
:symbols:`length` of the transmission line.
"""

propagation_constant = symbols.propagation_constant
"""
:symbols:`propagation_constant`.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _expression = propagation_constant * length

law = Eq(input_impedance,
    (cosh(_expression) * load_impedance + surge_impedance * sinh(_expression))
    / (load_impedance * sinh(_expression) / surge_impedance + cosh(_expression)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(surge_impedance_=surge_impedance,
    load_resistance_=load_impedance,
    constant_propagation_=propagation_constant,
    length_=length)
@validate_output(input_impedance)
def calculate_input_impedance(surge_impedance_: Quantity, load_resistance_: Quantity,
    constant_propagation_: Quantity, length_: Quantity) -> Quantity:
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        surge_impedance: surge_impedance_,
        load_impedance: load_resistance_,
        propagation_constant: constant_propagation_,
        length: length_,
    })
    return Quantity(result_expr)
