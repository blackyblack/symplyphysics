"""
Input impedance of lossy transmission line
==========================================

Knowing the length of the transmission line, the loss factor of the transmission line
and its characteristic resistance, as well as the propagation constant and load
resistance, it is possible to determine the input impedance of the transmission line.

..
    TODO: find link
"""

from sympy import Eq, solve, sinh, cosh, I, evaluate
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    SymbolNew,
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

propagation_constant = SymbolNew("b", 1 / units.length)
"""
The **propagation constant** is the inverse of the signal :symbols:`wavelength`.
"""

loss_factor = SymbolNew("a", 1 / units.length)
"""
The **loss factor** shows how many times the transmitted signal weakens per unit
:symbols:`length` of the transmission line.
"""

# the following block prevents the re-ordering of terms for the code printer
with evaluate(False):
    _expression = (loss_factor + I * propagation_constant) * length

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
    length_=length,
    loss_factor_=loss_factor)
@validate_output(input_impedance)
def calculate_input_impedance(surge_impedance_: Quantity, load_resistance_: Quantity,
    constant_propagation_: Quantity, length_: Quantity, loss_factor_: Quantity) -> Quantity:
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        surge_impedance: surge_impedance_,
        load_impedance: load_resistance_,
        propagation_constant: constant_propagation_,
        length: length_,
        loss_factor: loss_factor_,
    })
    return Quantity(result_expr)
