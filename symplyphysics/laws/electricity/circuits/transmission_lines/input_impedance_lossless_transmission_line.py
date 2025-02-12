"""
Input impedance of lossless transmission line
=============================================

Knowing the length of the transmission line and its surge impedance, as well as the
propagation constant and load impedance, it is possible to determine the input impedance
of the transmission line.

..
    TODO: find link
"""

from sympy import Eq, solve, tan, I
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
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

phase_constant = symbols.phase_constant
"""
:symbols:`phase_constant`.
"""

law = Eq(
    input_impedance,
    surge_impedance *
    (load_impedance + I * surge_impedance * tan(phase_constant * length)) /
    (surge_impedance + I * load_impedance * tan(phase_constant * length)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(characteristic_resistance_=surge_impedance,
    load_impedance_=load_impedance,
    constant_propagation_=phase_constant,
    length_=length)
@validate_output(input_impedance)
def calculate_input_impedance(characteristic_resistance_: Quantity, load_impedance_: Quantity,
    constant_propagation_: Quantity, length_: Quantity) -> Quantity:
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        surge_impedance: characteristic_resistance_,
        load_impedance: load_impedance_,
        phase_constant: constant_propagation_,
        length: length_,
    })
    return Quantity(result_expr)
