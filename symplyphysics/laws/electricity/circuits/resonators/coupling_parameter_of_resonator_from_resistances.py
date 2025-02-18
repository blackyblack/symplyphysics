"""
Coupling parameter of resonator from resistance
===============================================

There is a coupling parameter to describe the resonator and the load. The parameter is
equal to the ratio of the resonator's resistance to the load resistance.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

coupling_parameter = Symbol("g", dimensionless)
"""
Coupling parameter of the resonator.
"""

resonator_resistance = clone_as_symbol(symbols.electrical_resistance, subscript="0")
"""
:symbols:`electrical_resistance` of the resonator.
"""

load_resistance = clone_as_symbol(symbols.electrical_resistance, display_symbol="R_L", display_latex="R_\\text{L}")
"""
:symbols:`electrical_resistance` of the resonator.
"""

law = Eq(coupling_parameter, resonator_resistance / load_resistance)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resonator_resistance_=resonator_resistance, load_resistance_=load_resistance)
@validate_output(coupling_parameter)
def calculate_coupling_parameter(resonator_resistance_: Quantity,
    load_resistance_: Quantity) -> float:
    result_expr = solve(law, coupling_parameter, dict=True)[0][coupling_parameter]
    result_expr = result_expr.subs({
        resonator_resistance: resonator_resistance_,
        load_resistance: load_resistance_,
    })
    return convert_to_float(result_expr)
