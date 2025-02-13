"""
Quality factor of loaded resonator from circuit parameters
==========================================================

If the resonator is an oscillatory circuit to which an external circuit is connected,
then the loaded Q-factor of the resonator depends on the resistances of the resonator
and the external circuit, as well as on the inductance of the resonator and the
angular oscillation frequency.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    clone_as_symbol,
)

loaded_quality_factor = clone_as_symbol(symbols.quality_factor, subscript="1")
"""
:symbols:`quality_factor` of the loaded resonator.
"""

resistance = clone_as_symbol(symbols.electrical_resistance, subscript="0")
"""
:symbols:`electrical_resistance` of the oscillating circuit.
"""

inductance = symbols.inductance
"""
:symbols:`inductance` of the oscillating circuit.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the current.
"""

load_resistance = clone_as_symbol(symbols.electrical_resistance, display_symbol="R_L", display_latex="R_\\text{L}")
"""
:symbols:`electrical_resistance` of the load.
"""

law = Eq(loaded_quality_factor, (load_resistance * resistance) /
    (angular_frequency * inductance * (load_resistance + resistance)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resistance_=resistance,
    inductance_=inductance,
    frequency_=angular_frequency,
    load_resistance_=load_resistance)
@validate_output(loaded_quality_factor)
def calculate_quality_factor(resistance_: Quantity, inductance_: Quantity, frequency_: Quantity,
    load_resistance_: Quantity) -> float:
    result_expr = solve(law, loaded_quality_factor,
        dict=True)[0][loaded_quality_factor]
    result_expr = result_expr.subs({
        resistance: resistance_,
        inductance: inductance_,
        angular_frequency: frequency_,
        load_resistance: load_resistance_,
    })
    return convert_to_float(result_expr)
