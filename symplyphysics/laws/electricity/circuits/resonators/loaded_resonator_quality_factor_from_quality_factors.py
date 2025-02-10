"""
Quality resonator of loaded resonator from quality factors
==========================================================

The loaded quality factor of the resonator depends on the own quality factor of the
resonator and the quality factor of the external circuit connected to the resonator.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
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

quality_factor = clone_as_symbol(symbols.quality_factor, subscript="0")
"""
:symbols:`quality_factor` of the resonator without load.
"""

external_circuit_quality_factor = clone_as_symbol(symbols.quality_factor, display_symbol="Q_e", display_latex="Q_\\text{e}")
"""
:symbols:`quality_factor` of the external circuit.
"""

law = Eq(1 / loaded_quality_factor,
    1 / quality_factor + 1 / external_circuit_quality_factor)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(quality_factor_=quality_factor,
    external_circuit_quality_factor_=external_circuit_quality_factor)
@validate_output(loaded_quality_factor)
def calculate_quality_factor(quality_factor_: float,
    external_circuit_quality_factor_: float) -> float:
    result_expr = solve(law, loaded_quality_factor,
        dict=True)[0][loaded_quality_factor]
    result_expr = result_expr.subs({
        quality_factor: quality_factor_,
        external_circuit_quality_factor: external_circuit_quality_factor_,
    })
    return convert_to_float(result_expr)
