"""
Coupling parameter of resonator from quality factor
===================================================

There is a coupling parameter to describe the resonator and the load. The parameter is
equal to the ratio of the resonator's own quality factor to the quality factor of the
external circuit.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
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

resonator_quality_factor = clone_as_symbol(symbols.quality_factor, subscript="0")
"""
:symbols:`quality_factor` of the resonator.
"""

external_circuit_quality_factor = clone_as_symbol(symbols.quality_factor, display_symbol="Q_e", display_latex="Q_\\text{e}")
"""
:symbols:`quality_factor` of the external circuit.
"""

law = Eq(coupling_parameter, resonator_quality_factor / external_circuit_quality_factor)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(resonator_quality_factor_=resonator_quality_factor,
    external_circuit_quality_factor_=external_circuit_quality_factor)
@validate_output(coupling_parameter)
def calculate_coupling_parameter(resonator_quality_factor_: float,
    external_circuit_quality_factor_: float) -> float:
    result_expr = solve(law, coupling_parameter, dict=True)[0][coupling_parameter]
    result_expr = result_expr.subs({
        resonator_quality_factor: resonator_quality_factor_,
        external_circuit_quality_factor: external_circuit_quality_factor_,
    })
    return convert_to_float(result_expr)
