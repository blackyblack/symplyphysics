"""
Relative operating bandwidth of quarter-wave transformer
========================================================

The relative operating bandwidth of a quarter-wave transformer depends on the reflection
coefficient, the surge impendance of the transmission line to which the transformer is
connected and the load resistance.

..
    TODO: find link
"""

from sympy import Eq, solve, pi, acos, sqrt
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

relative_bandwidth = Symbol("b", dimensionless)
"""
Relative operating bandwidth is the ratio of the bandwidth to the center frequency.
"""

load_resistance = clone_as_symbol(symbols.electrical_resistance, display_symbol="R_L", display_latex="R_\\text{L}")
"""
:symbols:`electrical_resistance` of the load.
"""

surge_impedance = symbols.surge_impedance
"""
:symbols:`surge_impedance` of the transmission line.
"""

reflection_coefficient = symbols.reflection_coefficient
"""
:symbols:`reflection_coefficient`.
"""

law = Eq(
    relative_bandwidth, 2 -
    (4 / pi) * acos(2 * reflection_coefficient * sqrt(load_resistance * surge_impedance) /
    (sqrt(1 - reflection_coefficient**2) * abs(load_resistance - surge_impedance))))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(load_resistance_=load_resistance,
    characteristic_resistance_=surge_impedance,
    reflection_coefficient_=reflection_coefficient)
@validate_output(relative_bandwidth)
def calculate_relative_bandwidth(load_resistance_: Quantity, characteristic_resistance_: Quantity,
    reflection_coefficient_: float) -> float:
    result_expr = solve(law, relative_bandwidth, dict=True)[0][relative_bandwidth]
    result_expr = result_expr.subs({
        load_resistance: load_resistance_,
        surge_impedance: characteristic_resistance_,
        reflection_coefficient: reflection_coefficient_
    })
    return convert_to_float(result_expr)
