"""
Input impedance from transmission matrix
========================================

Knowing the transmission matrix of the device, it is possible to determine the input
impedance of this device.

**Notes:**

#. See :ref:`Transmission matrix`.

**Links:**

#. `Cadence System Analysis, derivable from here <https://resources.system-analysis.cadence.com/blog/msa2021-abcd-parameters-of-transmission-lines>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

input_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_in", display_latex="Z_\\text{in}")
"""
Input :symbols:`electrical_impedance` of the transmission line.
"""

load_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_L", display_latex="Z_\\text{L}")
"""
Load :symbols:`electrical_impedance`.
"""

voltage_voltage_parameter = SymbolNew("A", dimensionless)
"""
Ratio of input voltage to output voltage at idle at the output.
"""

voltage_current_parameter = SymbolNew("B", units.impedance)
"""
Ratio of input voltage to output current in case of a short circuit at the output.
"""

current_voltage_parameter = SymbolNew("C", units.conductance)
"""
Ratio of input current to output voltage at idle at the output.
"""

current_current_parameter = SymbolNew("D", dimensionless)
"""
Ratio of input current to output current in case of a short circuit at the output.
"""

law = Eq(input_impedance, (voltage_voltage_parameter * load_impedance + voltage_current_parameter) /
    (current_voltage_parameter * load_impedance + current_current_parameter))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(load_resistance_=load_impedance)
@validate_output(input_impedance)
def calculate_input_impedance(
        load_resistance_: Quantity, parameters_: tuple[tuple[float, Quantity], tuple[Quantity,
    float]]) -> Quantity:
    assert_equivalent_dimension(parameters_[0][1], "parameters_[0][1]", "calculate_input_impedance",
        units.impedance)
    assert_equivalent_dimension(parameters_[1][0], "parameters_[1][0]", "calculate_input_impedance",
        units.conductance)
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        load_impedance: load_resistance_,
        voltage_voltage_parameter: parameters_[0][0],
        voltage_current_parameter: parameters_[0][1],
        current_voltage_parameter: parameters_[1][0],
        current_current_parameter: parameters_[1][1],
    })
    return Quantity(result_expr)
