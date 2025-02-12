"""
Transmission matrix for serial load in line
===========================================

Knowing the impedance of the load connected in series to the transmission line, it is
possible to calculate the parameters :math:`A, B, C, D` of the transmission matrix of
the load.

**Notes:**

#. See :ref:`Transmission matrix`.

**Conditions:**

#. The load is connected to the transmission line *in series*.

..
    TODO: find link
"""

from sympy import Eq, solve, Matrix, S
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    dimensionless,
    convert_to,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.dimensions import assert_equivalent_dimension

voltage_voltage_parameter = SymbolNew("A", dimensionless)
"""
Ratio of input :symbols:`voltage` to output :symbols:`voltage` at idle at the output.
"""

voltage_current_parameter = SymbolNew("B", units.impedance)
"""
Ratio of input :symbols:`voltage` to output :symbols:`current` in case of a short
circuit at the output.
"""

current_voltage_parameter = SymbolNew("C", units.conductance)
"""
Ratio of input :symbols:`current` to output :symbols:`voltage` at idle at the output.
"""

current_current_parameter = SymbolNew("D", dimensionless)
"""
Ratio of input :symbols:`current` to output :symbols:`current` in case of a short
circuit at the output.
"""

load_impedance = clone_as_symbol(symbols.electrical_impedance, display_symbol="Z_L", display_latex="Z_\\text{L}")
"""
Load :symbols:`electrical_impedance`.
"""

law = Eq(
    Matrix([[voltage_voltage_parameter, voltage_current_parameter],
    [current_voltage_parameter, current_current_parameter]]), Matrix([[1, load_impedance], [0, 1]]))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(load_impedance_=load_impedance)
def calculate_transmission_matrix(
        load_impedance_: Quantity) -> tuple[tuple[float, Quantity], tuple[Quantity, float]]:
    result = solve(law, [
        voltage_voltage_parameter, voltage_current_parameter, current_voltage_parameter,
        current_current_parameter
    ],
        dict=True)[0]
    result_a = result[voltage_voltage_parameter]
    result_b = result[voltage_current_parameter]
    result_c = result[current_voltage_parameter]
    result_d = result[current_current_parameter]
    substitutions = {
        load_impedance: load_impedance_,
    }
    result_a = float(convert_to(Quantity(result_a.subs(substitutions)), S.One).evalf())
    result_b = Quantity(result_b.subs(substitutions))
    result_c = Quantity(result_c.subs(substitutions))
    result_d = float(convert_to(Quantity(result_d.subs(substitutions)), S.One).evalf())
    assert_equivalent_dimension(result_b, 'result_b', "calculate_transmission_matrix",
        units.impedance)
    return ((result_a, result_b), (result_c, result_d))
