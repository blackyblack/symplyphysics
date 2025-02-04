"""
Admittance of rectangular loop coupler
======================================

The rectangular loop coupler consists of four sections. The admittance of each section
can be calculated by calculating the admittance of the transmission line to which the
coupler is connected and the power ratio at the outputs.

..
    TODO: find link
    TODO: rename file
"""

from sympy import Eq, solve, Matrix, sqrt
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

first_admittance = clone_as_symbol(symbols.admittance, subscript="1")
"""
:symbols:`admittance` of the first section.
"""

second_admittance = clone_as_symbol(symbols.admittance, subscript="2")
"""
:symbols:`admittance` of the second section.
"""

third_admittance = clone_as_symbol(symbols.admittance, subscript="3")
"""
:symbols:`admittance` of the third section.
"""

fourth_admittance = clone_as_symbol(symbols.admittance, subscript="4")
"""
:symbols:`admittance` of the fourth section.
"""

transmission_line_admittance = clone_as_symbol(symbols.admittance, subscript="0")
"""
:symbols:`admittance` of the transmission line.
"""

power_ratio = SymbolNew("k", dimensionless)
"""
Ratio of the power at the outputs of the coupler.
"""

law = Eq(
    Matrix([first_admittance, second_admittance, third_admittance, fourth_admittance]),
    transmission_line_admittance * Matrix([
        1 / sqrt(power_ratio),
        sqrt((power_ratio + 1) / power_ratio),
        sqrt((power_ratio + 1) / power_ratio),
        1 / sqrt(power_ratio),
    ]))
r"""
..
    Code printers don't work with matrices yet.

:code:`[Y_1, Y_2, Y_3, Y_4] = Y_0 * [1 / sqrt(k), sqrt((k + 1) / k), sqrt((k + 1) / k), 1 / sqrt(k)]`

Latex:
    .. math::
        \begin{pmatrix} Y_1 \\ Y_2 \\ Y_3 \\ Y_4 \end{pmatrix}
        = Y_0 \begin{pmatrix}
            \frac{1}{\sqrt{k}} \\
            \sqrt{\frac{k + 1}{k}} \\
            \sqrt{\frac{k + 1}{k}} \\
            \frac{1}{\sqrt{k}} \\
        \end{pmatrix}
"""


@validate_input(transmission_line_admittance_=transmission_line_admittance,
    ratio_of_power_=power_ratio)
@validate_output(units.conductance)
def calculate_conductivities(
        transmission_line_admittance_: Quantity,
        ratio_of_power_: float) -> tuple[Quantity, Quantity, Quantity, Quantity]:
    result = solve(law,
        [first_admittance, second_admittance, third_admittance, fourth_admittance],
        dict=True)[0]
    result_y1 = result[first_admittance]
    result_y2 = result[second_admittance]
    result_y3 = result[third_admittance]
    result_y4 = result[fourth_admittance]
    substitutions = {
        transmission_line_admittance: transmission_line_admittance_,
        power_ratio: ratio_of_power_,
    }
    result_y1 = Quantity(result_y1.subs(substitutions))
    result_y2 = Quantity(result_y2.subs(substitutions))
    result_y3 = Quantity(result_y3.subs(substitutions))
    result_y4 = Quantity(result_y4.subs(substitutions))
    return (result_y1, result_y2, result_y3, result_y4)
