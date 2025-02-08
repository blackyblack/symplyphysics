"""
Impedance of Wilkinson microstrip divider
=========================================

The Wilkinson divider is a device designed to divide the power of a microwave signal
into two output ports. Different sections of the divider consist of a microstrip line of
different widths. There are four such sections in total and each has its own impedance.

.. image:: https://habrastorage.org/getpro/habr/upload_files/c24/031/52e/c2403152e2b320ab1c4c44f970dee1f2.gif
    :width: 400px
    :align: center

..
    TODO: find link
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

first_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="1")
"""
:symbols:`electrical_impedance` of the first section.
"""

second_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="2")
"""
:symbols:`electrical_impedance` of the second section.
"""

third_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="3")
"""
:symbols:`electrical_impedance` of the third section.
"""

fourth_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="4")
"""
:symbols:`electrical_impedance` of the fourth section.
"""

transmission_line_impedance = clone_as_symbol(symbols.electrical_impedance, subscript="0")
"""
:symbols:`electrical_impedance` of the transmission line to which the divider is
connected.
"""

power_ratio = SymbolNew("k", dimensionless)
"""
Ratio of the :symbols:`power` at the outputs of the divider.
"""

law = Eq(
    Matrix([first_impedance, second_impedance, third_impedance, fourth_impedance]),
    transmission_line_impedance * Matrix([
        sqrt(power_ratio * (1 + power_ratio**2)),
        sqrt((1 + power_ratio**2) / power_ratio**3),
        sqrt(power_ratio),
        1 / sqrt(power_ratio),
    ]))
r"""
..
    Code printers don't work with matrices yet.

    :code:`[Z_1, Z_2, Z_3, Z_4] = Z_0 * [sqrt(k * (1 + k^2)), sqrt((1 + k^2) / k^3), sqrt(k), 1 / sqrt(k)]`

    Latex:
        .. math::
            \begin{pmatrix} Z_1 \\ Z_2 \\ Z_3 \\ Z_4 \end{pmatrix}
            = Z_0 \begin{pmatrix}
                \sqrt{k (1 + k^2)} \\
                \sqrt{\frac{1 + k^2}{k^3}} \\
                \sqrt{k} \\
                \frac{1}{\sqrt{k}}
   
            \end{pmatrix}

:laws:symbol::

:laws:latex::
"""


@validate_input(characteristic_resistance_=transmission_line_impedance,
    ratio_of_power_=power_ratio)
@validate_output(units.impedance)
def calculate_impedances(characteristic_resistance_: Quantity,
    ratio_of_power_: float) -> tuple[Quantity, Quantity, Quantity, Quantity]:
    result = solve(law, [first_impedance, second_impedance, third_impedance, fourth_impedance],
        dict=True)[0]
    result_z1 = result[first_impedance]
    result_z2 = result[second_impedance]
    result_z3 = result[third_impedance]
    result_z4 = result[fourth_impedance]
    substitutions = {
        transmission_line_impedance: characteristic_resistance_,
        power_ratio: ratio_of_power_,
    }
    result_z1 = Quantity(result_z1.subs(substitutions))
    result_z2 = Quantity(result_z2.subs(substitutions))
    result_z3 = Quantity(result_z3.subs(substitutions))
    result_z4 = Quantity(result_z4.subs(substitutions))
    return (result_z1, result_z2, result_z3, result_z4)
