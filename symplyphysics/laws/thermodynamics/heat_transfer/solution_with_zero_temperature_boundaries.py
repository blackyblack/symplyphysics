r"""
Solution with zero temperature boundaries
=========================================

The heat equation coupled with a boundary condition can be solved to get a unique solution.
In this boundary-value problem the heat transfer within a thin rod is observed, the temperature
on both ends being zero. This restricts the number of functions :math:`f(x)` that can satisfy
the boundary equation :math:`f(x) = T(x, 0)`.

**Notes:**

#. :math:`f(x)` represents initial spatial distribution of temperature.

    .. _heat_transfer_zero_temperature_solution_coefficient_note:

#. Values :math:`B_n` are found using the boundary condition :math:`f(x) = \sum_n T_n(x, 0)`
   with the help of the Fourier method.
#. The total solution :math:`T(x, t) = \sum_n T_n(x, t)`.

**Conditions:**

#. Position :math:`x \in [0, L]`.
#. Temperature on both ends is zero: :math:`T_n(0, t) = 0`, :math:`T_n(L, t) = 0`

**Links:**

#. `Paul's Online Math Notes, Example 1 <https://tutorial.math.lamar.edu/classes/de/solvingheatequation.aspx>`__.
"""

from sympy import Eq, sin, exp, pi
from symplyphysics import (
    units,
    SymbolNew,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

temperature = clone_as_symbol(symbols.temperature, display_symbol="T_n(x, t)", display_latex="T_n")
"""
Solution to the heat equation corresponding to the :math:`n`:sup:`th` mode.
See :symbols:`temperature`.
"""

scaling_coefficient = SymbolNew("B_n", units.temperature)
"""
Scaling coefficient of the solution, see :ref:`Notes <heat_transfer_zero_temperature_solution_coefficient_note>`.
"""

thermal_diffusivity = symbols.thermal_diffusivity
"""
:symbols:`thermal_diffusivity`.
"""

mode_number = symbols.positive_number
"""
Number of the mode. See :symbols:`positive_number`.
"""

maximum_position = clone_as_symbol(symbols.position, subscript="\\text{max}")
"""
Maximum possible :symbols:`position`.
"""

position = symbols.position
"""
:symbols:`position`, or spatial variable.
"""

time = symbols.time
"""
:symbols:`time`.
"""

law = Eq(
    temperature,
    scaling_coefficient * sin(mode_number * pi * position / maximum_position) *
    exp(-1 * thermal_diffusivity * (mode_number * pi / maximum_position)**2 * time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    coefficient_=scaling_coefficient,
    thermal_diffusivity_=thermal_diffusivity,
    mode_number_=mode_number,
    maximum_position_=maximum_position,
    position_=position,
    time_=time,
)
@validate_output(temperature)
def calculate_temperature(
    coefficient_: Quantity,
    thermal_diffusivity_: Quantity,
    mode_number_: int,
    maximum_position_: Quantity,
    position_: Quantity,
    time_: Quantity,
) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    if maximum_position_.scale_factor <= 0:
        raise ValueError("maximum position must be positive")
    if position_.scale_factor < 0:
        raise ValueError("position must be non-negative")
    if position_.scale_factor > maximum_position_.scale_factor:
        raise ValueError("position must be not greater than the maximum position")

    result = law.rhs.subs({
        scaling_coefficient: coefficient_,
        thermal_diffusivity: thermal_diffusivity_,
        mode_number: mode_number_,
        maximum_position: maximum_position_,
        position: position_,
        time: time_,
    })
    return Quantity(result)
