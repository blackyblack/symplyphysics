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
"""

from sympy import Eq, sin, exp, pi
from symplyphysics import (
    dimensionless,
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T_n(x, t)",
    display_latex="T_n")
"""
Solution to the heat equation corresponding to the :math:`n`:sup:`th` mode.
See :symbols:`temperature`.
"""

scaling_coefficient = Symbol("scaling_coefficient", units.temperature)
"""
Scaling coefficient of the solution, see :ref:`Notes <heat_transfer_zero_temperature_solution_coefficient_note>`.

Symbol:
    :code:`B_n`

Latex:
    :math:`B_n`
"""

thermal_diffusivity = Symbol("thermal_diffusivity", units.area / units.time)
r"""
`Thermal diffusivity <https://en.wikipedia.org/wiki/Thermal_diffusivity>`_.

Symbol:
    :code:`chi`

Latex:
    :math:`\chi`
"""

mode_number = Symbol("mode_number", dimensionless, integer=True, positive=True)
"""
Number of the mode, which is a positive integer.

Symbol:
    :code:`n`
"""

maximum_position = Symbol("maximum_position", units.length, positive=True)
"""
Maximum possible position.

Symbol:
    :code:`L`
"""

position = Symbol("position", units.length)
"""
Position, or spatial variable.

Symbol: 
    :code:`x`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(
    temperature,
    scaling_coefficient * sin(mode_number * pi * position / maximum_position) *
    exp(-1 * thermal_diffusivity * (mode_number * pi / maximum_position)**2 * time))
r"""
:code:`T_n(x, t) = B_n * sin(n * pi * x / L) * exp(-1 * chi * (n * pi / L)^2 * t)`

Latex:
    .. math::
        T_n(x, t) = B_n \sin \left( \frac{n \pi x}{L} \right) \exp \left[ -\chi \left( \frac{n \pi}{L} \right)^2 t \right]
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
