r"""
Volumetric coefficient of thermal expansion
===========================================

The *coefficient of thermal expansion* describes how the size of an object changes with a change in temperature
at constant pressure.

**Conditions:**

#. Pressure must be constant during the expansion process.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    clone_as_function,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

volumetric_expansion_coefficient = clone_as_symbol(
    symbols.thermal_expansion_coefficient,
    display_symbol="alpha_V",
    display_latex="\\alpha_V",
)
"""
Volumetric :symbols:`thermal_expansion_coefficient`.
"""

volume = clone_as_function(symbols.volume, display_symbol="V(T, p)")
"""
:symbols:`volume` of the body as a function of temperature and pressure.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the body.
"""

definition = Eq(
    volumetric_expansion_coefficient,
    Derivative(volume(temperature), temperature) / volume(temperature),
)
r"""
:laws:symbol::

.. only:: comment

    Manual formula due to partial derivative not rendered by auto-generation

Latex:
    .. math::
        \alpha_V = \frac{1}{V(T, p)} \left( \frac{\partial V}{\partial T} \right)_p
"""


@validate_input(
    volume_before_=volume,
    volume_after_=volume,
    temperature_before_=temperature,
    temperature_after_=temperature,
)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(
    volume_before_: Quantity,
    volume_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
) -> Quantity:
    # The RHS of the equation is calculated at the temperature point after expansion (`temperature_after_`)

    volume_function = two_point_function(
        Point2D(temperature_before_, volume_before_),
        Point2D(temperature_after_, volume_after_),
        temperature,
    )
    result = ((definition.rhs).subs(volume(temperature),
        volume_function).doit().subs(temperature, temperature_after_))
    return Quantity(result)
