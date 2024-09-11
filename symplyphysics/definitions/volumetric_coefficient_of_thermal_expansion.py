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
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

volumetric_expansion_coefficient = Symbol("volumetric_expansion_coefficient", 1 / units.temperature)
r"""
Volumetric coefficient of thermal expansion.

Symbol:
    :code:`alpha_V`

Latex:
    :math:`\alpha_V`
"""

volume = Function("volume", units.volume)
"""
Volume of the body as a function of temperature and pressure.

Symbol:
    :code:`V(T, p)`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the body.
"""

definition = Eq(
    volumetric_expansion_coefficient,
    Derivative(volume(temperature), temperature) / volume(temperature),
)
r"""
:code:`alpha_V = 1 / V(T, p) * Derivative(V(T, p), T)`

Lambda:
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
