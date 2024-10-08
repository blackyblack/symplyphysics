r"""
Linear coefficient of thermal expansion
=======================================

The coefficient of thermal expansion describes how the size of an object changes with a change in temperature
at constant pressure.

**Notation:**

#. :math:`\left( \frac{\partial}{\partial T} \right)_p` is the derivative w.r.t. temperature 
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

linear_expansion_coefficient = clone_as_symbol(
    symbols.thermal_expansion_coefficient,
    display_symbol="alpha_l",
    display_latex="\\alpha_l",
)
"""
Linear :symbols:`thermal_expansion_coefficient` of the object.
"""

length = clone_as_function(symbols.length, display_symbol="l(T, p)")
"""
:symbols:`length` of the object as a function of temperature
and, indirectly, pressure :math:`p`.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the object.
"""

definition = Eq(
    linear_expansion_coefficient,
    Derivative(length(temperature), temperature) / length(temperature),
)
r"""
:laws:symbol::

.. only:: comment

    Derivatives with one parameter held constant are not very pretty when autogenerated.

Latex:
    .. math::
        \alpha_l = \frac{1}{l(T, p)} \left( \frac{\partial l(T, p)}{\partial T} \right)_p
"""


@validate_input(
    length_before_=length,
    length_after_=length,
    temperature_before_=temperature,
    temperature_after_=temperature,
)
@validate_output(linear_expansion_coefficient)
def calculate_linear_expansion_coefficient(
    length_before_: Quantity,
    length_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
) -> Quantity:
    # The RHS of the equation is calculated at the temperature point after expansion (`temperature_after_`)

    length_function = two_point_function(
        Point2D(temperature_before_, length_before_),
        Point2D(temperature_after_, length_after_),
        temperature,
    )
    result = ((definition.rhs).subs(length(temperature),
        length_function).doit().subs(temperature, temperature_after_))
    return Quantity(result)
