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
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

linear_expansion_coefficient = Symbol("linear_expansion_coefficient", 1 / units.temperature)
r"""
Linear coefficient of thermal expansion of the object.

Symbol:
    :code:`alpha_l`

Latex:
    :math:`\alpha_l`
"""

length = Function("length", units.length)
"""
Length of the object as a function of :attr:`~symplyphysics.symbols.thermodynamics.temperature`.

Symbol:
    :code:`l`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the object.

Symbol:
    :code:`T`
"""

definition = Eq(
    linear_expansion_coefficient,
    Derivative(length(temperature), temperature) / length(temperature),
)
r"""
:code:`alpha_l = 1 / l(T) * (dl(T)/dT)_p`

Latex:
    .. math::
        \alpha_l = \frac{1}{l(T)} \left( \frac{d l(T)}{d T} \right)_p
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
