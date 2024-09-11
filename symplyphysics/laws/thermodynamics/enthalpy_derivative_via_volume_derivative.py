"""
Enthalpy derivative via volume derivative
=========================================

Isothermal derivative of enthalpy w.r.t. pressure can be found using volume as a function
of temperature and pressure.

**Conditions:**

#. Works for an infinitesimal quasi-static isothermal process.
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
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.line import two_point_function, Point2D

enthalpy = Function("enthalpy", units.energy)
"""
Enthalpy of the system.

Symbol:
    :code:`H(T, p)`
"""

temperature = symbols.temperature
"""
:attr:`~symplyphysics.symbols.temperature` of the system.
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

volume = Function("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

law = Eq(
    Derivative(enthalpy(temperature, pressure), pressure),
    volume(temperature, pressure) -
    temperature * Derivative(volume(temperature, pressure), temperature))
r"""
:code:`Derivative(H(T, p), p) = V(T, p) - T * Derivative(V(T, p), T)`

Latex:
    .. math::
        \left( \frac{\partial H}{\partial p} \right)_T = V(T, p) - T \left( \frac{\partial V}{\partial T} \right)_p
"""

# TODO: Derive from differential of enthalpy and Maxwell relation


@validate_input(
    volume_before_=volume,
    volume_after_=volume,
    temperature_before_=temperature,
    temperature_after_=temperature,
)
@validate_output(units.energy / units.pressure)
def calculate_enthalpy_derivative(
    volume_before_: Quantity,
    volume_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
) -> Quantity:
    volume_ = two_point_function(
        Point2D(temperature_before_, volume_before_),
        Point2D(temperature_after_, volume_after_),
        temperature,
    )
    result = law.rhs.subs(volume(temperature, pressure), volume_).doit()

    # the result does not depend on temperature
    assert expr_equals(result.diff(temperature), 0)

    return Quantity(result)
