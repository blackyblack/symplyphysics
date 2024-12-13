"""
Enthalpy derivative via volume derivative
=========================================

Isothermal derivative of enthalpy w.r.t. pressure can be found using volume as a function
of temperature and pressure.

**Conditions:**

#. Works for an infinitesimal quasi-static isothermal process.

**Links:**

#. `Wikipedia, see third table <https://en.wikipedia.org/wiki/Table_of_thermodynamic_equations#Maxwell's_relations>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.line import two_point_function, Point2D

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

enthalpy = clone_as_function(symbols.enthalpy, [temperature, pressure])
"""
:symbols:`enthalpy` of the system.
"""

volume = clone_as_function(symbols.volume, [temperature, pressure])
"""
:symbols:`volume` of the system.
"""

law = Eq(
    Derivative(enthalpy(temperature, pressure), pressure),
    volume(temperature, pressure) -
    temperature * Derivative(volume(temperature, pressure), temperature))
"""
:laws:symbol::

:laws:latex::
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
