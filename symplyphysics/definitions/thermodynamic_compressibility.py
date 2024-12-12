"""
Thermodynamic compressibility
=============================

**Compressibility**, or the **coefficient of compressibility**, is a measure of the instantaneous
relative volume change of a fluid or solid as a response to pressure or mean stress change.

**Notes:**

#. This definition is incomplete in the sense that the value of the compressibility coefficient
   depends on whether the process is isentropic or isothermal, hence the partial derivative should
   be taken at either constant entropy or constant temperature.
#. For solids the difference between isentropic and isothermal compressibility coefficients
   is usually negligible.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Compressibility#>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    SymbolNew,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D
from symplyphysics.core.dimensions import any_dimension

compressibility = symbols.thermodynamic_compressibility
"""
:symbols:`thermodynamic_compressibility` of the gas.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` in the gas.
"""

parameters = SymbolNew("q", any_dimension)
"""
Parameters other than :attr:`~pressure` on which the volume function depends.
"""

volume = clone_as_function(symbols.volume, [pressure, parameters])
"""
:symbols:`volume` of the gas as a function of :attr:`~pressure` and :attr:`~parameters`.
"""

definition = Eq(compressibility, -1 * Derivative(volume(pressure, parameters), pressure) / volume(pressure, parameters))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    volume_before_=volume,
    volume_after_=volume,
    pressure_before_=pressure,
    pressure_after_=pressure,
)
@validate_output(compressibility)
def calculate_compressibility(
    volume_before_: Quantity,
    volume_after_: Quantity,
    pressure_before_: Quantity,
    pressure_after_: Quantity,
) -> Quantity:
    # The value of the volume is calculated in the `pressure_after_` point

    volume_function = two_point_function(
        Point2D(pressure_before_, volume_before_),
        Point2D(pressure_after_, volume_after_),
        pressure,
    )
    expr = definition.rhs.subs(volume(pressure, parameters), volume_function).doit()
    result = expr.subs(pressure, pressure_after_)
    return Quantity(result)
