r"""
Thermodynamic compressibility
=============================

*Compressibility*, or the *coefficient of compressibility*, is a measure of the instantaneous
relative volume change of a fluid or solid as a response to pressure or mean stress change.

**Notation:**

#. :math:`\frac{\partial}{\partial p}` denotes a partial derivative w.r.t. pressure.

**Notes:**

#. This definition is incomplete in the sense that the value of the compressibility coefficient
   depends on whether the process is isentropic or isothermal, hence the partial derivative should
   be taken at either constant entropy or constant temperature.
#. For solids the difference between isentropic and isothermal compressibility coefficients
   is usually negligible.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

compressibility = Symbol("compressibility", 1 / units.pressure)
r"""
Compressibility of the gas.

Symbol:
    :code:`beta`

Latex:
    :math:`\beta`
"""

volume = Function("volume", units.volume)
"""
Volume of the gas as a function of pressure and other parameters.

Symbol:
    :code:`V`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure in the gas.

Symbol:
    :code:`p`
"""

definition = Eq(compressibility, -1 * Derivative(volume(pressure), pressure) / volume(pressure))
r"""
:code:`beta = -1 / V * (dV/dp)`

Latex:
    .. math::
        \beta = - \frac{1}{V} \left( \frac{\partial V}{\partial p} \right)
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
    expr = definition.rhs.subs(volume(pressure), volume_function).doit()
    result = expr.subs(pressure, pressure_after_)
    return Quantity(result)
