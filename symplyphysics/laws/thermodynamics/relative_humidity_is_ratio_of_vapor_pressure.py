"""
Relative humidity is ratio of vapor pressure
============================================

*Relative humidity* is the ratio of actual water vapor pressure :math:`p` to saturation vapor
pressure :math:`p_0` at a given temperature. *Saturation vapor* is vapor in dynamic
equilibrium with a liquid or solid of the same composition in a closed system. In such
a system, the number of evaporating molecules is equal to the number of condensing
molecules per unit time.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless)

relative_humidity = Symbol("relative_humidity", dimensionless)
r"""
Relative humidity of air.

Symbol:
    :code:`phi`

Latex:
    :math:`\varphi`
"""

actual_vapor_pressure = Symbol("actual_vapor_pressure", units.pressure)
"""
*Partial pressure of water vapor* in the medium, representing how much water vapor is
actually in the air.

Symbol:
    :code:`p`
"""

saturation_vapor_pressure = Symbol("saturation_vapor_pressure", units.pressure)
r"""
*Equilibrium vapor pressure* of water above a flat surface of liquid water or solid ice,
representing how much water vapor the air could potentially contain at a given temperature.

Symbol:
    :code:`p_s`

Latex:
    :math:`p_\text{s}`
"""

law = Eq(relative_humidity, actual_vapor_pressure / saturation_vapor_pressure)
r"""
:code:`phi = p / p_s`

Latex:
    .. math::
        \varphi = \frac{p}{p_\text{s}}
"""


@validate_input(water_vapor_pressure_=actual_vapor_pressure,
    saturated_vapor_pressure_=saturation_vapor_pressure)
@validate_output(relative_humidity)
def calculate_relative_humidity(water_vapor_pressure_: Quantity,
    saturated_vapor_pressure_: Quantity) -> Quantity:
    result_expr = solve(law, relative_humidity, dict=True)[0][relative_humidity]
    result_expr = result_expr.subs({
        actual_vapor_pressure: water_vapor_pressure_,
        saturation_vapor_pressure: saturated_vapor_pressure_
    })
    return Quantity(result_expr)
