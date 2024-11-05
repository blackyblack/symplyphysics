"""
Relative humidity is ratio of vapor pressure
============================================

**Relative humidity** is the ratio of actual water vapor pressure :math:`p` to saturation vapor
pressure :math:`p_0` at a given temperature. **Saturation vapor** is vapor in dynamic
equilibrium with a liquid or solid of the same composition in a closed system. In such
a system, the number of evaporating molecules is equal to the number of condensing
molecules per unit time.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

relative_humidity = symbols.relative_humidity
"""
:symbols:`relative_humidity` of air.
"""

actual_vapor_pressure = symbols.pressure
"""
Partial :symbols:`pressure` of water vapor in the medium, representing how much water vapor is
actually in the air.
"""

saturation_vapor_pressure = clone_as_symbol(
    symbols.pressure,
    display_symbol="p_s",
    display_latex="p_\\text{s}",
)
"""
Equilibrium :symbols:`pressure` of water vapor above a flat surface of liquid water or solid ice,
representing how much water vapor the air could potentially contain at a given temperature.
"""

law = Eq(relative_humidity, actual_vapor_pressure / saturation_vapor_pressure)
"""
:laws:symbol::

:laws:latex::
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
