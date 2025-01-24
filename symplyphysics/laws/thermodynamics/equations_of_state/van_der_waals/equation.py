"""
Van der Waals equation
======================

To more accurately describe the behavior of real gases at low temperatures,
a Van der Waals gas model was created, taking into account the forces of intermolecular interaction.
In this model, internal energy becomes a function not only of temperature, but also of molar volume.

The Van der Waals equation is one of the well-known approximate equations of state describing
the properties of a real gas, having a compact form and taking
into account the main characteristics of a gas with intermolecular interaction.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Van_der_Waals_equation#>`__.
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, Quantity, validate_input, validate_output,
    quantities)

pressure = symbols.pressure
"""
:symbols:`pressure` in the van der Waals fluid.
"""

molar_volume = symbols.molar_volume
"""
:symbols:`molar_volume` of the van der Waals fluid.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the van der Waals fluid.
"""

attractive_forces_parameter = symbols.attractive_forces_parameter
"""
:symbols:`attractive_forces_parameter`.
"""

excluded_volume_parameter = symbols.excluded_volume_parameter
"""
:symbols:`excluded_volume_parameter`.
"""

law = Eq(
    (pressure + attractive_forces_parameter / molar_volume**2) *
    (molar_volume - excluded_volume_parameter),
    quantities.molar_gas_constant * temperature,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    molar_volume_=molar_volume,
    temperature_=temperature,
    bonding_forces_parameter_=attractive_forces_parameter,
    molecules_volume_parameter_=excluded_volume_parameter,
)
@validate_output(pressure)
def calculate_pressure(
    molar_volume_: Quantity,
    temperature_: Quantity,
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    solved = solve(law, pressure, dict=True)[0][pressure]
    result_expr = solved.subs({
        molar_volume: molar_volume_,
        temperature: temperature_,
        attractive_forces_parameter: bonding_forces_parameter_,
        excluded_volume_parameter: molecules_volume_parameter_
    })
    return Quantity(result_expr)
