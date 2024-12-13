"""
Chemical potential of ideal gas
===============================

The chemical potential of an ideal gas can be calculated from its temperature, concentration,
and thermal wavelength.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is ideal.

**Links:**

#. Formula on p. 394 of "Statistical Mechanics" by Terrent L. Hill (1987)
"""

from sympy import Eq, log
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

chemical_potential = symbols.chemical_potential
r"""
:symbols:`chemical_potential` of ideal gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

concentration = symbols.number_density
"""
Concentration of the gas, or :symbols:`number_density`.
"""

thermal_wavelength = symbols.thermal_wavelength
"""
:symbols:`thermal_wavelength` of the gas. Also see :doc:`Thermal de Broglie wavelength
<definitions.thermal_de_broglie_wavelength>`.
"""

law = Eq(chemical_potential,
    quantities.boltzmann_constant * temperature * log(concentration * thermal_wavelength**3))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    temperature_=temperature,
    concentration_=concentration,
    thermal_wavelength_=thermal_wavelength,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    temperature_: Quantity,
    concentration_: Quantity,
    thermal_wavelength_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        concentration: concentration_,
        thermal_wavelength: thermal_wavelength_,
    })
    return Quantity(result)
