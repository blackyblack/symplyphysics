"""
Chemical potential of ideal gas
===============================

The chemical potential of an ideal gas can be calculated from its temperature, concentration,
and thermal wavelength.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is ideal.
"""

from sympy import Eq, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
Chemical potential of ideal gas.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

concentration = Symbol("concentration", 1 / units.volume)
"""
Concentration of the gas, or particle count per unit volume.

Symbol:
    :code:`n`
"""

thermal_wavelength = Symbol("thermal_wavelength", units.length)
r"""
:doc:`Thermal de Broglie wavelength <definitions.thermal_de_broglie_wavelength>` of the gas.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

law = Eq(chemical_potential,
    quantities.boltzmann_constant * temperature * log(concentration * thermal_wavelength**3))
r"""
:code:`mu = k_B * T * log(n * lambda^3)`

Latex:
    .. math::
        \mu = k_\text{B} T \log \left( n \lambda^3 \right)
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
