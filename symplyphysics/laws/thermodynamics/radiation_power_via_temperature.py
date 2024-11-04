r"""
Radiation power via temperature
===============================

One of the methods an object can exchange energy with its environment is via thermal radiation by
emitting or absorbing energy in the form of electromagnetic waves. Also known as the *Stefan—Boltzmann
law* of radiation, it states that the rate of thermal radiation is proportional to the fourth power of
the radiating body's temperature.

**Notation:**

#. :quantity_notation:`stefan_boltzmann_constant`.
"""

from sympy import Eq
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)

radiation_power = symbols.power
"""
:symbols:`power` of radiation emitted or absorbed by the body.
"""

emissivity = symbols.emissivity
"""
:symbols:`emissivity` of the body's material.
"""

surface_area = symbols.area
"""
Surface :symbols:`area` of the body.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the body.
"""

law = Eq(
    radiation_power,
    quantities.stefan_boltzmann_constant * emissivity * surface_area * temperature**4,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    surface_emissivity_=emissivity,
    surface_area_=surface_area,
    temperature_=temperature,
)
@validate_output(radiation_power)
def calculate_energy_radiation_rate(
    surface_emissivity_: float,
    surface_area_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        emissivity: surface_emissivity_,
        surface_area: surface_area_,
        temperature: temperature_,
    })
    return Quantity(result)
