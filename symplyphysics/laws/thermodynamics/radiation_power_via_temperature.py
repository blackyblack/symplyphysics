r"""
Radiation power via temperature
===============================

One of the methods an object can exchange energy with its environment is via thermal radiation by
emitting or absorbing energy in the form of electromagnetic waves. Also known as the *Stefan—Boltzmann
law* of radiation, it states that the rate of thermal radiation is proportional to the fourth power of
the radiating body's temperature.

**Notation:**

#. :math:`\sigma` (:code:`sigma`) is the Stefan—Boltzmann constant.
"""

from sympy import Eq
from symplyphysics import (
    symbols,
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

radiation_power = Symbol("radiation_power", units.power)
"""
Power of radiation emitted or absorbed by the body.

Symbol:
    :code:`P`
"""

emissivity = Symbol("emissivity", dimensionless)
r"""
Emissivity of the body's material.

Symbol:
    :code:`epsilon`

Latex:
    :math:`\varepsilon`
"""

surface_area = Symbol("surface_area", units.area)
"""
Surface area of the body.

Symbol:
    :code:`A`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the body.
"""

law = Eq(
    radiation_power,
    units.stefan_boltzmann_constant * emissivity * surface_area * temperature**4,
)
r"""
:code:`P = sigma * epsilon * A * T^4`

Latex:
    .. math::
        P = \sigma \varepsilon A T^4
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
