r"""
Reduced temperature in Maxwell—Jüttner statistics
=================================================

Reduced temperature is a dimensionless quantity used in the
:doc:`Maxwell—Jüttner distribution <laws.thermodynamics.relativistic.maxwell_juettner_distribution>`
function. It depends on the temperature of the system, as well as on the mass of a particle of
the gas in question.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.
#. :quantity_notation:`speed_of_light`.
"""

from sympy import Eq
from symplyphysics import (
    dimensionless,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
    quantities,
)

reduced_temperature = Symbol("reduced_temperature", dimensionless)
r"""
Reduced temperature.

Symbol:
    :code:`theta`

Latex:
    :math:`\theta`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

particle_mass = symbols.mass
"""
:symbols:`mass` of a relativistic particle comprising the gas.
"""

law = Eq(
    reduced_temperature,
    (quantities.boltzmann_constant * temperature) / (particle_mass * quantities.speed_of_light**2),
)
r"""
:code:`theta = (k_B * T) / (m * c^2)`

Latex:
    .. math::
        \theta = \frac{k_\text{B} T}{m c^2}
"""


@validate_input(
    temperature_=temperature,
    mass_=particle_mass,
)
@validate_output(reduced_temperature)
def calculate_reduced_temperature(
    temperature_: Quantity,
    mass_: Quantity,
) -> float:
    result = law.rhs.subs({
        temperature: temperature_,
        particle_mass: mass_,
    })
    return convert_to_float(result)
