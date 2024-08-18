r"""
Reduced temperature in Maxwell—Jüttner statistics
=================================================

Reduced temperature is a dimensionless quantity used in the
:doc:`Maxwell—Jüttner distribution <laws.thermodynamics.relativistic.maxwell_juettner_distribution>`
function. It depends on the temperature of the system, as well as on the mass of a particle of
the gas in question.

**Notation:**

#. :math:`k_\text{B}` is the Boltzmann constant.
#. :math:`c` is the speed of light.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

reduced_temperature = Symbol("reduced_temperature", dimensionless)
r"""
Reduced temperature.

Symbol:
    :code:`theta`

Latex:
    :math:`\theta`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the gas.

Symbol:
    :code:`T`
"""

particle_mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of a relativistic particle comprising the gas.

Symbol:
    :code:`m`
"""

law = Eq(
    reduced_temperature,
    (units.boltzmann_constant * temperature) / (particle_mass * units.speed_of_light**2),
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
