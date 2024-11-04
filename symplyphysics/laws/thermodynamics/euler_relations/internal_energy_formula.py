"""
Internal energy formula
=======================

The formula of the internal energy differential can be integrated using the Euler's theorem on
homogeneous functions to get the following expression.

**Notes:**

#. This formula words for a single-component system. For a multi-component system replace the
   product of chemical potential and particle count with a sum over each type of components.
"""

from sympy import Eq
from symplyphysics import Quantity, validate_input, validate_output, symbols

internal_energy = symbols.internal_energy
"""
:symbols:`internal_energy` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` in the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

chemical_potential = symbols.chemical_potential
"""
:symbols:`chemical_potential` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

law = Eq(internal_energy,
    temperature * entropy - pressure * volume + chemical_potential * particle_count)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from Euler's theorem and internal energy differential


@validate_input(
    temperature_=temperature,
    entropy_=entropy,
    pressure_=pressure,
    volume_=volume,
    chemical_potential_=chemical_potential,
    particle_count_=particle_count,
)
@validate_output(internal_energy)
def calculate_internal_energy(
    temperature_: Quantity,
    entropy_: Quantity,
    pressure_: Quantity,
    volume_: Quantity,
    chemical_potential_: Quantity,
    particle_count_: int,
) -> Quantity:
    # pylint: disable=too-many-arguments, too-many-positional-arguments
    result = law.rhs.subs({
        temperature: temperature_,
        entropy: entropy_,
        pressure: pressure_,
        volume: volume_,
        chemical_potential: chemical_potential_,
        particle_count: particle_count_,
    })
    return Quantity(result)
