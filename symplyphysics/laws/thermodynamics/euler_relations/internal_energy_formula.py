"""
Internal energy formula
=======================

The formula of the internal energy differential can be integrated using the Euler's theorem on
homogeneous functions to get the following expression.

**Notes:**

#. This formula words for a single-component system. For multi-component system replace the
   product of chemical potential and particle count with a sum over each type of components.
"""

from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

internal_energy = Symbol("internal_energy", units.energy)
"""
Internal energy of the system.

Symbol:
    :code:`U`
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure in the system.

Symbol:
    :code:`p`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
Chemical potential of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

particle_count = Symbol("particle_count", dimensionless)
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

law = Eq(internal_energy,
    temperature * entropy - pressure * volume + chemical_potential * particle_count)
r"""
:code:`U = T * S - p * V + mu * N`

Latex:
    .. math::
        U = T S - p V + \mu N
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
# pylint: disable=too-many-arguments
def calculate_internal_energy(
    temperature_: Quantity,
    entropy_: Quantity,
    pressure_: Quantity,
    volume_: Quantity,
    chemical_potential_: Quantity,
    particle_count_: int,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy: entropy_,
        pressure: pressure_,
        volume: volume_,
        chemical_potential: chemical_potential_,
        particle_count: particle_count_,
    })
    return Quantity(result)
