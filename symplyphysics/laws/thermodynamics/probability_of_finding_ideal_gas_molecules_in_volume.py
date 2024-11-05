r"""
Probability of finding ideal gas molecules in volume
====================================================

For ideal gas molecules, we can calculate the probability of N molecules occupying a certain region
in space, say, a container with gas particles.

**Conditions:**

#. The gas is ideal.
#. :math:`V \le V_0`
#. There are no external fields acting on the particles.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    clone_as_symbol,
    symbols,
)
from symplyphysics.core.symbols.probability import Probability

probability = clone_as_symbol(symbols.probability)
"""
:symbols:`probability` of finding :attr:`~particle_count` gas particles in :attr:`~partial_volume`.
"""

total_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
Total :symbols:`volume` of the region.
"""

partial_volume = symbols.volume
"""
:symbols:`volume` of part of the container where we are looking for the particles.
"""

particle_count = clone_as_symbol(symbols.particle_count, integer=True)
"""
:symbols:`particle_count` in :attr:`~partial_volume`.
"""

law = Eq(probability, (partial_volume / total_volume)**particle_count)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    total_volume_=total_volume,
    partial_volume_=partial_volume,
    particle_count_=particle_count,
)
@validate_output(probability)
def calculate_probability(
    total_volume_: Quantity,
    partial_volume_: Quantity,
    particle_count_: int,
) -> Probability:
    if partial_volume_.scale_factor > total_volume_.scale_factor:
        raise ValueError("Partial volume must not exceed total volume")

    result = law.rhs.subs({
        total_volume: total_volume_,
        partial_volume: partial_volume_,
        particle_count: particle_count_,
    })

    return Probability(convert_to_float(result))
