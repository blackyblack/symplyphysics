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
    Symbol,
    dimensionless,
    units,
    validate_input,
    validate_output,
    convert_to_float,
)
from symplyphysics.core.symbols.probability import Probability

probability = Symbol("probability", dimensionless)
"""
Probability of finding :math:`N` gas particles in volume :math:`V`.

Symbol:
    :code:`P`
"""

total_volume = Symbol("total_volume", units.volume)
"""
Total volume of the region.

Symbol:
    :code:`V0`

Latex:
    :math:`V_0`
"""

partial_volume = Symbol("partial_volume", units.volume)
"""
Volume of part of the container where we are looking for the particles.

Symbol:
    :code:`V`
"""

particle_count = Symbol("particle_count", dimensionless, integer=True)
"""
Number of gas particles in :math:`V`.

Symbol:
    :code:`N`
"""

law = Eq(probability, (partial_volume / total_volume)**particle_count)
r"""
:code:`P = (V / V0)^N`

Latex:
    .. math::
        P = \left( \frac{V}{V_0} \right)^N
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
