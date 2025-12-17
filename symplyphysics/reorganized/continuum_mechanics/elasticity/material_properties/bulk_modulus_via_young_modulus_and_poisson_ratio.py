r"""
Bulk modulus via Young modulus and Poisson ratio
================================================

Suppose a uniform isotropic body is subjected to bulk compression, i.e. forces are applied
to it from all its sides. Then the bulk modulus is the proportionality coefficient between
relative volume change of the body and the pressure inside of it. It is proportional to
the Young modulus of the material and also depends on its Poisson ratio.

**Conditions:**

#. The body is isotropic and uniform.
#. The :ref:`Poisson ratio <poisson-ratio>` :math:`\nu < \frac{1}{2}`,
   since :doc:`elastic energy density <laws.dynamics.deformation.elastic_energy_density_of_bulk_compression_via_pressure>`
   cannot be negative.

**Links:**

#. `Wikipedia, derivable from second part of the equation <https://en.wikipedia.org/wiki/Young%27s_modulus#Usage>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

bulk_modulus = symbols.bulk_modulus
"""
:symbols:`bulk_modulus` of the material.
"""

young_modulus = symbols.young_modulus
"""
:symbols:`young_modulus` of the material.
"""

poisson_ratio = symbols.poisson_ratio
"""
.. _poisson-ratio:

:symbols:`poisson_ratio` of the material.
"""

law = Eq(bulk_modulus, young_modulus / (3 * (1 - 2 * poisson_ratio)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    young_modulus_=young_modulus,
    poisson_ratio_=poisson_ratio,
)
@validate_output(bulk_modulus)
def calculate_bulk_modulus(
    young_modulus_: Quantity,
    poisson_ratio_: float,
) -> Quantity:
    result = law.rhs.subs({
        young_modulus: young_modulus_,
        poisson_ratio: poisson_ratio_,
    })
    return Quantity(result)
