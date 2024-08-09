"""
Mean free path of random motion
===============================

The mean free path of a molecule in random motion is its average path length between collisions.

**Conditions:**

#. Spherical model of molecules is assumed.
"""

from sympy import Eq, sqrt, pi
from symplyphysics import (
    Quantity,
    Symbol,
    units,
    validate_input,
    validate_output,
)

mean_free_path = Symbol("mean_free_path", units.length, positive=True)
r"""
Mean free path estimate of molecules.

Symbol:
    :code:`lambda`

Latex:
    :math:`\lambda`
"""

molecular_diameter = Symbol("molecular_diameter", units.length, positive=True)
"""
Diameter of molecules.

Symbol:
    :code:`d`
"""

number_density = Symbol("number_density", 1 / units.volume)
"""
:doc:`Number density <definitions.number_density_is_number_of_objects_per_unit_volume>` of the system.

Symbol:
    :code:`n`
"""

law = Eq(mean_free_path, 1 / (sqrt(2) * pi * molecular_diameter**2 * number_density))
r"""
:code:`lambda = 1 / (sqrt(2) * pi * d^2 * n)`

Latex:
    .. math::
        \lambda = \frac{1}{\sqrt{2} \pi d^2 n}
"""


@validate_input(
    molecular_diameter_=molecular_diameter,
    number_density_=number_density,
)
@validate_output(mean_free_path)
def calculate_mean_free_path(
    molecular_diameter_: Quantity,
    number_density_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        molecular_diameter: molecular_diameter_,
        number_density: number_density_,
    })
    return Quantity(result)
