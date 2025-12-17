"""
Mean free path of random motion
===============================

The mean free path of a molecule in random motion is its average path length between collisions.

**Conditions:**

#. Spherical model of molecules is assumed.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Mean_free_path#Kinetic_theory_of_gases>`__.

..
    NOTE: replace `pi * d^2` with cross-section?
"""

from sympy import Eq, sqrt, pi
from symplyphysics import Quantity, validate_input, validate_output, symbols

mean_free_path = symbols.mean_free_path
"""
:symbols:`mean_free_path` estimate of molecules.
"""

molecular_diameter = symbols.diameter
"""
:symbols:`diameter` of molecules.
"""

number_density = symbols.number_density
"""
:symbols:`number_density` of the system.
"""

law = Eq(mean_free_path, 1 / (sqrt(2) * pi * molecular_diameter**2 * number_density))
"""
:laws:symbol::

:laws:latex::
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
