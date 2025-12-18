"""
Elastic potential energy from displacement
==========================================

Spring accumulates energy while being deformed. This law is known as the **Hooke's law**.

**Conditions:**

#. The deformation is elastic (reversible).

**Links:**

#. `Wikipedia, first formula <https://en.wikipedia.org/wiki/Elastic_energy#>`__.

..
    TODO Move law to ./deformations/
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

elastic_potential_energy = symbols.potential_energy
"""
The :symbols:`potential_energy` of the spring.
"""

stiffness = symbols.stiffness
"""
The spring's :symbols:`stiffness`, or spring constant.
"""

displacement = symbols.euclidean_distance
"""
The displacement of the spring, or the :symbols:`euclidean_distance` between the initial position
and the rest position.
"""

law = Eq(elastic_potential_energy, stiffness * displacement**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(stiffness_=stiffness, deformation_=displacement)
@validate_output(elastic_potential_energy)
def calculate_energy(stiffness_: Quantity, deformation_: Quantity) -> Quantity:
    result_energy_expr = solve(law, elastic_potential_energy,
        dict=True)[0][elastic_potential_energy]
    result_expr = result_energy_expr.subs({stiffness: stiffness_, displacement: deformation_})
    return Quantity(result_expr)
