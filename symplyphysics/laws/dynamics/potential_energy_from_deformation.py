"""
Elastic potential energy from displacement
==========================================

Spring accumulates energy while being deformed. This law is known as the *Hooke's law*.

**Conditions:**

#. The deformation is elastic (reversible).

..
    TODO Move law to ./deformations/
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

elastic_potential_energy = Symbol("elastic_potential_energy", units.energy)
"""
The potential energy of the spring.

Symbol:
    :code:`U`
"""

stiffness = Symbol("stiffness", units.force / units.length)
"""
The spring's stiffness, or spring constant.

Symbol:
    :code:`k`
"""

displacement = Symbol("displacement", units.length)
"""
The displacement of the spring.

Symbol:
    :code:`x`
"""

law = Eq(elastic_potential_energy, stiffness * displacement**2 / 2)
r"""
:code:`U = k * x^2 / 2`

Latex:
    .. math::
        U = \frac{1}{2} k x^2
"""


@validate_input(stiffness_=stiffness, deformation_=displacement)
@validate_output(elastic_potential_energy)
def calculate_energy(stiffness_: Quantity, deformation_: Quantity) -> Quantity:
    result_energy_expr = solve(law, elastic_potential_energy, dict=True)[0][elastic_potential_energy]
    result_expr = result_energy_expr.subs({
        stiffness: stiffness_,
        displacement: deformation_
    })
    return Quantity(result_expr)
