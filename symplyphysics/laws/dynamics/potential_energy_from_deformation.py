"""
Potential energy from deformation
=================================

Spring accumulates energy while being deformated. This law is known as the *Hooke's law*.

**Conditions:**

#. The deformation is elastic (reversible).
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

spring_energy = Symbol("spring_energy", units.energy)
"""
The potential energy of the spring.

Symbol:
    U
"""

elastic_koefficient = Symbol("elastic_koefficient", units.force / units.length)
"""
The spring's elasticity, or spring constant.

Symbol:
    k
"""

deformation = Symbol("deformation", units.length)
"""
The deformation of the spring.

Symbol:
    x
"""

law = Eq(spring_energy, elastic_koefficient * deformation**2 / 2)
r"""
U = k * x^2 / 2

Latex:
    .. math:
        E = \frac{1}{2} k x^2
"""


@validate_input(elastic_koefficient_=elastic_koefficient, deformation_=deformation)
@validate_output(spring_energy)
def calculate_energy(elastic_koefficient_: Quantity, deformation_: Quantity) -> Quantity:
    result_energy_expr = solve(law, spring_energy, dict=True)[0][spring_energy]
    result_expr = result_energy_expr.subs({
        elastic_koefficient: elastic_koefficient_,
        deformation: deformation_
    })
    return Quantity(result_expr)
