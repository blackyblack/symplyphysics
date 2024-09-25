"""
Elastic energy density of compression via strain
================================================

Volumetric density of the elastic energy of a body is proportional to
its Young's modulus and the square of its strain. The
:doc:`Hooke's law <laws.dynamics.deformation.tensile_stress_is_youngs_modulus_times_strain>`
can be used to obtain analogous forms of this law.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

elastic_energy_density = symbols.energy_density
"""
Elastic energy of the deformed body per unit of its volume. See :symbols:`energy_density`
"""

young_modulus = symbols.young_modulus
"""
:symbols:`young_modulus` of the body's material.
"""

engineering_normal_strain = symbols.engineering_normal_strain
"""
:symbols:`engineering_normal_strain` of the deformed body.
"""

law = Eq(elastic_energy_density, young_modulus * engineering_normal_strain**2 / 2)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    young_modulus_=young_modulus,
    engineering_normal_strain_=engineering_normal_strain,
)
@validate_output(elastic_energy_density)
def calculate_elastic_energy_density(
    young_modulus_: Quantity,
    engineering_normal_strain_: float,
) -> Quantity:
    result = law.rhs.subs({
        young_modulus: young_modulus_,
        engineering_normal_strain: engineering_normal_strain_,
    })
    return Quantity(result)
