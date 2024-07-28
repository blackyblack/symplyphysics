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
    units,
    dimensionless,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

elastic_energy_density = Symbol("elastic_energy_density", units.energy / units.volume)
"""
Elastic energy of the deformed body per unit of its volume.

Symbol:
    :code:`u`
"""

young_modulus = Symbol("young_modulus", units.pressure)
"""
Young's modulus of the body's material.

Symbol:
    :code:`E`
"""

engineering_normal_strain = Symbol("engineering_normal_strain", dimensionless)
"""
Engineering normal strain of the deformed body.

Symbol:
    :code:`e`
"""

law = Eq(elastic_energy_density, young_modulus * engineering_normal_strain**2 / 2)
r"""
:code:`u = E * e^2 / 2`

Latex:
    .. math::
        u = \frac{1}{2} E e^2
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
