"""
Speed of sound via bulk modulus and density
===========================================

Sound waves are longitudinal mechanical waves that can travel through any type of material,
be it solid, liquid, or gas. The phase velocity of the a sound wave depends on the bulk modulus
of the medium it is traveling in and its density.
"""

from sympy import Eq, sqrt
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
)

phase_velocity = Symbol("phase_velocity", units.velocity, positive=True)
"""
Phase velocity of sound wave.

Symbol:
    :code:`v`
"""

bulk_modulus = Symbol("bulk_modulus", units.pressure, real=True)
"""
Bulk modulus of the medium.

Symbol:
    :code:`B`

..
    TODO add reference to laws.hydro.hydraulic_stress_is_bulk_modulus_times_strain
"""

density = Symbol("density", units.mass / units.volume, positive=True)
r"""
Density of the medium.

Symbol:
    :code:`rho`

Latex:
    :math:`\rho`
"""

law = Eq(phase_velocity, sqrt(bulk_modulus / density))
r"""
:code:`v = sqrt(B / rho)`

Latex:
    .. math::
        v = \sqrt{\frac{B}{\rho}}
"""


@validate_input(
    bulk_modulus_=bulk_modulus,
    density_=density,
)
@validate_output(phase_velocity)
def calculate_speed_of_sound(bulk_modulus_: Quantity, density_: Quantity) -> Quantity:
    result = law.rhs.subs({
        bulk_modulus: bulk_modulus_,
        density: density_,
    })
    return Quantity(result)
