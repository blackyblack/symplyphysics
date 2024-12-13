"""
Speed of sound via bulk modulus and density
===========================================

Sound waves are longitudinal mechanical waves that can travel through any type of material,
be it solid, liquid, or gas. The phase velocity of the a sound wave depends on the bulk modulus
of the medium it is traveling in and its density.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Speed_of_sound#Equations>`__.
"""

from sympy import Eq, sqrt
from symplyphysics import Quantity, validate_input, validate_output, symbols

phase_velocity = symbols.phase_speed
"""
:symbols:`phase_speed` of sound wave.
"""

bulk_modulus = symbols.bulk_modulus
"""
:symbols:`bulk_modulus` of the medium. Also see :doc:`laws.hydro.hydraulic_stress_is_bulk_modulus_times_strain`.
"""

density = symbols.density
"""
:symbols:`density` of the medium.
"""

law = Eq(phase_velocity, sqrt(bulk_modulus / density))
"""
:laws:symbol::

:laws:latex::
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
