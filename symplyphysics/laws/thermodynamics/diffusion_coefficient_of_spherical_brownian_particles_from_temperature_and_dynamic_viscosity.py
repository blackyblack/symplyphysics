r"""
Diffusion coefficient of spherical Brownian particles from temperature and dynamic viscosity
============================================================================================

Brownian motion is the random motion of microscopic visible suspended particles of a solid
substance in a liquid or gas caused by the thermal motion of particles of a liquid or gas.
The *diffusion coefficient* is a quantitative characteristic of the diffusion rate, equal to
the amount of matter passing per unit time through a section of a unit area as a result of
the thermal motion of molecules with a concentration gradient equal to one (corresponding to
a change from :math:`1 \frac{\text{mol}}{\text{L}}` to :math:`0 \frac{\text{mol}}{\text{L}}`
per unit length). The diffusion coefficient is determined by the properties of the medium and
the type of diffusing particles. This law is also known as the *Stokes—Einstein—Sutherland relation*.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.
#. :quantity_notation:`avogadro_constant`.

**Conditions:**

#. Particle displacements are equally likely in any direction.
#. The inertia of a Brownian particle can be neglected compared to the influence of friction forces.
#. Particles are spherical.
#. Low Reynolds number, i.e. non-turbulent flow.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Einstein_relation_(kinetic_theory)>`__.
"""

from sympy import Eq, solve, pi
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    quantities,
)

diffusion_coefficient = symbols.diffusion_coefficient
"""
:symbols:`diffusion_coefficient` of the particles.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

particle_radius = symbols.radius
"""
:symbols:`radius` of the particles.
"""

dynamic_viscosity = symbols.dynamic_viscosity
"""
:symbols:`dynamic_viscosity` of the particles.
"""

law = Eq(
    diffusion_coefficient, quantities.molar_gas_constant * temperature /
    (6 * quantities.avogadro_constant * pi * particle_radius * dynamic_viscosity))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(temperature_=temperature,
    particle_radius_=particle_radius,
    dynamic_viscosity_=dynamic_viscosity)
@validate_output(diffusion_coefficient)
def calculate_diffusion_coefficient(temperature_: Quantity, particle_radius_: Quantity,
    dynamic_viscosity_: Quantity) -> Quantity:
    result_expr = solve(law, diffusion_coefficient, dict=True)[0][diffusion_coefficient]
    result_diffusion_coefficient = result_expr.subs({
        temperature: temperature_,
        particle_radius: particle_radius_,
        dynamic_viscosity: dynamic_viscosity_
    })
    return Quantity(result_diffusion_coefficient)
