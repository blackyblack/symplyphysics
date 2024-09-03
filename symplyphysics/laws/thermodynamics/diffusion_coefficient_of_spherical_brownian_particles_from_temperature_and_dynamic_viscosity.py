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

#. :math:`R` is the molar gas constant.
#. :math:`N_\text{A}` is the Avogadro constant.

**Conditions:**

#. Particle displacements are equally likely in any direction.
#. The inertia of a Brownian particle can be neglected compared to the influence of friction forces.
#. Particles are spherical.
#. Low Reynolds number, i.e. non-turbulent flow.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input, validate_output)

diffusion_coefficient = Symbol("diffusion_coefficient", units.area / units.time)
"""
Diffusion coefficient of the particles.

Symbol:
    :code:`D`
"""

temperature = symbols.thermodynamics.temperature
"""
Temperature of the system.
"""

particle_radius = Symbol("particle_radius", units.length)
"""
Radius of the particles.

Symbol:
    :code:`r`
"""

dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
r"""
Dynamic viscosity of the particles.

Symbol:
    :code:`eta`

Latex:
    :math:`\eta`
"""

law = Eq(
    diffusion_coefficient, units.molar_gas_constant * temperature /
    (6 * units.avogadro_constant * pi * particle_radius * dynamic_viscosity))
r"""
:code:`D = (R * T) / (6 * N_A * pi * r * eta)`

Latex:
    .. math::
        D = \frac{R T}{6 N_\text{A} \pi r \eta}
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
