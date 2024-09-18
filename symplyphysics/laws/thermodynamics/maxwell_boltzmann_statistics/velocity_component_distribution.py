r"""
Velocity component distribution
===============================

For a system containing a large number of identical non-interacting non-relativistic
classical particles in thermodynamic equilibrium, the velocity component distribution
is a function :math:`f(v_k)` such that :math:`f(v_k) dv_k` gives the fraction of particles with
speeds in the interval :math:`dv_k` around velocity component :math:`v_k`.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Notes:**

#. Applicable for any velocity component in Cartesian coordinates.

**Conditions:**

#. Number of particles is big enough that the laws of thermodynamics can be applied.
#. Particles are identical, non-interacting, non-relativistic, and classical.
#. The ensemble of particles is at thermodynamic equilibrium.
"""

from sympy import Eq, sqrt, pi, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    clone_as_symbol,
    symbols,
    quantities,
)

velocity_component_distribution = Symbol("velocity_component_distribution",
    1 / units.velocity,
    positive=True)
r"""
Distribution function of velocity component :math:`v_k`.

Symbol:
    :code:`f(v_k)`

Latex:
    :math:`f(v_k)`
"""

velocity_component = Symbol("velocity_component", units.velocity, positive=True)
"""
Velocity component in Cartesian coordinates, :math:`k = x, y, z`.

Symbol:
    :code:`v_k`

Latex:
    :math:`v_k`
"""

particle_mass = clone_as_symbol(symbols.mass, positive=True)
"""
:symbols:`mass` of a particle.
"""

equilibrium_temperature = clone_as_symbol(symbols.temperature, positive=True)
"""
Equilibrium :symbols:`temperature` of the ensemble.
"""

law = Eq(
    velocity_component_distribution,
    sqrt(particle_mass / (2 * pi * quantities.boltzmann_constant * equilibrium_temperature)) *
    exp(-1 * particle_mass * velocity_component**2 /
    (2 * quantities.boltzmann_constant * equilibrium_temperature)))
r"""
:code:`f(v_k) = sqrt(m / (2 * pi * k_B * T)) * exp(-1 * m * v_k^2 / (2 * k_B * T))`

Latex:
    .. math::
        f(v_k) = \sqrt{\frac{m}{2 \pi k_\text{B} T}} \exp \left( - \frac{m v_k^2}{2 k_\text{B} T} \right)
"""


@validate_input(
    velocity_component_=velocity_component,
    particle_mass_=particle_mass,
    ensemble_temperature_=equilibrium_temperature,
)
@validate_output(velocity_component_distribution)
def calculate_velocity_component_distribution(
    velocity_component_: Quantity,
    particle_mass_: Quantity,
    ensemble_temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        velocity_component: velocity_component_,
        particle_mass: particle_mass_,
        equilibrium_temperature: ensemble_temperature_,
    })
    return Quantity(result)
