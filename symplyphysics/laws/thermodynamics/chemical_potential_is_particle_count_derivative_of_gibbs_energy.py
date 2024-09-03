"""
Chemical potential is particle count derivative of Gibbs energy
===============================================================

The chemical potential of the system is the amount of energy the system absorbs or releases
due to the introduction of a particle into the system, i.e. when the particle count increases
by one.
"""

from sympy import Eq, Derivative, Point2D
from symplyphysics import (
    symbols,
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
Chemical potential of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

gibbs_energy = Function("gibbs_energy", units.energy)
"""
Gibbs energy of the system as a function of its natural variables.

Symbol:
    :code:`G(T, p, N)`
"""

particle_count = Symbol("particle_count", dimensionless)
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

temperature = symbols.thermodynamics.temperature
"""
Temperature of the system.
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

law = Eq(
    chemical_potential,
    Derivative(gibbs_energy(temperature, pressure, particle_count), particle_count),
)
r"""
:code:`mu = Derivative(G(T, p, N), N)`

Latex:
    .. math::
        \mu = \left( \frac{\partial G}{\partial N} \right)_{T, p}
"""


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    gibbs_energy_before_=gibbs_energy,
    gibbs_energy_after_=gibbs_energy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    gibbs_energy_before_: Quantity,
    gibbs_energy_after_: Quantity,
) -> Quantity:
    gibbs_energy_ = two_point_function(
        Point2D(particle_count_before_, gibbs_energy_before_),
        Point2D(particle_count_after_, gibbs_energy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        gibbs_energy(temperature, pressure, particle_count),
        gibbs_energy_,
    ).doit()

    return Quantity(result)
