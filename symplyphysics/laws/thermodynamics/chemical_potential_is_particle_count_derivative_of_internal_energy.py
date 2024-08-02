"""
Chemical potential is particle count derivative of internal energy
==================================================================

The chemical potential of a system is the amount of energy the system absorbs or releases
due to the introduction of a particle into the system, i.e. when the particle count increases
by one.
"""

from sympy import Eq, Derivative, Point2D
from symplyphysics import (
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

internal_energy = Function("internal_energy", units.energy)
"""
Internal energy of the system as a function of its natural variables.

Symbol:
    :code:`U`
"""

particle_count = Symbol("particle_count", dimensionless)
"""
Number of particles in the system.

Symbol:
    :code:`N`
"""

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

volume = Symbol("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

law = Eq(
    chemical_potential,
    Derivative(internal_energy(entropy, volume, particle_count), particle_count),
)
r"""
:code:`mu = Derivative(U(S, V, N), N)`

Latex:
    .. math::
        \mu = \left( \frac{\partial U}{\partial N} \right)_{S, V}
"""


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    internal_energy_before_=internal_energy,
    internal_energy_after_=internal_energy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    internal_energy_before_: Quantity,
    internal_energy_after_: Quantity,
) -> Quantity:
    internal_energy_ = two_point_function(
        Point2D(particle_count_before_, internal_energy_before_),
        Point2D(particle_count_after_, internal_energy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        internal_energy(entropy, volume, particle_count),
        internal_energy_,
    ).doit()

    return Quantity(result)
