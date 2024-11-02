"""
Chemical potential is particle count derivative of internal energy
==================================================================

The chemical potential of a system is the amount of energy the system absorbs or releases
due to the introduction of a particle into the system, i.e. when the particle count increases
by one.
"""

from sympy import Eq, Derivative, Point2D
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)
from symplyphysics.core.geometry.line import two_point_function

chemical_potential = symbols.chemical_potential
r"""
:symbols:`chemical_potential` of the system.
"""

entropy = symbols.entropy
"""
:symbols:`entropy` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

internal_energy = clone_as_function(
    symbols.internal_energy,
    [entropy, volume, particle_count],
)
"""
:symbols:`internal_energy` of the system as a function of its natural variables.
"""

law = Eq(
    chemical_potential,
    Derivative(internal_energy(entropy, volume, particle_count), particle_count),
)
"""
:laws:symbol::

:laws:latex::
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
