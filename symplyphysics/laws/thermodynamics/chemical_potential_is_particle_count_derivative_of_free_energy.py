"""
Chemical potential is particle count derivative of free energy
==============================================================

The chemical potential of the system is the amount of energy the system absorbs or releases
due to the introduction of a particle into the system, i.e. when the particle count increases
by one.
"""

from sympy import Eq, Derivative, Point2D
from symplyphysics import (
    symbols,
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
)
from symplyphysics.core.geometry.line import two_point_function

chemical_potential = symbols.chemical_potential
r"""
:symbols:`chemical_potential` of the system.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

free_energy = clone_as_function(
    symbols.helmholtz_free_energy,
    [temperature, volume, particle_count],
)
"""
Helmholtz free energy of the system as a function of its natural variables.

..
    TODO add link to definition file
"""

law = Eq(
    chemical_potential,
    Derivative(free_energy(temperature, volume, particle_count), particle_count),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    free_energy_before_=free_energy,
    free_energy_after_=free_energy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    free_energy_before_: Quantity,
    free_energy_after_: Quantity,
) -> Quantity:
    free_energy_ = two_point_function(
        Point2D(particle_count_before_, free_energy_before_),
        Point2D(particle_count_after_, free_energy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        free_energy(temperature, volume, particle_count),
        free_energy_,
    ).doit()

    return Quantity(result)
