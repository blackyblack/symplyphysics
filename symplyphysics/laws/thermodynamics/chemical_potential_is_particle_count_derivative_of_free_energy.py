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
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function

# Law: mu = (dF/dN)_(T, V)
## mu - chemical potential
## F = F(T, V, N) - [Helmholtz free energy](./helmholtz_free_energy_via_internal_energy.py)
## N - particle count
## T - temperature
## V - volume
## (d/dN)_(T, V) - derivative w.r.t. particle count at constant temperature and volume

chemical_potential = Symbol("chemical_potential", units.energy)
r"""
Chemical potential of the system.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

free_energy = Function("free_energy", units.energy)
"""
Helmholtz free energy of the system as a function of its natural variables.

..
    TODO add link to definition file

Symbol:
    :code:`F`
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

volume = Symbol("volume", units.volume)
"""
Volume of the system.

Symbol:
    :code:`V`
"""

law = Eq(
    chemical_potential,
    Derivative(free_energy(temperature, volume, particle_count), particle_count),
)
r"""
:code:`mu = Derivative(F(T, V, N), N)`

Latex:
    .. math::
        \mu = \left( \frac{\partial F}{\partial N} \right)_{T, V}
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
