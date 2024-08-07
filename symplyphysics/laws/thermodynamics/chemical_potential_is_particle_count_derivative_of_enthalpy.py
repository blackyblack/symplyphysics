"""
Chemical potential is particle count derivative of enthalpy
===========================================================

The chemical potential of the system is the amount of energy the system absorbs or releases
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

enthalpy = Function("enthalpy", units.energy)
"""
Enthalpy as a function of its natural variables.

Symbol:
    :code:`H(S, p, N)`
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

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the system.

Symbol:
    :code:`p`
"""

law = Eq(
    chemical_potential,
    Derivative(enthalpy(entropy, pressure, particle_count), particle_count),
)
r"""
:code:`mu = Derivative(H(S, p, N), N)`

Latex:
    .. math::
        \mu = \left( \frac{\partial H}{\partial N} \right)_{S, p}
"""


@validate_input(
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
    enthalpy_before_=enthalpy,
    enthalpy_after_=enthalpy,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    particle_count_before_: int,
    particle_count_after_: int,
    enthalpy_before_: Quantity,
    enthalpy_after_: Quantity,
) -> Quantity:
    enthalpy_ = two_point_function(
        Point2D(particle_count_before_, enthalpy_before_),
        Point2D(particle_count_after_, enthalpy_after_),
        particle_count,
    )

    result = law.rhs.subs(
        enthalpy(entropy, pressure, particle_count),
        enthalpy_,
    ).doit()

    return Quantity(result)
