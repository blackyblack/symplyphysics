"""
Chemical potential is particle count derivative of enthalpy
===========================================================

The chemical potential of the system is the amount of energy the system absorbs or releases
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

pressure = symbols.pressure
"""
:symbols:`pressure` inside the system.
"""

particle_count = symbols.particle_count
"""
:symbols:`particle_count` of the system.
"""

enthalpy = clone_as_function(symbols.enthalpy, [entropy, pressure, particle_count])
"""
:symbols:`enthalpy` as a function of its natural variables.
"""

law = Eq(
    chemical_potential,
    Derivative(enthalpy(entropy, pressure, particle_count), particle_count),
)
"""
:laws:symbol::

:laws:latex::
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
