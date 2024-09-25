"""
Period of torsion pendulum from rotational inertia
==================================================

A torsion pendulum is an angular version of a linear harmonic oscillator: a disk
oscillates in a horizontal plane; the reference line oscillates with some angular amplitude.
The element of elasticity is associated with the twisting of the suspension wire.
"""

from sympy import Eq, pi, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)

period = symbols.period
"""
The :symbols:`period` of pendulum's oscillations.
"""

rotational_inertia = symbols.rotational_inertia
"""
The :symbols:`rotational_inertia` of the disk.
"""

torsion_stiffness = symbols.torsion_stiffness
"""
The :symbols:`torsion_stiffness`, which depends on the properties of the suspension wire.
"""

law = Eq(period, 2 * pi * sqrt(rotational_inertia / torsion_stiffness))
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from relation between restoring torque and twist angle


@validate_input(rotational_inertia_=rotational_inertia, torsion_constant_=torsion_stiffness)
@validate_output(period)
def calculate_period(rotational_inertia_: Quantity, torsion_constant_: Quantity) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        torsion_stiffness: torsion_constant_,
    })
    return Quantity(result)
