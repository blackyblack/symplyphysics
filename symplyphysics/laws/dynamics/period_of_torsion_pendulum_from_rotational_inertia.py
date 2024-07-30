"""
Period of torsion pendulum from rotational inertia
==================================================

A torsion pendulum is an angular version of a linear harmonic oscillator: a disk
oscillates in a horizontal plane; the reference line oscillates with some angular amplitude.
The element of elasticity is associated with the twisting of the suspension wire.
"""

from sympy import Eq, pi, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

period = Symbol("period", units.time)
"""
The period of pendulum's oscillations.

Symbol:
    :code:`T`
"""

rotational_inertia = Symbol("rotational_inertia", units.mass * units.length**2)
"""
The rotational inertia of the disk.

Symbol:
    :code:`I`
"""

torsion_constant = Symbol("torsion_constant", units.force * units.length)
r"""
The torsion constant, which depends on the properties of the suspension wire.

Symbol:
    :code:`kappa`

Latex:
    :math:`\kappa`
"""

law = Eq(period, 2 * pi * sqrt(rotational_inertia / torsion_constant))
r"""
:code:`T = 2 * pi * sqrt(I / kappa)`

Latex:
    .. math::
        T = 2 \pi \sqrt{\frac{I}{\kappa}}
"""

# TODO: derive from relation between restoring torque and twist angle


@validate_input(rotational_inertia_=rotational_inertia, torsion_constant_=torsion_constant)
@validate_output(period)
def calculate_period(rotational_inertia_: Quantity, torsion_constant_: Quantity) -> Quantity:
    result = law.rhs.subs({
        rotational_inertia: rotational_inertia_,
        torsion_constant: torsion_constant_,
    })
    return Quantity(result)
