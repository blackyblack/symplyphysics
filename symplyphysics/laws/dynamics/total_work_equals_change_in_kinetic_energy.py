"""
Total work is change in kinetic energy
======================================

The work-energy principle states that the work done by all forces acting on a particle
(the work of the resultant force) equals the change in the kinetic energy of the particle.
"""

from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

work = Symbol("work", units.energy)
"""
The total work done on the body.

Symbol:
    :code:`W`
"""

kinetic_energy = Function("kinetic_energy", units.energy)
r"""
The kinetic energy of the body.

Symbol:
    :code:`K(t)`
"""

time_before = Symbol("time_before", units.time)
"""
The time before the work has been done.

Symbol:
    :code:`t0`

Latex:
    :math:`t_0`
"""

time_after = Symbol("time_after", units.time)
"""
The time after the work has been done.

Symbol:
    :code:`t1`

Latex:
    :math:`t_1`
"""

law = Eq(work, kinetic_energy(time_after) - kinetic_energy(time_before))
r"""
:code:`W = K(t1) - K(t0)`

Latex:
    .. math::
        W = K(t_1) - K(t_0)
"""

# TODO Derive the law in case of rectilinear motion with infinitesimal constant total force acting on particle and
# integrating the resulting infinitesimal work over time


@validate_input(kinetic_energy_before_=kinetic_energy, kinetic_energy_after_=kinetic_energy)
@validate_output(work)
def calculate_total_work(kinetic_energy_before_: Quantity,
    kinetic_energy_after_: Quantity) -> Quantity:
    result = law.rhs.subs({
        kinetic_energy(time_before): kinetic_energy_before_,
        kinetic_energy(time_after): kinetic_energy_after_,
    })
    return Quantity(result)
