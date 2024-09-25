"""
Total work is change in kinetic energy
======================================

The work-energy principle states that the work done by all forces acting on a particle
(the work of the resultant force) equals the change in the kinetic energy of the particle.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    clone_as_symbol,
)

work = symbols.work
"""
The total :symbols:`work` done on the body.
"""

kinetic_energy = clone_as_function(symbols.kinetic_energy, display_symbol="K(t)")
"""
The :symbols:`kinetic_energy` of the body.
"""

time_before = clone_as_symbol(symbols.time, display_symbol="t_0", display_latex="t_0")
"""
The :symbols:`time` before the work has been done.
"""

time_after = clone_as_symbol(symbols.time, display_symbol="t_1", display_latex="t_1")
"""
The :symbols:`time` after the work has been done.
"""

law = Eq(work, kinetic_energy(time_after) - kinetic_energy(time_before))
"""
.. only:: comment

    Custom function arguments are not yet supported.

:code:`W = K(t_1) - K(t_0)`

:laws:latex::
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
