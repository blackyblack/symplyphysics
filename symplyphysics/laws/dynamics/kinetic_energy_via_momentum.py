"""
Kinetic energy via momentum
===========================

Kinetic energy can be expressed as a function of momentum and mass.

**Notes:**

#. This relation also holds in Quantum Mechanics.

**Links:**

#. `Wikipedia, fourth formula <https://en.wikipedia.org/wiki/Kinetic_energy#Kinetic_energy_of_rigid_bodies>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.dynamics import (
    kinetic_energy_from_mass_and_speed as energy_law,)
from symplyphysics.definitions import (
    momentum_is_mass_times_speed as momentum_def,)

kinetic_energy = symbols.kinetic_energy
"""
The :symbols:`kinetic_energy` of the object.
"""

momentum = symbols.momentum
"""
The :symbols:`momentum` of the object.
"""

mass = symbols.mass
"""
The :symbols:`mass` of the object.
"""

law = Eq(kinetic_energy, momentum**2 / (2 * mass))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from kinetic energy and momentum expressions

_energy_eqn = energy_law.law.subs({
    energy_law.kinetic_energy: kinetic_energy,
    energy_law.speed: momentum_def.speed,
    energy_law.mass: mass,
})

_momentum_eqn = momentum_def.definition.subs({
    momentum_def.momentum: momentum,
    momentum_def.mass: mass,
})

_kinetic_energy_expr = solve(
    (_energy_eqn, _momentum_eqn),
    (kinetic_energy, momentum_def.speed),
    dict=True,
)[0][kinetic_energy]

assert expr_equals(_kinetic_energy_expr, law.rhs)


@validate_input(
    momentum_=momentum,
    mass_=mass,
)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(
    momentum_: Quantity,
    mass_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        momentum: momentum_,
        mass: mass_,
    })
    return Quantity(result)
