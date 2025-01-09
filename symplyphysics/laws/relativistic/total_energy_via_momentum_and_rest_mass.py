"""
Total energy via momentum and rest mass
=======================================

The energy—momentum relation, also called relativistic dispersion relation, is a relativistic
equation relating total energy to invariant mass and momentum. It is the extension of
mass-energy equivalence (:ref:`Total energy via relativistic mass`) for bodies or systems with non-zero momentum.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Energy%E2%80%93momentum_relation>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    symbols,
    quantities,
    Quantity,
    validate_input,
    validate_output,
)

relativistic_energy = symbols.energy
"""
Total, or relativistic, :symbols:`energy` of the body.
"""

relativistic_momentum = symbols.momentum
"""
Relativistic :symbols:`momentum` of the body.
"""

invariant_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the body.
"""

law = Eq(
    relativistic_energy**2,
    (relativistic_momentum * quantities.speed_of_light)**2 +
    (invariant_mass * quantities.speed_of_light**2)**2,
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    relativistic_momentum_=relativistic_momentum,
    invariant_mass_=invariant_mass,
)
@validate_output(relativistic_energy)
def calculate_relativistic_energy(
    relativistic_momentum_: Quantity,
    invariant_mass_: Quantity,
) -> Quantity:
    expr = solve(law, relativistic_energy)[1]
    result = expr.subs({
        relativistic_momentum: relativistic_momentum_,
        invariant_mass: invariant_mass_,
    })
    return Quantity(result)
