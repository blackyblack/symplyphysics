"""
Relativistic kinetic energy
===========================

The kinetic energy of an object is the form of energy that it possesses due to its motion.

**Notation:**

#. :quantity_notation:`speed_of_light`.

**Notes**

#. The work expended accelerating an object from rest approaches infinity as the velocity approaches the
   speed of light. Thus it is impossible to accelerate an object across this boundary.

**Links:**

#. `Wikipedia, last formula in paragraph <https://en.wikipedia.org/wiki/Kinetic_energy#Derivation_2>`__.
"""

from sympy import Eq
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)

kinetic_energy = symbols.kinetic_energy
"""
:symbols:`kinetic_energy` of the body.
"""

lorentz_factor = symbols.lorentz_factor
"""
:symbols:`lorentz_factor` of the body.
"""

rest_mass = symbols.rest_mass
"""
:symbols:`rest_mass` of the body.
"""

law = Eq(kinetic_energy, (lorentz_factor - 1) * rest_mass * quantities.speed_of_light**2)
"""
:laws:symbol::

:laws:latex::
"""

# TODO: derive from relativistic momentum and expression for kinetic energy


@validate_input(
    lorentz_factor_=lorentz_factor,
    rest_mass_=rest_mass,
)
@validate_output(kinetic_energy)
def calculate_kinetic_energy(
    lorentz_factor_: float,
    rest_mass_: Quantity,
) -> Quantity:
    if lorentz_factor_ < 1:
        raise ValueError("Lorentz factor must be greater or equal to 1")

    result = law.rhs.subs({
        lorentz_factor: lorentz_factor_,
        rest_mass: rest_mass_,
    })
    return Quantity(result)
