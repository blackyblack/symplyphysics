"""
Photoelectron energy from photon energy
=======================================

The photoelectric effect is the emission of electrons from a material caused by
electromagnetic radiation. Electrons emitted in this manner are called photoelectrons.

In the photoemission process, when an electron within some material absorbs the
energy of a photon and acquires more energy than its binding energy, it is likely
to be ejected. If the photon energy is too low, the electron is unable to escape
the material.

The theory predicts that the highest kinetic energy of emitted electrons is equal
to the difference between the absorbed photon energy and the work function of the
surface, which is the minimum energy required to remove an electron from the surface
of the material.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Photoelectric_effect#Theoretical_explanation>`__.

..
    TODO move to `quantum_mechanics`?
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)

maximum_kinetic_energy = clone_as_symbol(
    symbols.kinetic_energy,
    display_symbol="K_max",
    display_latex="K_\\text{max}",
)
"""
Maximum :symbols:`kinetic_energy` of emitted electron.
"""

photon_energy = symbols.energy
"""
:symbols:`energy` of absorbed photon.
"""

work_function = symbols.work_function
"""
:symbols:`work_function` of the surface.
"""

law = Eq(maximum_kinetic_energy, photon_energy - work_function)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(photon_energy_=photon_energy, work_function_=work_function)
@validate_output(maximum_kinetic_energy)
def calculate_max_kinetic_energy(photon_energy_: Quantity, work_function_: Quantity) -> Quantity:
    result_energy_expr = solve(law, maximum_kinetic_energy, dict=True)[0][maximum_kinetic_energy]
    result_expr = result_energy_expr.subs({
        photon_energy: photon_energy_,
        work_function: work_function_
    })
    return Quantity(result_expr)
