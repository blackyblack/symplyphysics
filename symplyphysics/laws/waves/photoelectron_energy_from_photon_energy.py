"""
Photoelectron energy from frequency
===================================

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
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)

maximum_kinetic_energy = Symbol("maximum_kinetic_energy", units.energy)
r"""
Maximum kinetic energy of emitted electron.

Symbol:
    :code:`K_max`

Latex:
    :math:`K_\text{max}`
"""

photon_energy = Symbol("photon_energy", units.energy)
r"""
Energy of absorbed photon.

Symbol:
    :code:`E_ph`

Latex:
    :math:`E_\text{ph}`
"""

work_function = Symbol("work_function", units.energy)
"""
Work function of the surface.

Symbol:
    :code:`W`
"""

law = Eq(maximum_kinetic_energy, photon_energy - work_function)
r"""
:code:`K_max = E_ph - W`

Latex:
    .. math::
        K_\text{max} = E_\text{ph} - W
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
