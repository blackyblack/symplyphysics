"""
Period of rotation of charged particle in magnetic field
========================================================

When a charged particle enters a magnetic field, it experiences an :doc:`electromagnetic
force <laws.electricity.vector.lorentz_force_via_electromagnetic_field>` upon itself.
In the absence of the electric field, the particle starts moving in a circular orbit. The
period of the particle's rotation is determined by the mass and charge of the particle
as well as by the magnetic field and it does not depend on the particle speed.

**Conditions:**

#. The particle's speed and the magnetic field are perpendicular to each other.
#. The magnetic field is uniform.
#. The electric field is zero.
"""

from sympy import (Eq, solve, pi)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, symbols)

period = Symbol("period", units.time)
"""
Period of the particle's rotation.

Symbol:
    :code:`T`
"""

mass = symbols.basic.mass
"""
:attr:`~symplyphysics.symbols.basic.mass` of the particle.

Symbol:
    :code:`m`
"""

charge = Symbol("charge", units.charge)
"""
Charge of the particle.

Symbol:
    :code:`q`
"""

magnetic_flux_density = Symbol("magnetic_flux_density", units.magnetic_density)
"""
Magnitude of magnetic flux density.

Symbol:
    :code:`B`
"""

law = Eq(period, 2 * pi * mass / (charge * magnetic_flux_density))
r"""
:code:`T = 2 * pi * m / (q * B)`

Latex:
    .. math::
        T = 2 \pi \frac{m}{q B}
"""


@validate_input(mass_=mass, charge_=charge, induction_=magnetic_flux_density)
@validate_output(period)
def calculate_period(mass_: Quantity, charge_: Quantity, induction_: Quantity) -> Quantity:
    result_period_expr = solve(law, period, dict=True)[0][period]
    result_expr = result_period_expr.subs({
        mass: mass_,
        charge: charge_,
        magnetic_flux_density: induction_
    })
    return Quantity(result_expr)
