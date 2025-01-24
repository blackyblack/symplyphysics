"""
Electron distribution function in gas plasma per Druyvestein
============================================================

In a gas discharge, electrons have a wide range of energies, which is described by the
electron energy distribution function. Electrons in a gas-discharge plasma acquire their
energy under the action of an electric field. Energy consumption occurs due to elastic
and, especially, inelastic collisions with atoms. In addition, energy exchange between
electrons is also possible in plasma. Depending on the relationship between all these
factors, different electron energy distributions are established. Under equilibrium
conditions, the Maxwell distribution is most common. But in the case of intense
ionization, the number of fast electrons decreases in the distribution function, and it
passes into the Druyvestein distribution function.

**Notation:**

#. :quantity_notation:`elementary_charge`.

**Links:**

#. `Comsol, possible similar formula here <https://www.comsol.com/blogs/electron-energy-distribution-function>`__.

..
    TODO: find a more suitable link
"""

from sympy import Eq, Rational, solve, exp
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
    symbols,
)
from symplyphysics.quantities import elementary_charge
from symplyphysics.core.symbols.probability import Probability

distribution_function = SymbolNew("f", dimensionless)
"""
Electron distribution function.
"""

voltage = symbols.voltage
"""
:symbols:`voltage` between electrodes.
"""

electron_energy = symbols.energy
"""
Electron :symbols:`energy`.
"""

energy_constant = Quantity(1.04 * units.electronvolt, display_symbol="E_0")
"""
Constant equal to :math:`1.04 \\, \\text{eV}`.
"""

law = Eq(
    distribution_function,
    energy_constant * (elementary_charge * voltage)**Rational(1, 2) / electron_energy**Rational(3, 2)
    * exp(-0.55 * (elementary_charge * voltage)**2 / electron_energy**2),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(voltage_between_electrodes_=voltage,
    electron_energy_=electron_energy)
@validate_output(distribution_function)
def calculate_value_of_distribution_function(voltage_between_electrodes_: Quantity,
    electron_energy_: Quantity) -> Probability:
    result_expr = solve(law, distribution_function,
        dict=True)[0][distribution_function]
    result_expr = result_expr.subs({
        voltage: voltage_between_electrodes_,
        electron_energy: electron_energy_,
    })
    return Probability(convert_to_float(result_expr))
