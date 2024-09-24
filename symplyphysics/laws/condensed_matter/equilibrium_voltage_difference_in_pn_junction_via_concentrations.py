r"""
Equilibrium voltage difference in p-n junction via concentrations
=================================================================

The p-n junction has a potential barrier preventing the movement of charge carriers.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. :math:`n_\text{a}` is constant on the p side and zero on the n side.
#. :math:`n_\text{d}` is constant on the n side and zero on the p side.
#. The semiconductor is non-degenerate, i.e. :math:`n_\text{i}` is independent of
   the Fermi energy.

**Links:**

#. `Size of depleting region <https://en.wikipedia.org/wiki/P%E2%80%93n_junction#Size_of_depletion_region>`_.
"""

from sympy import (Eq, solve)
from symplyphysics import (symbols, Quantity, validate_input,
    validate_output, quantities, clone_as_symbol)
from symplyphysics.core.functions import log

equilibrium_voltage_difference = clone_as_symbol(
    symbols.voltage,
    display_symbol="Delta(V)",
    display_latex="\\Delta V"
)
"""
Equilibrium :symbols:`voltage` difference corresponding to the size of the depletion region.
"""

donor_concentration = clone_as_symbol(
    symbols.number_density,
    display_symbol="n_d",
    display_latex="n_\\text{d}",
)
"""
:symbols:`number_density` of positively-charged donor atoms.
"""

acceptor_concentration = clone_as_symbol(
    symbols.number_density,
    display_symbol="n_a",
    display_latex="n_\\text{a}",
)
"""
:symbols:`number_density` of negatively-charged acceptor atoms.
"""

charge_carriers_concentration = symbols.number_density
"""
:symbols:`number_density` of intrinsic charge carriers.
"""

temperature = symbols.temperature
r"""
:symbols:`temperature` of the semiconductor.
"""

charge_electron = symbols.charge
"""
Magnitude of the electron :symbols:`charge`.
"""

law = Eq(equilibrium_voltage_difference, (quantities.boltzmann_constant * temperature / charge_electron) *
    log(donor_concentration * acceptor_concentration / charge_carriers_concentration**2))
r"""
:laws:symbol::

:laws:latex::
"""


@validate_input(donors_concentration_=donor_concentration,
    acceptors_concentration_=acceptor_concentration,
    charge_carriers_concentration_=charge_carriers_concentration,
    temperature_=temperature,
    charge_electron_=charge_electron)
@validate_output(equilibrium_voltage_difference)
def calculate_height_barrier(donors_concentration_: Quantity, acceptors_concentration_: Quantity,
    charge_carriers_concentration_: Quantity, temperature_: Quantity,
    charge_electron_: Quantity) -> Quantity:
    result_expr = solve(law, equilibrium_voltage_difference, dict=True)[0][equilibrium_voltage_difference]
    result_expr = result_expr.subs({
        donor_concentration: donors_concentration_,
        acceptor_concentration: acceptors_concentration_,
        charge_carriers_concentration: charge_carriers_concentration_,
        temperature: temperature_,
        charge_electron: charge_electron_
    })
    return Quantity(result_expr)
