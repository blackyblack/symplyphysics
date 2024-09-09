r"""
Equilibrium voltage difference in p-n junction via concentrations
=================================================================

The p-n junction has a potential barrier preventing the movement of charge carriers.

**Notation:**

#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.

**Conditions:**

#. :math:`n_\text{a}` is constant on the p side and zero on the n side.
#. :math:`n_\text{d}` is constant on the n side and zero on the p side.
#. The semiconductor is non-degenerate, i.e. :math:`n_\text{i}` is independent of
   the Fermi energy.

**Links:**

#. `Size of depleting region <https://en.wikipedia.org/wiki/P%E2%80%93n_junction#Size_of_depletion_region>`_.
"""

from sympy import (Eq, solve, log)
from sympy.physics.units import boltzmann
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input,
    validate_output)

equilibrium_voltage_difference = Symbol("equilibrium_voltage_difference", units.voltage)
r"""
Equilibrium voltage difference corresponding to the size of the depletion region.

Symbol:
    :code:`Delta(V)`

Latex:
    :math:`\Delta V`
"""

donor_concentration = Symbol("donor_concentration", 1 / units.volume)
r"""
Concentration of positively-charged donor atoms.

Symbol:
    :code:`n_d`

Latex:
    :math:`n_\text{d}`
"""

acceptor_concentration = Symbol("acceptor_concentration", 1 / units.volume)
r"""
Concentration of negatively-charged acceptor atoms.

Symbol:
    :code:`n_a`

Latex:
    :math:`n_\text{a}`
"""

charge_carriers_concentration = Symbol("charge_carriers_concentration", 1 / units.volume)
r"""
Concentration of intrinsic charge carriers.

Symbol:
    :code:`n`
"""

temperature = symbols.thermodynamics.temperature
r"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the semiconductor.
"""

charge_electron = Symbol("charge_electron", units.charge)
"""
Magnitude of the electron charge.

Symbol:
    :code:`q`
"""

law = Eq(equilibrium_voltage_difference, (boltzmann * temperature / charge_electron) *
    log(donor_concentration * acceptor_concentration / charge_carriers_concentration**2))
r"""
:code:`Delta(V) = (k_B * T / q) * log(n_a * n_d / n^2)`

Latex:
    .. math::
        \Delta V = \frac{k_\text{B} T}{q} \log \left( \frac{n_\text{a} n_\text{d}}{n^2} \right)
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
