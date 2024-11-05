r"""
Quantum isochoric molar heat capacity of solids
===============================================

To derive the heat capacity of a solid, one should account for quantum effects. Albert Einstein used the
same model as in the classical case, namely the atoms are harmonic oscillators with three degrees
of freedom, located in the nodes of the crystal lattice, performing thermal oscillations around the
equlibrium positions with the same frequency. But he used a more correct expression for the energy
of the oscillators, and although the result still only qualitatively describes the heat capacity of
solids, it is a big achievement and the result has correct asymptotic behaviour for :math:`T \to 0`.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.
"""

from sympy import Eq, exp
from symplyphysics import (
    dimensionless,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    quantities,
    clone_as_symbol,
    symbols,
)

isochoric_molar_heat_capacity = clone_as_symbol(
    symbols.molar_heat_capacity,
    display_symbol="c_Vm",
    display_latex="c_{V, m}",
)
"""
:symbols:`molar_heat_capacity` at constant :symbols:`volume`.
"""

reduced_photon_energy = SymbolNew("x", dimensionless)
r"""
Reduced photon energy, defined as the ratio of photon energy :math:`\hbar \omega` or :math:`h \nu` to
thermal energy :math:`k_\text{B} T`.
"""

law = Eq(isochoric_molar_heat_capacity, (3 * quantities.molar_gas_constant) *
    (reduced_photon_energy**2 * exp(reduced_photon_energy) / (exp(reduced_photon_energy) - 1)**2))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(reduced_photon_energy_=reduced_photon_energy)
@validate_output(isochoric_molar_heat_capacity)
def calculate_isochoric_molar_heat_capacity(reduced_photon_energy_: float,) -> Quantity:
    result = law.rhs.subs({
        reduced_photon_energy: reduced_photon_energy_,
    })
    return Quantity(result)
