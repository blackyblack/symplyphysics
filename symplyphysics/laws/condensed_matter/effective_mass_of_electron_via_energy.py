"""
Effective mass of electron via energy
=====================================

In solid state physics, a particle's effective mass is the mass that it seems to have when responding to forces,
or the mass that it seems to have when interacting with other identical particles in a thermal distribution.
The effective mass is a quantity that is used to simplify band structures by modeling the behavior
of a free particle with that mass.

**Notation:**

#. :math:`\hbar` (:code:`hbar`) is the reduced Planck constant.
"""

from sympy import (
    Derivative,
    Expr,
    Eq,
    solve,
)
from symplyphysics import (
    symbols,
    clone_symbol,
    SI,
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

effective_mass = clone_symbol(symbols.basic.mass, "effective_mass", display_symbol="m_eff", display_latex=r"m_\text{eff}")
r"""
Effective mass of the electron.
"""

energy = Function("energy", units.energy)
"""
Electron energy as a function of wavenumber.

Symbol:
    :code:`E(k)`
"""

wavenumber = Symbol("wavenumber", 1 / units.length)
r"""
Wavenumber.

Symbol:
    :code:`k`
"""

law = Eq(effective_mass, units.hbar**2 / Derivative(energy(wavenumber),
    (wavenumber, 2)))
r"""
:code:`m_eff = hbar^2 * Derivative(E(k), (k, 2))^(-1)`

Latex:
    .. math::
        m_\text{eff} = \hbar^2 \left( \frac{d^2 E}{d k^2} \right)^{-1}
"""


def _apply_energy_function(energy_function_: Expr) -> Expr:
    applied_law = law.subs(energy(wavenumber), energy_function_)
    return applied_law


@validate_input(wavenumber_=wavenumber)
@validate_output(effective_mass)
def calculate_mass(energy_function_: Expr, wavenumber_: Quantity) -> Quantity:
    energy_function_quantity = Quantity(energy_function_)
    assert SI.get_dimension_system().equivalent_dims(energy_function_quantity.dimension,
        energy.dimension)
    applied_law = _apply_energy_function(energy_function_)
    result_expr = applied_law.subs(wavenumber, wavenumber_)
    result = solve(result_expr, effective_mass, dict=True)[0][effective_mass]
    return Quantity(result)
