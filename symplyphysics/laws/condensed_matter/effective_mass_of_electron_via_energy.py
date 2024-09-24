r"""
Effective mass of electron via energy
=====================================

In solid state physics, a particle's effective mass is the mass that it seems to have when responding to forces,
or the mass that it seems to have when interacting with other identical particles in a thermal distribution.
The effective mass is a quantity that is used to simplify band structures by modeling the behavior
of a free particle with that mass.

**Notation:**

#. :quantity_notation:`hbar`.
"""

from sympy import (
    Derivative,
    Expr,
    Eq,
    solve,
)
from symplyphysics import (
    symbols,
    clone_as_symbol,
    SI,
    quantities,
    Quantity,
    validate_input,
    validate_output,
    clone_as_function,
)

effective_mass = clone_as_symbol(symbols.mass, display_symbol="m_eff", display_latex="m_\\text{eff}")
r"""
Effective :symbols:`mass` of the electron.
"""

energy = clone_as_function(symbols.energy, display_symbol="E(k)")
"""
Electron energy as a function of angular wavenumber.
"""

angular_wavenumber = symbols.angular_wavenumber
r"""
:symbols:`angular_wavenumber`.
"""

law = Eq(effective_mass, quantities.hbar**2 / Derivative(energy(angular_wavenumber),
    (angular_wavenumber, 2)))
r"""
:laws:symbol::

:laws:latex::
"""


@validate_input(wavenumber_=angular_wavenumber)
@validate_output(effective_mass)
def calculate_mass(energy_function_: Expr, wavenumber_: Quantity) -> Quantity:
    energy_function_quantity = Quantity(energy_function_)
    assert SI.get_dimension_system().equivalent_dims(energy_function_quantity.dimension,
        energy.dimension)
    applied_law = law.subs(energy(angular_wavenumber), energy_function_)
    result_expr = applied_law.subs(angular_wavenumber, wavenumber_)
    result = solve(result_expr, effective_mass, dict=True)[0][effective_mass]
    return Quantity(result)
