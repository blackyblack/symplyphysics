"""
Inductance is magnetic flux over current
========================================

Inductance is the tendency of an electrical component to oppose a change in the electric current
flowing through it. It can be defined as the ratio of the total magnetic flux through a circuit
due to the component over the current in it.
"""

from sympy import Eq
from symplyphysics import (
    SymbolNew,
    units,
    Quantity,
    validate_input,
    validate_output,
)

inductance = SymbolNew("L", units.inductance)
"""
Inductance.
"""

magnetic_flux = SymbolNew("Phi_B", units.magnetic_flux, display_latex="\\Phi_B")
"""
Magnetic flux due to the component.
"""

current = SymbolNew("I", units.current)
"""
Current flowing through the circuit.
"""

law = Eq(inductance, magnetic_flux / current)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    magnetic_flux_=magnetic_flux,
    current_=current,
)
@validate_output(inductance)
def calculate_inductance(
    magnetic_flux_: Quantity,
    current_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        magnetic_flux: magnetic_flux_,
        current: current_,
    })
    return Quantity(result)
