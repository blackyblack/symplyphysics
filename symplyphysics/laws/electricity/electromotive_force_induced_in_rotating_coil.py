"""
Electromotive force induced in rotating coil
============================================

Suppose a coil is being rotated around the axis that lies in the coil's cross section
(see `Figure <https://www.schoolphysics.co.uk/age16-19/Electricity%20and%20magnetism/Electromagnetic%20induction/text/Induced_emf_in_a_rotating_coil/index.html>`__)
in a magnetic field under the conditions described below. Then an electromotive will
be induced in the contour of the coil. It will be proportional to the number of turns
in the coil and the rate of change of magnetic flux through the coil.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (
    dimensionless,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.symbols.symbols import clone_function
from symplyphysics.core.geometry.line import two_point_function, Point2D

electromotive_force = symbols.electromotive_force
r"""
Electromotive force induced in the coil.
"""

coil_turn_count = SymbolNew("N", dimensionless)
"""
Number of turns in the coil.
"""

magnetic_flux = clone_function(symbols.magnetic_flux)
"""
Magnetic flux through the coil.
"""

time = symbols.time
"""
Time.
"""

law = Eq(electromotive_force, -1 * coil_turn_count * Derivative(magnetic_flux(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    coil_turn_count_=coil_turn_count,
    magnetic_flux_change_=magnetic_flux,
    time_change_=time,
)
@validate_output(electromotive_force)
def calculate_electromotive_force(
    coil_turn_count_: int,
    magnetic_flux_change_: Quantity,
    time_change_: Quantity,
) -> Quantity:
    magnetic_flux_ = two_point_function(
        Point2D(0, 0),
        Point2D(time_change_, magnetic_flux_change_),
        time,
    )
    result = law.rhs.subs({
        magnetic_flux(time): magnetic_flux_,
        coil_turn_count: coil_turn_count_,
    }).doit()
    return Quantity(result)
