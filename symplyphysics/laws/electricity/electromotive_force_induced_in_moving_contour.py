"""
Electromotive force induced in moving contour
=============================================

The **Faraday's law** states that the electromotive force around a closed path is equal to
the negative of the time rate of change of the magnetic flux enclosed by the path. In case
of the current making several turns around the contour, e.g. in a coil, the electromotive force
would also be proportional to the number of turn the current makes.
"""

from sympy import (Eq, Derivative)
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

electromotive_force = symbols.electromotive_force
r"""
:symbols:`electromotive_force` induced in the contour.
"""

current_turn_count = symbols.positive_number
"""
Number of turns the current makes around the contour. See :symbols:`positive_number`.
"""

time = symbols.time
"""
:symbols:`time`.
"""

magnetic_flux = clone_as_function(symbols.magnetic_flux, [time])
"""
:symbols:`magnetic_flux` through the contour.
"""

law = Eq(electromotive_force, -1 * current_turn_count * Derivative(magnetic_flux(time), time))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    current_turn_count_=current_turn_count,
    magnetic_flux_change_=magnetic_flux,
    time_change_=time,
)
@validate_output(electromotive_force)
def calculate_electromotive_force(
    current_turn_count_: int,
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
        current_turn_count: current_turn_count_,
    }).doit()
    return Quantity(result)
