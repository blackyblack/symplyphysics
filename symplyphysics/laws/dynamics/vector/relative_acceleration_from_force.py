"""
Relative acceleration from force
================================

Suppose reference frame :math:`S` is fixed to a moving object (e.g. Earth). For some body
:math:`B` we can write an equation of motion in coordinates of :math:`S'` akin to the Newton's
second law of motion for inertial frames, although we obtain two additional components to the
equation: one corresponding to the Coriolis force, and another to the fictitious force of
translation between inertial frame :math:`S` and non-inertial frame :math:`S'`.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Coriolis_force#Formula>`__.
"""

from sympy import Eq
from symplyphysics import symbols, Quantity, validate_input, validate_output

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import solve_for_vector

relative_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_rel",
    display_latex="{\\vec a}_\\text{rel}",
)
"""
Vector of relative :symbols:`acceleration` of body :math:`B` relative to :math:`S'`
"""

force = clone_as_vector_symbol(symbols.force)
"""
Vector of the net physical :symbols:`force` exerted on body :math:`B`.
"""

mass = symbols.mass
"""
:symbols:`mass` of body :math:`B`.
"""

coriolis_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_Cor",
    display_latex="{\\vec a}_\\text{Cor}",
)
"""
Vector of the Coriolis :symbols:`acceleration` of body :math:`B`.

..
    TODO: add link to vector law
"""

translation_acceleration = clone_as_vector_symbol(
    symbols.acceleration,
    display_symbol="a_tr",
    display_latex="{\\vec a}_\\text{tr}",
)
"""
Vector of translation :symbols:`acceleration` of body :math:`B`.

..
    TODO: add link to vector law
"""

law = Eq(relative_acceleration, force / mass + coriolis_acceleration - translation_acceleration)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    mass_=mass,
    force_=force,
    coriolis_acceleration_=coriolis_acceleration,
    translation_acceleration_=translation_acceleration,
)
@validate_output(relative_acceleration)
def calculate_relative_acceleration(
    mass_: Quantity,
    force_: QuantityCoordinateVector,
    coriolis_acceleration_: QuantityCoordinateVector,
    translation_acceleration_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    expr = solve_for_vector(law, relative_acceleration)
    value = expr.subs({
        mass: mass_,
        force: force_,
        coriolis_acceleration: coriolis_acceleration_,
        translation_acceleration: translation_acceleration_,
    })

    return QuantityCoordinateVector.from_expr(value)
