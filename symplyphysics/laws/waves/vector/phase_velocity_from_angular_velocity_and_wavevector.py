"""
Phase velocity from angular velocity and wavevector
===================================================

The phase velocity of a wave is the rate at which the wave propagates in a medium. It is the
velocity at which the phase of one frequency component of the wave travels. The phase velocity is
collinear with the wavevector.

**Notes:**

#. Angular wavevector is a vector used in describing a wave. Its magnitude is the angular
   wavenumber of the wave. Its direction is perpendicular to the wavefront, and in isotropic media
   it is also the direction of wave propagation.

#. Also see the :ref:`scalar counterpart <Phase speed from angular frequency and wavenumber>` of
   this law.

**Links:**

#. `Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Phase_velocity>`__.
"""

from sympy import Eq, solve
from symplyphysics import Quantity, validate_input, validate_output, symbols, clone_as_symbol

from symplyphysics.core.experimental.vectors import clone_as_vector_symbol, VectorNorm
from symplyphysics.core.experimental.coordinate_systems import QuantityCoordinateVector
from symplyphysics.core.experimental.solvers import vector_equals, apply, solve_for_vector

phase_velocity = clone_as_vector_symbol(symbols.phase_speed)
"""
Vector of the phase velocity of the wave. See :symbols:`phase_speed`.
"""

angular_frequency = clone_as_symbol(symbols.angular_frequency, positive=True)
"""
:symbols:`angular_frequency` of the wave.
"""

angular_wavevector = clone_as_vector_symbol(symbols.angular_wavenumber)
"""
Angular wavevector of the wave. See :symbols:`angular_wavenumber`.
"""

phase_velocity_law_ = Eq(
    phase_velocity,
    (angular_frequency / VectorNorm(angular_wavevector)**2) * angular_wavevector,
)
"""
:laws:symbol::

:laws:latex::
"""

angular_wavevector_law = Eq(
    angular_wavevector,
    (angular_frequency / VectorNorm(phase_velocity)**2) * phase_velocity,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive `angular_wavevector_law` from `phase_velocity_law_`

_phase_velocity_eqn_norm = apply(phase_velocity_law_, VectorNorm)

_angular_wavenumber_expr = solve(
    _phase_velocity_eqn_norm,
    VectorNorm(angular_wavevector),
)[0]

_angular_wavevector_expr = solve_for_vector(
    phase_velocity_law_,
    angular_wavevector,
).subs(
    VectorNorm(angular_wavevector),
    _angular_wavenumber_expr,
)

assert vector_equals(_angular_wavevector_expr, angular_wavevector_law.rhs)


@validate_input(
    angular_frequency_=angular_frequency,
    wavevector_=angular_wavevector,
)
@validate_output(phase_velocity)
def calculate_phase_velocity(
    angular_frequency_: Quantity,
    wavevector_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = phase_velocity_law_.rhs.subs({
        angular_frequency: angular_frequency_,
        angular_wavevector: wavevector_,
    })

    return QuantityCoordinateVector.from_expr(result)


@validate_input(
    angular_frequency_=angular_frequency,
    phase_velocity_=phase_velocity,
)
@validate_output(angular_wavevector)
def calculate_wavevector(
    angular_frequency_: Quantity,
    phase_velocity_: QuantityCoordinateVector,
) -> QuantityCoordinateVector:
    result = angular_wavevector_law.rhs.subs({
        angular_frequency: angular_frequency_,
        phase_velocity: phase_velocity_,
    })

    return QuantityCoordinateVector.from_expr(result)
