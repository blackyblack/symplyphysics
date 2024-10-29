"""
Group velocity from dispersion relation
=======================================

Waves can form a group, called wave packets. The speed with which a wave packet travels
is called group velocity. In other words, it is the speed with which the overall envelope
shape of the wave's amplitudes — called envelope of modulation of the wave — propagates
through space.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

group_velocity = symbols.group_speed
"""
:symbols:`group_speed` of the wave packet.
"""

angular_frequency = clone_as_function(
    symbols.angular_frequency,
    display_symbol="w(k)",
    real=True,
)
"""
:symbols:`angular_frequency` of the wave as a function of :symbols:`angular_wavenumber`, also called the
*dispersion relation* of the wave.
"""

angular_wavenumber = symbols.angular_wavenumber
"""
:symbols:`angular_wavenumber` of the wave.
"""

law = Eq(
    group_velocity,
    Derivative(angular_frequency(angular_wavenumber), angular_wavenumber),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    angular_frequency_before_=angular_frequency,
    angular_frequency_after_=angular_frequency,
    wavenumber_before_=angular_wavenumber,
    wavenumber_after_=angular_wavenumber,
)
@validate_output(group_velocity)
def calculate_group_velocity(
    angular_frequency_before_: Quantity,
    angular_frequency_after_: Quantity,
    wavenumber_before_: Quantity,
    wavenumber_after_: Quantity,
) -> Quantity:
    dispersion_relation_ = two_point_function(
        Point2D(wavenumber_before_, angular_frequency_before_),
        Point2D(wavenumber_after_, angular_frequency_after_),
        angular_wavenumber,
    )
    result = law.rhs.subs(
        angular_frequency(angular_wavenumber),
        dispersion_relation_,
    ).doit()
    return Quantity(result)
