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
    units,
    angle_type,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

group_velocity = Symbol("group_velocity", units.velocity, real=True)
r"""
Group velocity of the wave packet.

Symbol:
    :code:`v_group`

Latex:
    :math:`v_\text{group}`
"""

angular_frequency = Function("angular_frequency", angle_type / units.time, real=True)
r"""
Angular frequency of the wave as a function of angular wavenumber, also called the
*dispersion relation* of the wave.

Symbol:
    :code:`w(k)`

Latex:
    :math:`\omega(k)`
"""

angular_wavenumber = Symbol("angular_wavenumber", angle_type / units.length, positive=True)
"""
Angular wavenumber of the wave.

Symbol:
    :code:`k`
"""

law = Eq(
    group_velocity,
    Derivative(angular_frequency(angular_wavenumber), angular_wavenumber),
)
r"""
:code:`v_group = Derivative(w(k), k)`

Latex:
    .. math::
        v_\text{group} = \frac{\partial \omega}{\partial k}
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
