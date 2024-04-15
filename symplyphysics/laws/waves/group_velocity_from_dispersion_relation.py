from sympy import Eq, Derivative
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function, Point2D

# Description
## Waves can form a group, called wave packets. The velocity with which a wave packet travels
## is called group velocity. In other words, it is the velocity with which the overall envelope
## shape of the wave's amplitudes - called envelope of modulation of the wave -  propagates
## through space.

# Law: v_group = dw(k)/dk
## v_group - group velocity of a collection of waves
## w = w(k) - dispersion relation which describes how the wave's angular frequency depends on wavenumber
## k - wavenumber
## d/dk - partial derivative w.r.t. k

group_velocity = Symbol("group_velocity", units.velocity, real=True)
dispersion_relation = Function("dispersion_relation", angle_type / units.time, real=True)
wavenumber = Symbol("wavenumber", angle_type / units.length, positive=True)

law = Eq(
    group_velocity,
    Derivative(dispersion_relation(wavenumber), wavenumber),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    angular_frequency_before_=dispersion_relation,
    angular_frequency_after_=dispersion_relation,
    wavenumber_before_=wavenumber,
    wavenumber_after_=wavenumber,
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
        wavenumber,
    )
    result = law.rhs.subs(
        dispersion_relation(wavenumber),
        dispersion_relation_,
    ).doit()
    return Quantity(result)
