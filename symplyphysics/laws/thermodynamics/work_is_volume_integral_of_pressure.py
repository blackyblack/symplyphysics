from sympy import Eq, Integral, Point2D
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)
from symplyphysics.core.geometry.line import two_point_function

# Description
## The pressure-volume work occurs when the volume of a system changes.

# Law: W = Integral(p(V), V)
## W - work done by gas
## p - gas pressure
## V - gas volume

# Conditions
## - The process is reversible or quasi-static.
## - The system is closed.

work = Symbol("work", units.energy)
pressure = Function("pressure", units.pressure)
volume = Symbol("volume", units.volume)
volume_before = Symbol("volume_before", units.volume)
volume_after = Symbol("volume_after", units.volume)

law = Eq(work, Integral(pressure(volume), (volume, volume_before, volume_after)))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    volume_before_=volume_before,
    volume_after_=volume_after,
    pressure_before_=pressure,
    pressure_after_=pressure,
)
@validate_output(work)
def calculate_work(
    volume_before_: Quantity,
    volume_after_: Quantity,
    pressure_before_: Quantity,
    pressure_after_: Quantity,
) -> Quantity:
    pressure_ = two_point_function(
        Point2D(volume_before_, pressure_before_),
        Point2D(volume_after_, pressure_after_),
        volume,
    )
    result = law.rhs.subs({
        pressure(volume): pressure_,
        volume_before: volume_before_,
        volume_after: volume_after_,
    }).doit()
    return Quantity(result)
