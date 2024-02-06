from sympy import Eq, cos
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    angle_type,
)
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## Power is the rate at which the work is done by a force. The instantaneous power is the power
## at a specific instant.

# Law: P = F * v * cos(phi)
## P - instantaneous power of force F
## F - magnitude of force
## v - speed
## phi - angle between force and velocity vectors

power = Symbol("power", units.power)
force = Symbol("force", units.force)
speed = Symbol("speed", units.speed)
angle = Symbol("angle", angle_type)

law = Eq(power, force * speed * cos(angle))


def print_law() -> str:
    return print_expression(law)


@validate_input(force_=force, speed_=speed, angle_=angle)
@validate_output(power)
def calculate_power(force_: Quantity, speed_: Quantity, angle_: Quantity | float) -> Quantity:
    angle_ = scale_factor(angle_)
    result = law.rhs.subs({
        force: force_,
        speed: speed_,
        angle: angle_,
    })
    return Quantity(result)
