from sympy import (Eq, solve, sin)
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless, angle_type)

# Description
## With uniform rotation of a multi-turn frame with the same areas of turns in a uniform magnetic field,
## an induction EMF occurs in each turn of this frame, so that the total EMF is numerically equal to the
## product of the number of turns of the frame N, the induction of the magnetic field, the area of one
## turn of the frame, the angular velocity of rotation of the frame and the sine of the angle of rotation
## of the frame.
## At the initial moment of time, the contour plane is perpendicular to the field lines of force and the
## frame rotates around an axis perpendicular to the field lines and passing through the middle of the
## opposite sides of the frame.
## EMF - electromotive force.

## Law is: U = -N * B * S * w * sin(w * t), where
## U - voltage (electromotive force),
## N - number of turns,
## B - magnetic field induction,
## S - contour area,
## w - rotation frequency,
## t - time.

voltage = Symbol("voltage", units.voltage)

number_turns = Symbol("number_turns", dimensionless)
induction = Symbol("induction", units.magnetic_density)
contour_area = Symbol("contour_area", units.area)
rotation_frequency = Symbol("rotation_frequency", angle_type / units.time)
time = Symbol("time", units.time)

law = Eq(
    voltage,
    -number_turns * induction * contour_area * rotation_frequency * sin(rotation_frequency * time))


@validate_input(number_turns_=number_turns,
    induction_=induction,
    contour_area_=contour_area,
    rotation_frequency_=rotation_frequency,
    time_=time)
@validate_output(voltage)
def calculate_voltage(number_turns_: int, induction_: Quantity, contour_area_: Quantity,
    rotation_frequency_: Quantity, time_: Quantity) -> Quantity:
    result_expr = solve(law, voltage, dict=True)[0][voltage]
    result_expr = result_expr.subs({
        number_turns: number_turns_,
        induction: induction_,
        contour_area: contour_area_,
        rotation_frequency: rotation_frequency_,
        time: time_
    })
    return Quantity(result_expr)
