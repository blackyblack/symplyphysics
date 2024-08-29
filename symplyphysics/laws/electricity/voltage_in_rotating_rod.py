from sympy import (Eq, solve)
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    angle_type,
)

# Description
## Let the rod rotate in a uniform magnetic field. The plane of rotation is perpendicular to the magnetic
## field lines. The axis of rotation passes through one of the ends of the rod. Then the voltage generated
## at the ends of the rod depends on the magnitude of the magnetic induction, the rotation frequency and
## the length of the rod.

## Law is: V = (1 / 2) * B * w * L^2, where
## V - voltage,
## B - magnetic induction,
## w - rotation frequency,
## L - rod length.

voltage = Symbol("voltage", units.voltage)

magnetic_induction = Symbol("magnetic_induction", units.magnetic_density)
rotation_frequency = Symbol("rotation_frequency", angle_type / units.time)
rod_length = Symbol("rod_length", units.length)

law = Eq(voltage, (1 / 2) * magnetic_induction * rotation_frequency * rod_length**2)


@validate_input(magnetic_induction_=magnetic_induction,
    rotation_frequency_=rotation_frequency,
    rod_length_=rod_length)
@validate_output(voltage)
def calculate_voltage(magnetic_induction_: Quantity, rotation_frequency_: Quantity,
    rod_length_: Quantity) -> Quantity:
    result_expr = solve(law, voltage, dict=True)[0][voltage]
    result_expr = result_expr.subs({
        magnetic_induction: magnetic_induction_,
        rotation_frequency: rotation_frequency_,
        rod_length: rod_length_
    })
    return Quantity(result_expr)
