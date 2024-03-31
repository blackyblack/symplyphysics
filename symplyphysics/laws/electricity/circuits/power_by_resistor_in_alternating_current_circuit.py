from sympy import (Eq, solve, sin)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)

# Description
## The resistor is in an alternating current circuit. Then the thermal power dissipated on the resistor will depend on the amplitude of the current,
## the frequency of the alternating current, the time and the resistance of the resistor.

## Law is: P = I^2 * R * sin(w * t)^2, where
## P - instantaneous power,
## I - current amplitude,
## R - resistance of the resistor,
## w - frequency of the current,
## t - time.

power = Symbol("power", units.power)

current_amplitude = Symbol("current_amplitude", units.current)
resistance = Symbol("resistance", units.impedance)
current_frequency = Symbol("current_frequency", angle_type / units.time)
time = Symbol("time", units.time)

law = Eq(power, current_amplitude**2 * resistance * sin(current_frequency * time)**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(current_amplitude_=current_amplitude,
    resistance_=resistance,
    time_=time,
    current_frequency_=current_frequency)
@validate_output(power)
def calculate_power(current_amplitude_: Quantity, resistance_: Quantity, time_: Quantity,
    current_frequency_: Quantity) -> Quantity:
    result_expr = solve(law, power, dict=True)[0][power]
    result_expr = result_expr.subs({
        current_amplitude: current_amplitude_,
        resistance: resistance_,
        time: time_,
        current_frequency: current_frequency_,
    })
    return Quantity(result_expr)
