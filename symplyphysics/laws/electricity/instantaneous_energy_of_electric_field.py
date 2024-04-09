from sympy import (Eq, solve, cos)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)

# Description
## There is an oscillatory circuit with alternating current. Then the energy of the electric field will depend on the inductance,
## the maximum value of the current, the frequency of the current, the time and the initial phase.

## Law is: W = (L * I^2 / 2) * cos(w * t + phi)^2, where
## W - energy,
## L - inductance,
## I - maximum current,
## w - frequency,
## t - time,
## phi - initial phase.

energy = Symbol("energy", units.energy)

inductance = Symbol("inductance", units.inductance)
maximum_current = Symbol("maximum_current", units.current)
frequency = Symbol("frequency", angle_type / units.time)
time = Symbol("time", units.time)
initial_phase = Symbol("initial_phase", angle_type)

law = Eq(energy, (inductance * maximum_current**2 / 2) * cos(frequency * time + initial_phase)**2)


def print_law() -> str:
    return print_expression(law)


@validate_input(inductance_=inductance,
    maximum_current_=maximum_current,
    frequency_=frequency,
    time_=time,
    initial_phase_=initial_phase)
@validate_output(energy)
def calculate_energy(inductance_: Quantity, maximum_current_: Quantity, frequency_: Quantity,
    time_: Quantity, initial_phase_: float | Quantity) -> Quantity:
    result_expr = solve(law, energy, dict=True)[0][energy]
    result_expr = result_expr.subs({
        inductance: inductance_,
        maximum_current: maximum_current_,
        frequency: frequency_,
        time: time_,
        initial_phase: initial_phase_
    })
    return Quantity(result_expr)
