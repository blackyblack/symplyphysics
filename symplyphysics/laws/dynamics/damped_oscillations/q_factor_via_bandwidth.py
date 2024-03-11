from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Q factor (quality factor) is a dimensionless parameter that describes how underdamped 
## an oscillator or resonator is: the larger the Q factor is, the less damped it is. There 
## are two nearly equivalent definitions of it that become approximately equivalent as Q 
## becomes larger, meaning that the resonator becomes less damped.

# Law: Q = f_r / delta_f
## Q - Q factor
## f_r - oscillator's resonant frequency
## delta_f - resonance width, or full width at half maximum, i.e. the bandwidth over which
##           the power of vibration is greater than half the power at the resonant frequency

# Note
## - An equivalent definition uses angular frequencies instead of linear ones.

q_factor = Symbol("q_factor", dimensionless)
resonant_frequency = Symbol("resonant_frequency", units.frequency)
resonance_width = Symbol("resonance_width", units.frequency)

law = Eq(q_factor, resonant_frequency / resonance_width)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    resonant_frequency_=resonant_frequency,
    resonance_width_=resonance_width,
)
@validate_output(q_factor)
def calculate_q_factor(
    resonant_frequency_: Quantity,
    resonance_width_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        resonant_frequency: resonant_frequency_,
        resonance_width: resonance_width_,
    })
    return Quantity(result)
