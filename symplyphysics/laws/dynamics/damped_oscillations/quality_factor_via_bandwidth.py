from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
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
## delta_f - resonance width, or full width at half maximum, of the peak in the graph of the
##           dissipated power as a function of driving frequency, i.e. the difference between
##           the frequencies at which the dissipated power is half the peak dissipated power,
##           which happens ad the resonant frequency, see [this figure](http://spiff.rit.edu/classes/phys283/lectures/forced_ii/half_power.png)

# Note
## - An equivalent definition uses angular frequencies instead of linear ones.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Q_factor#Bandwidth_definition>

quality_factor = Symbol("quality_factor", dimensionless)
resonant_frequency = Symbol("resonant_frequency", units.frequency)
resonance_width = Symbol("resonance_width", units.frequency)

law = Eq(quality_factor, resonant_frequency / resonance_width)


@validate_input(
    resonant_frequency_=resonant_frequency,
    resonance_width_=resonance_width,
)
@validate_output(quality_factor)
def calculate_quality_factor(
    resonant_frequency_: Quantity,
    resonance_width_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        resonant_frequency: resonant_frequency_,
        resonance_width: resonance_width_,
    })
    return Quantity(result)
