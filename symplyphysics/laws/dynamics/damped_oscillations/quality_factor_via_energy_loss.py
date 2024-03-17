from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
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

# Law: Q = w_r * (E_stored / P_loss)
## Q - Q factor of oscillator
## w_r - resonant angular frequency of oscilaltor
## E_stored - energy stored in oscillator
## P_loss - power loss of oscillator, i.e. rate of energy dissipation from it

quality_factor = Symbol("quality_factor", dimensionless)
resonant_angular_frequency = Symbol("resonant_angular_frequency", angle_type / units.time)
energy_stored = Symbol("energy_stored", units.energy)
power_loss = Symbol("power_loss", units.power)

law = Eq(quality_factor, resonant_angular_frequency * (energy_stored / power_loss))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    resonant_angular_frequency_=resonant_angular_frequency,
    energy_stored_=energy_stored,
    power_loss_=power_loss,
)
@validate_output(quality_factor)
def calculate_quality_factor(
    resonant_angular_frequency_: Quantity,
    energy_stored_: Quantity,
    power_loss_: Quantity,
) -> float:
    result = law.rhs.subs({
        resonant_angular_frequency: resonant_angular_frequency_,
        energy_stored: energy_stored_,
        power_loss: power_loss_,
    })
    return Quantity(result).scale_factor
