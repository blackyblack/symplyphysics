from sympy import Eq, Derivative
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

# Description
## Maxwell relations are a set of equations that unite the most common thermodynamic quantities between
## one another. They are derived from the fundamental themodynamic relations featuring differentials
## of thermodynamic potentials, and this method of derivation is called the method of thermodynamic
## potentials.

# Law: (dS/dp)_T = -(dV/dT)_p
## S - entropy
## p - pressure
## V - volume
## T - absolute temperature

# Conditions
## - Changes in article count are not taken into account and it is assumed to be constant.

entropy = Function("entropy", units.energy / units.temperature)
pressure = Symbol("pressure", units.pressure)
volume = Function("volume", units.volume)
temperature = symbols.thermodynamics.temperature

law = Eq(
    Derivative(entropy(temperature, pressure), pressure),
    -1 * Derivative(volume(temperature, pressure), temperature)
)


@validate_input(
    volume_change_=volume,
    temperature_change_=temperature,
)
@validate_output(units.energy / (units.temperature * units.pressure))
def calculate_entropy_differential(
    volume_change_: Quantity,
    temperature_change_: Quantity,
) -> Quantity:
    volume_ = (volume_change_ / temperature_change_) * temperature
    result = law.rhs.subs(volume(temperature, pressure), volume_).doit()
    return Quantity(result)
