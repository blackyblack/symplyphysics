from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)

# Description
## Temperature of a thermodynamic system can be found if the function of the internal energy is known.

# Law: T = (dU/dS)_V
## T - absolute temperature
## U - internal energy
## S - entropy
## V - volume
## (d/dS)_V - derivative with respect to entropy at constant volume

temperature = symbols.thermodynamics.temperature
internal_energy = Function("internal_energy", units.energy)
entropy = Symbol("entropy", units.energy / units.temperature)
volume = Symbol("volume", units.volume)

law = Eq(temperature, Derivative(internal_energy(entropy, volume), entropy))

# TODO: derive from fundamental relation for internal energy


def print_law() -> str:
    return print_expression(law)


@validate_input(
    internal_energy_change_=internal_energy,
    entropy_change_=entropy,
)
@validate_output(temperature)
def calculate_temperature(
    internal_energy_change_: Quantity,
    entropy_change_: Quantity,
) -> Quantity:
    internal_energy_ = (internal_energy_change_ / entropy_change_) * entropy
    result = law.rhs.subs(
        internal_energy(entropy, volume), internal_energy_
    ).doit()
    return Quantity(result)
