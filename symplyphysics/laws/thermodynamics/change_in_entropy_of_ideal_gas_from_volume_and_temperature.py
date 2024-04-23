from sympy import (Eq, solve, log)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols, clone_symbol)

# Description
## Entropy is a property of a thermodynamic system that expresses the direction or outcome of spontaneous changes in the system.
## The term was introduced to explain the relationship of the internal energy that is available or unavailable for transformations in form of heat and work.
## Entropy predicts that certain processes are irreversible or impossible, despite not violating the conservation of energy.
## The definition of entropy is central to the establishment of the second law of thermodynamics, which states that the entropy of isolated systems cannot decrease with time,
## as they always tend to arrive at a state of thermodynamic equilibrium, where the entropy is highest.

## Law: S = (m / mu) * (C_v * ln(T_2 / T_1) + R * ln(V_2 / V_1))
## Where:
## S is the change in entropy during the transition from one state to another
## m is mass of gas
## mu is molar mass of gas
## C_v is molar heat capacity of a gas at constant volume
## T_2 is temperature in the final state
## T_1 is temperature in the initial state
## R is universal gas constant
## V_2 is volume in the final state
## V_1 is temperature in the initial state

## Conditions
## Gas is ideal

entropy_change = Symbol("entropy_change", units.energy / units.temperature)
molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
molar_heat_capacity = Symbol("molar_heat_capacity",
    units.energy / (units.temperature * units.amount_of_substance))
final_temperature = clone_symbol(symbols.thermodynamics.temperature, "final_temperature")
initial_temperature = clone_symbol(symbols.thermodynamics.temperature, "initial_temperature")
final_volume = Symbol("final_volume", units.volume)
start_volume = Symbol("start_volume", units.volume)
gas_mass = clone_symbol(symbols.basic.mass, "gas_mass")

law = Eq(entropy_change, (gas_mass / molar_mass) *
    ((molar_heat_capacity * log(final_temperature / initial_temperature)) +
    (units.molar_gas_constant * log(final_volume / start_volume))))


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=gas_mass,
    molar_mass_=molar_mass,
    molar_heat_capacity_=molar_heat_capacity,
    temperature_limits=initial_temperature,
    volume_limits=start_volume)
@validate_output(entropy_change)
def calculate_entropy_change(mass_: Quantity, molar_mass_: Quantity, molar_heat_capacity_: Quantity,
    temperature_limits: tuple[Quantity, Quantity], volume_limits: tuple[Quantity,
    Quantity]) -> Quantity:
    initial_temperature_, final_temperature_ = temperature_limits
    start_volume_, final_volume_ = volume_limits
    result_expr = solve(law, entropy_change, dict=True)[0][entropy_change]
    result_entropy_change = result_expr.subs({
        gas_mass: mass_,
        molar_mass: molar_mass_,
        molar_heat_capacity: molar_heat_capacity_,
        initial_temperature: initial_temperature_,
        final_temperature: final_temperature_,
        start_volume: start_volume_,
        final_volume: final_volume_,
    })
    return Quantity(result_entropy_change)
