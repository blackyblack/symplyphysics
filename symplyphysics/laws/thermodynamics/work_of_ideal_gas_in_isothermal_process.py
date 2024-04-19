from sympy import (Eq, solve, log)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols, clone_symbol)

## Description
## The isothermal process of expansion (or compression) of a gas can occur under conditions where heat exchange between the gas and the external environment is carried out at a constant temperature difference.
## To do this, the heat capacity of the external environment must be large enough, and the expansion (or compression) process must be slow enough.

## Law: A = (m / mu) * R * T * ln(V_2 / V_1)
## Where:
## A is work of gas
## m is mass
## mu is molar mass
## V_1 is start volume
## V_2 is final volume
## R is ideal gas constant,
## T is temperature

## Note
## Since the internal energy of an ideal gas does not change during the isothermal process, all the heat received from the source is converted into work.
## When gas is compressed, the work of external forces is positive, and when expanding, the gas does the positive work.

work = Symbol("work", units.energy)
molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
start_volume = Symbol("start_volume", units.volume)
final_volume = Symbol("final_volume", units.volume)
temperature = symbols.thermodynamics.temperature
gas_mass = clone_symbol(symbols.basic.mass, "gas_mass")

law = Eq(work, (gas_mass / molar_mass) * units.molar_gas_constant * temperature *
    log(final_volume / start_volume))


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_=gas_mass,
    molar_mass_=molar_mass,
    start_volume_=start_volume,
    final_volume_=final_volume,
    temperature_=temperature)
@validate_output(work)
def calculate_work(mass_: Quantity, molar_mass_: Quantity, start_volume_: Quantity,
    final_volume_: Quantity, temperature_: Quantity) -> Quantity:
    solved = solve(law, work, dict=True)[0][work]
    result_expr = solved.subs({
        gas_mass: mass_,
        molar_mass: molar_mass_,
        start_volume: start_volume_,
        final_volume: final_volume_,
        temperature: temperature_
    })
    return Quantity(result_expr)
