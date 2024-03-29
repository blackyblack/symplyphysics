from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, symbols)

# Description
## Inner energy of ideal gas is sum of kinetic energy of all it's molecules.
## Inner energy of Ideal gas law: U = 1.5 * m * R * T / M
## Where:
## U is inner energy,
## m is mass of a gas,
## R is ideal gas constant,
## T is absolute temperature in Kelvins,
## M is mass of one mole of this gas.

# Conditions
## Gas is ideal - no potential energy of any molecules interaction.

inner_energy = Symbol("energy", units.energy)
temperature = Symbol("temperature", units.temperature)
mole_mass = Symbol("mole_mass", units.mass / units.amount_of_substance)

law = Eq(inner_energy,
    1.5 * symbols.basic.mass * units.molar_gas_constant * temperature / mole_mass)


def print_law() -> str:
    return print_expression(law)


@validate_input(mass_of_gas_=symbols.basic.mass, temperature_=temperature, mole_mass_=mole_mass)
@validate_output(inner_energy)
def calculate_inner_energy(mass_of_gas_: Quantity, temperature_: Quantity,
    mole_mass_: Quantity) -> Quantity:
    solved = solve(law, inner_energy, dict=True)[0][inner_energy]
    result_expr = solved.subs({
        symbols.basic.mass: mass_of_gas_,
        temperature: temperature_,
        mole_mass: mole_mass_
    })
    return Quantity(result_expr)
