from sympy import (Eq, solve, exp)
from sympy.physics.units import molar_gas_constant
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    dimensionless
)

# Description
## The equilibrium constant is a value that determines for a given chemical reaction the ratio between
## thermodynamic activities (or, depending on the conditions of the reaction, partial pressures,
## concentrations or fugitives) of starting substances and products in a state of chemical equilibrium.

## Law is: K = exp(-delta_G / R * T), where
## K - equilibrium constant of reaction,
## delta_G - change of isobaric potential of reaction,
## R - molar gas constant,
## T - temperature.

# Conditions:
## - the process is isobaric-isothermal.

equilibrium_constant = Symbol("equilibrium_constant", dimensionless)

standart_change_isobaric_potential = Symbol("standart_change_isobaric_potential", units.energy / units.amount_of_substance)
temperature = Symbol("temperature", units.temperature)

celsius_to_kelvin = Quantity(273.15 * units.kelvin)

law = Eq(equilibrium_constant, exp(-standart_change_isobaric_potential / (molar_gas_constant * temperature)))


def print_law() -> str:
    return print_expression(law)


@validate_input(standart_change_isobaric_potential_=standart_change_isobaric_potential,
    temperature_=temperature)
@validate_output(equilibrium_constant)
def calculate_equilibrium_constant(standart_change_isobaric_potential_: Quantity,
    temperature_: Quantity) -> float:
    result_expr = solve(law, equilibrium_constant, dict=True)[0][equilibrium_constant]
    result_expr = result_expr.subs({
        standart_change_isobaric_potential: standart_change_isobaric_potential_,
        temperature: temperature_
    })
    return Quantity(result_expr)
