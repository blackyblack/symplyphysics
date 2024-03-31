from sympy import (Eq, solve, exp)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The barometric formula determines the dependence of the pressure or density of a gas on the height in the gravity field.

## Law: p = p_0 * exp((-g * m * (h - h_0)) / (R * T))
## Where:
## p is the gas pressure in the layer located at height h
## p_0 is the gas pressure at the initial height
## g is acceleration of free fall
## m is atomic weight of the gas
## h is final height
## h_0 is initial height
## R is ideal gas constant
## T is temperature

## Conditions
## Gas is ideal
## The gas is in a uniform gravity field

final_pressure = Symbol("final_pressure", units.pressure)
initial_pressure = Symbol("initial_pressure", units.pressure)
atomic_weight = Symbol("atomic_weight", units.mass / units.amount_of_substance)
final_height = Symbol("final_height", units.length)
initial_height = Symbol("initial_height", units.length)
temperature = Symbol("temperature", units.temperature)

law = Eq(
    final_pressure,
    initial_pressure * exp(-units.acceleration_due_to_gravity * atomic_weight *
    (final_height - initial_height) / (units.molar_gas_constant * temperature)))


def print_law() -> str:
    return print_expression(law)


@validate_input(initial_pressure_=initial_pressure,
    atomic_weight_=atomic_weight,
    final_height_=final_height,
    initial_height_=initial_height,
    temperature_=temperature)
@validate_output(final_pressure)
def calculate_final_pressure(initial_pressure_: Quantity, atomic_weight_: Quantity,
    initial_height_: Quantity, final_height_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, final_pressure, dict=True)[0][final_pressure]
    result_final_pressure = result_expr.subs({
        initial_pressure: initial_pressure_,
        atomic_weight: atomic_weight_,
        final_height: final_height_,
        initial_height: initial_height_,
        temperature: temperature_
    })
    return Quantity(result_final_pressure)
