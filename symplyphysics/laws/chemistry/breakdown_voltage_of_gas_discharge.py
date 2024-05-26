from sympy import Eq, solve, log
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless
)

# Description
## The potential difference between the electrodes, at which the discharge from an independent one passes into an independent one,
## is called a breakdown voltage. In the independent mode, the secondary emission of electrons by the cathode under the action of
## ions bombarding it is very significant. Due to the formation of ions in the volume and the knocking out of secondary electrons by
## them, the discharge ceases to depend on external ionizing influences, it switches to a self-sustaining mode and becomes independent.

## Law is: U = - B * p * L / ln(ln(1 / y + 1) / (A * p * L)), where
## U - breakdown voltage of the gas discharge,
## B - the second constant of gas,
## p - pressure
## L - distance between electrodes,
## y - the coefficient of secondary emission of the ion-electronic type,
## A - the first constant of gas.

breakdown_voltage = Symbol("breakdown_voltage", units.voltage)

first_constant_of_gas = Symbol("first_constant_of_gas", 1 / units.length / units.pressure)
second_constant_of_gas = Symbol("second_constant_of_gas", units.voltage / units.length / units.pressure)
pressure = Symbol("pressure", units.pressure)
distance_between_electrodes = Symbol("distance_between_electrodes", units.length)
secondary_emission_factor = Symbol("secondary_emission_factor", dimensionless)

expression_1 = second_constant_of_gas * pressure * distance_between_electrodes
expression_2 = log(1 / secondary_emission_factor + 1)
expression_3 = log(expression_2 / (first_constant_of_gas * pressure * distance_between_electrodes))
law = Eq(breakdown_voltage, - expression_1 / expression_3)


@validate_input(first_constant_of_gas_=first_constant_of_gas,
    second_constant_of_gas_=second_constant_of_gas,
    pressure_=pressure,
    distance_between_electrodes_=distance_between_electrodes,
    secondary_emission_factor_=secondary_emission_factor)
@validate_output(breakdown_voltage)
def calculate_breakdown_voltage(first_constant_of_gas_: Quantity,
    second_constant_of_gas_: Quantity, pressure_: Quantity,
    distance_between_electrodes_: Quantity, secondary_emission_factor_: float) -> Quantity:
    result_expr = solve(law, breakdown_voltage, dict=True)[0][breakdown_voltage]
    result_expr = result_expr.subs({
        first_constant_of_gas: first_constant_of_gas_,
        second_constant_of_gas: second_constant_of_gas_,
        pressure: pressure_,
        distance_between_electrodes: distance_between_electrodes_,
        secondary_emission_factor: secondary_emission_factor_,
    })
    return Quantity(result_expr)
