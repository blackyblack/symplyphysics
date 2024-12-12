from sympy import Eq, solve, exp
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## At a certain voltage, the gas discharge becomes independent. To find this voltage, it is necessary to know the volumetric ionization coefficient.
## And it, in turn, depends on the energy distribution of electrons and can be approximated by the expression below.

## Law is: a = A * p * exp(-B * p / E), where
## a - volumetric ionization coefficient,
## A - the first constant of gas,
## p - pressure,
## B - the second constant of gas,
## E - electric intensity.

#. Links: formula 1.12 <https://studfile.net/preview/3079348/page:3/>
# TODO: find English link

# TODO: move to `ionization` folder?

coefficient_of_volumetric_ionization = Symbol("coefficient_of_volumetric_ionization",
    1 / units.length)

first_constant_of_gas = Symbol("first_constant_of_gas", 1 / units.length / units.pressure)
second_constant_of_gas = Symbol("second_constant_of_gas",
    units.voltage / units.length / units.pressure)
pressure = Symbol("pressure", units.pressure)
electric_intensity = Symbol("electric_intensity", units.voltage / units.length)

law = Eq(
    coefficient_of_volumetric_ionization,
    first_constant_of_gas * pressure * exp(-second_constant_of_gas * pressure / electric_intensity))


@validate_input(first_constant_of_gas_=first_constant_of_gas,
    second_constant_of_gas_=second_constant_of_gas,
    pressure_=pressure,
    electric_intensity_=electric_intensity)
@validate_output(coefficient_of_volumetric_ionization)
def calculate_coefficient_of_volumetric_ionization(first_constant_of_gas_: Quantity,
    second_constant_of_gas_: Quantity, pressure_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, coefficient_of_volumetric_ionization,
        dict=True)[0][coefficient_of_volumetric_ionization]
    result_expr = result_expr.subs({
        first_constant_of_gas: first_constant_of_gas_,
        second_constant_of_gas: second_constant_of_gas_,
        pressure: pressure_,
        electric_intensity: electric_intensity_,
    })
    return Quantity(result_expr)
