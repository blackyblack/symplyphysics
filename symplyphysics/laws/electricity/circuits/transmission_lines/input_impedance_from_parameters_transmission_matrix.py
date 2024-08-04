from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)
from symplyphysics.core.dimensions import assert_equivalent_dimension

## Description
## Knowing the transfer matrix of the device, it is possible to determine the input impedance of this device.
## The transmission parameters matrix is one of the ways to describe a microwave device. The ABCD-parameters of the device act as elements
## of this matrix.

## Law is: Zin = (A * Zl + B) / (C * Zl + D), where
## Zin - input impedance of the transmission line,
## Zl - load resistance,
## A - ratio of input voltage to output voltage at idle at the output,
## B - ratio of input voltage to output current in case of a short circuit at the output,
## C - ratio of input current to output voltage at idle at the output,
## D - ratio of input current to output current in case of a short circuit at the output.

input_impedance = Symbol("input_impedance", units.impedance)

load_resistance = Symbol("load_resistance", units.impedance)
parameter_voltage_to_voltage = Symbol("parameter_voltage_to_voltage", dimensionless)
parameter_impedance = Symbol("parameter_impedance", units.impedance)
parameter_conductance = Symbol("parameter_conductance", units.conductance)
parameter_current_to_current = Symbol("parameter_current_to_current", dimensionless)

law = Eq(input_impedance, (parameter_voltage_to_voltage * load_resistance + parameter_impedance) /
    (parameter_conductance * load_resistance + parameter_current_to_current))


def print_law() -> str:
    return print_expression(law)


@validate_input(load_resistance_=load_resistance)
@validate_output(input_impedance)
def calculate_input_impedance(
        load_resistance_: Quantity, parameters_: tuple[tuple[float, Quantity], tuple[Quantity,
    float]]) -> Quantity:
    assert_equivalent_dimension(parameters_[0][1], "parameters_[0][1]", "calculate_input_impedance",
        units.impedance)
    assert_equivalent_dimension(parameters_[1][0], "parameters_[1][0]", "calculate_input_impedance",
        units.conductance)
    result_expr = solve(law, input_impedance, dict=True)[0][input_impedance]
    result_expr = result_expr.subs({
        load_resistance: load_resistance_,
        parameter_voltage_to_voltage: parameters_[0][0],
        parameter_impedance: parameters_[0][1],
        parameter_conductance: parameters_[1][0],
        parameter_current_to_current: parameters_[1][1],
    })
    return Quantity(result_expr)
