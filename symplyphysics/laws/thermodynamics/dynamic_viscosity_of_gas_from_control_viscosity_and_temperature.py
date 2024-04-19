from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, print_expression,
    validate_input, validate_output)

# Description
## Viscosity is one of the phenomena of transport, the property of fluid bodies (liquids and gases) to resist the movement of one part of them relative to another.
## As a result, the macroscopic work expended on this movement is dissipated as heat.
## Unlike liquids, the viscosity of gases increases with increasing temperature (for liquids, it decreases with increasing temperature).

## Law: mu = mu_0 * ((T_0 + C) / (T + C)) * (T / T_0)^1.5
## Where:
## T is set temperature
## T_0 is control temperature
## mu is dynamic viscosity at a set temperature T
## mu_0 is control viscosity at a certain control temperature T_0
## C is the Sutherland constant for the gas whose viscosity needs to be determined

## Conditions
## Gas is ideal

set_temperature = clone_symbol(symbols.thermodynamics.temperature, "set_temperature")
control_temperature = clone_symbol(symbols.thermodynamics.temperature, "control_temperature")
dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
control_viscosity = Symbol("control_viscosity", units.pressure * units.time)
sutherland_constant = Symbol("sutherland_constant", units.temperature)

law = Eq(
    dynamic_viscosity,
    control_viscosity * ((control_temperature + sutherland_constant) /
    (set_temperature + sutherland_constant)) * (set_temperature / control_temperature)**1.5)


def print_law() -> str:
    return print_expression(law)


@validate_input(control_viscosity_=control_viscosity,
    control_temperature_=control_temperature,
    sutherland_constant_=sutherland_constant,
    set_temperature_=set_temperature)
@validate_output(dynamic_viscosity)
def calculate_dynamic_viscosity(control_viscosity_: Quantity, control_temperature_: Quantity,
    sutherland_constant_: Quantity, set_temperature_: Quantity) -> Quantity:
    result_expr = solve(law, dynamic_viscosity, dict=True)[0][dynamic_viscosity]
    result_dynamic_viscosity = result_expr.subs({
        control_viscosity: control_viscosity_,
        control_temperature: control_temperature_,
        sutherland_constant: sutherland_constant_,
        set_temperature: set_temperature_
    })
    return Quantity(result_dynamic_viscosity)
