from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    Symbol,
    dimensionless,
    print_expression,
    units,
    validate_input,
    validate_output,
)
from symplyphysics.core.convert import convert_to_dimensionless

# Description
# The Prandtl number (Pr) is a dimensionless number, defined as the ratio of
# momentum diffusivity to thermal diffusivity.
# Law: Pr = Cp * mu / k, where
# Cp is heat capacity,
# mu is dynamic viscosity,
# k is thermal conductivity,
# Pr is Prandtl number.

heat_capacity = Symbol("heat_capacity", units.energy / units.mass / units.temperature)
dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
thermal_conductivity = Symbol("thermal_conductivity",
    units.power / units.length / units.temperature)
prandtl_number = Symbol("prandtl_number", dimensionless)

law = Eq(prandtl_number, heat_capacity * dynamic_viscosity / thermal_conductivity)


def print_law() -> str:
    return print_expression(law)


@validate_input(heat_capacity_=heat_capacity,
    dynamic_viscosity_=dynamic_viscosity,
    thermal_conductivity_=thermal_conductivity)
@validate_output(prandtl_number)
def calculate_prandtl_number(heat_capacity_: Quantity, dynamic_viscosity_: Quantity,
    thermal_conductivity_: Quantity) -> float:
    result_expr = solve(law, prandtl_number, dict=True)[0][prandtl_number]
    result_applied = result_expr.subs({
        heat_capacity: heat_capacity_,
        dynamic_viscosity: dynamic_viscosity_,
        thermal_conductivity: thermal_conductivity_
    })
    result = Quantity(result_applied)
    return convert_to_dimensionless(Quantity(result))
