from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
# heat flow calculation Ф = λ * (T1-T2)/ d * S
# where:
# λ = material_thermal_conductivity (Watts/(K*m))
# T1 = inside_temperature (C')
# T2 = outside_temperature (C')
# d = wall_thickness (m)
# S = wall_area (m2)
# Note: heat flow (Ф) is measured in Watts
# 
# This formula describes Fourier's law for heat transfer and is used to calculate heat flow (Ф) through a material or structure.
# The Fourier formula is used to determine the amount of heat that is transferred through a material from one side to the other due to the material's thermal conductivity, temperature gradient, material thickness, and cross-sectional area.

heat_flow = Symbol("heat_flow", units.power)
material_thermal_conductivity = Symbol("material_thermal_conductivity", units.power / (units.temperature * units.length))
inside_temperature = Symbol("inside_temperature", units.temperature)
outside_temperature = Symbol("outside_temperature", units.temperature)
wall_thickness = Symbol("wall_thickness", units.length)
wall_area = Symbol("wall_area", units.area)

law = Eq(heat_flow, material_thermal_conductivity * (inside_temperature - outside_temperature) / wall_thickness * wall_area)

def print_law() -> str:
    return print_expression(law)

@validate_input(material_thermal_conductivity_=material_thermal_conductivity,
    inside_temperature_=inside_temperature,
    outside_temperature_=outside_temperature,
    wall_thickness_=wall_thickness,
    wall_area_=wall_area)
@validate_output(heat_flow)
def calculate_heat_flow(material_thermal_conductivity_: Quantity, inside_temperature_: Quantity,
    outside_temperature_: Quantity, wall_thickness_: Quantity, wall_area_: Quantity) -> Quantity:

    result_heat_flow_expr = solve(law, heat_flow, dict=True)[0][heat_flow]
    result_expr = result_heat_flow_expr.subs({
        material_thermal_conductivity: material_thermal_conductivity_,
        inside_temperature: inside_temperature_,
        outside_temperature: outside_temperature_,
        wall_thickness: wall_thickness_,
        wall_area: wall_area_
    })
    return Quantity(result_expr)

