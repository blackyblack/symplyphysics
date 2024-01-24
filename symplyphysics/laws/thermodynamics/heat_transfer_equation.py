from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
# heat flow calculation Q = k * (T1-T2)/ d * S
# where:
# k = material thermal conductivity, ability of a material to transfer heat
# T1 = inside temperature
# T2 = outside temperature
# d = wall thickness
# S = wall area
# Note: heat flow is measured in Watts
# 
# This formula describes Fourier's law of heat transfer and is used to calculate the heat flow through a material or structure.
# Fourier's formula is used to determine the amount of heat that is transferred through a material from one side to the other due to the material's thermal conductivity, temperature gradient, material thickness and cross-sectional area. 
# The cross-sectional area can be any shape, including the surface area of ​​a circle.
# The heat flow sign indicates the direction of heat transfer and reflects the direction of heat flow from one temperature to another.
# If heat flow is positive, it means that thermal energy moves from a higher temperature to a lower temperature.
# If heat flow is negative, it means that thermal energy moves from a lower temperature to a higher temperature.

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

