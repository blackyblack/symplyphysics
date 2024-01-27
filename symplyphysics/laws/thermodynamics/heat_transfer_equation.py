from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
# Fourier's law Q = k * (T1-T2)/ d * S
# where:
# Q = heat flow 
# k = material thermal conductivity, ability of a material to transfer heat
# T1 = inside temperature
# T2 = outside temperature
# d = cross section thickness
# S = cross-sectional area
# 
# The Fourier formula is used to determine the amount of heat flow that is transferred through a material from one side to the other due to the material's thermal conductivity, temperature gradient, thickness and cross-sectional area. 
# The cross-sectional area can be any shape.
# And the formula is suitable for calculating the heat flow between the external and internal temperatures of the ball, indicating the thickness of the cross section of the ball and its surface area.
# The heat flow sign indicates the direction of heat transfer and reflects the direction of heat flow from one temperature to another.
# If heat flow is positive, it means that thermal energy moves from a higher temperature to a lower temperature.
# The minus sign means that the vectors of heat flow and temperature gradient are multidirectional. It should be understood that heat is transferred in the direction of decreasing temperature.

heat_flow = Symbol("heat_flow", units.power)
material_thermal_conductivity = Symbol("material_thermal_conductivity", units.power / (units.temperature * units.length))
inside_temperature = Symbol("inside_temperature", units.temperature)
outside_temperature = Symbol("outside_temperature", units.temperature)
cross_section_thickness = Symbol("cross_section_thickness", units.length)
cross_sectional_area = Symbol("cross_sectional_area", units.area)

law = Eq(heat_flow, material_thermal_conductivity * (inside_temperature - outside_temperature) / cross_section_thickness * cross_sectional_area)

def print_law() -> str:
    return print_expression(law)

@validate_input(material_thermal_conductivity_=material_thermal_conductivity,
    inside_temperature_=inside_temperature,
    outside_temperature_=outside_temperature,
    cross_section_thickness_=cross_section_thickness,
    cross_sectional_area_=cross_sectional_area)
@validate_output(heat_flow)
def calculate_heat_flow(material_thermal_conductivity_: Quantity, inside_temperature_: Quantity,
    outside_temperature_: Quantity, cross_section_thickness_: Quantity, cross_sectional_area_: Quantity) -> Quantity:

    result_heat_flow_expr = solve(law, heat_flow, dict=True)[0][heat_flow]
    result_expr = result_heat_flow_expr.subs({
        material_thermal_conductivity: material_thermal_conductivity_,
        inside_temperature: inside_temperature_,
        outside_temperature: outside_temperature_,
        cross_section_thickness: cross_section_thickness_,
        cross_sectional_area: cross_sectional_area_
    })
    return Quantity(result_expr)

