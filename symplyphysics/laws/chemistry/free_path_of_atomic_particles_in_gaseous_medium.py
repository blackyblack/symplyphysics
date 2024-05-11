from sympy import Eq, solve
from sympy.physics.units import boltzmann_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output,)

# Description
## The atoms of the target material evaporate and move towards the substrate inside the magnetron. At the same time,
## it collides with gas atoms. The free path length is the distance that a atomized atom travels between two collisions.

## Law is: l = k * T / (P * g), where
## l - free path length of atomized atom,
## k - boltzmann constant,
## T - temperature,
## P - pressure,
## g - the cross section of the interaction of a atomized atom and a gas atom.

free_path_length = Symbol("free_path_length", units.length)

pressure = Symbol("pressure", units.pressure)
temperature = Symbol("temperature", units.temperature)
cross_sectional_area_of_interaction = Symbol("cross_sectional_area_of_interaction", units.area)

law = Eq(free_path_length, boltzmann_constant * temperature / (pressure * cross_sectional_area_of_interaction))


@validate_input(pressure_=pressure,
    temperature_=temperature,
    cross_sectional_area_of_interaction_=cross_sectional_area_of_interaction)
@validate_output(free_path_length)
def calculate_free_path_length(pressure_: Quantity, temperature_: Quantity,
                                            cross_sectional_area_of_interaction_: Quantity) -> Quantity:
    result_expr = solve(law, free_path_length, dict=True)[0][free_path_length]
    result_expr = result_expr.subs({
        pressure: pressure_,
        temperature: temperature_,
        cross_sectional_area_of_interaction: cross_sectional_area_of_interaction_,
    })
    return Quantity(result_expr)
