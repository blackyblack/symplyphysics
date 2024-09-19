from sympy import Eq, solve, pi
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

# Description
## The effective cross section is a physical quantity characterizing the probability of transition of a system of
## two interacting particles to a certain final state, a quantitative characteristic of the acts of collision of
## particles of a stream hitting a target with target particles. The effective cross-section has the dimension of the area.

## Law is: g =  pi * d^2 * (1 + C / T), where
## g - the cross-sectional area of the interaction of particles,
## d - the diameter of the atom,
## C - Sutherland's constant,
## T - temperature.

cross_sectional_area_of_interaction = Symbol("cross_sectional_area_of_interaction", units.area)

diameter_of_atom = Symbol("diameter_of_atom", units.length)
sutherland_constant = Symbol("sutherland_constant", units.temperature)
temperature = symbols.temperature

law = Eq(cross_sectional_area_of_interaction,
    pi * diameter_of_atom**2 * (1 + sutherland_constant / temperature))


@validate_input(diameter_of_atom_=diameter_of_atom,
    sutherland_constant_=sutherland_constant,
    temperature_=temperature)
@validate_output(cross_sectional_area_of_interaction)
def calculate_cross_sectional_area_of_interaction(diameter_of_atom_: Quantity,
    sutherland_constant_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, cross_sectional_area_of_interaction,
        dict=True)[0][cross_sectional_area_of_interaction]
    result_expr = result_expr.subs({
        diameter_of_atom: diameter_of_atom_,
        sutherland_constant: sutherland_constant_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
