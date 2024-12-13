from sympy import Eq, solve, sqrt
from sympy.physics.units import boltzmann_constant
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, symbols)

# Description
## The atoms of the target material evaporate and move towards the substrate inside the magnetron. At the same time,
## it collides with gas atoms. The free path length is the distance that a traveling atom travels between two collisions.

## Law is: l = k * T / (sqrt(2) * P * sigma), where
## l - free path length of traveling atom,
## k - boltzmann constant,
## T - temperature,
## P - pressure,
## sigma - the cross section of the interaction of a traveling atom and a gas atom.

# Notes
## Assuming the model of spherical gas molecules, `sigma = pi * d**2` where `d` is the molecule diameter.

# Links
## Wikipedia, the fourth formula <https://en.wikipedia.org/wiki/Mean_free_path#Kinetic_theory_of_gases>
## Chemistry LibreTexts, "27.6.4. Mean Free Path" <https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Physical_Chemistry_(LibreTexts)/27%3A_The_Kinetic_Theory_of_Gases/27.06%3A_Mean_Free_Path>

# NOTE: remove the mention of a magnetron from the description?

free_path_length = Symbol("free_path_length", units.length)

pressure = Symbol("pressure", units.pressure)
temperature = symbols.temperature
cross_sectional_area_of_interaction = Symbol("cross_sectional_area_of_interaction", units.area)

law = Eq(free_path_length,
    boltzmann_constant * temperature / (sqrt(2) * pressure * cross_sectional_area_of_interaction))


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
