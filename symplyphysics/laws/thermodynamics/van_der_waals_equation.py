from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The Van der Waals equation more accurately describes the behavior of real gases at low temperatures and takes into account the main characteristics of a gas with intermolecular interaction, in comparison with the equation of an ideal gas.


## Law: (p + nu^2 * a / V^2) * (V - nu * b) = nu * R * T
## Where:
## p is pressure
## nu is amount of substance
## a is a constant (different for different substances) characterizing the mutual attraction of molecules
## V is volume
## b is a constant (different for different substances) related to the size of molecules, characterizing the mutual repulsion of molecules.
## R is universal gas constant
## T is temperature


pressure = Symbol("pressure", units.pressure)
amount_of_substance = Symbol("amount_of_substance", units.amount_of_substance)
mutual_attraction_constant = Symbol("mutual_attraction_constant", units.pressure * units.volume**2 / units.amount_of_substance**2)
volume = Symbol("volume", units.volume)
mutual_repulsion_constant = Symbol("mutual_repulsion_constant", units.volume / units.amount_of_substance)
temperature = Symbol("temperature", units.temperature)


law = Eq((pressure + (amount_of_substance**2 * mutual_attraction_constant / volume**2)) * (volume - (amount_of_substance * mutual_repulsion_constant)), amount_of_substance * units.molar_gas_constant * temperature)

def print_law() -> str:
    return print_expression(law)


@validate_input(amount_of_substance_=amount_of_substance, mutual_attraction_constant_=mutual_attraction_constant, volume_=volume, mutual_repulsion_constant_=mutual_repulsion_constant, pressure_=pressure)
@validate_output(temperature)
def calculate_temperature(amount_of_substance_: Quantity, mutual_attraction_constant_: Quantity,
    volume_: Quantity, mutual_repulsion_constant_: Quantity, pressure_: Quantity) -> Quantity:
    result_expr = solve(law, temperature, dict=True)[0][temperature]
    result_temperature = result_expr.subs({
        amount_of_substance: amount_of_substance_,
        mutual_attraction_constant: mutual_attraction_constant_,
        volume: volume_,
        mutual_repulsion_constant: mutual_repulsion_constant_,
        pressure: pressure_
    })
    return Quantity(result_temperature)
