from sympy import Eq, solve
from sympy.physics.units import elementary_charge, avogadro_constant
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output,
                           dimensionless, clone_symbol, symbols)

# Description
## Electrolysis is a physico-chemical process consisting in the release of components of dissolved substances or other substances on the electrodes,
## which are the result of secondary reactions on the electrodes, which occurs when an electric current passes through a solution or melt of an electrolyte.
## Current output is the proportion of electric current spent on the passage of the target electrochemical reaction.

## Law is: m =  J * M * B * t / (Z * q * Na), where
## m - mass of film,
## J - the current in the electrolyte,
## M - molar mass of the metall atom,
## B - current output,
## t - time,
## Z - valence of the metal,
## q - elementary charge,
## Na - avogadro constant.

mass_of_film = clone_symbol(symbols.basic.mass, "mass_of_film")

current = Symbol("current", units.current)
molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
current_output = Symbol("current_output", dimensionless)
valence = Symbol("valence", dimensionless)
time = Symbol("time", units.time)

law = Eq(mass_of_film, current * molar_mass * current_output * time / (valence * elementary_charge * avogadro_constant))


@validate_input(current_=current,
                molar_mass_=molar_mass,
                current_output_=current_output,
                valence_=valence,
                time_=time)
@validate_output(mass_of_film)
def calculate_mass_of_film(current_: Quantity, molar_mass_: Quantity,
                                                  current_output_: float, valence_: int, time_: Quantity) -> Quantity:
    result_expr = solve(law, mass_of_film, dict=True)[0][mass_of_film]
    result_expr = result_expr.subs({
        current: current_,
        molar_mass: molar_mass_,
        current_output: current_output_,
        valence: valence_,
        time: time_,
    })
    return Quantity(result_expr)
