from sympy import (Eq, solve)
from sympy.physics.units import faraday_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

# Description
## Faraday's second law of electrolysis. The equivalent mass of a substance in general in chemistry is its molar mass divided by an integer
## depending on the chemical reaction in which the substance participates; in this case, the equivalent
## is the molar mass of the substance formed during ion discharge divided by the sum of the ion charges
## (measured in elementary units), resulting in a molecule or atom of the substance.

## Law is: k = M / (F * n), where
## k - electrochemical equivalent,
## M - molar mass,
## F - faraday constant
## n - valence.

equivalent = Symbol("equivalent", units.mass / units.charge)

molar_mass = Symbol("molar_mass", units.mass / units.amount_of_substance)
valence = Symbol("valence", dimensionless)

law = Eq(equivalent, molar_mass / (faraday_constant * valence))


def print_law() -> str:
    return print_expression(law)


@validate_input(molar_mass_=molar_mass, valence_=valence)
@validate_output(equivalent)
def calculate_equivalent(molar_mass_: Quantity, valence_: int) -> Quantity:
    if not isinstance(valence_, int):
        raise ValueError("valence_ must be an integer.")
    if valence_ <= 0:
        raise ValueError("valence_ must be greater than 0.")

    result_expr = solve(law, equivalent, dict=True)[0][equivalent]
    result_expr = result_expr.subs({
        molar_mass: molar_mass_,
        valence: valence_,
    })
    return Quantity(result_expr)
