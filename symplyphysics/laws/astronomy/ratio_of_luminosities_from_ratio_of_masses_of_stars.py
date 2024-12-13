from sympy import (
    Eq,
    solve,
)
from symplyphysics import (
    clone_as_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
)

# Description
## Comparisons of masses and luminosities for most stars revealed the following relationship:
## luminosity is approximately proportional to the fourth power of mass.

## Law is: L2 / L1 = (M2 / M1)^4, where
## M2 - mass of second object,
## M1 - mass of first object,
## L2 - illuminance of second object,
## L1 - illuminance of first object.

# Conditions
## The masses must obey the inequality `0.43 * M1 < M2 < 2 * M1`, see link for more information.

# Links: Wikipedia <https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation>

mass_first = clone_as_symbol(symbols.mass, display_symbol="m_1", display_latex="m_1")
mass_second = clone_as_symbol(symbols.mass, display_symbol="m_2", display_latex="m_2")
illuminance_first = Symbol("illuminance_first", units.energy / units.area)
illuminance_second = Symbol("illuminance_second", units.energy / units.area)

law = Eq(illuminance_second / illuminance_first, (mass_second / mass_first)**4)


@validate_input(mass_first_=mass_first,
    mass_second_=mass_second,
    illuminance_first_=illuminance_first)
@validate_output(illuminance_second)
def calculate_illuminance_second(mass_first_: Quantity, mass_second_: Quantity,
    illuminance_first_: Quantity) -> Quantity:
    result_expr = solve(law, illuminance_second, dict=True)[0][illuminance_second]
    result_expr = result_expr.subs({
        mass_first: mass_first_,
        mass_second: mass_second_,
        illuminance_first: illuminance_first_
    })
    return Quantity(result_expr)
