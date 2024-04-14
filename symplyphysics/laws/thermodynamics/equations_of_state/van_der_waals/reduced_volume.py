from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to_float,
)

# Description
## Reduced units are used in the dimensionless van der Waals equation of state.

# Law: V* = V / V_c
## V* - reduced volume
## V - volume
## V_c - [critical volume](./critical_volume.py)

# Note: one can also use specific or molar volumes in the right-hand side

reduced_volume = Symbol("reduced_volume", dimensionless)
volume = Symbol("volume", units.volume)
critical_volume = Symbol("critical_volume", units.volume)

law = Eq(reduced_volume, volume / critical_volume)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    volume_=volume,
    critical_volume_=critical_volume,
)
@validate_output(reduced_volume)
def calculate_reduced_volume(
    volume_: Quantity,
    critical_volume_: Quantity,
) -> float:
    result = law.rhs.subs({
        volume: volume_,
        critical_volume: critical_volume_,
    })
    return convert_to_float(result)
