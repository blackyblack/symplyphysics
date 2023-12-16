from sympy import (I, Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, angle_type)

# Description
## The impedance of ideal coil depends on its inductivity and frequency. While having zero resistivity, the real part of
## coil impedance is zero.
## Law: Zl = jwL, where
## Zl is coil impedance,
## j is imaginary number,
## w is circular frequency,
## L is coil inductivity.

coil_impedance = Symbol("coil_impedance", units.impedance)
circular_frequency = Symbol("circular_frequency", angle_type / units.time)
coil_inductivity = Symbol("coil_inductivity", units.inductance)

law = Eq(coil_impedance, I * circular_frequency * coil_inductivity)


def print_law() -> str:
    return print_expression(law)


@validate_input(inductivity_=coil_inductivity, circular_frequency_=circular_frequency)
@validate_output(coil_impedance)
def calculate_impedance(inductivity_: Quantity, circular_frequency_: Quantity) -> Quantity:
    result_impedance_expr = solve(law, coil_impedance, dict=True)[0][coil_impedance]
    result_expr = result_impedance_expr.subs({
        coil_inductivity: inductivity_,
        circular_frequency: circular_frequency_
    })
    return Quantity(result_expr)
