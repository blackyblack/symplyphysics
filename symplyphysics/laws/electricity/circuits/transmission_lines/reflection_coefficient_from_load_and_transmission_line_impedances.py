from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

# Description
## The reflection coefficient is equal to the ratio of the reflected wave power to the incident wave power.
## Knowing the load impedance and the characteristic impedance of the transmission line, it is possible to calculate the reflection coefficient.

## Law is: G = (Zl - Z0) / (Zl + Z0), where
## G - reflection coefficient,
## Zl - load impedance,
## Z0 - characteristic impedance of the transmission line.

reflection_coefficient = Symbol("reflection_coefficient", dimensionless)
load_impedance = Symbol("load_impedance", units.impedance)
characteristic_impedance = Symbol("characteristic_impedance", units.impedance)

law = Eq(reflection_coefficient,
    (load_impedance - characteristic_impedance) / (load_impedance + characteristic_impedance))


@validate_input(load_impedance_=load_impedance, characteristic_impedance_=characteristic_impedance)
@validate_output(reflection_coefficient)
def calculate_reflection_coefficient(load_impedance_: Quantity,
    characteristic_impedance_: Quantity) -> float:
    result_expr = solve(law, reflection_coefficient, dict=True)[0][reflection_coefficient]
    result_expr = result_expr.subs({
        load_impedance: load_impedance_,
        characteristic_impedance: characteristic_impedance_,
    })
    return convert_to_float(result_expr)
