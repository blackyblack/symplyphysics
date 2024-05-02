from sympy import Eq, solve
from symplyphysics import (
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float
)

## Description
## The coplanar transmission line is a dielectric substrate on the surface of which 3 electrodes are located.
## When a wave propagates along a coplanar line, part of the field goes out, since the coplanar line does
## not have metal borders on all sides, unlike, for example, rectangular waveguides. Then imagine an environment
## in which the field will have the same magnitude as the field of a coplanar line. The permittivity of such a
## medium will be called the effective permittivity of the line.

## Law is: ef = 1 + er / 2, where
## ef - effective permittivity of the coplanar line,
## er - relative permittivity of the dielectric substrate of the coplanar line.

# Conditions:
# h - thickness of substrate,
# 2 * b - distance between the extreme electrodes.
# - h >= b / 2.


effective_permittivity = Symbol("effective_permittivity", dimensionless)

relative_permittivity = Symbol("relative_permittivity", dimensionless)


law = Eq(effective_permittivity, (1 + relative_permittivity) / 2)


@validate_input(relative_permittivity_=relative_permittivity)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float) -> float:
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
    })
    return convert_to_float(result_expr)
