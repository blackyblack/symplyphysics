from sympy import Eq, solve, sqrt, pi, log, sinh
from symplyphysics import (
    units,
    Quantity,
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

## Law is: ef = 1 + ((er - 1) / 2) * (ln(2 * (1 + sqrt(sqrt(1 - k))) / (1 - sqrt(sqrt(1 - k)))) / pi) * (pi / ln(2 * (1 + sqrt(sqrt(1 - k1))) / (1 - sqrt(sqrt(1 - k1))))), where
## k = a / b,
## k1 = sinh(pi * a / (2 * h)) / sinh(pi * b / (2 * h)),
## ef - effective permittivity of the coplanar line,
## er - relative permittivity of the dielectric substrate of the coplanar line,
## 2 * a - width of the central electrode of the coplanar line,
## h - thickness of substrate,
## 2 * b - distance between the extreme electrodes.

# Conditions:
# - h < b / 2;
# - 0 < k^2 <= 0.5;
# - 0 < k1^2 <= 0.5.


effective_permittivity = Symbol("effective_permittivity", dimensionless)

relative_permittivity = Symbol("relative_permittivity", dimensionless)
distance_between_electrodes = Symbol("distance_between_electrodes", units.length)
thickness_of_substrate = Symbol("thickness_of_substrate", units.length)
central_electrode_width = Symbol("central_electrode_width", units.length)

expression_1 = (relative_permittivity - 1) / 2
expression_2 = sinh(pi * central_electrode_width / (2 * thickness_of_substrate)) / sinh(pi * distance_between_electrodes / (2 * thickness_of_substrate))
expression_3 = sqrt(1 - expression_2**2)
expression_4 =  pi / log(2 * (1 + sqrt(expression_3)) / (1 - sqrt(expression_3)))
expression_5 = sqrt(1 - (central_electrode_width / distance_between_electrodes)**2)
expression_6 =  log(2 * (1 + sqrt(expression_5)) / (1 - sqrt(expression_5))) / pi


law = Eq(effective_permittivity, 1 + expression_1 * expression_6 * expression_4)


@validate_input(relative_permittivity_=relative_permittivity,
    distance_between_electrodes_=distance_between_electrodes,
    thickness_of_substrate_=thickness_of_substrate,
    central_electrode_width_=central_electrode_width)
@validate_output(effective_permittivity)
def calculate_effective_permittivity(relative_permittivity_: float,
    distance_between_electrodes_: Quantity, thickness_of_substrate_: Quantity,
    central_electrode_width_: Quantity) -> float:
    if thickness_of_substrate_.scale_factor >= distance_between_electrodes_.scale_factor / 2:
        raise ValueError("The thickness of substrate must be less than the distance between electrodes divided in half")
    if (central_electrode_width_.scale_factor / distance_between_electrodes_.scale_factor)**2 > 0.5:
        raise ValueError("k^2 must be less than or equal to the 0.5")
    if (sinh(pi * central_electrode_width_.scale_factor / 2 / thickness_of_substrate_.scale_factor)
        / sinh(pi* distance_between_electrodes_.scale_factor / 2 / thickness_of_substrate_.scale_factor))**2 > 0.5:
        raise ValueError("k1^2 must be less than or equal to the 0.5")
    result_expr = solve(law, effective_permittivity, dict=True)[0][effective_permittivity]
    result_expr = result_expr.subs({
        relative_permittivity: relative_permittivity_,
        distance_between_electrodes: distance_between_electrodes_,
        thickness_of_substrate: thickness_of_substrate_,
        central_electrode_width: central_electrode_width_
    })
    return convert_to_float(result_expr)
