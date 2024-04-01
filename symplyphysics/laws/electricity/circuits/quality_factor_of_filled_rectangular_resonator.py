from sympy import (
    Eq,
    solve,
)
from symplyphysics import (Symbol, print_expression, validate_input, validate_output,
    dimensionless)

# Description
## A rectangular resonator consists of metal walls and a material filling it.
## In the case when the resonator is empty, its quality factor depends only on the losses in the metal walls of the resonator.
## The quality factor shows the ratio of the energy stored in the resonator to the energy loss during one oscillation period.
## In the case of a filled resonator, the Q factor also depends on the losses in the dielectric.

## Law is: Q = 1 / ((1 / Qm) + tan(d)), where
## Q - quality factor of the filled resonator,
## Qm - quality factor of the empty resonator,
## tan(d) - tangent of the dielectric loss angle of the material filling the resonator.

quality_factor = Symbol("quality_factor", dimensionless)

empty_resonator_quality_factor = Symbol("empty_resonator_quality_factor", dimensionless)
tangent_dielectric_loss_angle = Symbol("tangent_dielectric_loss_angle", dimensionless)

law = Eq(quality_factor, 1 / ((1 / empty_resonator_quality_factor) + tangent_dielectric_loss_angle))


def print_law() -> str:
    return print_expression(law)


@validate_input(empty_resonator_quality_factor_=empty_resonator_quality_factor,
    tangent_dielectric_loss_angle_=tangent_dielectric_loss_angle)
@validate_output(quality_factor)
def calculate_quality_factor(empty_resonator_quality_factor_: float,
    tangent_dielectric_loss_angle_: float) -> float:
    result_expr = solve(law, quality_factor, dict=True)[0][quality_factor]
    result_expr = result_expr.subs({
        empty_resonator_quality_factor: empty_resonator_quality_factor_,
        tangent_dielectric_loss_angle: tangent_dielectric_loss_angle_,
    })
    return float(result_expr)
