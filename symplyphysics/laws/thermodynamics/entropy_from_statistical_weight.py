from sympy import (Eq, solve, log)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless)

# Description
## The entropy of a system depends on the statistical weight of the state of the system.
## Statistical weight is the average number of microstates of a system that implement its macrostate.

## Law: S = k * ln(W), where
## S - entropy,
## k - boltzmann constant,
## W - statistical weight of state of the system.

entropy = Symbol("entropy", units.energy / units.temperature)
statistical_weight = Symbol("statistical_weight", dimensionless)

law = Eq(entropy, units.boltzmann_constant * log(statistical_weight))


def print_law() -> str:
    return print_expression(law)


@validate_input(statistical_weight_=statistical_weight)
@validate_output(entropy)
def calculate_entropy(statistical_weight_: float) -> Quantity:
    result_entropy_expr = solve(law, entropy, dict=True)[0][entropy]
    result_expr = result_entropy_expr.subs({statistical_weight: statistical_weight_})
    return Quantity(result_expr)
