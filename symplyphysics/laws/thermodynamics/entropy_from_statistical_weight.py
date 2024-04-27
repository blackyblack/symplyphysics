from sympy import Eq, solve, log, symbols, Function as SymFunction, dsolve
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, assert_equal)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count as avogadro_law
from symplyphysics.laws.thermodynamics import (
    change_in_entropy_of_ideal_gas_from_volume_and_temperature as entropy_change_law,
)
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,
)

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

# Derive the law from the properties of entropy

# TODO: describe the systems

_time = symbols("time", real=True)
_first_probability = symbols("first_probability", cls=SymFunction, positive=True)(_time)
_second_probability = symbols("second_probability", cls=SymFunction, positive=True)(_time)
_entropy = symbols("entropy", cls=SymFunction, real=True)

# TODO: describe additive property (make a law?)

_additive_property = Eq(
    _entropy(_first_probability * _second_probability),
    _entropy(_first_probability) + _entropy(_second_probability)
)

_constant_condition = (_first_probability * _second_probability).diff(_time)

_first_probability_diff_expr = solve(
    _constant_condition,
    _first_probability.diff(_time)
)[0]

_additive_property_diff = Eq(
    _additive_property.lhs.diff(_time),
    _additive_property.rhs.diff(_time)
).subs(
    _first_probability.diff(_time),
    _first_probability_diff_expr
)

_first_probability_expr = solve(
    _additive_property_diff,
    _first_probability
)[0]

_differential_property = Eq(
    _first_probability * _entropy(_first_probability).diff(_first_probability),
    _first_probability_expr * _entropy(_first_probability).diff(_first_probability),
)

assert expr_equals(_differential_property.lhs.diff(_second_probability), 0)
assert expr_equals(_differential_property.rhs.diff(_first_probability), 0)

# TODO: write what this means for the proof

_probability = symbols("probability", positive=True)
_equation_constant = symbols("equation_constant", real=True)

_entropy_eqn = Eq(
    _differential_property.lhs.subs(_first_probability, _probability),
    _equation_constant,
)

_integration_constant = symbols("integration_constant", real=True)

_entropy_via_probability = dsolve(
    _entropy_eqn,
    _entropy(_probability)
).rhs.subs(
    "C1",
    _integration_constant,
)

_additive_property_updated = _additive_property.replace(
    _entropy,
    lambda probability: _entropy_via_probability.subs(_probability, probability)
)

_integration_constant_expr = solve(
    _additive_property_updated,
    _integration_constant,
)[0]

assert expr_equals(_integration_constant_expr, 0)

_entropy_via_probability = _entropy_via_probability.subs(
    _integration_constant,
    _integration_constant_expr,
)

# TODO: explain the following expression of probability

_volume, _total_volume, _particle_count = symbols(
    "volume total_volume particle_count",
    positive=True,
)

_probability_via_volume = (_volume / _total_volume)**_particle_count

_entropy_via_volume = _entropy_via_probability.subs(_probability, _probability_via_volume)

_first_volume, _second_volume = symbols(
    "first_volume second_volume",
    positive=True,
)

_entropy_difference_derived = (
    _entropy_via_volume.subs(_volume, _second_volume)
    - _entropy_via_volume.subs(_volume, _first_volume)
).simplify()

_entropy_difference_from_law = entropy_change_law.law.rhs.subs({
    entropy_change_law.gas_mass: molar_qty_law.law.rhs.subs(
        molar_qty_law.molar_quantity,
        entropy_change_law.molar_mass,
    ),
    entropy_change_law.final_temperature: entropy_change_law.initial_temperature,
    entropy_change_law.start_volume: _first_volume,
    entropy_change_law.final_volume: _second_volume,
})

_equation_constant_expr = solve(
    (
        Eq(_entropy_difference_derived, _entropy_difference_from_law),
        avogadro_law.law.subs({
            avogadro_law.particles_count: _particle_count,
            avogadro_law.mole_count: molar_qty_law.amount_of_substance,
        })
    ),
    (
        _equation_constant,
        molar_qty_law.amount_of_substance,
    ),
    dict=True,
)[0][_equation_constant]

# This in, in fact, the Boltzmann constant
assert_equal(_equation_constant_expr, units.boltzmann_constant)

_entropy_via_probability = _entropy_via_probability.subs(
    _equation_constant,
    units.boltzmann_constant
)

# TODO: explain the following formula

_probability_via_statistical_weight = statistical_weight * _probability_via_volume

_entropy_via_statistical_weight = _entropy_via_probability.subs(
    _probability,
    _probability_via_statistical_weight,
).expand()

# TODO: explain that we can get rid of the free term

_entropy_via_statistical_weight = (
    _entropy_via_statistical_weight
    - _entropy_via_statistical_weight.subs(statistical_weight, 1)
)

assert expr_equals(_entropy_via_statistical_weight, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(statistical_weight_=statistical_weight)
@validate_output(entropy)
def calculate_entropy(statistical_weight_: float) -> Quantity:
    result_entropy_expr = solve(law, entropy, dict=True)[0][entropy]
    result_expr = result_entropy_expr.subs({statistical_weight: statistical_weight_})
    return Quantity(result_expr)
