r"""
Entropy from statistical weight
===============================

Entropy of a system depends on the statistical weight of the system's state. Statistical weight
is the average number of microstates of a system that implement its macrostate.

**Notation:**

#. :math:`k_\text{B}` (:code:`k_B`) is the Boltzmann constant.
"""

from sympy import Eq, solve, log, symbols, Function as SymFunction, dsolve
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, dimensionless,
    assert_equal)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count as avogadro_law
from symplyphysics.laws.thermodynamics import (
    change_in_entropy_of_ideal_gas_from_volume_and_temperature as entropy_change_law,)
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,)

entropy = Symbol("entropy", units.energy / units.temperature)
"""
Entropy of the system.

Symbol:
    :code:`S`
"""

statistical_weight = Symbol("statistical_weight", dimensionless)
"""
Statistical weight of the system's state.

Symbol:
    :code:`W`
"""

law = Eq(entropy, units.boltzmann_constant * log(statistical_weight))
r"""
:code:`S = k_B * log(W)`

Latex:
    .. math::
        S = k_text{B} \log W
"""

# Derive the law from the properties of entropy

# Let us find entropy as a function of probability of the system's state, i.e. `S = f(P)`,
# where `S` is entropy, `P` is probability, and `f` is the functional dependency to be found.
# This dependency must be universal for all physical systems and should only depend on the system's
# probability.

# Suppose there are two independent subsystems in states with probabilities `P_1` and `P_2`
# respectively. If we join the two subsystems into one, the probability of the total system
# will be `P_12`, and using their independence we have `P_12 = P_1 * P_2`.

_time = symbols("time", real=True)
_first_probability = symbols("first_probability", cls=SymFunction, positive=True)(_time)
_second_probability = symbols("second_probability", cls=SymFunction, positive=True)(_time)
_entropy = symbols("entropy", cls=SymFunction, real=True)

# FIXME: use law for this step
# According to thermodynamics, the entropy of a complex system should equal the sum of entropies
# of its independent subsystems. This is the additive property of entropy:

_additive_property = Eq(_entropy(_first_probability * _second_probability),
    _entropy(_first_probability) + _entropy(_second_probability))

# To solve this functional equation, we can use the following method of variating
# the two probabilities.

_constant_condition = Eq((_first_probability * _second_probability).diff(_time), 0)

_first_probability_diff_expr = solve(_constant_condition, _first_probability.diff(_time))[0]

_additive_property_diff = Eq(_additive_property.lhs.diff(_time),
    _additive_property.rhs.diff(_time)).subs(_first_probability.diff(_time),
    _first_probability_diff_expr)

_first_probability_expr = solve(_additive_property_diff, _first_probability)[0]

_differential_property = Eq(
    _first_probability * _entropy(_first_probability).diff(_first_probability),
    _first_probability_expr * _entropy(_first_probability).diff(_first_probability),
)

# Now the left-hand side only depends on `P_1` and the right-hand side only depends on `P_2`

assert expr_equals(_differential_property.lhs.diff(_second_probability), 0)
assert expr_equals(_differential_property.rhs.diff(_first_probability), 0)

# 1. Since both `P_1` and `P_2` can take arbitrary values in the range of their definition,
# the function described by `_differential_property` must be constant for all values of its argument.

# 2. Since `f(P)` is universal for all bodies, the constant is universal for all bodies as well.

_probability = symbols("probability", positive=True)
_equation_constant = symbols("equation_constant", real=True)

_entropy_eqn = Eq(
    _differential_property.lhs.subs(_first_probability, _probability),
    _equation_constant,
)

_integration_constant = symbols("integration_constant", real=True)

# Now we can integrate the above differential property and show that the integration constant is
# actually zero.

_entropy_via_probability = dsolve(_entropy_eqn, _entropy(_probability)).rhs.subs(
    "C1",
    _integration_constant,
)

_additive_property_updated = _additive_property.replace(
    _entropy, lambda probability: _entropy_via_probability.subs(_probability, probability))

_integration_constant_expr = solve(
    _additive_property_updated,
    _integration_constant,
)[0]

assert expr_equals(_integration_constant_expr, 0)

_entropy_via_probability = _entropy_via_probability.subs(
    _integration_constant,
    _integration_constant_expr,
)

# Now, the question is to find the numeric value of the equation constant.
# To do that, we can compare the entropy difference in two arbitrary states and the
# logarithm of the probability ratio of these two states.

# Let `V` be some volume containing `N` ideal gas particles and `V_0` be a part of it.
# Then the probability of finding one molecule in `V_0` is `V_0 / V`, and the probability
# of finding all `N` particles is `(V_0 / V)**N`

_volume, _total_volume, _particle_count = symbols(
    "volume total_volume particle_count",
    positive=True,
)

# FIXME: use law for this step
_probability_via_volume = (_volume / _total_volume)**_particle_count

_entropy_via_volume = _entropy_via_probability.subs(_probability, _probability_via_volume)

# Let us assume that during a process the volume of the gas changes from `V_1` to `V_2` at
# constant temperature and particle count. Then we can find the entropy of the two states
# via the above formula.

_first_volume, _second_volume = symbols(
    "first_volume second_volume",
    positive=True,
)

_entropy_difference_derived = (_entropy_via_volume.subs(_volume, _second_volume) -
    _entropy_via_volume.subs(_volume, _first_volume)).simplify()

# We also utilise another law derived from completely different thermodynamic principles:

_entropy_difference_from_law = entropy_change_law.law.rhs.subs({
    entropy_change_law.mass:
    molar_qty_law.law.rhs.subs(
    molar_qty_law.molar_quantity,
    entropy_change_law.molar_mass,
    ),
    entropy_change_law.final_temperature:
    entropy_change_law.initial_temperature,
    entropy_change_law.initial_volume:
        _first_volume,
    entropy_change_law.final_volume:
        _second_volume,
})

_equation_constant_expr = solve(
    (Eq(_entropy_difference_derived, _entropy_difference_from_law),
    avogadro_law.law.subs({
    avogadro_law.particles_count: _particle_count,
    avogadro_law.mole_count: molar_qty_law.amount_of_substance,
    })),
    (
    _equation_constant,
    molar_qty_law.amount_of_substance,
    ),
    dict=True,
)[0][_equation_constant]

# Thus, we obtain that this is in, in fact, the Boltzmann constant
assert_equal(_equation_constant_expr, units.boltzmann_constant)

_entropy_via_probability = _entropy_via_probability.subs(_equation_constant,
    units.boltzmann_constant)

# The derivation of the following formula is too long and it can be found in "General Course of Physics"
# by Sivukhin D.V, Chapter 80, formula (80.8):
# FIXME: use law in this step
_probability_via_statistical_weight = statistical_weight * _probability_via_volume

_entropy_via_statistical_weight = _entropy_via_probability.subs(
    _probability,
    _probability_via_statistical_weight,
).expand()

# Entropy is defined up to a constant term, therefore we can get rid of it in the formula

_entropy_via_statistical_weight = (_entropy_via_statistical_weight -
    _entropy_via_statistical_weight.subs(statistical_weight, 1))

assert expr_equals(_entropy_via_statistical_weight, law.rhs)


@validate_input(statistical_weight_=statistical_weight)
@validate_output(entropy)
def calculate_entropy(statistical_weight_: float) -> Quantity:
    result_entropy_expr = solve(law, entropy, dict=True)[0][entropy]
    result_expr = result_entropy_expr.subs({statistical_weight: statistical_weight_})
    return Quantity(result_expr)
