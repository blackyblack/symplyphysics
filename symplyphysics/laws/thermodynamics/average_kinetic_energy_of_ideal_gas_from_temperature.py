r"""
Average kinetic energy of ideal gas from temperature
====================================================

The kinetic theory of ideal gases allows us to determine the *average kinetic energy* for an ideal gas.
It states that the average kinetic energy for all ideal gases is directly proportional to the absolute
temperature of the gas and only depends on the temperature.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is ideal.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Kinetic_theory_of_gases#Temperature_and_kinetic_energy>`__.
"""

from sympy import Eq, solve, Rational, stats, Interval, S
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import energy_distribution

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy, positive=True)
r"""
Average kinetic energy of the ideal gas.

Symbol:
    :code:`avg(K)`

Latex:
    :math:`\langle K \rangle`
"""

equilibrium_temperature = clone_as_symbol(
    symbols.temperature,
    positive=True,
)
"""
Equilibrium :symbols:`temperature` of the gas.
"""

law = Eq(average_kinetic_energy,
    Rational(3, 2) * quantities.boltzmann_constant * equilibrium_temperature)
r"""
:code:`avg(K) = 3/2 * k_B * T`

Latex:
    :math:`\langle K \rangle = \frac{3}{2} k_\text{B} T`
"""

# Derive from Maxwell-Boltzmann energy distribution

_distribution = energy_distribution.law.rhs.subs(
    energy_distribution.equilibrium_temperature,
    equilibrium_temperature,
)

_random_energy_variable = stats.ContinuousRV(
    energy_distribution.energy,
    _distribution,
    set=Interval(0, S.Infinity),
)

_average_energy_derived = stats.E(_random_energy_variable)

assert expr_equals(_average_energy_derived, law.rhs)


@validate_input(temperature_=symbols.temperature)
@validate_output(average_kinetic_energy)
def calculate_average_kinetic_energy(temperature_: Quantity) -> Quantity:
    result_expr = solve(law, average_kinetic_energy, dict=True)[0][average_kinetic_energy]
    result_average_kinetic_energy = result_expr.subs(equilibrium_temperature, temperature_)
    return Quantity(result_average_kinetic_energy)
