"""
Average square speed in Maxwell—Boltzmann statistics
====================================================

For an ideal gas, the average square of speed is directly proportional to its temperature
and inversely proportional to the mass of the gas.

**Notation:**

#. :quantity_notation:`boltzmann_constant`.

**Conditions:**

#. The gas is in thermal equilibrium with the environment.
#. The gas particles are distributed according to Maxwell—Boltzmann statistics.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution#Typical_speeds>`__.

..
    TODO: move to `./maxwell_boltzmann_statistics.py`
"""

from sympy import (Eq, solve, S, stats, Interval)
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol,
    quantities)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.maxwell_boltzmann_statistics import speed_distribution
from symplyphysics.core.operations.symbolic import Average

average_square_speed = Average(symbols.speed**2)
"""
Average square of the :symbols:`speed` of gas molecules.
"""

equilibrium_temperature = clone_as_symbol(symbols.temperature, positive=True)
"""
Equilibrium :symbols:`temperature` of the gas.
"""

molecular_mass = clone_as_symbol(symbols.mass, positive=True)
"""
:symbols:`mass` of a gas molecule.
"""

law = Eq(
    average_square_speed,
    3 * quantities.boltzmann_constant * equilibrium_temperature / molecular_mass,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from Maxwell-Boltzmann distribution function

_distribution = speed_distribution.law.rhs.subs({
    speed_distribution.particle_mass: molecular_mass,
    speed_distribution.equilibrium_temperature: equilibrium_temperature,
})

_speed_random_variable = stats.ContinuousRV(
    speed_distribution.particle_speed,
    _distribution,
    set=Interval(0, S.Infinity),
)

_average_square_of_speed_derived = stats.E(_speed_random_variable**2)

assert expr_equals(_average_square_of_speed_derived, law.rhs)


@validate_input(temperature_in_gas_=equilibrium_temperature, mass_of_molecule_=molecular_mass)
@validate_output(average_square_speed)
def calculate_average_square_velocity(temperature_in_gas_: Quantity,
    mass_of_molecule_: Quantity) -> Quantity:
    result_average_square_velocity = solve(law, average_square_speed,
        dict=True)[0][average_square_speed]
    result_expr = result_average_square_velocity.subs({
        equilibrium_temperature: temperature_in_gas_,
        molecular_mass: mass_of_molecule_
    })
    return Quantity(result_expr)
