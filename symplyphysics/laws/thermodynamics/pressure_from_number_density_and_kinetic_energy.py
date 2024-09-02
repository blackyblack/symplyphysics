"""
Pressure from number density and kinetic energy
===============================================

The ideal gas equation can be written in an equivalent way using the number density
of the gas and the average kinetic energy of gas particles.

**Condition:**

#. The gas is ideal.
"""

from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation as ideal_gas_law
from symplyphysics.laws.thermodynamics import average_kinetic_energy_of_ideal_gas_from_temperature as kinetic_energy
from symplyphysics.laws.chemistry import avogadro_number_from_mole_count as avogadro_number
from symplyphysics.definitions import number_density_is_number_of_objects_per_unit_volume as number_density_def

average_kinetic_energy = Symbol("average_kinetic_energy", units.energy)
r"""
Average kinetic energy of gas particles.

Symbol:
    :code:`avg(K)`

Latex:
    :math:`\langle K \rangle`
"""

number_density = Symbol("number_density", 1 / units.volume)
"""
Number density of the gas.

Symbol:
    :code:`n`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure inside the gas.

Symbol:
    :code:`p`
"""

law = Eq(pressure, 2 / 3 * number_density * average_kinetic_energy)
r"""
:code:`p = 2/3 * n * avg(K)`

Latex:
    .. math::
        p = \frac{2}{3} n \langle K \rangle
"""

## Proof

_temperature_eq = kinetic_energy.law.subs({
    kinetic_energy.average_kinetic_energy: average_kinetic_energy,
})

_derived_temperature = solve(_temperature_eq, kinetic_energy.equilibrium_temperature)[0]

_particles_number_eq = number_density_def.definition.subs({
    number_density_def.volume: ideal_gas_law.volume,
    number_density_def.number_density: number_density
})

_derived_particles_number = solve(_particles_number_eq, number_density_def.number_of_objects)[0]

_mole_count_eq = avogadro_number.law.subs(
    {avogadro_number.particles_count: _derived_particles_number})

_derived_mole_count = solve(_mole_count_eq, avogadro_number.mole_count)[0]

_ideal_gas_law_subs = ideal_gas_law.law.subs({
    ideal_gas_law.temperature: _derived_temperature,
    units.molar_gas_constant: units.boltzmann * units.avogadro,
    ideal_gas_law.amount_of_substance: _derived_mole_count
})

_derived_pressure = solve(_ideal_gas_law_subs, ideal_gas_law.pressure)[0]

assert expr_equals(_derived_pressure, law.rhs)


@validate_input(molecules_concentration_=number_density,
    average_kinetic_energy_=average_kinetic_energy)
@validate_output(pressure)
def calculate_pressure(molecules_concentration_: Quantity,
    average_kinetic_energy_: Quantity) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_pressure = result_expr.subs({
        number_density: molecules_concentration_,
        average_kinetic_energy: average_kinetic_energy_,
    })
    return Quantity(result_pressure)
