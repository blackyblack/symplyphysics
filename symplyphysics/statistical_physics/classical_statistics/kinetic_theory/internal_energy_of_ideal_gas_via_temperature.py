"""
Internal energy of ideal gas via temperature
============================================

Internal energy of an ideal gas is the sum of the kinetic energy of all of its molecules.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Conditions:**

#. The gas is ideal.
#. The gas is monatomic.

**Links:**

#. `Physics LibreTexts <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/University_Physics_II_-_Thermodynamics_Electricity_and_Magnetism_(OpenStax)/03%3A_The_First_Law_of_Thermodynamics/3.03%3A_Work_Heat_and_Internal_Energy>`__.

#. `Wikipedia, see text <https://en.wikipedia.org/wiki/Ideal_gas#Internal_energy>`__.

..
    TODO replace `mass/molar_mass` with `amount_of_substance`
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.statistical_physics.classical_statistics.kinetic_theory import average_kinetic_energy_of_ideal_gas_from_temperature as average_energy_law
from symplyphysics.chemistry.molecular_properties import avogadro_constant_is_particle_count_over_amount_of_substance as avogadro_law
from symplyphysics.quantity_relations import quantity_is_molar_quantity_times_amount_of_substance as molar_quantity_law

internal_energy = symbols.internal_energy
"""
:symbols:`internal_energy` of the gas.
"""

mass = symbols.mass
"""
:symbols:`mass` of the gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass` of the gas.
"""

law = Eq(internal_energy, 3 * mass * quantities.molar_gas_constant * temperature / (2 * molar_mass))
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the average kinetic energy of ideal gas molecules.

## Express the amount of substance of the gas via its mass and molar mass.
_amount_of_substance = solve(
    molar_quantity_law.law.subs({
    molar_quantity_law.extensive_quantity: mass,
    molar_quantity_law.molar_quantity: molar_mass,
    }),
    molar_quantity_law.amount_of_substance,
)[0]

## Express the number of gas molecules via the amount of substance.
_particle_count = solve(avogadro_law.law, avogadro_law.particle_count)[0].subs(
    avogadro_law.amount_of_substance, _amount_of_substance)

## The average kinetic energy of a molecule at the equilibrium temperature of the gas.
_average_energy = average_energy_law.law.rhs.subs(average_energy_law.equilibrium_temperature,
    temperature)

## The gas is ideal, so the molecules do not interact and there is no potential energy, and
## monatomic, so all of the kinetic energy is translational. Hence the internal energy of
## the gas is the average kinetic energy of a molecule times the number of molecules.
_internal_energy_derived = _particle_count * _average_energy

## The molar gas constant is the Boltzmann constant times the Avogadro constant.
_internal_energy_derived = _internal_energy_derived.subs(
    quantities.boltzmann_constant * quantities.avogadro_constant,
    quantities.molar_gas_constant,
)

assert expr_equals(_internal_energy_derived, law.rhs)


@validate_input(mass_of_gas_=mass, temperature_=temperature, mole_mass_=molar_mass)
@validate_output(internal_energy)
def calculate_inner_energy(mass_of_gas_: Quantity, temperature_: Quantity,
    mole_mass_: Quantity) -> Quantity:
    solved = solve(law, internal_energy, dict=True)[0][internal_energy]
    result_expr = solved.subs({
        mass: mass_of_gas_,
        temperature: temperature_,
        molar_mass: mole_mass_
    })
    return Quantity(result_expr)
