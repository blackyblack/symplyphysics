"""
Molar mass via molecular mass
=============================

The Avogadro constant is also the constant of proportionality between molar mass and
molecular mass.

**Notation:**

#. :quantity_notation:`avogadro_constant`

**Links:**

#. `Wikipedia, formula in the second paragraph <https://en.wikipedia.org/wiki/Avogadro_constant#>`__.
"""

from sympy import Eq, solve, Symbol as SymSymbol, Idx
from symplyphysics import (
    clone_as_symbol,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
    global_index,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.chemistry import (
    avogadro_constant_is_particle_count_over_amount_of_substance as avogadro_law,)
from symplyphysics.laws.conservation import (
    mixture_mass_equal_sum_of_components_masses as mass_sum_law,)
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,)

molar_mass = Symbol("molas_mass", units.mass / units.amount_of_substance)
"""
Mass per unit amount of substance.

Symbol:
    :code:`M`
"""

molecular_mass = clone_as_symbol(symbols.mass, display_symbol="m_0", display_latex="m_0")
"""
:symbols:`mass` of a single molecule.
"""

law = Eq(molar_mass, molecular_mass * quantities.avogadro_constant)
r"""
:code:`M = m_0 * N_A`

Latex:
    .. math::
        M = m_0 N_\text{A}
"""

# Derive from another definition of molar mass

_number_of_particles = SymSymbol("number_of_particles", integer=True)

_amount_of_substance = solve(avogadro_law.law,
    avogadro_law.amount_of_substance)[0].subs(avogadro_law.particle_count, _number_of_particles)

_local_index = Idx("local_index", (1, _number_of_particles))

_total_mass = mass_sum_law.law.rhs.subs(global_index,
    _local_index).subs(mass_sum_law.component_mass[_local_index], molecular_mass).doit()

_molar_mass_derived = solve(
    molar_qty_law.law,
    molar_qty_law.molar_quantity,
)[0].subs({
    molar_qty_law.extensive_quantity: _total_mass,
    molar_qty_law.amount_of_substance: _amount_of_substance,
})

assert expr_equals(_molar_mass_derived, law.rhs)


@validate_input(particle_mass_=molecular_mass)
@validate_output(molar_mass)
def calculate_molar_mass(particle_mass_: Quantity) -> Quantity:
    result = law.rhs.subs(molecular_mass, particle_mass_)
    return Quantity(result)
