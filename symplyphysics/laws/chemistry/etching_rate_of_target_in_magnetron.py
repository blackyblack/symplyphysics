"""
Etching rate of target in magnetron
===================================

The ions of the gas-discharge plasma in the magnetron fall on the target and knock the
atoms out of it. The etching rate is how many nanometers of the target
substance are etched per unit of time. In other words, this is how much thinner a target
becomes per unit of time.

**Notation:**

#. :quantity_notation:`elementary_charge`.
#. :quantity_notation:`avogadro_constant`.

..
    TODO: find link
    TODO: move to `magnetron` folder?
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    symbols,
)
from symplyphysics.quantities import elementary_charge, avogadro_constant

etching_rate = symbols.speed
"""
Target etching rate. See :symbols:`speed`.
"""

ion_current_density = symbols.current_density
"""
Ion flux :symbols:`current_density` incident on the target.
"""

target_molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass` of target atoms.
"""

sputtering_coefficient = Symbol("Y", dimensionless)
"""
Sputtering coefficient. Shows how many target atoms are knocked out of the target by a single ion.
"""

target_density = symbols.density
"""
Target :symbols:`density`.
"""

law = Eq(
    etching_rate,
    ion_current_density * target_molar_mass * sputtering_coefficient /
    (elementary_charge * target_density * avogadro_constant),
)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(ion_current_density_=ion_current_density,
    molar_mass_of_target_atom_=target_molar_mass,
    sputtering_coefficient_=sputtering_coefficient,
    target_density_=target_density)
@validate_output(etching_rate)
def calculate_etching_rate(ion_current_density_: Quantity, molar_mass_of_target_atom_: Quantity,
    sputtering_coefficient_: float, target_density_: Quantity) -> Quantity:
    result_expr = solve(law, etching_rate, dict=True)[0][etching_rate]
    result_expr = result_expr.subs({
        ion_current_density: ion_current_density_,
        target_molar_mass: molar_mass_of_target_atom_,
        sputtering_coefficient: sputtering_coefficient_,
        target_density: target_density_,
    })
    return Quantity(result_expr)
