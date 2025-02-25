"""
Current density from mobility
=============================

The current density can be calculated from the concentration and mobilities of holes and
electrons and the electric field.

**Notation:**

#. :quantity_notation:`elementary_charge`.

**Links:**

#. `Wikipedia, last formula <https://en.wikipedia.org/wiki/Electron_mobility#Relation_to_current_density>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.condensed_matter import current_density_via_number_density_and_drift_velocity as density_velocity_law
from symplyphysics.laws.condensed_matter import drift_velocity_of_charge_carriers as velocity_law

current_density = symbols.current_density
"""
:symbols:`current_density` of charge carriers.
"""

electrons_concentration = clone_as_symbol(symbols.number_density, display_symbol="n_e", display_latex="n_\\text{e}")
"""
Concentration (:symbols:`number_density`) of electrons.
"""

holes_concentration = clone_as_symbol(symbols.number_density, display_symbol="n_h", display_latex="n_\\text{h}")
"""
Concentration (:symbols:`number_density`) of holes.
"""

electrons_mobility = clone_as_symbol(symbols.mobility, display_symbol="mu_e", display_latex="\\mu_\\text{e}")
"""
:symbols:`mobility` of electrons.
"""

holes_mobility = clone_as_symbol(symbols.mobility, display_symbol="mu_h", display_latex="\\mu_\\text{h}")
"""
:symbols:`mobility` of holes.
"""

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength`.
"""

law = Eq(
    current_density,
    quantities.elementary_charge *
    (-electrons_concentration * electrons_mobility + holes_concentration * holes_mobility) *
    electric_field_strength)
"""
:laws:symbol::

:laws:latex::
"""

## This law might be derived via law for current density in metals.
_velocity_law_electrons = velocity_law.law.subs({
    velocity_law.mobility: electrons_mobility,
    velocity_law.electric_field_strength: electric_field_strength,
})
_velocity_electrons = solve(_velocity_law_electrons, velocity_law.drift_velocity,
    dict=True)[0][velocity_law.drift_velocity]
_density_velocity_law_electrons = density_velocity_law.law.subs({
    density_velocity_law.charge: -quantities.elementary_charge,
    density_velocity_law.number_density: electrons_concentration,
    density_velocity_law.drift_velocity: _velocity_electrons,
})

_velocity_law_holes = velocity_law.law.subs({
    velocity_law.mobility: holes_mobility,
    velocity_law.electric_field_strength: electric_field_strength,
})
_velocity_holes = solve(_velocity_law_holes, velocity_law.drift_velocity,
    dict=True)[0][velocity_law.drift_velocity]
_density_velocity_law_holes = density_velocity_law.law.subs({
    density_velocity_law.charge: quantities.elementary_charge,
    density_velocity_law.number_density: holes_concentration,
    density_velocity_law.drift_velocity: _velocity_holes,
})
_density_current_electrons_derived = solve(_density_velocity_law_electrons,
    density_velocity_law.current_density,
    dict=True)[0][density_velocity_law.current_density]
_density_current_holes_derived = solve(_density_velocity_law_holes,
    density_velocity_law.current_density,
    dict=True)[0][density_velocity_law.current_density]
_density_current_derived = _density_current_electrons_derived + _density_current_holes_derived

# Check if derived density current is same as declared.
assert expr_equals(_density_current_derived, law.rhs)


@validate_input(electrons_concentration_=electrons_concentration,
    holes_concentration_=holes_concentration,
    electrons_mobility_=electrons_mobility,
    holes_mobility_=holes_mobility,
    electric_intensity_=electric_field_strength)
@validate_output(current_density)
def calculate_current_density(electrons_concentration_: Quantity, holes_concentration_: Quantity,
    electrons_mobility_: Quantity, holes_mobility_: Quantity,
    electric_intensity_: Quantity) -> Quantity:
    result_expr = solve(law, current_density, dict=True)[0][current_density]
    result_expr = result_expr.subs({
        electrons_concentration: electrons_concentration_,
        holes_concentration: holes_concentration_,
        electrons_mobility: electrons_mobility_,
        holes_mobility: holes_mobility_,
        electric_field_strength: electric_intensity_
    })
    return Quantity(result_expr)
