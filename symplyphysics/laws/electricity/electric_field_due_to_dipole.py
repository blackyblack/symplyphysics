r"""
Electric field due to dipole
============================

The strength of the electric field set up by a dipole at a distant point on the dipole axis,
which runs through both point charges comprising the dipole, is proportional to the inverse
cube of the distance to the dipole and the value of the electric dipole moment.

**Notation:**

#. :quantity_notation:`vacuum_permittivity`.

**Conditions:**

#. :math:`r / d \gg 1` where :math:`d` is the distance between point charges of the
   dipole. This means that the point where electric field is measured must be far enough
   from the dipole itself.
"""

from sympy import Eq, solve, series, pi, Symbol as SymSymbol
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    quantities,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as point_field
from symplyphysics.laws.electricity import electric_dipole_moment_is_charge_times_distance as dipole_law

electric_field_strength = symbols.electric_field_strength
"""
:symbols:`electric_field_strength` due to dipole.
"""

electric_dipole_moment = symbols.electric_dipole_moment
"""
See :doc:`laws.electricity.electric_dipole_moment_is_charge_times_distance` and
:symbols:`electric_dipole_moment`.
"""

distance = symbols.distance_to_origin
"""
:symbols:`distance` to dipole.
"""

law = Eq(electric_field_strength, 1 / (2 * pi * quantities.vacuum_permittivity) * electric_dipole_moment / distance**3)
"""
:laws:symbol::

:laws:latex::
"""

# Derive the law from the expression for the electric field of point charges.
# Assuming the dipole is made up of two point charges q and -q (q > 0) with distance d in between.
# Assume the z-axis running through both charges, and let its origin be at their middle point,
# the negative charge located below and the positive one above the origin.

_charge = SymSymbol("_charge")
_distance_between_charges = SymSymbol("_distance_between_charges")
_distance_to_origin = SymSymbol("_distance_to_origin")

_positive_charge_field = point_field.law.rhs.subs({
    point_field.charge: _charge,
    point_field.distance: _distance_to_origin - _distance_between_charges / 2
})
_negative_charge_field = point_field.law.rhs.subs({
    point_field.charge: -_charge,
    point_field.distance: _distance_to_origin + _distance_between_charges / 2
})

# The net electric field is the sum of electric fields due to both charges
_net_field = _positive_charge_field + _negative_charge_field

# The condition that _distance_to_origin/_distance_between_charges >> 1 can be analyzed as such:
# Let _distance_between_charges = _factor * _distance_to_origin, where _factor -> 0.
_factor = SymSymbol("_factor")

# Use the above definition of _factor to substitute _distance_to_origin in the _net_field formula
_net_field_sub = _net_field.subs(_distance_between_charges, _factor * _distance_to_origin)

# Expand _net_field_sub with respect to _factor around 0 up to the first power
# in order to find the first approximation of the current formula.
_net_field_approx = series(_net_field_sub, _factor, 0, 2).removeO()

# Substitute _distance_to_origin back into the _net_field formula
_net_field_approx_sub = _net_field_approx.subs(_factor, _distance_between_charges / _distance_to_origin)

_dipole_eqn = dipole_law.law.subs({
    dipole_law.electric_dipole_moment: electric_dipole_moment,
    dipole_law.charge: _charge,
    dipole_law.distance: _distance_between_charges,
})

# Replace _charge*_distance_between_charges back with electric_dipole_moment
_net_field_derived = solve(
    [
        Eq(electric_field_strength, _net_field_approx_sub),
        _dipole_eqn
    ],
    (_charge, electric_field_strength),
    dict=True,
)[0][electric_field_strength]

_net_field_from_law = law.rhs.subs(distance, _distance_to_origin)

assert expr_equals(_net_field_from_law, _net_field_derived)


@validate_input(dipole_moment_=electric_dipole_moment, distance_to_dipole_=distance)
@validate_output(electric_field_strength)
def calculate_electric_field(dipole_moment_: Quantity, distance_to_dipole_: Quantity) -> Quantity:
    result = solve(law, electric_field_strength)[0]
    result_field = result.subs({
        electric_dipole_moment: dipole_moment_,
        distance: distance_to_dipole_,
    })
    return Quantity(result_field)
