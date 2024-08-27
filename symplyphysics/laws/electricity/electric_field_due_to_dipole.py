r"""
Electric field due to dipole
============================

The value of the electric field set up by a dipole at a distant point on the dipole axis,
which runs through both point charges comprising the dipole, is proportional to the inverse
cube of the distance to the dipole and the value of the dipole moment.

**Notation:**

#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity.

**Conditions:**

#. :math:`r / d \gg 1` where :math:`d` is the distance between point charges of the
   dipole. This means that the point where electric field is measured must be far enough
   from the dipole itself.
"""

from sympy import Eq, solve, series, pi
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.electricity import electric_field_due_to_point_charge as point_field
from symplyphysics.laws.electricity import electric_dipole_moment_is_charge_times_distance as dipole_law

electric_field = Symbol("electric_field", units.force / units.charge)
"""
Electric field due to dipole.

Symbol:
    :code:`E`
"""

electric_dipole_moment = Symbol("electric_dipole_moment", units.charge * units.length)
"""
:doc:`laws.electricity.electric_dipole_moment_is_charge_times_distance`.

Symbol:
    :code:`p`
"""

distance = Symbol("distance", units.length)
"""
Distance to dipole.

Symbol:
    :code:`r`
"""

law = Eq(electric_field, 1 / (2 * pi * units.vacuum_permittivity) * electric_dipole_moment / distance**3)
r"""
:code:`E = 1 / (2 * pi * epsilon_0) * p / r^3`

Latex:
    .. math::
        E = \frac{1}{2 \pi \varepsilon_0} \frac{p}{r^3}
"""

# Derive the law from the expression for the electric field of point charges.
# Assuming the dipole is made up of two point charges q and -q (q > 0) with distance d in between.
# Assume the z-axis running through both charges, and let its origin be at their middle point,
# the negative charge located below and the positive one above the origin.

_charge = Symbol("_charge", units.charge)
_distance_between_charges = Symbol("_distance_between_charges", units.length)
_distance_to_origin = Symbol("_distance_to_origin", units.length)

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
_factor = Symbol("_factor", dimensionless)

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
        Eq(electric_field, _net_field_approx_sub),
        _dipole_eqn
    ],
    (_charge, electric_field),
    dict=True,
)[0][electric_field]

_net_field_from_law = law.rhs.subs(distance, _distance_to_origin)

assert expr_equals(_net_field_from_law, _net_field_derived)


@validate_input(dipole_moment_=electric_dipole_moment, distance_to_dipole_=distance)
@validate_output(electric_field)
def calculate_electric_field(dipole_moment_: Quantity, distance_to_dipole_: Quantity) -> Quantity:
    result = solve(law, electric_field)[0]
    result_field = result.subs({
        electric_dipole_moment: dipole_moment_,
        distance: distance_to_dipole_,
    })
    return Quantity(result_field)
